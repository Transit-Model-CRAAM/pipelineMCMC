from scipy import interpolate
import emcee
from Planet.Planeta import Planeta
from Star.Estrela import Estrela #estrela e eclipse:: extensões de programas auxiliares que realizam o cálculo da curva de luz.
from Planet.Eclipse import Eclipse
from Misc.Verify import converte
import numpy

class Ajuste:
    
    def __init__(self,tratamento, time, flux, nwalkers, niter, burnin, rsun = 1, periodo = 1):
        """
        Attributes
        ----------
        gc_trigger_counter : int
            Contador de iterações usado para verificar quando chamar o garbage collector.
        """
        self.tratamento = tratamento

        self.u1_p0 = self.tratamento.u1
        self.u2_p0 = self.tratamento.u2
        self.a_p0 = self.tratamento.AU
        self.inc_p0 = self.tratamento.inc
        self.rp_p0 = self.tratamento.modelo.R_jup
        self.rsun = rsun
        self.periodo = periodo

        self.time = time
        self.flux = flux

        self.flux_err = numpy.var(self.flux)
        self.data = (self.time, self.flux, self.flux_err)

        self.nwalkers = nwalkers
        self.niter = niter
        self.burnin = burnin

        self.initial = numpy.array([self.u1_p0, self.u2_p0, self.a_p0, self.inc_p0, self.rp_p0])
        self.ndim = len(self.initial)

        variations = numpy.array([0.001, 0.001, 0.001, 0.5, 0.01])

        self.p0 = [numpy.array(self.initial) + variations * numpy.random.randn(self.ndim) for i in range(self.nwalkers)]

        self.gc_trigger_counter = 0

    #--------------------------------------------------#
    #----------------------MCMC------------------------#
    #--------------------------------------------------#
    def eclipse_mcmc(self, time, theta):
        u1, u2, semiEixoUA, anguloInclinacao, raioPlanJup = theta

        raioStar, raioPlanetaRstar, semiEixoRaioStar = converte(self.rsun,raioPlanJup,semiEixoUA)
        
        estrela_ = Estrela(373, self.rsun, 240., u1, u2, 856)
        Nx = estrela_.getNx()
        Ny = estrela_.getNy()
        raioEstrelaPixel = estrela_.getRaioStar()
        
        # semiEixoUA, raioPlanJup, periodo, anguloInclinacao, ecc, anom, raioStar,mass): 
        planeta_ = Planeta(semiEixoUA, raioPlanJup, self.periodo, anguloInclinacao, 0, 0, estrela_.getRaioSun(), self.tratamento.mass, self.tratamento.planet_name)
        
        # Nx, Ny, raio_estrela_pixel, estrela_manchada: Estrela, planeta_: Planeta
        eclipse = Eclipse(Nx,Ny,raioEstrelaPixel,estrela_, planeta_)
        
        eclipse.setTempoHoras(1.)

        self.gc_trigger_counter += 1
        # Coleta o lixo apenas a cada ciclo rodado 
        # O garbage collector custa cerca de 0.08 segundo pra rodar, o que é pouco se você rodar apenas uma vez,
        # mas como o MCMC pode rodar milhares de iterações, pode acabar aumentando bastante o custo de execução.
        if self.gc_trigger_counter % self.niter == 0:
            eclipse.criarEclipse(anim = False, plot= False)
        else:
            eclipse.criarEclipse(anim = False, plot= False, collect_garbage=False)

        lc0 = numpy.array(eclipse.getCurvaLuz(), copy=True) 
        ts0 = numpy.array(eclipse.getTempoHoras(), copy=True) 

        del eclipse, planeta_, estrela_

        return interpolate.interp1d(ts0,lc0,fill_value="extrapolate")(time)
        
    #--------------------------------------------------#
    def lnlike(self, theta, time, flux, flux_err):
        return -0.5 * numpy.sum(((flux - self.eclipse_mcmc(time, theta))/flux_err) ** 2)
    #--------------------------------------------------#
    def lnprior(self, theta):
        u1, u2, semiEixoUA, anguloInclinacao, rp = theta
        if u1 > 0.0:
            if u1 + 2*u2 > 0.0 and u1 + u2 < 1 and 0.001 < semiEixoUA < 1 and 80. < anguloInclinacao < 90 and 0.01 < rp < 5:
                return 0.0
        elif u1 < 0.0:
            if u1 + 2*u2 < 0.0 and u1 + u2 > -1 and 0.001 < semiEixoUA < 1 and 80. < anguloInclinacao < 90 and 0.01 < rp < 5:
                return 0.0
        return -numpy.inf

    #--------------------------------------------------#
    def lnprob(self, theta, time, flux, flux_err):
        lp = self.lnprior(theta)
        if not numpy.isfinite(lp):
            return -numpy.inf
        return lp + self.lnlike(theta, time, flux, flux_err)
    #--------------------------------------------------#
    def main(self):
        self.sampler = emcee.EnsembleSampler(self.nwalkers, self.ndim, self.lnprob, args=self.data)

        print("Running burn-in...")
        self.p0, _, _ = self.sampler.run_mcmc(self.p0, self.burnin, progress=True)
        self.sampler.reset()

        print("Running production...")    
        self.pos, self.prob, self.state = self.sampler.run_mcmc(self.p0, self.niter, progress=True)

        return self.sampler, self.pos, self.prob, self.state
    #--------------------------------------------------#

class AjusteManchado: 
    def __init__(self,tratamento, time, flux, nwalkers, niter, burnin, ndim, eclipse: Eclipse, rsun = 1, periodo = 1):
        """
        Attributes
        ----------
        gc_trigger_counter : int
            Contador de iterações usado para verificar quando chamar o garbage collector.
        """
        
        self.manchas: Estrela.Mancha = eclipse.estrela_.manchas

        self.u1 = eclipse.estrela_.coeficienteHum
        self.u2 = eclipse.estrela_.coeficienteDois
        self.semiEixoUA = eclipse.planeta_.semiEixoUA
        self.anguloInclinacao = eclipse.planeta_.anguloInclinacao
        self.raioPlanJup = eclipse.planeta_.raioPlanJup

        self.rsun = rsun
        self.periodo = periodo

        self.time = time
        self.flux = flux

        self.flux_err = numpy.var(self.flux)
        self.data = (self.time, self.flux, self.flux_err)

        self.nwalkers = nwalkers
        self.niter = niter
        self.burnin = burnin


        self.initial = numpy.array([])

        # limitacao do numero de manchas
        if(ndim > len(self.manchas)):
            ndim = len(self.manchas)
        elif(ndim < 1):
            ndim = 1
            
        for i in range(ndim):
            self.initial = numpy.append(self.initial, [self.manchas[i].longitude, self.manchas[i].latitude, self.manchas[i].raio, self.manchas[i].intensidade])

        variations = numpy.array([0.8, 0.8, 0.001, 0.01])

        ndim_variations = numpy.tile(variations, ndim)
        
        self.ndim = len(self.initial)
        self.p0 = [numpy.array(self.initial) + ndim_variations * numpy.random.randn(self.ndim) for i in range(self.nwalkers)]
        self.tratamento = tratamento

        self.gc_trigger_counter = 0

    #--------------------------------------------------#
    #----------------------MCMC------------------------#
    #--------------------------------------------------#
    def eclipse_mcmc(self, time, theta):
        raioStar, raioPlanetaRstar, semiEixoRaioStar = converte(self.rsun,self.raioPlanJup,self.semiEixoUA)
        
        estrela_ = Estrela(373, self.rsun, 240., self.u1, self.u2, 856)
        Nx = estrela_.getNx()
        Ny = estrela_.getNy()
        raioEstrelaPixel = estrela_.getRaioStar()
        
        
        for i in range(len(theta)//4):
            # TO-DO: Mudar essa função que esta sendo chamada aqui 
            # intensidade, raio, latitude, longitude
            raioRStar = theta[(i*4)+2]
            intensidade = theta[(i*4)+3]
            lat = theta[(i*4)+1]
            long = theta[i*4]

            estrela_.addMancha(Estrela.Mancha(intensidade, raioRStar, lat, long))
        
        estrela_.criaEstrelaManchada()
        # semiEixoUA, raioPlanJup, periodo, anguloInclinacao, ecc, anom, raioStar,mass): 
        planeta_ = Planeta(self.semiEixoUA, self.raioPlanJup, self.periodo, self.anguloInclinacao, 0, 0, estrela_.getRaioSun(), self.tratamento.mass, self.tratamento.planet_name)
        
        # Nx, Ny, raio_estrela_pixel, estrela_manchada: Estrela, planeta_: Planeta
        eclipse = Eclipse(Nx,Ny,raioEstrelaPixel,estrela_, planeta_)

        eclipse.setTempoHoras(1.)

        self.gc_trigger_counter += 1
        # Coleta o lixo apenas a cada ciclo rodado 
        # O garbage collector custa cerca de 0.08 segundo pra rodar, o que é pouco se você rodar apenas uma vez,
        # mas como o MCMC pode rodar milhares de iterações, pode acabar aumentando bastante o custo de execução.
        if self.gc_trigger_counter % self.niter == 0:
            eclipse.criarEclipse(anim = False, plot= False)
        else:
            eclipse.criarEclipse(anim = False, plot= False, collect_garbage=False)
        
        lc0 = numpy.array(eclipse.getCurvaLuz())
        ts0 = numpy.array(eclipse.getTempoHoras()) 
        
        del eclipse, planeta_, estrela_

        return interpolate.interp1d(ts0,lc0,fill_value="extrapolate")(time)
    #--------------------------------------------------#
    def lnlike(self, theta, time, flux, flux_err):
        return -0.5 * numpy.sum(((flux - self.eclipse_mcmc(time, theta))/flux_err) ** 2)
    #--------------------------------------------------#
    def lnprior(self, theta):
        for i in range(len(theta)//4):
            #if (0.0 < lat) and (0.0 < long) and (0.0 < raioRstar < 0.5) and (0.0 < intensidade <= 1):
            if (-70 <= theta[i*4] <= 70) and (-70 <= theta[(i*4)+1] <= 70) and (0.0 < theta[(i*4)+2] < 0.5) and (0.0 < theta[(i*4)+3] <= 1):
                continue
            return -numpy.inf
        return 0.0
    #--------------------------------------------------#
    def lnprob(self, theta, time, flux, flux_err):
        lp = self.lnprior(theta)
        if not numpy.isfinite(lp):
            return -numpy.inf
        return lp + self.lnlike(theta, time, flux, flux_err)
    #--------------------------------------------------#
    def main(self):
        self.sampler = emcee.EnsembleSampler(self.nwalkers, self.ndim, self.lnprob, args=self.data)

        print("Running burn-in...")
        self.p0, _, _ = self.sampler.run_mcmc(self.p0, self.burnin, progress=True)
        self.sampler.reset()

        print("Running production...")    
        self.pos, self.prob, self.state = self.sampler.run_mcmc(self.p0, self.niter, progress=True)

        return self.sampler, self.pos, self.prob, self.state
    #--------------------------------------------------#