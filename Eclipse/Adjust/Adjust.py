from scipy import interpolate
import emcee
from Star.Estrela import Estrela #estrela e eclipse:: extensões de programas auxiliares que realizam o cálculo da curva de luz.
from Planet.Eclipse import Eclipse
from Misc.Verify import converte
import numpy

class Ajuste:
    
    def __init__(self,tratamento, time, flux, nwalkers, niter, burnin, rsun = 1, periodo = 1):

        self.u1_p0 = 0.5
        self.u2_p0 = 0.1
        self.a_p0 = 0.05
        self.inc_p0 = 88.
        self.rp_p0 = 1
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

        self.p0 = [numpy.array(self.initial) + 1e-4 * numpy.random.randn(self.ndim) for i in range(self.nwalkers)]

        self.tratamento = tratamento

    #--------------------------------------------------#
    #----------------------MCMC------------------------#
    #--------------------------------------------------#
    def eclipse_mcmc(self, time, theta):
        u1, u2, semiEixoUA, anguloInclinacao, raioPlanJup = theta
        raioStar, raioPlanetaRstar, semiEixoRaioStar = converte(self.rsun,raioPlanJup,semiEixoUA)
        
        estrela_ = Estrela(373, raioStar, 240., u1, u2, 856)
        Nx = estrela_.getNx()
        Ny = estrela_.getNy()
        raioEstrelaPixel = estrela_.getRaioStar()
        estrelaManchada = estrela_.getEstrela()
        
        eclipse = Eclipse(Nx,Ny,raioEstrelaPixel,estrelaManchada)
        eclipse.setTempoHoras(1.)
        eclipse.criarEclipse(semiEixoRaioStar, semiEixoUA, raioPlanetaRstar, raioPlanJup, self.periodo, anguloInclinacao, 0,0 , 0, False, False)
        lc0 = numpy.array(eclipse.getCurvaLuz()) 
        ts0 = numpy.array(eclipse.getTempoHoras()) 
        return interpolate.interp1d(ts0,lc0,fill_value="extrapolate")(time)
        
    #--------------------------------------------------#
    def lnlike(self, theta, time, flux, flux_err):
        return -0.5 * numpy.sum(((flux - self.eclipse_mcmc(time, theta))/flux_err) ** 2)
    #--------------------------------------------------#
    def lnprior(self, theta):
        u1, u2, semiEixoUA, anguloInclinacao, rp = theta
        if 0.0 < u1 < 1.0 and 0.0 < u2 < 1.0 and 0.001 < semiEixoUA < 1 and 80. < anguloInclinacao < 90 and 0.01 < rp < 5:
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
    def __init__(self,tratamento, time, flux, nwalkers, niter, burnin, ndim, u1, u2, semiEixoUA, anguloInclinacao, raioPlanJup, rsun = 1, periodo = 1):
        
        # parametros das 4 manchas (limite) 
        # TODO: Arrumar multiplas manchas
        self.lat = [-10, -10, 0, 0]
        self.long = [10, 40, 0, 0] 
        self.raioRStar = [0.1, 0.1, 0.1, 0.1]  # em relacao ao raio da estrela 
        self.intensidade = [0.3, 0.3, 0.5, 0.5]  # em relacao ao raio da estrela

        self.u1 = u1
        self.u2 = u2
        self.semiEixoUA = semiEixoUA
        self.anguloInclinacao = anguloInclinacao
        self.raioPlanJup = raioPlanJup
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
        if(ndim > 4):
            ndim = 4
        elif(ndim < 1):
            ndim = 1
            
        for i in range(ndim):
            self.initial = numpy.append(self.initial, [self.long[i], self.lat[i], self.raioRStar[i], self.intensidade[i]])
        
        self.ndim = len(self.initial)
        self.p0 = [numpy.array(self.initial) + 1e-4 * numpy.random.randn(self.ndim) for i in range(self.nwalkers)]
        self.tratamento = tratamento

    #--------------------------------------------------#
    #----------------------MCMC------------------------#
    #--------------------------------------------------#
    def eclipse_mcmc(self, time, theta):
        raioStar, raioPlanetaRstar, semiEixoRaioStar = converte(self.rsun,self.raioPlanJup,self.semiEixoUA)
        
        estrela_ = Estrela(373, raioStar, 240., self.u1, self.u2, 856)
        Nx = estrela_.getNx()
        Ny = estrela_.getNy()
        raioEstrelaPixel = estrela_.getRaioStar()
        
        
        for i in range(len(theta)//4):
            # manchas(raioRStar, intensidade, lat, long)
            estrela=estrela_.manchas(theta[(i*4)+2],theta[(i*4)+3],theta[(i*4)+1],theta[i*4])
        # estrela=estrela_.manchas(r,intensidadeMancha,lat,longt) #recebe a escolha de se irá receber manchas ou não
        # count+=1
        
        estrelaManchada = estrela_.getEstrela()
        eclipse = Eclipse(Nx,Ny,raioEstrelaPixel,estrelaManchada)
        eclipse.setTempoHoras(1.)
        eclipse.criarEclipse(semiEixoRaioStar, self.semiEixoUA, raioPlanetaRstar, self.raioPlanJup, self.periodo, self.anguloInclinacao, 0,0 , 0, False, False)
        lc0 = numpy.array(eclipse.getCurvaLuz())
        ts0 = numpy.array(eclipse.getTempoHoras()) 
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
