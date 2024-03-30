__author__ = "Yuri Neto, Felipe Pinho, Beatriz Duque, Adriana Valio (Orientadora)"
__copyright__ = "..."
__credits__ = ["Universidade Presbiteriana Mackenzie, CRAAM"]
__license__ = ""
__version__ = ""
__maintainer__ = ""
__email__ = "starsandexoplanets@gmail.com"
__status__ = "Production"

#%matplotlib nbagg
#%matplotlib inline

import scipy
from Star.Estrela import Estrela #estrela e eclipse:: extensões de programas auxiliares que realizam o cálculo da curva de luz.
from Planet.Eclipse import Eclipse
import numpy
import matplotlib.pyplot as plt
from lightkurve import search_lightcurve


class Modelo:

    '''
    parâmetro estrela :: classe estrela 
     parâmetro eclipse :: classe Eclipse
    parâmetro missão :: Missão selecionada para a coleta da estrela (KEPLER, K2 OU TESS)
    '''
    def __init__(self,estrela, eclipse, mission):
       
        #coletando objetos de estrela 
        self.u1 = estrela.getu1()
        self.u2 = estrela.getu2()
        self.n = estrela.getTamanhoMatriz()
        self.r = estrela.getRaioStar() #raio da estrela em pixel 
        self.r_Sun = estrela.getRaioSun() #raio da estrela em RSun
        self.mx = estrela.getIntensidadeMaxima()
        self.star_name = estrela.getStarName() #nome da estrela
        self.cadence = estrela.getCadence() #cadencia da estrela (short ou long)
        self.mission = mission

        #coletando objetos de Eclipse
        self.raioPlan = eclipse.getRaioPlan()
        self.R_jup = eclipse.getRplanJup()
        self.AU = eclipse.getSemiEixo() #semieixo em UA
        self.semiEixoRaioStar = eclipse.getsemiEixoRaioStar() #semieixo em relacao ao raio da estrela
        self.porb = eclipse.getPeriodo()
        self.inc = eclipse.getInc()
        self.ecc,self.anom = eclipse.getEccAnom()
        self.lua = eclipse.getLua()

        #variaveis que serao retornadas a partir da função rd_data
        self.time = 0.
        self.flux = 0.
        self.flux_err = 0.


    def rd_data(self,plot,save_data):
                
        ''' 
        Funcao criada para acessar os dados de curva de luz de estrelas e extrair o tempo e fluxo. 
        
        lightkurve packeage:: utilizado para abrir as curvas de Luz do Kepler e
        também abrir curva de luz do TESS/ curvas de luz em geral
        Documentação:: https://docs.lightkurve.org/api/index.html documentação do uso do lightkurve
        
        outputs esperados:: time, flux 
        
        plot = 1 (plota a curva de luz, para não plotar basta digitar qualquer valor)
        save_data = 1 plota a curva de luz(para não plotar, basta digitar qualquer valor)
        '''
    ##--------------------------------------------------------------------------------------------------------------------------------------------------##
    
        # utiiza-se o PDCSAP_FLUX porque será realizado a análise no trânsito.
        lc = search_lightcurve(self.star_name, cadence = self.cadence, mission=self.mission).download_all()
        time = [] # time = array com os dados de tempo
        flux = [] # flux = array com os dados de fluxo
        flux_err = [] # flux_err = array com os dados de erro do fluxo
        time_temp = [] #_variavel temporaria
        flux_temp = []
        flux_err_temp = []

        for i in range(0, len(lc)):
            try:
                flux_temp.append(lc[i].sap_flux)
                flux_err_temp.append(lc[i].sap_flux_err)
                time_temp.append(lc[i].time.to_value('bkjd','float'))
            except:
                pass

        for i in range(0, len(flux_temp)):   
            
            # Limita os índices aos valores de tempo lidos
            flux_temp[i] = flux_temp[i][0:time_temp[i].size]
            flux_err_temp[i] = flux_err_temp[i][0:time_temp[i].size]
        
            # Elimina todos os valores NaN (Not a Number) dos dados
            time_temp[i] = time_temp[i][~numpy.isnan(flux_err_temp[i])]
            flux_temp[i] = flux_temp[i][~numpy.isnan(flux_err_temp[i])]
            flux_err_temp[i] = flux_err_temp[i][~numpy.isnan(flux_err_temp[i])]
        
            # Normaliza cada quarter
            flux_err_temp[i] = flux_err_temp[i]/ abs(numpy.median(flux_temp[i]))
            flux_temp[i] = flux_temp[i]/ abs(numpy.median(flux_temp[i]))

        for i in range(0, len(flux_temp)):
            flux = numpy.append(flux, flux_temp[i])
            time = numpy.append(time, time_temp[i])
            flux_err = numpy.append(flux_err,flux_err_temp[i])
        
        #plot da curva de luz completa
        if plot == (1):
            plt.rcParams['figure.figsize'] = 10,4
            graf1, ax = plt.subplots()

            ax.set_xlabel('Time (BJD - 2454833)')
            ax.set_ylabel('Normalized Flux')

            ax.set_title("Light Curve - " + self.star_name)
            ax.set_xlim(min(time), max(time))
            ax.set_ylim(min(flux), max(flux))
                    
            ax.errorbar(time, flux, yerr = flux_err, fmt = '.k', capsize = 0,alpha = 0.5)
            
        #salva os dados em um arquivo .dat
        if save_data == 1:
            numpy.savetxt('%s_LC.dat'%self.star_name, numpy.c_[(time, flux, flux_err)])

        self.time = time 
        self.flux = flux
        self.flux_err = flux_err

        return self.time, self.flux, self.flux_err

    def det_x0(self, plot):
        
        '''
        Função para obter o centro do primeiro transito. 
        porb é utilizado para criar uma curva de luz em fase é aplicado um smooth na curva de luz em fase.
        parâmetro time::
        parâmetro flux:: 
        parâmetro porb:: periodo orbital do planeta (per em dias)
        parâmetro plot:: 
        
        returns
        x0 = valor do centro do primeiro transito 
        nt = numero de transitos possiveis 
        plot = 1 (plota a curva de luz, para nao plotar basta digitar qualquer valor)
        '''
        time = self.time 
        porb = self.porb
        flux = self.flux

        phase = (time % porb)/ porb
        jj = numpy.argsort(phase)
        ff = phase[jj]

        smoothed_LC = scipy.ndimage.filters.uniform_filter(flux[jj], size = 100) # equivalente ao smooth do idl com edge_truncade
        smoothed_LC[0:200] = 1
        smoothed_LC[len(flux[jj])-200:len(flux[jj])] = 1
        x = phase[jj]
        y = 1 - smoothed_LC
        yh = 0.002

        kk = numpy.where(y >= yh)

        x1 = min(x[kk])
        x2 = max(x[kk])
        fa0 = (x1 + x2)/ 2 # valor central dos transitos em fase

        self.x0 = (numpy.fix(time[0] / porb) + fa0) * porb # tempo central do primeiro transito
        
        self.nt = numpy.fix(max(time - time[0])/ porb) + 1 # numero de transitos possiveis

        #plot da curva de luz completa
        if plot == 1:
        
            plt.rcParams['figure.figsize'] = 10,8
            graf2, ax = plt.subplots(2, 1, sharex = False, gridspec_kw = {'height_ratios': [1, 1]})
            graf2.subplots_adjust(hspace = 0.5)
        
            ax[0].set_ylim(0.9, 1.1)
            ax[0].set_title("Phased LC")
        
            ax[0].set_xlabel('Time')
            ax[0].set_ylabel('Normalized Flux')

            ax[0].plot(phase, flux, "k.", ms = 2)
            ax[0].plot(phase[jj], smoothed_LC, "r.", ms = 2)        
        
            ax[1].set_ylim(0.9, 1.1)
            ax[1].set_title("x0 = tempo central do primeiro transito")
        
            ax[1].set_xlabel('Time (BJD - 2454833)')
            ax[1].set_ylabel('Normalized Flux')
        
            ax[1].set_xlim(min(time) - 2, min(time) + 2)
            ax[1].set_ylim(0.95, 1.01)
        
            ax[1].plot(time, flux, "k.", ms=2)

            plt.axvline(x = self.x0,linewidth = 2, color='r')
        
        return self.x0, self.nt

    #--------------------------------------------------#
    def limb(self, plot):
        '''
        Funcao que gera uma estrela sintetizada com 
        obscurecimento de limbo dado por 1-u1*(1-cos(mu))-u2*(1-cos(mu))^2), onde mu=angulo heliocentrico
        
        -- coeficiente de escurecimento de limbo --
        parâmetro u1 :: coeficiente linear 
        parâmetro u2 :: coeficiente do termo quadratico
        
        returns 
        parâmetro wl:: matriz com a intensidade da estrela
        se plot = 1, plota o perfil da estrela (para não plotar, basta digitar qualquer valor)

        '''
        #PERGUNTAR

        if self.u1 == 999:
            self.u1 = 0.59
        if self.u2 == 999:
            self.u2 = 0.0
        
        self.wl = numpy.zeros((self.n, self.n))
        x = numpy.arange(self.n/2-self.n, self.n - self.n/2, 1)
        
        for j in range (0, self.n):     
            z = numpy.sqrt(x**2 + (j - self.n/2.)**2)
            kk = numpy.where(z <= self.r)
            if kk[0].size > 0:
                    m = numpy.cos(numpy.arcsin(z[kk]/self.r))
                    self.wl[kk, j] = self.mx*(1 - self.u1*(1 - m) - self.u2*(1 - m)**2)
        
        data = self.wl[:,int(self.n/2)]
        smoothed = scipy.ndimage.filters.uniform_filter(data, size = 5)

        if plot == 1:
            plt.rcParams['figure.figsize']= 12, 5
            graf1,ax = plt.subplots()
            ax.plot(x,smoothed)
            ax.plot(x, self.wl[:, 428],'k')
            ax.set_title('limb')


        return self.wl, self.u1, self.u2 
    #--------------------------------------------------#

    def eclipse_model(self):
        '''
        Chamada de programas auxiliares para a criacao do modelo da curva de luz, podendo ela conter:
        - Um planeta ou mais 
        - Uma mancha ou mais 
        - Uma lua ou mais 
        
        parâmetro u1 :: coeficiente de escurecimento de limbo 1
        parâmetro u2 :: coeficiente de escurecimento de limbo 1
        parâmetro per :: periodo do transito em dias 
        parâmetro a :: semieixo em UA 
        parâmetro inc :: ângulo de inclinacao em graus 
        parâmetro rp :: raio do planeta em relacao ao raio da estrela 
        
        returns 
        parâmetro lc_model :: curva de luz 
        parâmetro ts_model :: tempo do trânsito em Horas
        
        '''
        estrela_1 = Estrela(self.r,self.r_Sun, self.mx , self.u1, self.u2, self.n)  #cria o objeto estrela 
        Nx1 = estrela_1.getNx() #coleta parametros da matriz estrela 
        Ny1 = estrela_1.getNy()
        raioEstrelaPixel1 = estrela_1.getRaioStar() #coleta raio da estrela em pixel 
        estrelaManchada1 = estrela_1.getEstrela() #coleta estrela manchada 

        eclipse1 = Eclipse(Nx1, Ny1, raioEstrelaPixel1, estrelaManchada1)  #cria o objeto eclipse

        eclipse1.setTempoHoras(1.)

        eclipse1.criarEclipse(self.semiEixoRaioStar, self.AU, self.raioPlan,self.R_jup, self.porb,self.inc,self.lua,self.ecc,self.anom, False)

        self.lc_model = numpy.array(eclipse1.getCurvaLuz())
        self.ts_model = numpy.array(eclipse1.getTempoHoras())
        
        return self.lc_model, self.ts_model

    #--------------------------------------------------#
    def retornaParametros(self):
        return self.u1,self.u2,self.porb,self.time,self.flux,self.flux_err,self.raioPlan,self.AU,self.inc, self.x0, self.nt,self.ts_model
    
    def setTime(self,time):
        self.time = time

    def setFlux(self,flux):
        self.flux = flux
    
    def setFluxErr(self,flux_err):
        self.flux_err = flux_err
