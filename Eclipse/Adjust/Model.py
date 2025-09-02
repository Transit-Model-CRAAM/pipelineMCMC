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

from typing import List
import scipy
from Planet.Planeta import Planeta
from Star.Estrela import Estrela #estrela e eclipse:: extensões de programas auxiliares que realizam o cálculo da curva de luz.
from Planet.Eclipse import Eclipse
import numpy as np
import matplotlib.pyplot as plt
from lightkurve import search_lightcurve, KeplerLightCurve
import pandas as pd 
import os
import pyvo


class Modelo:
    '''
    parâmetro estrela :: classe estrela 
    parâmetro eclipse :: classe Eclipse
    parâmetro missão :: Missão selecionada para a coleta da estrela (KEPLER, K2 OU TESS)
    '''
    def __init__(self,estrela: Estrela, eclipse: Eclipse, mission):
       
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

        self.manchas: List[Estrela.Mancha] = estrela.manchas

        #coletando objetos de Eclipse
        self.planet_name = eclipse.planeta_.getPlanetName()
        self.raioPlan = eclipse.planeta_.getRaioPlan()
        self.R_jup = eclipse.planeta_.getRplanJup()
        self.AU = eclipse.planeta_.getSemiEixo() #semieixo em UA
        self.semiEixoRaioStar = eclipse.planeta_.getsemiEixoRaioStar() #semieixo em relacao ao raio da estrela
        self.porb = eclipse.planeta_.getPeriodo()
        self.inc = eclipse.planeta_.getInc()
        self.mass = eclipse.planeta_.mass
        self.ecc,self.anom = eclipse.planeta_.getEccAnom()

        if eclipse.planeta_.hasMoons(): 
            self.lua = eclipse.planeta_.luas[0] # por enquanto apenas uma lua 

        #variaveis que serao retornadas a partir da função rd_data
        self.lightcurve_collection = []

    def rd_data(self):
        '''
        Function responsible for accessing the data from the star's light curve and extract the total light curve.

        Returns a lightkurve.lightcurve.KeplerLightCurve object in case a curve is found; None if an error occurs.
        -----------------------------------------------------------------------------------------------------------
        Função responsável por acessar os dados da curva de luz da estrela e extrair a curva de luz total.

        Retorna um objeto lightkurve.lightcurve.KeplerLightCurve caso ache a curva; None caso ocorra algum erro.
        '''
        try:
            lc_collection = search_lightcurve(self.star_name, cadence = self.cadence, mission = self.mission).download_all()

            # **Filter data with flag = 0**
            # Iterate through the collection and filter each light curve
            filtered_lcs = []
            for light_curve in lc_collection:
                filtered_lcs.append(light_curve[light_curve.quality == 0])

            # Normalize each filtered light curve before stitching
            normalized_lcs = []
            for light_curve in filtered_lcs:
                # Normalize each light curve individually
                normalized_lcs.append(light_curve.normalize())

            # Stitch the normalized light curves into a single LightCurve object
            if normalized_lcs:
                stitched_lc = normalized_lcs[0]
                for i in range(1, len(normalized_lcs)):
                    # Append the normalized light curves
                    stitched_lc = stitched_lc.append(normalized_lcs[i])

            self.lightcurve_collection = stitched_lc
            return stitched_lc
        except Exception as e:
            print(f"Error: {e}")
            return None
        
    def get_transit_parameters(self):
        '''
        Função responsável por obter os parâmetros necessários para cortar os trânsitos:

        Period: Periodo do trânsito em dias
        Epoch: Época do primeiro trânsito em BJD
        Duration: Duração do trânsito em dias
        '''
        # Access to ExoplanetArchive API
        service = pyvo.dal.TAPService("https://exoplanetarchive.ipac.caltech.edu/TAP")

        # Se quiser procurar por um planeta:
        target_planet = self.planet_name
        query = f"SELECT * FROM pscomppars WHERE LOWER(pl_name) LIKE '%{target_planet.lower()}%'"

        results = service.search(query)

        # Verifying if there are results
        if len(results) == 0:
            return None
        
        self.period = results[0].get("pl_orbper")         # dias
        self.epoch = results[0].get("pl_tranmid")         # BJD
        self.duration = results[0].get("pl_trandur")/24   # dias

        return self.period, self.epoch, self.duration

    def get_transits(self):
        '''
        Função responsável por dividir os trânsitos.

        Retorna uma lista dos trânsitos divididos.
        '''
        self.get_transit_parameters()

        end_time = 0
        nn = 0
        self.transit_list = []
        while(end_time < self.lightcurve_collection[-1]["time"].jd):
            # Plot the light curve around the nth transit
            start_time = self.epoch - 2*self.duration + nn * self.period
            end_time = self.epoch + 2* self.duration + nn * self.period

            # Select out-of-transit data for normalization
            out_of_transit_flux = self.lightcurve_collection.flux[np.where((self.lightcurve_collection.time.jd  > start_time) & (self.lightcurve_collection.time.jd < start_time + 0.05))]

            # Normalize the flux using median of out-of-transit data
            normalized_lc = self.lightcurve_collection.copy() # creating a copy of the lc to avoid modifying the original
            normalized_lc.flux = self.lightcurve_collection.flux / np.nanmedian(out_of_transit_flux)

            # Filter the data to focus on the transit
            transit_indices = ((normalized_lc.time.jd >= start_time) & (normalized_lc.time.jd <= end_time))

            if normalized_lc[transit_indices]["flux"].size != 0:
                self.transit_list.append(normalized_lc[transit_indices])

            nn += 1

        return self.transit_list

    def get_transits_csv(): 
        self.transit_list = []
        
        return

    def rd_data_csv(self, path):
        '''
        Função que le um csv e retorna no frame da lib KeplerLightCurve
        '''
        lc = pd.read_csv(path) 

        lc_lc = KeplerLightCurve(
            time=lc["time"].values,
            flux=lc["normal_flux"].values
        )

        self.transit_list=[lc_lc]

        return self.transit_list

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
        jj = np.argsort(phase)
        ff = phase[jj]

        smoothed_LC = scipy.ndimage.filters.uniform_filter(flux[jj], size = 100) # equivalente ao smooth do idl com edge_truncade
        smoothed_LC[0:200] = 1
        smoothed_LC[len(flux[jj])-200:len(flux[jj])] = 1
        x = phase[jj]
        y = 1 - smoothed_LC
        yh = 0.002

        kk = np.where(y >= yh)

        x1 = min(x[kk])
        x2 = max(x[kk])
        fa0 = (x1 + x2)/ 2 # valor central dos transitos em fase

        self.x0 = (np.fix(time[0] / porb) + fa0) * porb # tempo central do primeiro transito
        
        self.nt = np.fix(max(time - time[0])/ porb) + 1 # numero de transitos possiveis

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
        
        self.wl = np.zeros((self.n, self.n))
        x = np.arange(self.n/2-self.n, self.n - self.n/2, 1)
        
        for j in range (0, self.n):     
            z = np.sqrt(x**2 + (j - self.n/2.)**2)
            kk = np.where(z <= self.r)
            if kk[0].size > 0:
                    m = np.cos(np.arcsin(z[kk]/self.r))
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
        
        if self.manchas: 
            for mancha in self.manchas: 
                estrela_1.addMancha(mancha)
            estrela_1.criaEstrelaManchada()
            
        Nx1 = estrela_1.getNx() #coleta parametros da matriz estrela 
        Ny1 = estrela_1.getNy()
        raioEstrelaPixel1 = estrela_1.getRaioStar() #coleta raio da estrela em pixel 
        # elf, semiEixoUA, raioPlanJup, periodo, anguloInclinacao, ecc, anom, raioStar,mass
        planeta_1 = Planeta(self.AU, self.R_jup, self.porb, self.inc, self.ecc, self.anom, self.r_Sun, self.mass, self.planet_name)

        # self, Nx, Ny, raio_estrela_pixel, estrela_manchada: Estrela, planeta_: Planeta):
        eclipse1 = Eclipse(Nx1, Ny1, raioEstrelaPixel1, estrela_1, planeta_1)  #cria o objeto eclipse

        eclipse1.setTempoHoras(1.)

        eclipse1.criarEclipse(anim = False, plot = False)

        self.lc_model = np.array(eclipse1.getCurvaLuz())
        self.ts_model = np.array(eclipse1.getTempoHoras())
        
        return self.lc_model, self.ts_model

    #--------------------------------------------------#
    def retornaParametros(self):
        return self.u1,self.u2,self.porb,self.raioPlan,self.AU,self.inc,self.ts_model, self.mass
    
    def setTime(self,time):
        self.time = time

    def setFlux(self,flux):
        self.flux = flux
    
    def setFluxErr(self,flux_err):
        self.flux_err = flux_err
