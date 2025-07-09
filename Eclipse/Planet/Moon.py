import string
import numpy as np
import math


############ adição de luas ###########
class Moon:
    '''
    Classe Moon (lua), criada de acordo com a adição de planetas.
    '''
    pos = np.random.choice([-1, 1])

    def __init__(self, raioM, massM ,periodoM,tempoHoras, planetaAnguloInclinacao, massPlaneta, raioPlanetaPixel, raioEstrelaPixel, raioStar, perPlan, dMoon = None):

        '''
        ; PURPOSE:
        ; Calculate the light curve of a planetary and its moon eclipse of its parent
        ; star.
        ;
        INPUT:
        ; ap = semi-major axis of planetary orbit (in Rstar)
        ; P  = Period of planetary orbit (in days), assumed circular orbit
        ; Rplan = planetary radius in Rstar (1 Rjup=6.9911e4 km=0.100447 Rsun)
        ; inc  = inclination angle of orbital plane (in deg)
        ; Pm = moon period
        ; Rmoon = moon radius (in Rstar)
        ; dmoon = moon-planet distance (in Rstar), if >0 moon transits first
        ; tetam0 = inicial position of moon in orbit (in deg)
        ; incm = inclinaiton angle of moon orbital plane (default = inc) in deg
        ; dt = time interval between light curve points (in minutes)
        ; wl = white light image of star (or Sun)
        ; 
        OUTPUT:
        ; ts  = time array (in hours w.r.t. transit central time)
        ; lc = light curve (dL/L)

        :parâmetro raioM:: raio da lua em unidades de raio da Terra
        :parâmetro massM:: massa da Lua em unidades de massa da Terra
        :parâmetro anguloInclinacao:: angulo de inclinação do planeta em graus
        :parâmetro periodo:: periodo da órbita da lua em dias 
        :parâmetro raioPlanetaPixel:: raio do planeta em pixel
        :parâmetro tempoHoras:: tempo do trânsito do planeta em horas
        :parâmetro distancia:: distância lua-planeta em km
        '''
        
        tm0 = 0 # moon first transit time
        self.raioM = raioM*6371 #r moon em relacao ao raio da Terra, #multiplicando pelo R da terra em Km
        self.massM = massM*(5.972*(10**24)) #em relacao a massa da Terra
        self.periodo = periodoM #em dias
        self.tm0 = tm0 #default
        self.tempoHoras = tempoHoras

        # planeta
        self.planetaAnguloInclinacao = planetaAnguloInclinacao
        self.massPlaneta = massPlaneta
        self.raioPlanetaPixel = raioPlanetaPixel
        self.perPlan = perPlan
        # estrela 
        self.raioEstrelaPixel = raioEstrelaPixel
        self.raioStar = raioStar # raio da estrela em relacao ao raio do sol

        if (dMoon is not None): 
            self.distancia = dMoon * 100 # em relacao ao planeta
        else: 
            self.distancia = (self.getDistancia())
        
    # moon orbit in equatorial plane of planet
    def moonOrbit(self):
        '''
        funcao que calcula a orbita da lua, necessario apenas passar o raio da estrela como raioStar em km
        '''
        raio_estrela_terra = self.raioStar * 109.076
        self.Rmoon = self.raioM / self.raioStar # raio da lua em relacao ao raio da estrela 
        self.RmoonPixel = self.Rmoon * self.raioEstrelaPixel #raio da lua calculado em pixel 
        
        self.dmoon = self.getDMoon() #calculo da distancia em pixel
        
        self.theta_m = 2*np.pi * self.tempoHoras / (self.perPlan*24.) - self.tm0
        self.xm = self.dmoon * np.cos(self.theta_m)
        self.ym = self.dmoon * np.sin(self.theta_m) * np.cos(self.planetaAnguloInclinacao) 

    def getppMoon(self,tamanhoMatriz):
        #calculando a orbita projetada da lua
        dtor = np.pi/180.
        xlua = self.xm + tamanhoMatriz/2
        ylua = self.ym + tamanhoMatriz/2
        if(self.planetaAnguloInclinacao > 90.): 
            ylua = -self.dmoon*np.sin(self.theta_m)*np.cos(self.planetaAnguloInclinacao*dtor) + tamanhoMatriz/2

        #orbita projetada da Lua
        ppMoon, = np.where((xlua >= 0) & (xlua < tamanhoMatriz) & (ylua >= 0) & (ylua < tamanhoMatriz)) 
        self.xl = xlua[ppMoon]
        self.yl = ylua[ppMoon]
        return ppMoon   

    def getDistancia(self): 
        G = (6.674184*(10**(-11)))
        return ((((self.periodo*24.*3600./2./np.pi)**2)*G*((self.massPlaneta*1e-3)+self.massM))**(1./3))/self.raioStar
    
    def getxl(self):
        return self.xl
    def getyl(self):
        return self.yl

    def getRmoon(self):
        return self.RmoonPixel

    def getDMoon(self):
        return self.distancia * self.raioPlanetaPixel

    def getxm(self):
        return self.xm

    def getym(self):
        return self.ym

    def setMoonName(self, name: string): 
        self.name = name