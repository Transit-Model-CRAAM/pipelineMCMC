import numpy as np
import math


############ adição de luas ###########
class Moon:
    '''
    Classe Moon (lua), criada de acordo com a adição de planetas.
    '''
    pos = np.random.choice([-1, 1])

    def __init__(self, raioM, massM, raioEstrelaPixel, anguloInclinacao ,periodoM, raioPlanetaPixel,tempoHoras,distancia):

        '''
        :parâmetro raioM:: raio da lua em unidades de raio da Terra
        :parâmetro massM:: massa da Lua em unidades de massa da Terra
        :parâmetro anguloInclinacao:: angulo de inclinação do planeta em graus
        :parâmetro periodo:: periodo da órbita da lua em dias 
        :parâmetro raioPlanetaPixel:: raio do planeta em pixel
        :parâmetro tempoHoras:: tempo do trânsito do planeta em horas
        :parâmetro distancia:: distância lua-planeta em km
        '''
        
        tm0 = 0 # moon first transit time
        self.raioM = raioM
        self.massM = massM
        self.raioEstrelaPixel = raioEstrelaPixel
        self.anguloInclinacao = anguloInclinacao #em graus
        self.periodo = periodoM #em dias
        self.tm0 = tm0 #default
        self.raioPlanetaPixel = raioPlanetaPixel
        self.tempoHoras = tempoHoras
        self.distancia = distancia
        
        
    # moon orbit in equatorial plane of planet
    def moonOrbit(self, raioStar):
        '''
        funcao que calcula a orbita da lua, necessario apenas passar o raio da estrela como raioStar em km
        '''
        self.Rmoon = self.raioM / raioStar #raio da lua em relacao ao raio da estrela 
        self.RmoonPixel = self.Rmoon * self.raioEstrelaPixel #raio da lua calculado em pixel 
        
        self.dmoon = self.distancia * self.raioEstrelaPixel #calculo da distancia em pixel
        
        self.theta_m = 2*np.pi * self.tempoHoras / (self.periodo*24.) - self.tm0
        self.xm = self.dmoon * np.cos(self.theta_m)
        self.ym = self.dmoon * np.sin(self.theta_m) * np.cos(self.anguloInclinacao) 

    def getppMoon(self,tamanhoMatriz):
        #calculando a orbita projetada da lua
        dtor = np.pi/180.
        xlua = self.xm + tamanhoMatriz/2
        ylua = self.ym + tamanhoMatriz/2
        if(self.anguloInclinacao > 90.): 
            ylua = -self.dmoon*np.sin(self.theta_m)*np.cos(self.anguloInclinacao*dtor) + tamanhoMatriz/2

        #orbita projetada da Lua
        ppMoon, = np.where((xlua >= 0) & (xlua < tamanhoMatriz) & (ylua >= 0) & (ylua < tamanhoMatriz)) 
        self.xl = xlua[ppMoon]
        self.yl = ylua[ppMoon]
        return ppMoon   
    
    def getxl(self):
        return self.xl
    def getyl(self):
        return self.yl

    def getRmoon(self):
        return self.RmoonPixel

    def dMoon(self):
        return self.distancia * self.raioPlanetaPixel

    def getxm(self):
        return self.xm

    def getym(self):
        return self.ym