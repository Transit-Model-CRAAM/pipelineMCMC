__author__ = "Adriana Valio, Beatriz Duque, Felipe Pereira Pinho"
__copyright__ = "..."
__credits__ = ["Universidade Presbiteriana Mackenzie, CRAAM"]
__license__ = ""
__version__ = ""
__maintainer__ = ""
__email__ = "biaduque7@hotmail.com"
__status__ = "Production"

'''
Programa que simula o eclipse e a curva de luz de um planeta ao transitar 
sua host star.
Nesse programa é calculada a curva de luz da estrela em relação aos parâmetros do planeta adicionados
pelo usuário.
***Bibliotecas importadas***
numpy:
matplotlib:
estrela: arquivo de programa onde são calculados os parâmetros da estrela, dado os inputs do usuário (raio, intensidade,etc)
verify:função criada para validar entradas, por exemplo numeros nao float/int ou negativos
'''

import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import pyplot
from Star.Estrela import Estrela
from Planet.Moon import Moon
from Misc.Verify import Validar
#from keplerAux import keplerfunc  #biblioteca auxiliar caso a biblioteca kepler nao funcione
import matplotlib.animation as animation
from kepler._core import solve #para o calculo de orbitas excentricas (pip install kepler)
import os
from ctypes import *
from numpy.ctypeslib import ndpointer
import time
import gc
import sys
import platform

class Eclipse:

   
    def __init__(self,Nx,Ny,raioEstrelaPixel,estrelaManchada):
        

        '''
        :parâmetro Nx e Ny: tamanho da matriz estrela 
        :parâmetro raioEstrelaPixel: raio da estrela em pixel 
        :parâmetro estrelaManchada: objeto ESTRELA passado como estrelaManchada apos a inserção de manchas
        '''
        self.Nx = Nx
        self.Ny = Ny
        self.raioEstrelaPixel = raioEstrelaPixel
        self.estrelaManchada = estrelaManchada
        
        # OUTPUT
        curvaLuz =[ 1.0 for i in range(self.Nx)]
        self.curvaLuz = curvaLuz

    def geraTempoHoras(self):
        '''
        Função chamada na Main para o cálculo do tempo de Trânsito em Horas
        '''
        x=int(input("Intervalo de tempo=1. Deseja alterar? 1. SIM | 2. NÃO: "))
        if x ==1:
            self.intervaloTempo=float(input('Digite o intervalo de tempo em minutos: '))
        elif x==2:
            self.intervaloTempo = 1.   # em minutos


        self.tamanhoMatriz= self.Nx #Nx ou Ny
        tempoHoras = (np.arange(self.tamanhoMatriz)-self.tamanhoMatriz/2)*self.intervaloTempo/60.   # em horas
        self.tempoHoras= tempoHoras

    def setTempoHoras(self,intervalo):
        self.intervaloTempo = intervalo   # em minutos
        self.tamanhoMatriz= self.Nx #Nx ou Ny
        tempoHoras = (np.arange(self.tamanhoMatriz)-self.tamanhoMatriz/2)*self.intervaloTempo/60.   # em horas
        self.tempoHoras= tempoHoras

    #a partir do momento em que a lua é instanciada na main, esses objetos se tornam objetos da classe com self.
    def criarLua(self, raioM, massM, raioPlanetaPixel, raioStar,tempoHoras,anguloInclinacao,periodo,distancia):
        moon = Moon(raioM, massM, self.raioEstrelaPixel,anguloInclinacao ,periodo, raioPlanetaPixel, self.tempoHoras,distancia)
        moon.moonOrbit(raioStar)
        Rmoon = moon.getRmoon()

        #coleta de dados necessarias para a plotagem do eclipse
        self.xxm = moon.getxm()
        self.yym = moon.getym()
        self.Rmoon = Rmoon #em pixel 
        self.massM = massM
        self.tamanhoMatriz= self.Nx
        #coletando dados da lua
        self.ppMoon = moon.getppMoon(self.tamanhoMatriz)
        self.xl = moon.getxl()
        self.yl = moon.getyl()
        return moon
        


    def criarEclipse(self,semiEixoRaioStar,semiEixoUA, raioPlanetaRstar, raioPlanJup ,periodo,anguloInclinacao,lua,ecc,anom,anim=True,plot=True):

        '''
        Criação da classe eclipse, que retornará a curva de luz do trânsito do planeta ao redor da estrela


        ****parâmetros atribuidos ao planeta****
        :parâmetro periodo: período de rotação do planeta
        :parâmetro semiEixoRaioStar: semi eixo do planeta em relação ao raio da estrela
        :parâmetro semiEixoUA: semi eixo do planeta em UA
        :parâmetro anguloInclinacao: angulo de inclinação do planeta
        :parâmetro raioPlanetaRstar: raio do planeta em relacao ao raio da estrela 
        :parâmetro raioPlanJup: raio do planeta em relacao ao raio de Jupiter
        :parâmetro lua: lua que orbita o planeta (entra como True or False)
        :parâmetro ecc: excêntricidade da órbita do planeta
        :parâmetro anom: anomalia da órbita do planeta
        :parâmetro anim: verifica se a animação será mostrada para o usuário (True por default)
        :parâmetro plot: verifica se o gráfico da curva de luz será mostrado para o usuário (True por default)
        '''

        intervaloTempo = self.intervaloTempo
        tamanhoMatriz = self.tamanhoMatriz
        self.semiEixoRaioStar = semiEixoRaioStar
        self.semiEixoUA = semiEixoUA
        self.raioPlanetaRstar = raioPlanetaRstar
        self.raioPlanJup = raioPlanJup
        self.periodo = periodo
        self.anguloInclinacao = anguloInclinacao
        self.lua = lua
        self.ecc = ecc
        self.anom = anom 

        dtor = np.pi/180.
        semiEixoPixel = self.semiEixoRaioStar * self.raioEstrelaPixel

        '''Inicio do calculo do TEMPO TOTAL de trânsito através dos parâmetros passados ao planeta.'''

        #ecc = 0. #default
        #anom = 0.  #default

        #calculando obliquidade

        '''
        Parâmetros de órbita
        :parâmetro xplaneta: x na matriz que projetará o planeta
        :parâmetro yplaneta: y na matriz que projetará o planeta
        '''

        nk=2*np.pi/(self.periodo*24)    # em horas^(-1)
        Tp=self.periodo*anom/360.*24. # tempo do pericentro (em horas)
        m = nk*(self.tempoHoras-Tp)     # em radianos

        # calculando a anomalia excentrica em radianos
        
        eccanom = solve(m,ecc)  # subrotina em anexo
        xs=semiEixoPixel*(np.cos(eccanom)-ecc)
        ys=semiEixoPixel*(math.sqrt(1-(ecc**2))*np.sin(eccanom))

        ang=anom*dtor-(np.pi/2)
        xp=xs*np.cos(ang)-ys*np.sin(ang)
        yp=xs*np.sin(ang)+ys*np.cos(ang)

        ie, = np.where(self.tempoHoras == min(abs(self.tempoHoras)))

        xplaneta=xp-xp[ie[0]]
        yplaneta=yp*np.cos(self.anguloInclinacao*dtor)

        pp, = np.where((abs(xplaneta) < 1.2 * tamanhoMatriz/2) & (abs(yplaneta) < tamanhoMatriz/2)) #rearranja o vetor apenas com os pontos necessários para a análise da curva de luz
        xplan = xplaneta[pp] + tamanhoMatriz/2
        yplan = yplaneta[pp] + tamanhoMatriz/2

        raioPlanetaPixel = self.raioPlanetaRstar * self.raioEstrelaPixel

       

        '''
        Inicio do calculo do tempo em Horas e da curva de Luz na matriz
        :parâmetro nn: calculo do numero de pontos na curva de luz
        :parâmetro tamanhoMatriz: recebe a estrela manchada para depois plotar o planeta
        :parâmetro tempoHoras: calcula o tempo do transito em horas, transformando-o em objeto da classe Eclipse
        :parâmetro curvaLuz: calcula a curva de luz do transito do planeta ao eclipsar a estrela, também se torna 
        objeto de Eclipse       
        '''
        latitudeTransito = -np.arcsin(self.semiEixoRaioStar*np.cos(self.anguloInclinacao*dtor))/dtor # latitude Sul (arbitraria)
        # duracao do transito em horas
        duracaoTransito=2 * (90.-np.arccos((np.cos(latitudeTransito*dtor))/self.semiEixoRaioStar)/dtor)*self.periodo/360*24. 
        tempoTotal = 3 * duracaoTransito
        self.tempoTotal= tempoTotal

        
        # calculo do numero de pontos na curva de luz
        nn=np.fix(tempoTotal*60./intervaloTempo)

        #seleciona a maior orbita para que a curva de luz seja plotada de maneira correta (observando ela inteira)
        if(lua == True):
            if (len(pp)>len(self.ppMoon)):
                rangeloop = pp
            else: 
                rangeloop = self.ppMoon
                xplan = xplaneta[self.ppMoon] + tamanhoMatriz/2 #x plan e y plan se alteram caso haja o acrescimo de luas 
                yplan = yplaneta[self.ppMoon] + tamanhoMatriz/2
        else:
            rangeloop = pp


        ''''
        Curva de Luz e normalização da intensidade
        '''
        # maximo da curva de luz, usado na normalizacao da intensidade
        maxCurvaLuz = np.sum(self.estrelaManchada)

        # definição de variaveis para utilizacao da função de calculo da curva de luz em C
        tamanho = self.tamanhoMatriz*self.tamanhoMatriz

        # Matriz em auxiliar para ser passada como parametro para o script em C
        em = (c_double*tamanho)()
        for j in range (self.tamanhoMatriz):
            for i in range(self.tamanhoMatriz):
                index = i*self.tamanhoMatriz + j
                num = self.estrelaManchada[i][j]
                num = (c_double)(num)
                em[index] = num

        # Matriz plan auxiliar para ser passada como parametro para o script em C
        plan = (c_double*tamanho)()
        for i in range(tamanho):
            num = 1.
            num = (c_double)(num)
            plan[i] = num

        kk=np.arange(tamanhoMatriz*tamanhoMatriz)

        # Matriz kk auxiliar para ser passada como parametro para o script em C
        kk2 = (c_double * len(kk))(*kk)

        # Obter o caminho absoluto do diretório atual
        dir_atual = os.path.dirname(os.path.abspath(__file__))

        # Voltar um diretório para chegar ao diretório pai
        dir_pai = os.path.dirname(dir_atual)
        
        # Verifica o SO e se o Python é 32 ou 64 bit
        if(platform.system() == "Windows"):
            if(platform.architecture()[0] == "32bit"):
                script_path = os.path.join(dir_pai, 'scripts', 'func32.dll')
                my_func = WinDLL("a", winmode = 0x8)
            elif(platform.architecture()[0] == "64bit"):
                script_path = os.path.join(dir_pai, 'scripts', 'func64.dll')
                my_func = WinDLL(script_path, winmode = 0x8)
        elif(platform.system() == "Darwin"):
            script_path = os.path.join(dir_pai, 'scripts', 'func64.dylib')
            my_func = cdll.LoadLibrary(script_path)
        else:
            script_path = os.path.join(dir_pai, 'scripts', 'func64.so')
            my_func = CDLL(script_path)

        # Prepara os tipos de cada variável dos argumentos e do retorno da função do calculo da curva de luz
        my_func.curvaLuz.restype = c_double
        my_func.curvaLuz.argtypes = c_double,c_double,c_int,c_double,POINTER(c_double),POINTER(c_double),c_double
        my_func.curvaLuzLua.restype = c_double
        my_func.curvaLuzLua.argtypes = c_double,c_double,c_double,c_double,c_double,c_int,c_double,POINTER(c_double),POINTER(c_double),c_double

        raioPlanetaPixel = float(raioPlanetaPixel)

        '''
        Criação da matriz para plotagem:
        '''
        if(anim):
            #criacao de variaveis para plotagem da animacao 
            fig, (ax1, ax2) = plt.subplots(2,1)
            ims = []
            plota = True #variavel FLAG que indica quando armazenar a imagem do PLOT 
            numAux = 0 #variavel FLAG que indica quantidade de imagens no vetor de PLOT

            print("\nAguarde um momento, a animacao do trânsito está sendo gerada.\n")
            #Inicio dos loops para a plotagem e calculo do trânsito
            #start = time.time()
            intervalo = math.ceil(len(rangeloop)/400)
            if (lua == False):
                for i in range(0,len(rangeloop)):

                                x0 = xplan[i]
                                y0 = yplan[i]

                                self.curvaLuz[rangeloop[i]]=my_func.curvaLuz(x0,y0,self.tamanhoMatriz,raioPlanetaPixel,em,kk2,maxCurvaLuz)

                                if(plota and self.curvaLuz[rangeloop[i]] != 1 and numAux<200):
                                    plan = np.zeros(tamanhoMatriz*tamanhoMatriz)+1.
                                    ii = np.where(((kk/tamanhoMatriz-y0)**2+(kk-tamanhoMatriz*np.fix(kk/tamanhoMatriz)-x0)**2 <= raioPlanetaPixel**2))
                                    plan[ii]=0.
                                    plan = plan.reshape(self.tamanhoMatriz, self.tamanhoMatriz) #posicao adicionada na matriz
                                    plt.axis([0,self.Nx,0,self.Ny])
                                    im = ax1.imshow(self.estrelaManchada*plan,cmap="hot", animated = True)
                                    ims.append([im]) #armazena na animação os pontos do grafico (em imagem)
                                    numAux+=1
                                plota = not(plota) #variavel auxiliar que seleciona o intervalo correto para plotagem
            else:
                for i in range(0,len(rangeloop)):

                                x0 = xplan[i] 
                                y0 = yplan[i]

                                ### adicionando luas ###
                                xm = x0-self.xxm[i]         
                                ym = y0-self.yym[i]  
                                
                                self.curvaLuz[rangeloop[i]]=my_func.curvaLuzLua(x0,y0,xm,ym,self.Rmoon,self.tamanhoMatriz,raioPlanetaPixel,em,kk2,maxCurvaLuz)

                                if(plota and self.curvaLuz[rangeloop[i]] != 1 and numAux<200):
                                    plan = np.zeros(tamanhoMatriz*tamanhoMatriz)+1.
                                    ii = np.where(((kk/tamanhoMatriz-y0)**2+(kk-tamanhoMatriz*np.fix(kk/tamanhoMatriz)-x0)**2 <= raioPlanetaPixel**2))
                                    ll = np.where((kk/tamanhoMatriz-ym)**2+(kk-tamanhoMatriz*np.fix(kk/tamanhoMatriz)-xm)**2 <= self.Rmoon**2)
                                    plan[ii]=0.
                                    plan[ll]=0.
                                    plan = plan.reshape(self.tamanhoMatriz, self.tamanhoMatriz) #posicao adicionada na matriz
                                    plt.axis([0,self.Nx,0,self.Ny])
                                    im = ax1.imshow(self.estrelaManchada*plan,cmap="hot", animated = True)
                                    ims.append([im]) #armazena na animação os pontos do grafico (em imagem)
                                    numAux+=1
                                plota = not(plota) #variavel auxiliar que seleciona o intervalo correto para plotagem

            #end = time.time()
            #print(end-start)
            ax2.plot(self.tempoHoras,self.curvaLuz)
            ax2.axis([-self.tempoTotal/2,self.tempoTotal/2,min(self.curvaLuz)-0.001,1.001])
            ani =animation.ArtistAnimation(fig, ims, interval=50, blit=True,repeat_delay=0.1)
            
            plt.show()
            #ani.save('animacao_transito.gif',writer="PillowWriter") #salva o gif gerado na raiz do arquivo, para utilizacao do usuario
        else:
            #Inicio dos loops para a plotagem e calculo do trânsito
            start = time.time()
            if (lua == False):
                for i in range(0,len(rangeloop)):

                                x0 = xplan[i]
                                y0 = yplan[i]

                                self.curvaLuz[rangeloop[i]]=my_func.curvaLuz(x0,y0,self.tamanhoMatriz,raioPlanetaPixel,em,kk2,maxCurvaLuz)
            else:
                for i in range(0,len(rangeloop)):

                                x0 = xplan[i] 
                                y0 = yplan[i]

                                ### adicionando luas ###
                                xm = x0-self.xxm[i]         
                                ym = y0-self.yym[i]  
                                
                                self.curvaLuz[rangeloop[i]]=my_func.curvaLuzLua(x0,y0,xm,ym,self.Rmoon,self.tamanhoMatriz,raioPlanetaPixel,em,kk2,maxCurvaLuz)
            if(plot):
                end = time.time()
                print(end-start)
                plt.plot(self.tempoHoras,self.curvaLuz)
                plt.axis([-self.tempoTotal/2,self.tempoTotal/2,min(self.curvaLuz)-0.001,1.001])
                plt.show()

        locals().clear # Limpa qualquer possível sujeira de memória
        del my_func

        error=0
        self.error=error

    '''Chamada dos objetos atribuídos à classe Eclipse.'''
    def getTempoTransito(self):
        '''Retorna o parâmetro tempoTotal, representando o tempo de trânsito do planeta em sua host star.'''
        return self.tempoTotal
    def getTempoHoras(self):
        '''Retorna o parâmetro tempoHoras, representando o tempo de trânsito do planeta em sua host star em Horas.'''
        return self.tempoHoras
    def getCurvaLuz(self):
        '''Retorna o parâmetro curvaLuz, representando a curva de luz da estrela que possui um planeta a orbitar nela.'''
        return self.curvaLuz

    def getRaioPlan(self):
        return self.raioPlanetaRstar

    def getRplanJup(self):
        return self.raioPlanJup

    def getSemiEixo(self):
        return self.semiEixoUA
    
    def getsemiEixoRaioStar(self):
        return self.semiEixoRaioStar

    def getPeriodo(self):
        return self.periodo

    def getInc(self):
        return self.anguloInclinacao

    def getEccAnom(self):
        return self.ecc, self.anom 

    def getLua(self):
        return self.lua

    def getError(self):
        '''
        Retorna o valor de erro, ocorrendo ou não algum. Se não houver erro, recebe 0. Se houver, a variável terá
        seu valor de inicio (que é -1)
        '''
        return self.error
    def setEstrela(self,estrela):
        '''
        com essa funcao, é possivel passar a estrela atualizada para o eclipse que esta se formando, caso sejam adicionadas mais manchas.
        '''
        self.estrelaManchada = estrela