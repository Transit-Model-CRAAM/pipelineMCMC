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
from IPython.display import HTML
import cv2 as cv
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import pyplot
from Star.Estrela import Estrela
from Planet.Moon import Moon
from Planet.Planeta import Planeta
from Misc.Verify import calculaLat
# from keplerAux import keplerfunc  #biblioteca auxiliar caso a biblioteca kepler nao funcione
import matplotlib.animation as animation
import matplotlib.colors as mcolors
# para o calculo de orbitas excentricas (pip install kepler)
from kepler._core import solve
import os
from ctypes import *
from numpy.ctypeslib import ndpointer
import time
import platform
import gc


class Eclipse:

    def __init__(self, Nx, Ny, raio_estrela_pixel, estrela_manchada: Estrela, planeta_: Planeta):
        '''
        :parâmetro Nx e Ny: tamanho da matriz estrela 
        :parâmetro raioEstrelaPixel: raio da estrela em pixel 
        :parâmetro estrelaManchada: objeto ESTRELA passado como estrelaManchada apos a inserção de manchas
        '''
        self.Nx = Nx
        self.Ny = Ny
        self.intervaloTempo = 1
        self.tamanhoMatriz = self.Nx  # Nx ou Ny

        # Estrela
        self.raio_estrela_pixel = raio_estrela_pixel
        self.estrela_ = estrela_manchada
        self.estrela_matriz = estrela_manchada.getMatrizEstrela()

        # Planeta
        self.planeta_ = planeta_

        # Ruídos
        self.ejecaoMassa = self.getMatrizTransformada(self.estrela_matriz)

        # OUTPUT
        self.curvaLuz = [1.0 for i in range(self.Nx)]

    def geraTempoHoras(self, intervalo):
        '''
        Função chamada na Main para o cálculo do tempo de Trânsito em Horas
        '''
        #x=int(input("Intervalo de tempo=1. Deseja alterar? 1. SIM | 2. NÃO: "))

        self.intervaloTempo = intervalo

        tempoHoras = (np.arange(self.tamanhoMatriz) -
                      self.tamanhoMatriz/2)*self.intervaloTempo/60.   # em horas
        self.tempoHoras = tempoHoras

    def setTempoHoras(self, intervalo):
        self.intervaloTempo = intervalo  # em minutos
        self.tamanhoMatriz = self.Nx  # Nx ou Ny
        tempoHoras = (np.arange(self.tamanhoMatriz) -
                      self.tamanhoMatriz/2)*self.intervaloTempo/60.   # em horas
        self.tempoHoras = tempoHoras

    # a partir do momento em que a lua é instanciada na main, esses objetos se tornam objetos da classe com self.
    def criarLua(self, moon: Moon):
        moon.moonOrbit()
        self.planeta_.addLua(moon)
        
    def criarEclipse(self, anim=True, plot=True, collect_garbage=True):
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
        :parâmetro self.planeta_.ecc: excêntricidade da órbita do planeta
        :parâmetro self.planeta_.anom: self.planeta_.anomalia da órbita do planeta
        :parâmetro anim: verifica se a animação será mostrada para o usuário (True por default)
        :parâmetro plot: verifica se o gráfico da curva de luz será mostrado para o usuário (True por default)
        '''

        intervaloTempo = self.intervaloTempo
        tamanhoMatriz = self.tamanhoMatriz
        dtor = np.pi/180.
        semiEixoPixel = self.planeta_.semiEixoRaioStar * self.raio_estrela_pixel
        self.geraTempoHoras(1)

        '''Inicio do calculo do TEMPO TOTAL de trânsito através dos parâmetros passados ao planeta.'''

        # calculando obliquidade

        '''
        Parâmetros de órbita
        :parâmetro xplaneta: x na matriz que projetará o planeta
        :parâmetro yplaneta: y na matriz que projetará o planeta
        '''

        nk = 2*np.pi/(self.planeta_.periodo * 24)    # em horas^(-1)
        Tp = self.planeta_.periodo*self.planeta_.anom / \
            360. * 24.  # tempo do pericentro (em horas)
        m = nk*(self.tempoHoras-Tp)     # em radianos

        # calculando a self.planeta_.anomalia excentrica em radianos

        eccanom = solve(m, self.planeta_.ecc)  # subrotina em anexo
        xs = semiEixoPixel*(np.cos(eccanom)-self.planeta_.ecc)
        ys = semiEixoPixel * \
            (math.sqrt(1-(self.planeta_.ecc**2))*np.sin(eccanom))

        ang = self.planeta_.anom*dtor-(np.pi/2)
        xp = xs*np.cos(ang)-ys*np.sin(ang)
        yp = xs*np.sin(ang)+ys*np.cos(ang)

        ie, = np.where(self.tempoHoras == min(abs(self.tempoHoras)))

        xplaneta = xp-xp[ie[0]]
        yplaneta = yp*np.cos(self.planeta_.anguloInclinacao*dtor)

        # Intervalo para calculo do transito
        # rearranja o vetor apenas com os pontos necessários para a análise da curva de luz
        pp, = np.where((abs(xplaneta) < 1.2 * tamanhoMatriz/2)
                       & (abs(yplaneta) < tamanhoMatriz/2))
        xplan = xplaneta[pp] + tamanhoMatriz/2
        yplan = yplaneta[pp] + tamanhoMatriz/2

        raioPlanetaPixel = self.planeta_.raioPlanetaRstar * self.raio_estrela_pixel

        '''
        Inicio do calculo do tempo em Horas e da curva de Luz na matriz
        :parâmetro nn: calculo do numero de pontos na curva de luz
        :parâmetro tamanhoMatriz: recebe a estrela manchada para depois plotar o planeta
        :parâmetro tempoHoras: calcula o tempo do transito em horas, transformando-o em objeto da classe Eclipse
        :parâmetro curvaLuz: calcula a curva de luz do transito do planeta ao eclipsar a estrela, também se torna 
        objeto de Eclipse       
        '''
        latitudeTransito = -np.arcsin(self.planeta_.semiEixoRaioStar*np.cos(
            self.planeta_.anguloInclinacao*dtor))/dtor  # latitude Sul (arbitraria)

        # duracao do transito em horas
        duracaoTransito = 2 * (90.-np.arccos((np.cos(latitudeTransito*dtor)) /
                               self.planeta_.semiEixoRaioStar)/dtor)*self.planeta_.periodo/360*24.
        tempoTotal = 3*duracaoTransito
        self.tempoTotal = tempoTotal

        # calculo do numero de pontos na curva de luz
        nn = np.fix(tempoTotal*60./intervaloTempo)

        # seleciona a maior orbita para que a curva de luz seja plotada de maneira correta (observando ela inteira)
        rangeloop = pp

        ''''
        Curva de Luz e normalização da intensidade
        '''

        # definição de variaveis para utilizacao da função de calculo da curva de luz em C
        tamanho = self.tamanhoMatriz*self.tamanhoMatriz

        # Matriz plan auxiliar para ser passada como parametro para o script em C
        # Create the numpy array filled with 1.0
        plan_np = np.ones(tamanho, dtype=np.float64)

        # Get the ctypes pointer
        plan = plan_np.ctypes.data_as(POINTER(c_double))
        kk = np.arange(tamanhoMatriz*tamanhoMatriz)

        # Matriz kk auxiliar para ser passada como parametro para o script em C
        kk2 = (c_double * len(kk))(*kk)

        # Obter o caminho absoluto do diretório atual
        dir_atual = os.path.dirname(os.path.abspath(__file__))

        # Voltar um diretório para chegar ao diretório pai
        dir_pai = os.path.dirname(dir_atual)

        # Verifica o SO e se o Python é 32 ou 64 bit
        if (platform.system() == "Windows"):
            if (platform.architecture()[0] == "32bit"):
                script_path = os.path.join(dir_pai, 'scripts', 'func32.dll')
                my_func = WinDLL("a", winmode=0x8)
            elif (platform.architecture()[0] == "64bit"):
                script_path = os.path.join(dir_pai, 'scripts', 'func64.dll')
                my_func = WinDLL(script_path, winmode=0x8)
        elif (platform.system() == "Darwin"):
            script_path = os.path.join(dir_pai, 'scripts', 'func64.dylib')
            my_func = cdll.LoadLibrary(script_path)
        else:
            script_path = os.path.join(dir_pai, 'scripts', 'func64.so')
            my_func = CDLL(script_path)

        # Prepara os tipos de cada variável dos argumentos e do retorno da função do calculo da curva de luz
        my_func.curvaLuz.restype = c_double
        my_func.curvaLuz.argtypes = c_double, c_double, c_int, c_double, POINTER(
            c_double), POINTER(c_double), c_double
        my_func.curvaLuzLua.restype = c_double
        my_func.curvaLuzLua.argtypes = c_double, c_double, c_double, c_double, c_double, c_int, c_double, POINTER(
            c_double), POINTER(c_double), c_double

        raioPlanetaPixel = float(raioPlanetaPixel)

        '''
        Criação da matriz para plotagem:
        '''
        # Inicio dos loops para a plotagem e calculo do trânsito
        # maximo da curva de luz, usado na normalizacao da intensidade
        maxCurvaLuz = np.sum(self.estrela_matriz)

        em = self.getMatrizTransformada(self.estrela_matriz)

        # Animacao de trânsito
        if (anim):
            # criacao de variaveis para plotagem da animacao
            fig, (ax1, ax2) = plt.subplots(2, 1)
            ims = []
            plota = True  # variavel FLAG que indica quando armazenar a imagem do PLOT
            numAux = 0  # variavel FLAG que indica quantidade de imagens no vetor de PLOT

            print("\nAguarde um momento, a animacao do trânsito está sendo gerada...\n")

            if (self.planeta_.hasMoons()):
                tamPp = len(pp)
                tamLoop = len(rangeloop)

                for lua in self.planeta_.luas:
                    ppMoon = lua.getppMoon(self.tamanhoMatriz)
                    tamMoon = len(ppMoon)
                    if (tamPp > tamMoon):
                        rangeloop = pp
                    else:
                        if (tamLoop < tamMoon):
                            rangeloop = ppMoon

                        # x plan e y plan se alteram caso haja o acrescimo de luas
                        xplan = xplaneta[ppMoon] + self.tamanhoMatriz/2
                        yplan = yplaneta[ppMoon] + self.tamanhoMatriz/2

                    # Adiciona lua ao eclipse
                    self.addLua(rangeloop, xplan, yplan, raioPlanetaPixel, kk2,
                                maxCurvaLuz, numAux, ims, ax1, ax2, fig, kk, my_func, plota, lua)

            else:
                for i in range(0, len(rangeloop)):

                    x0 = xplan[i]
                    y0 = yplan[i]

                    self.curvaLuz[rangeloop[i]] = my_func.curvaLuz(
                        x0, y0, self.tamanhoMatriz, raioPlanetaPixel, em, kk2, maxCurvaLuz)

                    if (plota and self.curvaLuz[rangeloop[i]] != 1 and numAux < 200):
                        plan = np.zeros(tamanhoMatriz*tamanhoMatriz)+1.
                        ii = np.where(((kk/tamanhoMatriz-y0)**2+(kk-tamanhoMatriz *
                                      np.fix(kk/tamanhoMatriz)-x0)**2 <= raioPlanetaPixel**2))
                        plan[ii] = 0.
                        # posicao adicionada na matriz
                        plan = plan.reshape(
                            self.tamanhoMatriz, self.tamanhoMatriz)

                        plt.axis([0, self.Nx, 0, self.Ny])

                        cmap = plt.cm.get_cmap("hot").copy()
                        cmap.set_over('white')  # valores acima de vmax ficam brancos

                        # Normalização com vmax = intensidade máxima
                        norm = mcolors.Normalize(vmin=0, vmax=self.estrela_.intensidadeMaxima * 1.2, clip=False)

                        # Multiplica a máscara e aplica normalização
                        image_to_plot = self.estrela_matriz * plan
                        image_to_plot = np.where(image_to_plot > self.estrela_.intensidadeMaxima * 1.2, self.estrela_.intensidadeMaxima * 1.2+1, image_to_plot)

                        im = ax1.imshow(image_to_plot, cmap=cmap, norm=norm, animated=True, aspect='equal')

                        # armazena na animação os pontos do grafico (em imagem)
                        ims.append([im])
                        numAux += 1
                    # variavel auxiliar que seleciona o intervalo correto para plotagem
                    plota = not (plota)

                self.plot_anim(ax1, ax2, plan, fig, ims)
            
        # Sem animação de trânsito
        else:
            # Inicio dos loops para a plotagem e calculo do trânsito
            start = time.time()

            if (self.planeta_.hasMoons() == False):
                for i in range(0, len(rangeloop)):

                    x0 = xplan[i]
                    y0 = yplan[i]

                    self.curvaLuz[rangeloop[i]] = my_func.curvaLuz(
                        x0, y0, self.tamanhoMatriz, raioPlanetaPixel, em, kk2, maxCurvaLuz)
            else:
                for i in range(0, len(rangeloop)):
                    ### adicionando luas ###
                    for lua in self.planeta_.luas:
                        x0 = xplan[i]
                        y0 = yplan[i]

                        xxm = lua.getxm()
                        yym = lua.getxym()

                        xm = x0-xxm[i]
                        ym = y0-yym[i]

                        self.curvaLuz[rangeloop[i]] = my_func.curvaLuzLua(
                            x0, y0, xm, ym, lua.Rmoon, self.tamanhoMatriz, raioPlanetaPixel, em, kk2, maxCurvaLuz)
            if (plot):
                end = time.time()
                print(end-start)
                # Altere os valores conforme necessário
                plt.figure(figsize=(10, 5))
                plt.plot(self.tempoHoras, self.curvaLuz)
                plt.axis([-self.tempoTotal/2, self.tempoTotal /
                         2, min(self.curvaLuz)-0.001, 1.001])

                # Adicionando título e rótulos dos eixos
                plt.title('Modelo de Eclipse na estrela hospedeira')
                plt.xlabel('Tempo (Horas)')
                plt.ylabel('Brilho da estrela (fluxo normalizado)')

                plt.show()
                
        # Libera memória de lixo que ficou alocado
        if collect_garbage:
            gc.collect()

    def addLua(self, rangeloop, xplan, yplan, raioPlanetaPixel, kk2, maxCurvaLuz, numAux, ims, ax1, ax2, fig, kk, my_func, plota, lua):
        em = self.getMatrizTransformada(self.estrela_matriz)

        for i in range(0, len(rangeloop)):
            x0 = xplan[i]
            y0 = yplan[i]
            ### adicionando luas ###
            xxm = lua.getxm()
            yym = lua.getym()

            xm = x0-xxm[i]
            ym = y0-yym[i]

            self.curvaLuz[rangeloop[i]] = my_func.curvaLuzLua(x0, y0, xm, ym, lua.getRmoon(
            ), self.tamanhoMatriz, raioPlanetaPixel, em, kk2, maxCurvaLuz)

            if (plota and self.curvaLuz[rangeloop[i]] != 1 and numAux < 200):
                plan = np.zeros(self.tamanhoMatriz*self.tamanhoMatriz)+1.
                ii = np.where(((kk/self.tamanhoMatriz-y0)**2+(kk-self.tamanhoMatriz *
                              np.fix(kk/self.tamanhoMatriz)-x0)**2 <= raioPlanetaPixel**2))
                ll = np.where((kk/self.tamanhoMatriz-ym)**2+(kk-self.tamanhoMatriz *
                              np.fix(kk/self.tamanhoMatriz)-xm)**2 <= lua.getRmoon()**2)
                plan[ii] = 0.
                plan[ll] = 0.
                # posicao adicionada na matriz
                plan = plan.reshape(self.tamanhoMatriz, self.tamanhoMatriz)
                plt.axis([0, self.Nx, 0, self.Ny])

                cmap = plt.cm.get_cmap("hot").copy()
                cmap.set_over('white')  # valores acima de vmax ficam brancos

                # Normalização com vmax = intensidade máxima
                norm = mcolors.Normalize(vmin=0, vmax=self.estrela_.intensidadeMaxima * 1.2, clip=False)

                # Multiplica a máscara e aplica normalização
                image_to_plot = self.estrela_matriz * plan
                image_to_plot = np.where(image_to_plot > self.estrela_.intensidadeMaxima * 1.2, self.estrela_.intensidadeMaxima * 1.2+1, image_to_plot)

                im = ax1.imshow(image_to_plot, cmap=cmap, norm=norm, animated=True, aspect='equal')

                # armazena na animação os pontos do grafico (em imagem)
                ims.append([im])
                numAux += 1
            # variavel auxiliar que seleciona o intervalo correto para plotagem
            plota = not (plota)

        self.plot_anim(ax1, ax2, plan, fig, ims)


    def calculaLatMancha(self):
        latsugerida = calculaLat(
            self.planeta_.semiEixoRaioStar, self.planeta_.anguloInclinacao)
        print("A latitude sugerida para que a mancha influencie na curva de luz da estrela é:", latsugerida)
        return latsugerida

    def calculaLongMancha(self, a, time, lat):
        # Converte latitude para radianos
        latitude_rad = math.radians(lat)
        # Calcula o ângulo interno
        angle = math.radians(90) - (math.radians(360) *
                                    time) / (24 * self.planeta_.periodo)
        # Calcula a longitude da mancha
        longitude = math.degrees(
            math.asin((a * math.cos(angle)) / math.cos(abs(latitude_rad))))

        print("A longitude sugerida para que a mancha influencie na curva de luz da estrela é:", longitude)
        return longitude

    def getMatrizTransformada(self, estrela):
        matriz_aux = np.array(estrela, dtype=np.float64)
        matriz_estrela_manchada = matriz_aux.ctypes.data_as(POINTER(c_double))
        return matriz_estrela_manchada

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

    def getLuas(self):
        return self.planeta_.luas

    def setEstrela(self, estrela):
        '''
        com essa funcao, é possivel passar a estrela atualizada para o eclipse que esta se formando, caso sejam adicionadas mais manchas.
        '''
        self.estrela_matriz = estrela

    def getError(self):
        '''
        Retorna o valor de erro, ocorrendo ou não algum. Se não houver erro, recebe 0. Se houver, a variável terá
        seu valor de inicio (que é -1)
        '''
        return self.error

    def plot_anim(self, ax1, ax2, plan, fig, ims):
        cmap = plt.cm.get_cmap("hot").copy()
        cmap.set_over('white')  # values above vmax will be white

        norm = mcolors.Normalize(vmin=0, vmax=self.estrela_.intensidadeMaxima * 1.2, clip=False)

        image_to_plot = self.estrela_matriz * plan
        image_to_plot = np.where(
            image_to_plot > self.estrela_.intensidadeMaxima * 1.2,
            self.estrela_.intensidadeMaxima * 1.2 + 1,
            image_to_plot
        )

        im = ax1.imshow(
            image_to_plot,
            cmap=cmap,
            norm=norm,
            animated=True,
            origin='lower',
            aspect='equal'
        )

        cbar = plt.colorbar(im, ax=ax1, extend='max')
        cbar.set_label('Intensidade')

        pos = ax1.get_position()

        new_width = 0.27                # Fração da largura da figura (Ajustar caso preciso)
        new_left = 0.5 - new_width / 2  # Centralizando horizontalmente

        # Aplicando a nova posição
        ax1.set_position([new_left, pos.y0, new_width, pos.height])

        box = cbar.ax.get_position()
        cbar.ax.set_position([
            box.x0 - 0.1,   # Aumentar = mover pra esquerda / Diminuir = mover pra direita
            box.y0,         # Manter posição vertical
            box.width,      # Manter largura
            box.height      # Manter altura
        ])

        ax2.plot(self.tempoHoras, self.curvaLuz, label='Curva de Luz')

        # Definindo os limites dos eixos
        ax2.axis([-self.tempoTotal/2, self.tempoTotal /
                  2, min(self.curvaLuz)-0.001, 1.001])

        # Adicionando título e legendas dos eixos
        ax1.set_title('Modelo de Eclipse na estrela hospedeira')
        ax2.set_xlabel('Tempo (Horas)')
        ax2.set_ylabel('Brilho da estrela (fluxo normalizado)')

        # Exibindo a legenda
        # Ajuste a posição conforme necessário
        ax2.legend(loc='lower right')

        ani = animation.ArtistAnimation(fig, ims, interval=50, blit=True, repeat_delay=0.1)

        plt.show(block=True)
