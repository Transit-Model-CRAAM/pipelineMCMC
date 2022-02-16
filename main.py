__author__ = "Adriana Valio, Beatriz Duque"
__copyright__ = "..."
__credits__ = ["Universidade Presbiteriana Mackenzie, CRAAM"]
__license__ = ""
__version__ = ""
__maintainer__ = ""
__email__ = "biaduque7@hotmail.com"
__status__ = "Production"



import numpy as np
from matplotlib import pyplot
from estrela_nv1 import estrela
from eclipse_nv1 import Eclipse
from verify import Validar,calSemiEixo,calculaLat


############CRIAÇÃO DA ESTRELA#############


'''
No main, as funções estrela e eclipse sao chamadas. 
Primeiramente a estrela é criada, enquanto o usuário adiciona os parâmetros que desejar
como por exemplo, manchas, flares, fáculas, etc.
Após a conclusão dessa etapa, o usuário deve adicionar parâmetros do planeta
em relação a estrela para que assim seja plotada a curva de luz característica deste.

importações 
verify:função criada para validar entradas, por exemplo numeros nao float/int ou negativos
estrela: classe Estrela que possui como objetos manchas, faculas e flares
eclipse: classe eclipse que recebe estrela e calcula sua curva de luz 
'''

print('''\033[1;36m
╔═════════════»»» ESTRELA EXEMPLO «««════════════╗
║-RAIO EM PIXEL: 373                             ║      
║-INTENSIDADE DO CENTRO DA ESTRELA: 240          ║
║ COEFICIENTE DE ESCURECIMENTO DE LIMBO          ║
║-u1: 0.5                                        ║
║-u2: 0.3                                        ║
╚════════════════════════════════════════════════╝\033[m\n\n''')

#default
#raio da estrela em pixels
raio = 373.
#intensidade do centro da estrela (ou intensidade máxima)
intensidadeMaxima=240.
#iniciando a coleta de dados dos parametros da estrela
error=-1
while error==-1:
    try:
        x=int(input('1. Plotar a estrela exemplo|2. Adicionar parâmetros a estrela:'))
        if x==1:
            raioStar=1.
            raioSun=1
        #parametros dos coeficientes de escurecimento de limbo
            coeficienteHum = 0.5
            coeficienteDois = 0.3
        elif x==2:
            raioSun= Validar('Raio da estrela (em relação ao raio do sol):')
            raioStar=raioSun*696340 #multiplicando pelo raio solar em Km 
            coeficienteHum = float(input('Coeficiente de escurecimento de limbo 1 (u1):'))
            coeficienteDois = float(input('Coeficiente de escurecimento de limbo 2 (u2):'))

        #fazer um if de config recomendada
        x=int(input('Tamanho da matriz definido como 856, deseja alterar? 1.Sim 2.Não :'))
        if x== 1:
            tamanhoMatriz=int(input('Digite o tamanho da Matriz:'))
        elif x==2:
            tamanhoMatriz = 856
        # Construcao da estrela com escurecimento de limbo 
        #Criação do objeto estrela dentro de estrela.
        error=0 #se passar de todas as etapas, erro=0 e sai do laço
    except Exception as erro:
        print(f'\033[0;31mO valor digitado é inválido. Por favor, digite novamente. O tipo de problema encontrado foi{erro.__class__}\n\n\033[m') #retorna o tipo de erro


error=-1 #trata o erro
#ESTRELA
while error==-1: #enquanto o erro for -1, o while irá rodar, ou seja, o usuario terá que refazer o processo se estiver apresentando algum erro
    #se passar por todos os inputs sem erro, a função retornará 0, e assim, sairá do laço e passará para a próxima etapa
    try:
        estrela_ = estrela(raio,raioSun,intensidadeMaxima,coeficienteHum,coeficienteDois,tamanhoMatriz)
        estrela = estrela_.getEstrela()
        estrela_.Plotar(tamanhoMatriz,estrela)
        error=estrela_.getError() #se não houver nenhum erro, ele receberá 0, o que o fará sair do laço.
    except Exception as erro:
        print(f'\033[0;31mO valor digitado é inválido. Por favor, digite novamente. O tipo de problema encontrado foi{erro.__class__}\n\n\033[m') #retorna o tipo de erro


#CHAMADA DE VARIÁVEIS NECESSÁRIAS PARA O CÁLCULO DA CURVA DE LUZ 
Nx= estrela_.getNx() #Nx e Ny necessarios para a plotagem do eclipse
Ny= estrela_.getNy()
raioEstrelaPixel = estrela_.getRaioStar()
estrelaManchada = estrela_.getEstrela() #retorna na verdade a estrela sem manchas para plotar no final
#atribuição de variáveis dadas pelos parâmetros para que sejam plotadas.


#instanciando o ECLIPSE 
eclipse= Eclipse(Nx,Ny,raioEstrelaPixel,estrelaManchada) #instancia o Eclipse 

##################### Main que proporciona a criação de mais de um planeta, cada um com sua lua #####################
def criandoPlanetas(estrela_,planetas,qtd):

    '''
    Função utilizada para a adição de mais de um planeta que orbita a estrela. Essa função
    roda dentro de um while de acordo com a quantidade de planetas desejada pelo usuario
    parametro_ :: objeto estrela
    parametro planetas :: quantidade de planetas a ser adicionado
    parametro qtd :: variavel inicial = 0 que é incrementada a cada laço (utilizada como auxiliar para verificacoes 
                     dentro da funcao)
    '''

    ##################### INICIANDO O ECLIPSE ##########################

    #planeta,transito,curvaLuz 
    print('''\033[1;33m
            ╔══════════════════════════»»» ECLIPSE EXEMPLO «««══════════════════════════════╗
            ║- PERÍODO: 10 (em dias)                                                        ║           
            ║- SEMI EIXO (em relação ao raio da estrela): 15                               	║
            ║- ÂNGULO DE INCLINAÇÃO (em graus): 88                                       	║
            ║- RAIO DO PLANETA (em relação ao raio da estrela): 0.1                        	║
            ╚═══════════════════════════════════════════════════════════════════════════════╝\033[m'''
            )
    dtor = np.pi/180.  
    aux= True
    while aux == True: 
        try:
            x=int(input('1.Eclipse Exemplo. 2.Alterar:')) #escolha de como passar os parametros do eclipse
            aux= False
        except Exception as erro:
            print(f'\033[0;31mO valor digitado é inválido. Por favor, digite novamente. O tipo de problema encontrado foi{erro.__class__}\n\n\033[m')
                    
    if x==1:
    #entradas default
        periodo = 10.  # em dias
        semiEixoRaioStar = 15   # em unidades de Rstar
        anguloInclinacao = 88.  # em graus
        raioPlanetaRstar = 0.1   # em unidades de Rstar      
                
    elif x==2:
        aux= True
        while aux == True:
            try:
                periodo = Validar("Período:")
                anguloInclinacao = float(input('Angulo de inclinação:')) 
                raioPlanJup = Validar('Raio do planeta (em relação ao raio de Júpiter):')
                raioPlanetaRstar = (raioPlanJup*69911)/raioStar #multiplicando pelo raio de jupiter em km 
                anom = float(input("Digite e Anomalia: (Default = 0)"))
                ecc = float(input("Digite a Excentricidade: (Default = 0)"))

                dec=int(input("Deseja calular o semieixo Orbital do planeta através da 3a LEI DE KEPLER? 1. Sim 2.Não |"))
                if dec==1:
                    mass= Validar("Digite a massa da estrela em unidades de MassSun:")
                    semieixoorbital = calSemiEixo(periodo,mass)
                    print("resultado = ", semieixoorbital)
                    semiEixoRaioStar = ((semieixoorbital/1000)/raioStar)
                    #transforma em km para fazer em relação ao raio da estrela
                else:
                    semiEixoUA = Validar('Semi eixo (em UA:)')
                    # em unidades de Rstar
                    semiEixoRaioStar = ((1.469*(10**8))*semiEixoUA)/raioStar
                    #multiplicando pelas UA (transformando em Km) e convertendo em relacao ao raio da estrela
                
                while semiEixoRaioStar*np.cos(anguloInclinacao*dtor) >= 1: 
                    print('O Planeta não está eclipsando a estrela. Por favor, mude o angulo de inclinação.')
                    anguloInclinacao = float(input('Angulo de inclinação:')) 
                
                aux= False
            
            except Exception as erro:
                print(f'\033[0;31mO valor digitado é inválido. Por favor, digite novamente. O tipo de problema encontrado foi{erro.__class__}\n\n\033[m')
                    

    #MANCHA DA ESTRELA
    error=-1
    while error==-1:
        try:
            escolha= Validar('\033[1;35mDeseja adicionar manchas em sua estrela? 1. Sim 2. Não |\033[m')  #x define a escolha, se haverá mancha ou não.
            if escolha==1:
                latsugerida = calculaLat(semiEixoRaioStar,anguloInclinacao)
                print("A latitude sugerida para que a mancha influencie na curva de luz da estrela é:", latsugerida)
                quantidade= Validar('\033[1;35mDigite a quantidade de manchas a serem adicionadas:\033[m')
                quantidade=int(quantidade)
                count=0
                #cria vetores do tamanho de quantidade para colocar os parametros das manchas
                fa = [0.]*quantidade #vetor area manchas
                fi = [0.]*quantidade #vetor intensidade manchas
                li = [0.]*quantidade #vetor longitude manchas
                while count!=quantidade: #o laço ira rodar a quantidade de manchas selecionada pelo usuario
                    print('\033[1;35m\n\n══════════════════ Parâmetros da mancha ',count+1,'═══════════════════\n\n\033[m')
                    r = Validar('Digite o raio da mancha em função do raio da estrela: ')
                    
                    intensidadeMancha= Validar('Digite a intensidade da mancha em função da intensidade máxima da estrela:')
                    fi[count]=intensidadeMancha
                    lat=float(input('Latitude da mancha:'))

                    longt=float(input('Longitude da mancha:'))
                    li[count]=longt

                    raioMancha= r*raioStar
                    area = np.pi *(raioMancha**2)
                    fa[count] = area

                    estrela=estrela_.manchas(r,intensidadeMancha,lat,longt) #recebe a escolha de se irá receber manchas ou não
                    error=estrela_.getError()
                    count+=1
                print("Intensidades:",fi)
                print("Longitudes:",li)
                print("Areas:",fa)
                estrela = estrela_.getEstrela()
                eclipse.setEstrela(estrela) #passando a estrela para o eclipse
            else:
                estrela = estrela_.getEstrela()  # armazena o que esta guardado no objeto estrela para sobrescrever
                error=0    
        except Exception as erro:
            print(f'\033[0;31mO valor digitado é inválido. Por favor, digite novamente. O tipo de problema encontrado foi{erro.__class__}\n\n\033[m')


    estrela_.Plotar(tamanhoMatriz,estrela)
    eclipse.geraTempoHoras()
    tempoHoras=eclipse.getTempoHoras() #gerar o tempo do Eclipse em Horas

    ##################### INICIANDO A LUA ##########################

    error=-1 #tratando os erros em Lua
    lua = False
    while error==-1:
        try:
            escolha= Validar('\033[1;35mDeseja adicionar luas neste Planeta? 1. Sim 2. Não |\033[m')  #x define a escolha, se haverá luas ou não.
            if escolha==1:
                lua = True
                quantidade= Validar('\033[1;35mDigite a quantidade de LUAS a serem adicionadas:\033[m')
                quantidade=int(quantidade)
                count=0
                while count!=quantidade: #o laço ira rodar a quantidade de luas selecionada pelo usuario
                    print('\033[1;35m\n\n══════════════════ Parâmetros da LUA ',count+1,'═══════════════════\n\n\033[m')
                    # radius, mass,distance, raioPlanetaPixel, raioStar,tempoHoras
                    rmoon = Validar('Digite o raio da Lua em função do raio da Terra: ')
                    rmoon = rmoon *6371 #multiplicando pelo R da terra em Km
                    mass = Validar('Digite a massa da Lua (em unidades de massa Terra): ')  
                    mass = mass * (5.972*(10**24))
                    massPlaneta = Validar('Digite a massa do Planeta (em unidades de massa de Júpiter ): ') #deve ser convertido para Kg
                    massPlaneta = massPlaneta * (1.898 *(10**27)) #passar para gramas por conta da constante G
                    G = (6.674184*(10**(-11)))
                    perLua = Validar('Digite o período da órbita da Lua:') #em dias 
                    distancia=((((perLua*24.*3600./2./np.pi)**2)*G*(massPlaneta+mass))**(1./3))/raioStar
                    distancia = distancia/100

                    #x1000**(1/3)/10ˆ5 para a conversao por conta do parametro G
                    print("distancia = ",distancia)
                    print("rmoon = ",rmoon/raioStar)
                    raioPlanetaPixel = (raioPlanetaRstar)*(raioEstrelaPixel)
                    #instanciando LUA
                    moon = eclipse.criarLua(rmoon,mass,raioPlanetaPixel,raioStar,tempoHoras,anguloInclinacao,periodo,distancia)
                    estrela = estrela_.getEstrela()
                    count+=1
                break      
            else:
                estrela = estrela_.getEstrela()
                error=0    
        except Exception as erro:
            print(f'\033[0;31mO valor digitado é inválido. Por favor, digite novamente. O tipo de problema encontrado foi{erro.__class__}\n\n\033[m')

    error = -1
    #inicia o calculo do eclipse com planetas, luas, manchas e o que o usuario desejou adicionar 
    while error==-1:
        try:
            eclipse.criarEclipse(semiEixoRaioStar, semiEixoUA, raioPlanetaRstar,raioPlanJup,periodo,anguloInclinacao,lua,ecc,anom)
            error=eclipse.getError()
        except Exception as erro:
            print(f'\033[0;31mO valor digitado é inválido. Por favor, digite novamente. O tipo de problema encontrado foi{erro.__class__}\n\n\033[m')

    #deseja printar a curva de luz agora? adicionar opcao
    decisao = Validar("Deseja printar a curva de luz agora? 1.SIM | 2.NÃO: ")
    if (decisao==2 and qtd!=planetas-1):
        print("\033[1;33mO proximo planeta pode ser adicionado :) \033[m")
    else:
        print ("Tempo Total (Trânsito):",eclipse.getTempoTransito()) 
        tempoTransito=eclipse.getTempoTransito()
        #print("Curva Luz:",eclipse.getCurvaLuz())
        curvaLuz=eclipse.getCurvaLuz()
        #Plotagem da curva de luz sozinha
        pyplot.plot(tempoHoras,curvaLuz)
        pyplot.axis([-tempoTransito/2,tempoTransito/2,min(curvaLuz)-0.001,1.001])                       
        pyplot.show()

#############################################################################################################
#cria objeto do grafico
#passa pra main esse objeto 
planetas = Validar("Digite a quantidade de planeta(s) desejado(s): ")
qtd = 0
while (qtd!= planetas):
    criandoPlanetas(estrela_,planetas,qtd)
    qtd+=1 
#plota a orbita dos planetas adicionados 
