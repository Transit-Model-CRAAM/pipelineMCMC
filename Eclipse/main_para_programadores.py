import numpy as np
from matplotlib import pyplot
from Star.Estrela import Estrela
from Planet.Eclipse import Eclipse
from Misc.Verify  import Validar,calSemiEixo,calculaLat

'''
main programado para profissionais e estudantes familiarizados com a área 
--- estrela ---
parâmetro raio:: raio da estrela em pixel
parâmetro intensidadeMaxima:: intensidade da estrela que sera plotada 
parâmetro tamanhoMatriz:: tamanho em pixels da matriz estrela
parâmetro raioStar:: raio da estrela em relação ao raio do sol
parâmetro coeficienteHum:: coeficiente de escurecimento de limbo 1 (u1)p
parâmetro coeficienteDois:: coeficiente de escurecimento de limbo 2 (u2)
objeto estrela_ :: é o objeto estrela onde é guardada a matriz estrela de acordo com os parâmetros. Chamadas das funções da classe
estrela são feitas através dele.
parâmetro estrela :: variavel que recebe o objeto estrela 
--- planeta ---
parâmetro periodo:: periodo de órbita do planeta em dias 
parâmetro anguloInclinacao:: ângulo de inclinação do planeta em graus
parâmetro semieixoorbital:: semi-eixo orbital do planeta
parâmetro semiEixoRaioStar:: conversão do semi-eixo orbital em relação ao raio da estrela 
parâmetro raioPlanetaRstar:: conversão do raio do planeta em relação ao raio de Júpiter para em relação ao raio da estrela
---  mancha --- 
parâmetro latsugerida:: latitude sugerida para a mancha
parâmetro fa:: vetor com a área de cada mancha
parâmetro fi:: vetor com a intensidade de cada mancha
parâmetro li:: vetor com a longitude de cada mancha 
parâmetro quantidade:: variavel que armazena a quantidade de manchas
parâmetro r:: raio da mancha em relação ao raio da estrela
parâmetro intensidadeMancha:: intensidade da mancha em relação a intensidade da estrela
parâmetro lat:: latitude da mancha 
parâmetro longt:: longitude da mancha 
parâmetro raioMancha:: raio real da mancha
parâmetro area::  area da mancha 
--- eclipse ---
parâmetro eclipse:: variavel que guarda o objeto da classe eclipse que gera a curva de luz. Chamadas das funções da classe 
Eclipse () são feitas através dele. 
parâmetro tempoTransito:: tempo do transito do planeta 
parâmetro curvaLuz:: matriz curva de luz que sera plotada em forma de grafico 
parâmetro tempoHoras:: tempo do transito em horas
'''

raio= 373. #default (pixel)
intensidadeMaxima=240 #default
tamanhoMatriz = 856 #default
raioSun=0.805 #raio da estrela em relacao ao raio do sol
raioStar=raioSun*696340 #multiplicando pelo raio solar em Km 
coeficienteHum=0.377
coeficienteDois=0.024


#cria estrela
estrela_ = Estrela(raio,raioSun,intensidadeMaxima,coeficienteHum,coeficienteDois,tamanhoMatriz)

Nx= estrela_.getNx() #Nx e Ny necessarios para a plotagem do eclipse
Ny= estrela_.getNy()
dtor = np.pi/180.  

periodo = 2.219 # em dias
anguloInclinacao = 85.51  # em graus
ecc = 0
anom = 0 
raioPlanJup = 1.138 #em relação ao raio de jupiter
raioPlanetaRstar = (raioPlanJup*69911)/raioStar #multiplicando pelo raio de jupiter em km 
semiEixoUA = 0.031

dec=int(input("Deseja calular o semieixo Orbital do planeta através da 3a LEI DE KEPLER? 1. Sim 2.Não | "))
if dec==1:
    mass=0. #colocar massa da estrela em relação a massa do sol
    semieixoorbital = calSemiEixo(periodo,mass)
    semiEixoRaioStar = ((semieixoorbital/1000)/raioStar)
    #transforma em km para fazer em relação ao raio da estrela
else:
    #semiEixoUA = Validar('Semi eixo (em UA:)') #descomentar essa linha caso queira adicionar o valor em runtime
    semiEixoUA = 0.028 # Adicionar apenas valores maiores que 0 

    # em unidades de Rstar
    semiEixoRaioStar = ((1.469*(10**8))*semiEixoUA)/raioStar
    #multiplicando pelas UA (transformando em Km) e convertendo em relacao ao raio da estrela 

latsugerida = calculaLat(semiEixoRaioStar,anguloInclinacao)
print("A latitude sugerida para que a mancha influencie na curva de luz da estrela é: ", latsugerida)

#manchas
count = 0
quantidade = 0 #quantidade de manchas desejadas, se quiser acrescentar, mude essa variavel
#cria vetores do tamanho de quantidade para colocar os parametros das manchas
fa = [0.]*quantidade #vetor area manchas
fi = [0.]*quantidade #vetor intensidade manchas
li = [0.]*quantidade #vetor longitude manchas

while count!=quantidade: #o laço ira rodar a quantidade de manchas selecionada pelo usuario
    print('\033[1;35m\n\n══════════════════ Parâmetros da mancha ',count+1,'═══════════════════\n\n\033[m')
    r = Validar('Digite o raio da mancha em função do raio da estrela em pixels:')
                
    intensidadeMancha= float(input('Digite a intensidade da mancha em função da intensidade máxima da estrela:'))
    fi[count]=intensidadeMancha
    lat=float(input('Latitude da mancha:'))
    longt=float(input('Longitude da mancha:'))
    li[count]=longt

    raioMancha= r*raioStar
    area = np.pi *(raioMancha**2)
    fa[count]= area

    estrela=estrela_.manchas(r,intensidadeMancha,lat,longt) #recebe a escolha de se irá receber manchas ou não
    count+=1

#print vetor de intensidade, longitude e area da mancha para testes
print("Intensidades:",fi)
print("Longitudes:",li)
print("Areas:",fa)

estrela = estrela_.getEstrela()
#para plotar a estrela 
#caso nao queira plotar a estrela, comentar linhas abaixo
if (quantidade>0): #se manchas foram adicionadas. plotar
    estrela_.Plotar(tamanhoMatriz,estrela)

#criando lua 
lua = True #se nao quiser luas, mudar para False
eclipse= Eclipse(Nx,Ny,raio,estrela)
estrela_.Plotar(tamanhoMatriz,estrela)
eclipse.geraTempoHoras()
tempoHoras=eclipse.getTempoHoras()
#instanciando LUA
rmoon = 0.5 #em relacao ao raio da Terra
rmoon = rmoon *6371 #multiplicando pelo R da terra em Km
mass = 0.001 #em relacao a massa da Terra
mass = mass * (5.972*(10**24))
massPlaneta = 0.002 #em relacao ao R de jupiter
massPlaneta = massPlaneta * (1.898 *(10**27)) #passar para gramas por conta da constante G
G = (6.674184*(10**(-11)))
perLua = 0.1 #em dias 
distancia=((((perLua*24.*3600./2./np.pi)**2)*G*(massPlaneta+mass))**(1./3))/raioStar
distancia = distancia/100
moon = eclipse.criarLua(rmoon,mass,raio,raioStar,tempoHoras,anguloInclinacao,periodo,distancia)
estrela = estrela_.getEstrela()



#eclipse
eclipse.criarEclipse(semiEixoRaioStar, semiEixoUA, raioPlanetaRstar,raioPlanJup,periodo,anguloInclinacao,lua,ecc,anom)


print ("Tempo Total (Trânsito):",eclipse.getTempoTransito()) 
tempoTransito=eclipse.getTempoTransito()
curvaLuz=eclipse.getCurvaLuz()
tempoHoras=eclipse.getTempoHoras()

#Plotagem da curva de luz 
pyplot.plot(tempoHoras,curvaLuz)
pyplot.axis([-tempoTransito/2,tempoTransito/2,min(curvaLuz)-0.001,1.001])                       
pyplot.show()
