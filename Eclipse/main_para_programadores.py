# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2020-2026 The ECLIPSE authors (CRAAM, Universidade Presbiteriana Mackenzie)
'''
main programado para profissionais e estudantes familiarizados com a área 
'''

import numpy as np
from matplotlib import pyplot
from Star.Estrela import Estrela
from Planet.Eclipse import Eclipse
from Misc.Verify  import Validar,calSemiEixo,calculaLat
from Planet.Planeta import Planeta
from Planet.Moon import Moon

''''
--- ⭐️ Estrela ---
parâmetro raio:: raio da estrela em pixel
parâmetro intensidadeMaxima:: intensidade da estrela que sera plotada 
parâmetro tamanhoMatriz:: tamanho em pixels da matriz estrela
parâmetro raioStar:: raio da estrela em relação ao raio do sol
parâmetro coeficienteHum:: coeficiente de escurecimento de limbo 1 (u1)p
parâmetro coeficienteDois:: coeficiente de escurecimento de limbo 2 (u2)
objeto estrela_ :: é o objeto estrela 
estrela são feitas através dele.
parâmetro estrela_matriz :: variavel que recebe a matriz da estrela
'''

raio_estrela_pixel = 373. #default (pixel)
intensidadeMaxima = 240 #default
tamanhoMatriz = 856 #default
raioSun = 0.805 #raio da estrela em relacao ao raio do sol
raioStar = raioSun*696340 #multiplicando pelo raio solar em Km 
coeficienteHum = 0.377
coeficienteDois = 0.024

# cria estrela
estrela_ = Estrela(raio_estrela_pixel,raioSun,intensidadeMaxima,coeficienteHum,coeficienteDois,tamanhoMatriz)

Nx = estrela_.getNx() #Nx e Ny necessarios para a plotagem do eclipse
Ny = estrela_.getNy()
dtor = np.pi/180.  

'''
--- 🪐 Planeta ---
parâmetro periodo:: periodo de órbita do planeta em dias 
parâmetro anguloInclinacao:: ângulo de inclinação do planeta em graus
parâmetro semieixoorbital:: semi-eixo orbital do planeta
parâmetro semiEixoRaioStar:: conversão do semi-eixo orbital em relação ao raio da estrela 
parâmetro raioPlanetaRstar:: conversão do raio do planeta em relação ao raio de Júpiter para em relação ao raio da estrela
'''
periodo = 2.219 # em dias
anguloInclinacao = 85.51  # em graus
ecc = 0
anom = 0 
raioPlanJup = 1.138 #em relação ao raio de jupiter
semiEixoUA = 0.031
massPlaneta = 0.002 #em relacao ao R de jupiter

# cria planeta
planeta_ = Planeta(semiEixoUA, raioPlanJup, periodo, anguloInclinacao, ecc, anom, estrela_.getRaioSun(), massPlaneta, "modelo_planeta")

dec = int(input("Deseja calular o semieixo Orbital do planeta através da 3a LEI DE KEPLER? 1. Sim 2.Não | "))
if dec==1:
    mass = 0. #colocar massa da estrela em relação a massa do sol
    semieixoorbital = calSemiEixo(periodo,mass)
    semiEixoRaioStar = ((semieixoorbital/1000)/raioStar)
    #transforma em km para fazer em relação ao raio da estrela
else:
    #semiEixoUA = Validar('Semi eixo (em UA:)') #descomentar essa linha caso queira adicionar o valor em runtime
    semiEixoUA = 0.028 # Adicionar apenas valores maiores que 0 

    # em unidades de Rstar
    semiEixoRaioStar = ((1.469*(10**8))*semiEixoUA)/raioStar
    #multiplicando pelas UA (transformando em Km) e convertendo em relacao ao raio da estrela 

'''
--- 🌖 Eclipse ---
parâmetro eclipse:: variavel que guarda o objeto da classe eclipse que gera a curva de luz. Chamadas das funções da classe 
Eclipse () são feitas através dele. 
parâmetro tempoTransito:: tempo do transito do planeta 
parâmetro curvaLuz:: matriz curva de luz que sera plotada em forma de grafico 
parâmetro tempoHoras:: tempo do transito em horas
'''
estrela_matriz = estrela_.getMatrizEstrela()
eclipse_ = Eclipse(Nx,Ny,raio_estrela_pixel, estrela_, planeta_)

tempoHoras = 1
eclipse_.geraTempoHoras(tempoHoras)
eclipse_.criarEclipse(False, anim=True)

print ("Tempo Total (Trânsito):", eclipse_.getTempoTransito()) 
tempoTransito = eclipse_.getTempoTransito()
curvaLuz = eclipse_.getCurvaLuz()
tempoHoras = eclipse_.getTempoHoras()

#Plotagem da curva de luz 
pyplot.plot(tempoHoras, curvaLuz)
pyplot.axis([-tempoTransito/2, tempoTransito/2, min(curvaLuz)-0.001, 1.001])                       
pyplot.show()

#### Adicionando interferências na curva de luz #### 
latsugerida = eclipse_.calculaLatMancha()
print("A latitude sugerida para que a mancha influencie na curva de luz da estrela é: ", latsugerida)


'''
---  🎯 Mancha --- 
parâmetro latsugerida:: latitude sugerida para a mancha
parâmetro quantidade:: variavel que armazena a quantidade de manchas
parâmetro raio_mancha:: raio da mancha em relação ao raio da estrela
parâmetro intensidade:: intensidade da mancha em relação a intensidade da estrela
parâmetro latitude:: latitude da mancha 
parâmetro longitude:: longitude da mancha 
parâmetro raioMancha:: raio real da mancha
parâmetro area::  area da mancha 
'''
count = 0
quantidade = 1 #quantidade de manchas desejadas, se quiser acrescentar, mude essa variavel

while count!=quantidade: #o laço ira rodar a quantidade de manchas selecionada pelo usuario
    print('\033[1;35m\n\n══════════════════ Parâmetros da mancha ',count+1,'═══════════════════\n\n\033[m')
    raio_mancha = Validar('Digite o raio da mancha em função do raio da estrela em pixels:')   
    intensidade = float(input('Digite a intensidade da mancha em função da intensidade máxima da estrela:'))
    latitude = float(input('Latitude da mancha:'))
    longitude = float(input('Longitude da mancha:'))

    mancha = Estrela.Mancha(intensidade, raio_mancha, latitude, longitude)

    raioMancha = raio_mancha*raioStar
    area = np.pi *(raioMancha**2)
    mancha.area = area

    estrela_.addMancha(mancha)
    estrela_.criaEstrelaManchada()
    count+=1

estrela_matriz = estrela_.getMatrizEstrela()

# para plotar a estrela manchada
plota_manchada = True # caso nao queira plotar a estrela manchada, mudar para False
if (quantidade>0 and plota_manchada): #se manchas foram adicionadas. plotar
    eclipse_.setEstrela(estrela_matriz)

    estrela_.Plotar(tamanhoMatriz, estrela_matriz)
    eclipse_.criarEclipse(False, anim=True)

    tempoTransito = eclipse_.getTempoTransito()
    curvaLuz = eclipse_.getCurvaLuz()
    tempoHoras = eclipse_.getTempoHoras()

    #Plotagem da curva de luz 
    pyplot.plot(tempoHoras,curvaLuz)
    pyplot.axis([-tempoTransito/2,tempoTransito/2,min(curvaLuz)-0.001,1.001])                       
    pyplot.show()

'''
---  🌙 Lua --- 
parâmetro rmoon:: raio da lua em relacao ao raio da Terra
parâmetro mass:: em relacao a massa da Terra
parâmetro perLua:: periodo da lua em dias 
'''
addMoon = False #mudar para True se quiser adicionar 
if addMoon: 
    rmoon = 50 
    mass = 0.5 
    perLua = 0.03 

    moon = Moon(rmoon, mass, perLua, tempoHoras, planeta_.anguloInclinacao, planeta_.mass, estrela_.raio, planeta_.getRaioPlanPixel(estrela_.raio, estrela_.raioSun), estrela_.getRaioSun(), planeta_.periodo)
    eclipse_.criarLua(moon) #adiciona lua no planeta que esta no eclipse

    moon.setMoonName("Moon name")
    estrela_.Plotar(tamanhoMatriz, estrela_matriz)
    eclipse_.criarEclipse(False, anim=True)

    # Criando planeta com lua 
    estrela_matriz = estrela_.getMatrizEstrela()

    # Passa para o eclipse a estrela atualizada
    eclipse_.setEstrela(estrela_matriz)

    # Análise final 
    print ("Tempo Total (Trânsito):",eclipse_.getTempoTransito()) 
    tempoTransito = eclipse_.getTempoTransito()
    curvaLuz = eclipse_.getCurvaLuz()
    tempoHoras = eclipse_.getTempoHoras()

    #Plotagem da curva de luz 
    pyplot.plot(tempoHoras,curvaLuz)
    pyplot.axis([-tempoTransito/2,tempoTransito/2,min(curvaLuz)-0.001,1.001])                       
    pyplot.show()