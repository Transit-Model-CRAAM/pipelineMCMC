# Snippet para acesso à API dos dados do exoplanet.eu

# https://exoplanet.eu/API/

# pip install pyvo

import pyvo

# Acesso à API
service = pyvo.dal.TAPService("http://voparis-tap-planeto.obspm.fr/tap")

# ======================================================================
# Criando a query para procurar por planetas de uma estrela
# específica ou por um planeta específico. As tabelas para serem
# acessadas são as do link abaixo. O exoplanet.eu usa a tabela
# exoplanet.epn_core
# http://voparis-tap-planeto.obspm.fr/__system__/dc_tables/list/form

# Se quiser procurar por uma estrela:
target_star = "Kepler-468"
query = f"SELECT * FROM exoplanet.epn_core WHERE star_name ILIKE '%{target_star}%'"

# Se quiser procurar por um planeta:
target_planet = "Kepler-468 b"
#query = f"SELECT * FROM exoplanet.epn_core WHERE target_name ILIKE '%{target_planet}%'"
# ======================================================================

# Aplicando a query.
# É importante notar que results e uma lista com todos os valores
# possíveis da pesquisa. Caso seja uma pesquisa específica por apenas
# 1 planeta, o resultado esperado virá em results[0].
results = service.search(query)

if len(results) == 0:
    print("Não existe nenhuma estrela/planeta com o nome requisitado")

# Pegando os parâmetros de acordo com o site
# As colunas para serem acessadas via get podem ser vistas no link abaixo
# http://voparis-tap-planeto.obspm.fr/__system__/dc_tables/show/tableinfo/exoplanet.epn_core

# Pegando todos os resultados:
for result in results:
    # Parâmetros da Estrela
    estrela = result.get('star_name')
    raioSun = result.get('star_radius')
    # Não possui u1 nem u2


    # Parâmetros do Planeta
    planeta = result.get('target_name')
    semiEixo = result.get("semi_major_axis")
    periodo = result.get("period")
    anguloInclinacao = result.get("inclination")
    raioPlanJup = result.get("radius")
    ecc = result.get("eccentricity")

    # Link para acesso da informação do planeta
    # Útil já que no final tem artigos pra poder procurar o coeficiente de limbo
    link = result.get('external_link')

    # Printando os dados
    print(f"Estrela: {estrela}")
    print(f"Raio da Estrela: {raioSun}")

    print("")

    print(f"Planeta: {planeta}")
    print(f"SemiEixo: {semiEixo}")
    print(f"Periodo: {periodo}")
    print(f"Ângulo de Inclinação: {anguloInclinacao}")
    print(f"Raio: {raioPlanJup}")
    print(f"Excentricidade: {ecc}")

    print("")

    print(f"Link: {link}")

    print("\n=====================================================\n")

print("\n=================== Result[0] =======================\n")

# Ou pegando só o primeiro mesmo:
# Parâmetros da Estrela
estrela = results[0].get('star_name')
raioSun = results[0].get('star_radius')
# Não possui u1 nem u2


# Parâmetros do Planeta
planeta = results[0].get('target_name')
semiEixo = results[0].get("semi_major_axis")
periodo = results[0].get("period")
anguloInclinacao = results[0].get("inclination")
raioPlanJup = results[0].get("radius")
ecc = results[0].get("eccentricity")

# Link para acesso da informação do planeta
# Útil já que no final tem artigos pra poder procurar o coeficiente de limbo
link = results[0].get('external_link')

# Printando os dados
print(f"Estrela: {estrela}")
print(f"Raio da Estrela: {raioSun}")

print("")

print(f"Planeta: {planeta}")
print(f"SemiEixo: {semiEixo}")
print(f"Periodo: {periodo}")
print(f"Ângulo de Inclinação: {anguloInclinacao}")
print(f"Raio: {raioPlanJup}")
print(f"Excentricidade: {ecc}")

print("")

print(f"Link: {link}")
