# Base branch: ```main```
>[!IMPORTANT]
>```main``` represents more stable version of the code
> 
> For a better development experience, we recommend using [VSCode](https://code.visualstudio.com/download), with the extension for [jupyter notebook](https://marketplace.visualstudio.com/items?itemName=ms-toolsai.jupyter), [c and gcc](https://marketplace.visualstudio.com/items?itemName=ms-vscode.cpptools)

# Pipeline MCMC

[ 🇧🇷 ] Repositório com o projeto de programação para calcular a curva de luz de uma estrela com possíveis exoplanetas em sua órbita através de programação Python híbrida para programação C usando Pipeline em tempo de execução.

[ 🇺🇸 ] Repository with the programming project in order to calculate the light curve of a Star with possible exoplanets in its orbit through hybrid Python programming to C programming using Pipeline at runtime. 

>[!TIP]
> 🔭 [ 🇧🇷 ] Conheça também o nosso site de divulgação científica, **[Stars and Exoplanets](https://www.starsandexoplanets.org/)**, onde você pode simular curvas de luz de estrelas com planetas em órbita!
>
> 🔭 [ 🇺🇸 ] Also check out our science outreach website, **[Stars and Exoplanets](https://www.starsandexoplanets.org/)**, where you can simulate light curves of stars with orbiting planets!

### Origin of the project

[ 🇧🇷 ] Este repositório é uma continuação do trabalho iniciado pela Profa. Dra. Adriana Valio, que desenvolveu originalmente em IDL o modelo de trânsito planetário em frente a estrelas com manchas, apresentado no artigo:

[ 🇺🇸 ] This repository is a continuation of the work started by Prof. Dr. Adriana Valio, who originally developed in IDL the model of planetary transits across spotted stars, presented in the paper:

- Silva, A. V. R. (2003). **Method for Spot Detection on Solar-like Stars**. *The Astrophysical Journal*, 585(2), L147–L150. [doi:10.1086/374324](https://doi.org/10.1086/374324)

*** 
## Preparing the enviroment 

[ 🇧🇷 ] Bibliotecas necessárias para rodar o programa e como instalá-las:

Com o Python instalado na sua máquina, abra o prompt de comando (cmd) e insira os comandos a seguir para instalar cada biblioteca. Pode ser necessário adicionar "python" ou "py3" antes de "pip install" (Ex.: py3 pip install numpy) dependendo de como o Python estiver instalado no seu computador.

[ 🇺🇸 ] Required libraries to run the software and how to install them:

With Python installed in your machine, open the command prompt (cmd) and type the following commands to install each library. It may be required to add "python" or "py3" before the "pip install" (Ex.: py3 pip install numpy) depending on how Python is installed in your computer.

## Preparing multitrheading environment 

[ 🇧🇷 ] Para habilitar o paralelismo de código (multithreading) no ECLIPSE, é necessário seguir o tutorial disponível em: 

https://github.com/Transit-Model-CRAAM/pipelineMCMC/tree/feat/new-transit-method/Eclipse/scripts#readme

Lembre-se de possuir instalado o ambiente para programação em C (C e gcc)

[ 🇺🇸 ] To enable code parallelism (multithreading) in ECLIPSE, you must follow the tutorial available at:

https://github.com/Transit-Model-CRAAM/pipelineMCMC/blob/feat/new-transit-method/Eclipse/scripts/README-EN.md

Remember to have the C programming environment installed (C and gcc)

*** 

## How to install 

[ 🇧🇷 ] Para instalar todas as bibliotecas de uma vez, abra o cmd dentro da pasta onde os arquivos foram extraídos e digite o comando ``` pip install -r requirements.txt --user ```.

Caso a instalação dê algum problema, também é possível instalar individualmente cada biblioteca através dos comandos abaixo.

[ 🇺🇸 ] To install all libraries at once, open cmd inside the folder where files were extracted and type in the command ``` pip install -r requirements.txt --user ```.

In case the installation goes wrong, it's also possible to install each library individually using the commands below.


MatPlotLib:  ``` pip install matplotlib ```

NumPy: ``` pip install numpy ```

SciPy: ``` pip install scipy ```

PyAstronomy: ``` pip install pyastronomy ```

Pandas: ``` pip install pandas ```

LmFit: ``` pip install lmfit ```

Emcee: ``` pip install emcee ```

Corner: ``` pip install corner ```

Requests: ``` pip install requests ```

Pytz: ``` pip install pytz ```

IPython: ``` pip install ipython ```

LightKurve: ``` pip install lightkurve ```

Kepler.py: ``` pip install kepler.py ```

*** 
## How to run

[ 🇧🇷 ] Para o uso da versão main.ipynb, recomenda-se o uso do ```Jupyter Notebook```. Para instalá-lo, basta digitar no prompt de comando o seguinte comando:

``` pip install notebook ```

Para utilizá-lo, basta ir até a pasta onde está instalado o programa e digitar no clipboard "jupyter notebook".

[ 🇺🇸 ] To use the main.ipynb version, it is recommended using the ```Jupyter Notebook```. To install it, just type in the command prompt the following command:

``` pip install notebook ```

To use it, go to the folder where the software is installed and type in the clipboard "jupyter notebook".

***
## Publications

[ 🇧🇷 ] Artigos produzidos com o uso deste código:

[ 🇺🇸 ] Papers produced using this code:

- Duque Estrada, B., Pinho, F., & Valio, A. (2026). **Modelagem em Python do trânsito de um planeta em frente a uma estrela com manchas: um estudo aplicado a KEPLER-423**. *Revista Brasileira de Ensino de Física*, 48. [doi:10.1590/1806-9126-RBEF-2025-0333](https://doi.org/10.1590/1806-9126-RBEF-2025-0333)

- Duque Estrada, B., Pinho, F., & Valio, A. (2025). **Simulations of stellar activity in Python and its impacts on the star's light curve**. *Proceedings of the International Astronomical Union*, 21(S400), 188–191. [doi:10.1017/S1743921326107170](https://doi.org/10.1017/S1743921326107170)

### Undergraduate research (Iniciação Científica)

[ 🇧🇷 ] Projetos de iniciação científica que deram forma a este repositório:

[ 🇺🇸 ] Undergraduate research projects that shaped this repository:

- **Beatriz Duque** (orientação/advisor: Adriana Valio)
  - [ 🇧🇷 ] Conversão do modelo original em IDL para Python, adição de documentação e reestruturação do código com orientação a objetos.
  - [ 🇺🇸 ] Conversion of the original IDL model to Python, addition of documentation and restructuring of the code with object-oriented programming.
  - Silva, B. D. E. T., & Valio, A. B. M. (2020). **Atividade estelar e o raio de exoplanetas**. *XVI Jornada de Iniciação Científica e X Mostra de Iniciação Tecnológica*, Universidade Presbiteriana Mackenzie. [Link](https://eventoscopq.mackenzie.br/jornada/pt_BR/article/view/2986)

- **Felipe Pinho** (orientação/advisor: Adriana Valio)
  - [ 🇧🇷 ] Paralelismo do ajuste MCMC, com o suporte de Beatriz Duque, que iniciou essa implementação.
  - [ 🇺🇸 ] Parallelization of the MCMC fitting, with the support of Beatriz Duque, who started this implementation.
  - Pinho, F. P., & Valio, A. B. M. (2022). **Determinação de parâmetros de exoplanetas por ajuste do trânsito da curva de luz pelo algoritmo MCMC**. *XVIII Jornada de Iniciação Científica*, Universidade Presbiteriana Mackenzie. [Link](https://eventoscopq.mackenzie.br/jornada/pt_BR/article/view/452)

***
## Authors

[ 🇧🇷 ] O modelo de trânsito foi originalmente proposto e desenvolvido por Adriana Valio, em IDL ([Silva 2003](https://doi.org/10.1086/374324)). A versão em Python deste repositório, bem como a continuidade do projeto, foram desenvolvidas a partir dessa base.

[ 🇺🇸 ] The transit model was originally proposed and developed by Adriana Valio, in IDL ([Silva 2003](https://doi.org/10.1086/374324)). The Python version in this repository, as well as the continuation of the project, were developed on that foundation.

 <table>
    <tr>
     <td align="center">
       <a href="https://br.linkedin.com/in/adriana-valio-a5964347">
         <img src="https://www.mackenzie.br/fileadmin/_processed_/c/6/csm_Docente_mackenzista_%C3%A9_nomeada_diretora_do_Laborat%C3%B3rio_Nacional_de_Astrof%C3%ADsica_991e9d41db.png" height="100" width="auto" style="corner-radius:50%" alt="Foto da Profa. Dra. Adriana Valio"/><br>
         <sub>
           <b> Adriana Valio</b>
         </sub>
       </a>
       <br>
       <sub>Original model & IDL code<br>Advisor</sub>
     </td>
     <td align="center">
       <a href="https://github.com/biaduque">
         <img src="https://avatars.githubusercontent.com/u/53840501?v=4" height="auto" width="100" style="corner-radius:50%" alt="Foto de perfil do GitHub da Beatriz Duque"/><br>
         <sub>
           <b> Beatriz Duque</b>
         </sub>
       </a>
       <br>
       <sub>Python version, OOP<br>& documentation</sub>
     </td>
     <td align="center">
       <a href="https://github.com/felipexpinho">
         <img src="https://avatars.githubusercontent.com/u/70541163?v=4" height="auto" width="100" style="corner-radius:50%" alt="Foto de perfil do GitHub do Pinho"/><br>
         <sub>
           <b> Felipe Pinho</b>
         </sub>
       </a>
       <br>
       <sub>MCMC parallelization</sub>
     </td>
  </tr>
  </table>
