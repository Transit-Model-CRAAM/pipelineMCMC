# Base branch: ```main```
>[!IMPORTANT]
>```main``` represents more stable version of the code
> 
> For a better development experience, we recommend using [VSCode](https://code.visualstudio.com/download), with the extension for [jupyter notebook](https://marketplace.visualstudio.com/items?itemName=ms-toolsai.jupyter), [c and gcc](https://marketplace.visualstudio.com/items?itemName=ms-vscode.cpptools)

# Pipeline MCMC

[ 🇧🇷 ] Repositório com o projeto de programação para calcular a curva de luz de uma estrela com possíveis exoplanetas em sua órbita através de programação Python híbrida para programação C usando Pipeline em tempo de execução.

[ 🇺🇸 ] Repository with the programming project in order to calculate the light curve of a Star with possible exoplanets in its orbit through hybrid Python programming to C programming using Pipeline at runtime. 

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
## Authors

 <table>
    <tr>
      <td align="center">
       <a href="https://github.com/biaduque">
         <img src="https://avatars.githubusercontent.com/u/53840501?v=4" height="auto" width="100" style="corner-radius:50%" alt="Foto de perfil do GitHub da Beatriz Duque"/>          <br>
         <sub>
           <b> Beatriz Duque</b>
         </sub>
       </a>
     </td>
     <td align="center">
       <a href="https://github.com/felipexpinho">
         <img src="https://avatars.githubusercontent.com/u/70541163?v=4" height="auto" width="100" style="corner-radius:50%" alt="Foto de perfil do GitHub do Pinho"/><br>
         <sub>
           <b> Felipe Pinho</b>
         </sub>
       </a>
     </td>
  </tr>
  </table>
