## Compile

    git clone --recursive https://github.com/matband/PEF5762_trabalhos/

Go to the folder **Trabalho 2**

     cd Trabalho\ 2/

Run the command g++ to compile

    g++ -I ./eigen -g ./finitedif.cpp -o finitedif.exe

This should create the executable `finitedif.exe`

Running `finitedif.exe` will create `prandtl.csv`, `shear.csv`, `warp.csv` and `convergence.csv` (optional) 

Then, run plot.py with the command

    python3 plot.py
    
It'll generate 3d plots based on csv files:

## Plots 

<img src="https://github.com/matband/PEF5762_trabalhos/blob/main/Trabalho%202/plots/ShearStress.png?raw=true" width = '40%'> <img src="https://github.com/matband/PEF5762_trabalhos/blob/main/Trabalho%202/plots/convergence3DScatter.png?raw=true" width = '40%'>
<img src="https://github.com/matband/PEF5762_trabalhos/blob/main/Trabalho%202/plots/warpColormap.png?raw=true" width = '40%'> 
<img src="https://github.com/matband/PEF5762_trabalhos/blob/main/Trabalho%202/plots/warpSurface.png?raw=true" width = '40%'> 
<img src="https://github.com/matband/PEF5762_trabalhos/blob/main/Trabalho%202/plots/prandtlColormap.png?raw=true" width = '45%'> 
<img src="https://github.com/matband/PEF5762_trabalhos/blob/main/Trabalho%202/plots/membrane3DPlot.png?raw=true" width = '40%'> 
