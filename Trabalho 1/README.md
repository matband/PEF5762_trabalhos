<h1> Trabalho 1 </h1> 

<h2>Torsional vibration of a prismatic beam with a circular cross section fixed on both sides </h2>

<h3>Numerical Values:</h3> 

<strong>Degrees of freedom:</strong><br> $DoF_i = 2^{i+2} \ , \ i=0,1,...,n$ 
<br><strong>Linear mass density:</strong><br> $\rho = 0.5$
<br><strong>Shear Modulus:</strong><br> $G = 200$ 
<br><strong>Bar length:</strong><br> $l = 10$ 
<br><strong>Added mass inertia at $\dfrac{l}{2}$ :</strong><br> $M_j = 2^{j} \ , \ j=0,1,...,9$  
<h3>Convergence: Degrees of Freedom x Ratio</h3>

$Ratio(i) = \dfrac{DiscreteFreq_m}{TheoreticalFreq_m}$ , $m=1,...,modes$

<h2>Results:</h2>
<h3>Frequencies:</h3>
<div>
<img src="https://github.com/matband/PEF5762_trabalhos/blob/main/Trabalho%201/plots/freqtable.png?raw=true" alt="Frequency table" 
style=" width:40%;">
<img src="https://github.com/matband/PEF5762_trabalhos/blob/main/Trabalho%201/plots/convergence.png?raw=true" alt="Frequency convergence" 
style=" width:35%;">  
</div>

<h3>Modes:</h3>
<div>
<img src="https://github.com/matband/PEF5762_trabalhos/blob/main/Trabalho%201/plots/normalmodes/modes_4_dof.png?raw=true" alt="modes for 4 dof" 
style=" width:37%;">
<img src="https://github.com/matband/PEF5762_trabalhos/blob/main/Trabalho%201/plots/normalmodes/modes_1024_dof.png?raw=true" alt="modes for 1024 dof" 
style=" width:37%;">  
</div>
<h3>Added mass inertia problem:</h3>
<div>
<img src="https://github.com/matband/PEF5762_trabalhos/blob/main/Trabalho%201/plots/freqtable_addedmass.png?raw=true" alt="frequency table added mass inertia" 
style=" width:40%;">
<img src="https://github.com/matband/PEF5762_trabalhos/blob/main/Trabalho%201/plots/normalmodes_addedmass/modes_am_16.png?raw=true" alt="modes for 512 dof 16 added mass inertia" 
style=" width:37%;">
</div>
