import os
from math import sqrt
from scipy import linalg
import numpy as np
import matplotlib.pyplot as plt

#######
# Barra bi engastada de seção circular prismática homogênea
# Para n graus de liberdade 
# Sistema massa mola com n+1 molas e n massas
# Ki = K/n ; i = 0,1,2,...,n
# Mi = m/n ; i = 0,1,2,...,n-1
# G*Ip = Módulo de elasticidade transversal * Momento polar de inércia = Módulo de rigidez à torção- Dado
# rho = Massa linear da barra - Dado
# l = Comprimento da barra - Dado

# K = matriz tridiagonal

def modal(j,freq_conv, problemadois = False):

    n = 2**(j+2) 
    pi = 3.14159265


    # Inserir dados 
    #
    gip = 200 # G Ip
    rho = 0.5
    l = 10
    

    # Montagem das matrizes
    # K = g Ip /l 
    k = gip/l *(n+1)
    # m = rho * l -> massa discreta
    m = (rho*l) / n
    
    # Matriz diagonal de massas
    md = m*np.ones(n)

    # Matriz tridiagonal simétrica de rigidez

    #Elementos da diagonal da matriz tridiagonal simétrica de rigidez K
    kd = 2*k*np.ones(n)
    inerciapontual = gip
    if(problemadois):
        kd[int(n/2)] += inerciapontual
        kd[int(n/2)-1] += inerciapontual
    
    #Elementos fora da diagonal da matriz tridiagonal simétrica de rigidez K
    ke = -1*k*np.ones(n-1)
    if(problemadois):
        
        kd[int(n/2)] += inerciapontual
        kd[int(n/2)-1] += inerciapontual


    ######
    # Matriz de rigide de ordem n+1
    # [2k -k  0  0 ... 0]
    # [-k 2k -k  0 ... 0]
    # [ 0 -k 2k -k ... 0]
    # [       ...       ]
    # [0  ...  0 0 -k 2k]


    # Matriz A da resolução do problema de autovalor clássico
    # A {phi} = omega^2 {phi} ; A = [M]^-1 * [K]
    ad = np.divide(kd,md)
    ae = np.divide(ke,md[1:])

    # Resolução do problema de autovalor clássico para matrizes tridiagonais simétricas pela biblioteca scipy linalg
    eigvals, eigvecs = linalg.eigh_tridiagonal(ad, ae, eigvals_only=False)

    # Adicionando phi(0) e phi(l) na discretização
    eigvecs = np.vstack([eigvecs, np.zeros((1,n))])
    eigvecs = np.vstack([np.zeros((1,n)),eigvecs])

    # frequências naturais = raiz quadrada dos autovalores
    freq = np.sqrt(eigvals)

    # Guardando os valores de convergência do modelo discreto em relação ao modelo contínuo para os três primeiros modos de Vibração 
    conv = np.zeros(4)
    analitico = np.zeros(4)
    for i in range(0,3):
        if(problemadois):
            conv[i] = freq[i]
        else:
            analitico[i] = (((i+1)*pi/l)* sqrt(gip/rho))
            conv[i] = freq[i]/analitico[i]

    # Impressão dos valores numéricos e analíticos (opcional)
    
    print(f"resultado numérico: para {n} graus de liberdade \n")
    for i in range(0,4):
        print(f"omega {i+1} = ", freq[i])


    # Plotando os modos naturais de vibração
    plt.clf()
    abcissas = np.arange(0,n+2)
    abcissas = l * abcissas / n 
    modo = np.ones(n)
    for i in range(n):
        modo[i] = np.sin((pi*i)/l)
    for i in range(3):
        plt.plot(abcissas,eigvecs.T[i], label = f"modo {i}")
    plt.legend(loc="best")
    plt.xlabel("x")
    plt.ylabel(r'$\theta (x)$')
    plt.title(f"Modos naturais de vibração para {n} graus de liberdade")
    plt.savefig(absolute_path+f"/modo_natural_{n}", bbox_inches="tight")
    freq_conv[j] = conv

def main():
    # n: quantidade de repetições do refinamento da discretização 
    # Graus de liberdade = 2*(j+2); j = 0,1,2,...,n-1; 
    # para j = 0 => GLs = 4
    # para j = 7 => GLs = 2^9 = 512 
    n = 7
    
    freq_conv = np.zeros((n+1,4))
    for i in range(n+1):
        # Para adição de inércia pontual, Chamar a função com a flag True
        modal(i,freq_conv,False)
    abcissas = np.arange(0,n+1) 
    abcissas = [2**(j+2) for j in range(n+1)]
    plt.clf()
    for i in range(3):
            plt.plot(abcissas, freq_conv.T[i], label = f"modo {i}")
    plt.xscale('log',basex=2)
    plt.title("Convergência da razão entre as soluções numérica e analítica")
    plt.legend(loc="lower right")
    plt.xlabel("Graus de liberdade")
    plt.savefig(absolute_path+"/Convergência", bbox_inches="tight")


####
#caminho absoluto para salvar figuras
absolute_path = os.path.dirname(__file__) + "/plots"
if not os.path.exists(absolute_path):
    os.mkdir(absolute_path)

main()