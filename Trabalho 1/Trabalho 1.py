import os
from scipy import linalg
from scipy.sparse import diags
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#######
# Numerical modal analysis 
# Torsional vibration of a prismatic beam with a circular cross section fixed on both sides
# Multi degrees of freedom - DOF = n
# N coupled undamped oscillators
# Ki = K/n ; i = 0,1,2,...,n
# Mi = m/n ; i = 0,1,2,...,n-1
#
# |----ki----Mi----ki--...--ki----Mi----ki----|
#
# 
# G: shear modulus / Ip: polar moment of inertia
# rho = linear mass density
# l = length of the bar
# G*Ip = torsional stiffness
# rho*Ip = masspolar moment of inertia per unit length
# K = tridiagonal matrix


def modal(j,df,freq_conv,modes, pmConv = 0):
    """
    Parameters
    j: iteration of DoF from main function
    df:  pd dataFrame storing frequencies
    freq_conv: ratio of frequencies discrete/theoretical for given DoF
    modes: maximum mode number for the analysis
    pmConv: Point mass inertia problem convergence -> added mass inertia = 2^i , i=0,...,pmConv
    """
    n = 2**(j+2) 
    pi = 3.14159265

    # Data to set: g,rho,l
    g = 200 
    rho = 0.5
    l = 10
    #---------------------#

    
    # Matrices construction
    # K = G /l -> Discrete stiffness
    k = g/l *(n+1)
    # m = rho * l -> Discrete mass 
    m = (rho*l) / n
    
    # Mass matrix -> diagonal size (n x n)
    md = m*np.ones(n)

    # Tridiagonal Stiffness matrix size (n x n)

    # Main diagonal elements of stiffness matrix
    kd = 2*k*np.ones(n)
    
    # Subdiagonal elements of stiffness matrix
    ke = -1*k*np.ones(n-1)
    
       ######
    # Stiffnes Matrix  size(n x n)
    # [2k -k  0  0 ... 0]
    # [-k 2k -k  0 ... 0]
    # [ 0 -k 2k -k ... 0]
    # [       ...       ]
    # [0  ...  0 0 -k 2k]
    
    # pointMassConvergence:
    # added mass inertia at l/2 
    # convergence of frequencies -> mass = 2^i , i=0,1,...,9

    if(pmConv != 0):
        for ai in range(pmConv):
            ad = 2**ai
            coefs = [ke,kd,ke]
            offset = [-1,0,1]
            K = diags(coefs,offset).toarray()
            
            md[int(n/2)-1] += ad
            M = np.diag(md)
            eigvals,eigvecs = linalg.eigh(K,M)
            eigvecs = np.vstack([eigvecs, np.zeros((1,n))])
            eigvecs = np.vstack([np.zeros((1,n)),eigvecs])
            freq = np.sqrt(eigvals)
            freq_s = list(map(lambda i: f"{i/pi:.2f} pi" if (i==i) else "-" ,list(np.resize(freq,modes))))

            # print(f'\n Values for added mass inertia of {ad}')
            # for i in range(modes):
            #     print(f'mode {i} {ad};',freq[i])
        
            if(ai == 0):
                theo = [f'{4*((mode+1)//2)} pi' for mode in range(modes)]
                df.loc[0] = ['Theoretical:\n'+r'$mass-> \infty $'] + theo

            df.loc[ai+1] = [str(ad)] + freq_s

            # Plotting normal modes of vibration
            # plt.clf()
            abcissa = np.arange(0,n+2)
            abcissa = l * abcissa / n 
            for i in range(modes):
                plt.plot(abcissa,eigvecs.T[i], label = f"mode {i}")
            plt.xlabel("x")
            plt.ylabel(r'$\theta (x)$')
            plt.legend(loc="lower right")
            plt.title(f"Normal modes: {n} DoF - {ad} added mass inertia")
            plt.savefig(absolute_path+f"/normalmodes_addedmass/modes_am_{ad}", bbox_inches="tight")
            plt.clf()
                    
            
    else:
        # Matrix A of Standard Eigenvalue Problem
        # A {phi} = omega^2 {phi} ; A = [M]^-1 * [K]
        ad = np.divide(kd,md)
        ae = np.divide(ke,md[1:])

        # Solving Standard Eigenvalue Problem for tridiagonal symmetric matrices by scipy linalg lib
        eigvals, eigvecs = linalg.eigh_tridiagonal(ad, ae, eigvals_only=False)

        # Adding phi(0)=0 e phi(l)=0 
        eigvecs = np.vstack([eigvecs, np.zeros((1,n))])
        eigvecs = np.vstack([np.zeros((1,n)),eigvecs])

        # natural frequencies = square root of eigenvalues
        freq = np.sqrt(eigvals)
        nullFreqs = modes-freq.size
        if (nullFreqs < 0 ):
            nullFreqs = 0
        for nullValues in range (nullFreqs):
            freq = np.append(freq,[np.nan])
        
        # Storing convergence values:
        # Discrete/Continuous for given normal modes of vibration 
        conv = np.zeros(modes)
        theoretical = np.zeros(modes)
        for i in range(modes):
            theoretical[i] = (((i+1)*pi/l)* np.sqrt(g/rho))
            conv[i] = freq[i]/theoretical[i]

        # Values to string
        theoretical_s = list(map(lambda i: f"{i/pi:.2f} pi",theoretical))
        freq_s = list(map(lambda i: f"{i/pi:.2f} pi" if (i==i) else "-" ,list(np.resize(freq,modes))))

        # Setting up dataframe
        if(j == 0):
            df.loc[0] = ['Theoretical'] + theoretical_s
        df.loc[j+1] = [str(n)] + freq_s

        # Print of numerical and theoretical values (optional)
        # if(j == 0):
        #     print(f"theoretical values:")
        #     for i in range(0,modes):
        #         print(f"omega {i+1} = {theoretical[i]/pi:.3f} * pi")
        #     print("numerical values:")
        # print(f"{n} DoF")
        # for i in range(0,modes):
        #     print(f"omega {i+1} = {freq[i]/pi:.3f} * pi")

        # Plotting natural modes of torsional vibration
        plt.clf()
        abcissa = np.arange(0,n+2)
        abcissa = l * abcissa / (n+1) 
        for i in range(modes-nullFreqs):
            plt.plot(abcissa,eigvecs.T[i], label = f"mode {i}")

        plt.legend(loc="lower right")
        plt.xlabel("x")
        plt.ylabel(r'$\theta (x)$')
        plt.title(f"Normal modes of vibration: {n} degrees of freedom")
        plt.savefig(absolute_path+f"/normalmodes/modes_{n}_dof", bbox_inches="tight")
        if(j<len(freq_conv)):
            freq_conv[j] = conv

def main():
    # n: DoF refinement
    # DoF = 2*(j+2); for j = 0,1,2,...,n-1; 
    # e.g.: 
    # j = 0 => DoF = 4
    # j = 7 => DoF = 2^9 = 512 
    
    modes = int(input("Insert maximum mode number for the analysis: \n"))
    print('Maxima degrees of freedom = 2^(N+2); N<9')
    n = int(input("Insert N: \n"))


    # Setting up the table to export data
    normalModeCols = ['DoF']
    for i in range(modes):
        normalModeCols.append(f'omega {i+1}')
        
    df = pd.DataFrame(columns=normalModeCols)

    freq_conv = np.zeros((n+1,modes))
    for i in range(0,n+1):
        modal(i,df,freq_conv,modes)

    # Generating Table
    fig, ax = plt.subplots(1, 1)
    ax.table(cellText=df.values, colLabels=df.keys(), loc='center')
    plt.axis('off')
    plt.savefig(absolute_path+"/freqtable", bbox_inches="tight")
    # Generating plots
    abcissa = np.arange(0,n+1) 
    abcissa = [2**(j+2) for j in range(n+1)]
    plt.clf()
    for i in range(modes):
        plt.plot(abcissa, freq_conv.T[i], label = f"mode {i}")
    plt.xscale('log',base=2)
    plt.title("Convergence of theoretical and numerical frequencies ratio")
    plt.legend(loc="lower right")
    plt.xlabel("Degrees of Freedom")
    plt.savefig(absolute_path+"/convergence", bbox_inches="tight")
    plt.clf()
    #############
    # Same problem with point mass inertia at l/2
    print(f"Point mass inertia at l/2 for {2**(n+2)} DoF")
    print('Convergence for exponentially progressive values of mass inertia from 2^0 to 2^9')
    option = int(input("Enter 1 to save the results ; 0 to close \n"))
    if (option == 1):
        modal(n,df,freq_conv,modes,10)
        fig, ax = plt.subplots(1, 1)
        ytable = ax.table(cellText=df.values, colLabels=df.keys(), loc='center')
        for r in range(0, len(df.keys())):
            cell = ytable[1, r]
            cell.set_height(0.1)
        plt.axis('off')
        plt.savefig(absolute_path+"/freqtable_addedmass", bbox_inches="tight")
    print("Data stored on "+absolute_path)
    input("Press enter to close...")

####
#Absolute plots path
absolute_path = os.path.dirname(__file__) + "/plots"
if not os.path.exists(absolute_path):
    os.mkdir(absolute_path)
    if not os.path.exists(absolute_path+'/normalmodes'):
        os.mkdir(absolute_path+'/normalmodes')
    if not os.path.exists(absolute_path+'/normalmodes_addedmass'):
        os.mkdir(absolute_path+'/normalmodes_addedmass')


main()