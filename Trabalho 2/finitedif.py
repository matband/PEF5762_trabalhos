import numpy as np
#Parâmetros do problema
Ly = 1
Lz = 1
G = 10
theta_prime = 1

#Valores de refinamento
n_values = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50]

#Dicionário para armazenar os resultados
results = {}

#Loop sobre os valores de refinamento
for n in n_values:

    # Malha
    hx = Ly/n
    hy = Lz/n
    x = np.linspace(-Ly/2, Ly/2, n+1)
    y = np.linspace(-Lz/2, Lz/2, n+1)
    X, Y = np.meshgrid(x, y)
    # Matriz de coeficientes
    A = np.zeros(((n+1)*(n+1), (n+1)*(n+1)))
    for i in range(n+1):
        for j in range(n+1):
            k = i*(n+1) + j
            A[k, k] = 1
            if i == 0 or i == n:
                B = 2*G*theta_prime*hy/hx
                if j < n:
                    A[k, k+1] = -1
                if j > 0:
                    A[k, k-1] = -1
                A[k, k] = 3
            elif j == 0 or j == n:
                B = 2*G*theta_prime*hx/hy
                if i < n:
                    A[k, k+n+1] = -1
                if i > 0:
                    A[k, k-n-1] = -1
                A[k, k] = 3
            else:
                A[k, k+n+1] = -1
                A[k, k-n-1] = -1
                A[k, k+1] = -1
                A[k, k-1] = -1
                A[k, k] = 4
    # Vetor de termos independentes
    B = np.zeros((n+1)*(n+1))
    for i in range(n+1):
        for j in range(n+1):
            k = i*(n+1) + j
            if i == 0:
                B[k] = 2*G*theta_prime*hy/hx*X[i,j]
            elif i == n:
                B[k] = -2*G*theta_prime*hy/hx*X[i,j]
            elif j == 0:
                B[k] = -2*G*theta_prime*hx/hy*Y[i,j]
            elif j == n:
                B[k] = 2*G*theta_prime*hx/hy*Y[i,j]

    # Solução do sistema linear
    U = np.linalg.solve(A, B)

    # Momento de inércia a torção
    I = 0
    for i in range(n):
        for j in range(n):
            k = i*(n+1) + j
            I += hy*(U[k+1]-U[k])*(X[i,j+1]-X[i,j]) - hx*(U[k+n+1]-U[k])*(Y[i+1,j]-Y[i,j])

    # Armazenando os resultados
    results[n] = I
    #Imprimindo os resultados
for n, I in results.items():
    print(f'n = {n}: I = {I}')