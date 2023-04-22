#include <iostream>
#include <Eigen/Dense>
#include <omp.h>
#include <fstream>
#include <chrono>
using namespace Eigen;


 
MatrixXd recover_psi(const MatrixXd& dphi_dy, const MatrixXd& dphi_dz, const double dy, const double dz, const int ny, const int nz,double G, double theta_prime)
{
    //Relating phi to the partial derivatives of psi
    MatrixXd dpsi_dy = MatrixXd::Zero(ny, nz);
    MatrixXd dpsi_dz = MatrixXd::Zero(ny, nz);
    
    for (int i = 1; i < ny ; i++) {
        for (int j = 1; j < nz ; j++) {
            dpsi_dy(i,j) = (dphi_dz(i,j)/G/theta_prime + (j-nz/2)*dz); // ∂φ/∂z=Gθ' (∂ψ/∂y - z)
            dpsi_dz(i,j) = (-1 * dphi_dy(i,j)/G/theta_prime - (i-ny/2)*dy); // ∂φ/∂y=-Gθ' (∂ψ/∂z + y)
        }
    }

    // psi1: integral of dpsi_dy dy
    // psi2: integral of dpsi_dz dz

    MatrixXd psi1 = MatrixXd::Zero(ny, nz);
    MatrixXd psi2 = MatrixXd::Zero(ny, nz);

    for (int i = 1; i < ny; ++i)
    {
        for (int j = 1; j < nz; ++j)
        {
            // Integral em y
            psi1(i, j) = psi1(i-1, j) + dy * dpsi_dy(i, j);
            // Integral em z
            psi2(i, j) = psi2(i, j-1) + dz * dpsi_dz(i, j);
        }
    }
    MatrixXd psi = (psi1+psi2)/2;

    return psi;
}

void grad(MatrixXd &dfunc_dy, MatrixXd &dfunc_dz, MatrixXd &func, const int ny, const int nz, const double dy, const double dz){
    for (int i = 1; i < ny - 1; i++) {
        for (int j = 1; j < nz - 1; j++) {
            
            dfunc_dy(i,j) = (func(i+1,j) - func(i-1,j)) / (2*dy);
            dfunc_dz(i,j) = (func(i,j+1) - func(i,j-1)) / (2*dz);
        }
    }
}

MatrixXd getphi(double Ly, double Lz,double gt, int ny, int nz, double dy, double dz){

    
    MatrixXd A = MatrixXd::Zero(ny*nz, ny*nz);
    VectorXd b = VectorXd::Zero(ny*nz);

    // Second order central difference 

    for (int i = 1; i < ny - 1; i++) {
        for (int j = 1; j < nz - 1; j++) {
            int k = i*nz + j;  // Index of point (i,j) in matrix
            A(k,k-nz) = 1.0 / (dy*dy);  // (i-1,j)
            A(k,k-1) = 1.0 / (dz*dz);  // (i,j-1)
            A(k,k) = -2.0 / (dy*dy) - 2.0 / (dz*dz);  // (i,j)
            A(k,k+1) = 1.0 / (dz*dz);  // (i,j+1)
            A(k,k+nz) = 1.0 / (dy*dy);  // (i+1,j)
            b(k) = -2.0*gt;
        }
    }

    // Boundary Conditions
    
    for (int i = 0; i < ny; i++) {
        int k1 = i*nz;  // Index of point (i,0) in matrix
        int k2 = i*nz + nz - 1;  // Index of point (i,nz-1) in matrix

        // Dirichlet Boundary conditions on Linear System
        A(k1,k1) = 1.0; 
        A(k2,k2) = 1.0;
        b(k1) = 0.0;
        b(k2) = 0.0;
    }
    
    for (int j = 0; j < nz; j++) {
        int k1 = j;  // Index of (0,j) in matrix
        int k2 = (ny-1)*nz + j;  // Index of (ny-1,j) in matrix

        // Dirichlet Boundary conditions on Linear System
        A(k1,k1) = 1.0;
        A(k2,k2) = 1.0;
        b(k1) = 0.0;
        b(k2) = 0.0;
    }
    
    // Phi from linear system
    VectorXd phi = A.colPivHouseholderQr().solve(b);
    MatrixXd phiM = phi.reshaped(nz,ny).transpose();

    return phiM;

}

int main()
{
    // Parameters
    double Ly, Lz, dy, dz, G, theta_prime, dA, gt;
    int ny, nz, op;
    
    //////////////
    
    std::cout << "Enter '1' to use default values (Ly=1, Lz=1.2, ny=40, nz=40, G=10, theta_prime=1): "<< std::endl;
    std::cout << "Or enter any other key to set custom values." << std::endl;
    std::cin >> op;
    if (op == 1){
        Ly = 1;
        Lz = 1.2;
        ny = 40;
        nz = 40;
        G = 10;
        theta_prime = 1;
    }else{
        std::cout << "Ly: ";
        std::cin >> Ly;
        std::cout << "Lz: ";
        std::cin >> Lz;
        std::cout << "ny: ";
        std::cin >> ny;
        std::cout << "nz: ";
        std::cin >> nz;
        std::cout << "G: ";
        std::cin >> G;
        std::cout << "theta_prime: ";
        std::cin >> theta_prime;
    }

    //////////////
    dy = Ly / (ny - 1);
    dz = Lz / (nz - 1);
    dA = dz*dy;
    gt = G*theta_prime;
    // Calculating Prandtl stress function by finite difference method
    MatrixXd phiM = getphi(Ly, Lz, gt, ny, nz, dy, dz);
    
    // Rectangle method
    double T = 0.0;
    for (int i = 0; i<ny;i++){
        for (int j = 0; j<nz ;j++) 
            T += phiM(i,j)*dA;
    }
    T*=2;
    double Ip = T /(G*theta_prime);

    printf("IP = %.4f \n", Ip);
    printf("T = %.4f \n", T);
    
    std::ofstream myFile4("config.txt");
    myFile4 << Ly << std::endl << Lz << std::endl << ny << std::endl << nz << std::endl;
    myFile4 << G << std::endl << theta_prime << std::endl << T << std::endl << Ip;
    myFile4.close();

    //Partial first order derivatives 
    MatrixXd dphi_dy = MatrixXd::Zero(ny, nz);
    MatrixXd dphi_dz = MatrixXd::Zero(ny, nz);
    //Grad phi
    grad(dphi_dy,dphi_dz,phiM,ny,nz,dy,dz);

    // Calculating shear stress at cross section
    MatrixXd shearStress = dphi_dy.array()*dphi_dy.array()+dphi_dz.array()*dphi_dz.array();
    shearStress = shearStress.cwiseSqrt();

    MatrixXd psi = MatrixXd::Zero(ny, nz);
    psi = recover_psi(dphi_dy,dphi_dz,dy,dz,ny,nz, G, theta_prime);

    // Convergence
    int conv;
    std::cout << "Enter 1 for convergence analysis..." << std::endl;
    std::cin >> conv;
    if (conv == 1){
        int n = 7; // n x n convergence iterations
        MatrixXd convergence = MatrixXd::Zero(n+1,n+1);
        for (int i = 0; i < n;i++){
            ny = 5*(i+1);
            dy = Ly / (ny - 1);
            for(int j = 0; j < n;j++){
                nz = 5*(j+1);
                dz = Lz / (nz - 1);
                dA = dz*dy;
                MatrixXd phiC = MatrixXd::Zero(ny, nz);
                phiC = getphi(Ly, Lz, gt, ny, nz, dy, dz);
                // Rectangle method
                double Tsum = 0.0;
                for (int i = 0; i<ny;i++){
                    for(int j = 0; j< nz;j++){
                        Tsum += phiC(i,j)*dA;
                    }
                }
                convergence(0,j+1) = nz;
                convergence(i+1,j+1) = 2*Tsum/(gt);
            }
            convergence(i+1,0) = ny;
        }    
        std::ofstream myFilec("convergence.csv");
        myFilec << convergence;
        myFilec.close();
    }
    

    //Generating shear and warp csv
    std::ofstream myFile("shear.csv");
    myFile << shearStress;
    myFile.close();
    std::ofstream myFile2("warp.csv");
    myFile2 << psi;
    myFile2.close();
    std::ofstream myFile3("prandtl.csv");
    myFile3 << phiM;
    myFile3.close();
    return 0;

}

