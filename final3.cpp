#include <iostream> // cout
#include <cmath> // math lib
#include <fstream> // file I/O
#include <stdlib.h> // malloc
#include <unistd.h> // sleep
#include <string> // string
using namespace std;




void visualize(double **var, int nx, int ny){
  // visualization
  for(int i = 0; i <= nx-1; i++){
    for(int j = 0; j <= ny-1; j++){    
      cout << var[i][j] << "\t";          // \n  <= newline, \t <= tab
    }
    cout << "\n";
  }
  cout << "\n";
}



void update(double **ux, double **ux_new,double **uy, double **uy_new, int nx, int ny){
  for(int i = 0; i <= nx-1; i++){
    for(int j = 0; j <= ny-1; j++){    
      ux[i][j] = ux_new[i][j];
      uy[i][j] = uy_new[i][j];
    }
  }
}

void update_two(double **P, double **P_new, int nx, int ny){
  for(int i = 0; i <= nx-1; i++){
    for(int j = 0; j <= ny-1; j++){
      P[i][j] = P_new[i][j];
    }
  }
}


void navier_one(double **ux, double **ux_star,double **uy, double **uy_star
		, int nx, int ny, double c, double dx, double dy, double dt,double re){
  //double RHS;
  // simulation
 
  for(int i = 1; i <= nx-2; i++){
    for(int j = 1; j <= ny-2; j++){
      double na_1_1;
      double na_1_2;
      double na_1_3;
      double na_1_4;
      double na_1_5;
      double na_1_6;
      double na_1_7;
      double na_1_8;
      
      
      na_1_1 = (ux[i+1][j]-(2*ux[i][j])+ux[i-1][j])/(dx*dx);
      na_1_2 = (ux[i][j+1]-(2*ux[i][j])+ux[i][j-1])/(dy*dy);
      na_1_3 = (1/dx)*(((ux[i][j]+ux[i+1][j])/2)*((ux[i][j]+ux[i+1][j])/2)-((ux[i-1][j]+ux[i][j])/2)*((ux[i-1][j]+ux[i][j])/2));
      na_1_4 = ((((uy[i][j]+uy[i+1][j])/2)*((ux[i][j]+ux[i][j+1])/2))-(((uy[i][j-1]+uy[i+1][j-1])/2)*((ux[i][j-1]+ux[i][j])/2)))/dy;
      na_1_5 = (uy[i+1][j]-(2*uy[i][j])+uy[i-1][j])/(dx*dx);
      na_1_6 = (uy[i][j+1]-(2*uy[i][j])+uy[i][j-1])/(dy*dy);
      na_1_7 = ((((ux[i][j]+ux[i][j+1])/2)*((uy[i][j]+uy[i+1][j])/2))-(((ux[i-1][j]+ux[i-1][j+1])/2)*((uy[i-1][j]+uy[i][j])/2)))/dx;
      na_1_8 = ((((uy[i][j]+uy[i][j+1])/2)*((uy[i][j]+uy[i][j+1])/2))-(((uy[i][j-1]+uy[i][j])/2)*((uy[i][j-1]+uy[i][j]))))/dy;
	 
      ux_star[i][j] = ux[i][j]+ dt*(((na_1_1+na_1_2)/re)-na_1_3-na_1_4);
      uy_star[i][j] = uy[i][j]+ dt*(((na_1_5+na_1_6)/re)-na_1_7-na_1_8);
    }
  }

}

void initialize(double **ux, double **ux_new,
		double **uy, double **uy_new,
		double **P,int nx, int ny){
  // initialize var
  for(int i = 0; i <= nx-1; i++){
    for(int j = 0; j <= ny-1; j++){    
      ux[i][j] = 0;
      ux_new[i][j] = 0.;
      uy[i][j] = 0.;
      uy_new[i][j] = 0.;
      P[i][j] = 1;
    }
  }
}

void boundary_condition (double **ux,double **uy,double **P, int nx , int ny, double **F, double **G ){
 for(int i = 0; i <= nx-1; i++){
    for(int j =0; j <= ny-1; j++){
    // ux inlet = 1
    if(i ==0 & j>= 1 & j<= ny-2 ){
      ux[i][j] = 1;
    }
    // no slip cond. @ lower wall
    if(i>=1 & i <= nx-2 & j == 0){
      uy[i][j] = 0;
      ux[i][j] = -ux[i][j+1];
    }
    // no slip cond. @ upper wall for v
    if(i>=1 & i <= nx-2 & j == ny-2){
      uy[i][j] = 0;

    }
    // no slip cond.c@ upper wall for u
    if(i>= 1 & i <= nx-2 & j == ny-1){
       ux[i][j] = -ux[i][j-1];

    }
    // outflow for ux
    if( j >= 1 & j <= ny-2 & i == nx-1){
       ux[i][j] = ux[i-1][j];
    }
    // outflow for uy
    if(j>= 0  & j <= ny-2   & i == nx-1){
      uy[i][j] = uy[i-1][j];
    }
    // pressure for inlet
    if(i == 0 & j>= 1 & j <= ny-2){
      P[i][j] = P[i+1][j];
    }
    // pressure for lower
    if(j == 0 & i>=1 & i <= nx-2){
      P[i][j] = P[i][j+1];
    }
    // pressure for outlet
    if(i == nx-1 & j>= 1 & j <= ny-2){
      P[i][j] = P[i-1][j];
    }
    // pressure for upper
    if(j == ny-1 & i>=1 & i <= nx-2){
      P[i][j] = P[i][j-1];
    }
    // F inlet
    if(i == 0 & j>=1 & j <= ny-2){
      F[i][j] = ux[i][j];
    }
    // F outlet
    if(i == nx-1 & j>= 1 & j <= ny-2){
      F[i][j] = ux[i][j];
    }
    // G lower
    if(i >=1 & i <= nx-2 & j == 0){
      G[i][j] = uy[i][j];
    }
    // G upper
    if(i >=1 & i <= nx-2 & j == ny-1){
      G[i][j] = uy[i][j];
    }

       }    
   }
}

void navier_two(double **F,double **G,double **P,double **P_new,double dx,double dy,double dt,int nx,int ny){
        int itmax=999;
        int it=0;
        double rit=1.;
        double eps=0.00001;
        double w=1.7;
        int ew=0;int ee=0;int es=0;int en=0;
        double RHS=0.;
        double pnorm=0.;
        //set boundary
        
        while ((it<itmax) && (rit>eps*pnorm )){
            rit =0.;  pnorm=0.;
            //set boundary
            for (int j =1; j<= ny ;j++){
              P_new[0][j]=P[1][j]; 
              P_new[nx+1][j]=P[nx][j];
            }
            for (int i =1; i<= nx ;i++){
              P_new[i][0]=P[i][1]; 
              P_new[i][ny+1]=P[i][ny];
            }
            
            
            for (int i=1 ;i<=nx;i++){
                for (int j=1;j<=ny;j++){
                    if (i==1){
                        ew=0;}
                    else{
                        ew=1;}
                    if (i==nx){
                        ee=0;}
                    else{
                        ee=1;}
                    if (j==1){
                        es=0;}
                    else{
                        es=1;}
                    if (j==ny){
                        en=0;}
                    else{
                        en=1;}
                   //
                      
                      
                    RHS=((F[i][j]-F[i-1][j])/dx + (G[i][j]-G[i][j-1])/dy)/dt;
                    P_new[i][j]=(1-w)*P[i][j] +  (w/((ee+ew)/pow(dx,2) + (en+es)/pow(dy,2)))*((ee*P[i+1][j]+ew*P_new[i-1][j])/pow(dx,2) + (en*P[i][j+1]+es*P_new[i][j-1])/pow(dy,2) - RHS);
                    rit=rit+pow( (ee*(P[i+1][j]-P[i][j])-ew*(P[i][j]-P[i-1][j]))/pow(dx,2) + (en*(P[i][j+1]-P[i][j])-es*(P[i][j]-P[i][j-1]))/pow(dy,2) - RHS,2);
                    pnorm=pnorm+pow(P_new[i][j],2);
                    //cout<<"RHS = "<<RHS;
                    
                    }           
            } 
            
            rit=sqrt(rit/(nx*ny));
            pnorm=sqrt(pnorm/(nx*ny));
            //cout<<"\n"<<it<<"\t"<<rit<<"\n"<<"pnorm "<<pnorm;
            //cout <<"<<<<<<<< P SOR>>>>>>>>>>>>>>>>>>>";
            //test(P,nx,ny);
            //cout <<"<<<<<<<< P NEW SOR>>>>>>>>>>>>>>>>>>>";
            //test(P_new,nx,ny);
            update_two(P,P_new,nx,ny);
            
            
           
            it=it+1;
        }
}

void navier_two(double **P, double **P_new, int nx, int ny, double **F, double **G,  double dx, double dy, double dt){
  double RHS=0.;
  double M =1;
  double Mmax = 0.001;
  double W;
  double E;
  double S;
  double N;
  double w = 1.7;
  double R = 0; 

  //for(int i = 1; i <= nx-2; i++){
  //for(int j =1; j <= ny-2; j++){
  //P[i][j] = 1;
  //}}
  while(M > Mmax && R <= 1000 ){
    double  Mloop = 0;
    R = R+1;
  for(int i = 1; i <= nx-2; i++){
    for(int j = 1; j <= ny-2; j++){
      if(i == 1){
        W = 0;
      }
      if(i > 1){
        W = 1;
      }
      if(i < nx-2){
        E = 1;
      }
      if(i == nx-2){
        E = 0;
      }
      if(j == 1){
	S = 0;
      }
      if(j>1){
	S = 1;
      }
      if(j < ny-2){
	N = 1;
      }
      if(j == ny-2){
        N = 0;
      }
      RHS = (1/dt)*((F[i][j]-F[i-1][j])/(dx)+(G[i][j]-G[i][j-1])/(dy));
      P_new[i][j] = (1-w)*P[i][j]+(w/((E+W)/(dx*dx)+(N+S)/(dy*dy)))*(((E*P[i+1][j]+W*P_new[i-1][j])/(dx*dx))+((N*P[i][j+1]+S*P_new[i][j-1])/(dy*dy))-RHS);
      M = abs(P_new[i][j]-P[i][j])/P[i][j];
      if( M > Mloop){
	Mloop = M;
      }
      // P[i][j] = P_new[i][j];
    }
  }
  M=Mloop;
  for(int i = 1; i <= nx-2; i++){
   for(int j = 1; j <= ny-2; j++){
     P[i][j] = P_new[i][j];
   }}
   }
}
   


void navier_three(double** P,double** F,double** G,double** ux_new,double** uy_new,double dx,double dy,double dt,int nx,int ny){
  double RHS;
  for(int i =1;i<=nx-2;i++){
    for (int j=1;j<=ny-2;j++){                                                         
      ux_new[i][j]=F[i][j]-dt*((P[i+1][j]-2*P[i][j]+P[i-1][j])/(dx)); //from explicit euler
      uy_new[i][j]=G[i][j]-dt*((P[i][j+1]-2*P[i][j]+P[i][j-1])/(dy)); 
    }
  }}


void paraview(string fileName, double **var, int nx, int ny, double dx, double dy){
  ofstream myfile;
  myfile.open(fileName);
  //------------------------------------------------------------//                                                                                                                                                                                                                                               
    // Paraview header                                                                                                                                                                                                                                                                                           
  myfile << "# vtk DataFile Version 2.0\n";
  myfile << "FlowField\n";
  myfile << "ASCII\n";

    // Grid                                                                                                                                                                                                                                                                                                      
  myfile << "DATASET STRUCTURED_GRID\n";
  myfile << "DIMENSIONS " << nx << " " << 1 << " " << ny << "\n";
  myfile << "POINTS " << nx*1*ny << " float\n";
  for(int j = 0; j <= ny-1; j++){
    for(int i = 0; i <= nx-1; i++){
      myfile << dx*i << " " << dy*j << " 0\n";
    }
  }

  // Data                                                                                                                                                                                                                                                                                                        
  myfile << "\n";
  myfile << "POINT_DATA";
  myfile << " " << nx*ny << "\n";

  myfile << "\n";
  myfile << "SCALARS PHI float 1\n";
  myfile << "LOOKUP_TABLE default\n";
  for(int j = 0; j <= ny-1; j++){
    for(int i = 0; i <= nx-1; i++){
      myfile << var[i][j] << "\n";
    }
  }
  myfile.close();
}



int main(){

  int nx = 151;
  int ny = 11;  
  double k = 1.;
  double dx = 0.05;
  double dy = 0.01;  
  double dt = 0.01;
  string fileName; 
  double re = 300;
  //-----------------------------------
  double **ux;
  ux = (double **) malloc (nx * sizeof(double));
  for(int i = 0; i < nx; i++){
    ux[i] = (double *) malloc (ny * sizeof(double));
  }
  double **ux_new;
  ux_new = (double **) malloc (nx * sizeof(double));
  for(int i = 0; i < nx; i++){
    ux_new[i] = (double *) malloc (ny * sizeof(double));
  }
  
  double **F;
  F = (double **) malloc (nx * sizeof(double));
  for(int i = 0; i < nx; i++){
    F[i] = (double *) malloc (ny * sizeof(double));
  }

  double **uy;
  uy = (double **) malloc (nx * sizeof(double));
  for(int i = 0; i < nx; i++){
    uy[i] = (double *) malloc (ny * sizeof(double));
  }
  double **uy_new;
  uy_new = (double **) malloc (nx * sizeof(double));
  for(int i = 0; i < nx; i++){
    uy_new[i] = (double *) malloc (ny * sizeof(double));
  }
  
  double **G;
  G = (double **) malloc (nx * sizeof(double));
  for(int i = 0; i < nx; i++){
    G[i] = (double *) malloc (ny * sizeof(double));
  }
   double **P;
  P = (double**) malloc (nx*sizeof(double));
  for(int i =0;i<=nx;i++){
    P[i]=(double*) malloc (ny*sizeof(double));
  }
  double **P_new;

  P_new = (double**) malloc (nx*sizeof(double));
  for(int i =0;i<=nx;i++){
    P_new[i]=(double*) malloc (ny*sizeof(double));
  }

  //----------------------------------

  initialize(ux, ux_new,
	     uy, uy_new,
	     P, nx, ny);
  //
  for(int n = 1; n<= 10000  ; n++){
    boundary_condition(ux_new,uy_new,P, nx , ny, F, G);

    navier_one(ux, F,uy, G , nx, ny,  k,  dx,  dy, dt, re);

    // navier_two(F, G, P, P_new, dx, dy, dt);
    navier_two(F, G, P, P_new, dx, dy, dt, nx, ny);
    navier_three(P_new,F, G, ux_new, uy_new, dx, dy, dt,nx, ny);
  
    update(ux,ux_new,uy, uy_new,  nx,ny);

    visualize(ux, nx, ny);

    cout << n << "\n";
        if( n%50 == 0) {
      string fileName = "ux_" + to_string(n) + ".vtk";
      paraview(fileName, ux, nx, ny, dx, dy);
      }
  }
  //paraview(phi, nx, ny, dx, dy)

  
  // begin time loop
  //for(int n = 1; n <= 100000; n++){
  //  simulation(phi, phi_new, nx, ny, k, dx, dy, dt);
    // update(phi, phi_new, nx, ny);
    // visualize(phi, nx, ny);
    // sleep(1);
    //cout << "n = " << n << "\n";
    // if( n%100 == 0) {
    //fileName = "phi_" + to_string(n) + ".vtk";
    //paraview(fileName, phi, nx, ny, dx, dy);
    // }
    // }
  //fileName = "var.vtk";
  //paraview(fileName, phi, nx, ny, dx, dy);

}



