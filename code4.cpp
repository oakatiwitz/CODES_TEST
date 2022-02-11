#include <iostream> // cout
#include <cmath> // use pow&sqrt
#include <fstream> // file I/O
using namespace std; // std = standard


//void initialize(){
// for(i = 0; i <= nx-1; i++){
//   for(j = 0; j <= ny-1; j++){
//     phi[i][j] = 0.0;
//   }
// }
//}

//void visualize(){
// for(i = 0; i <= nx-1; i++){
//   for(j = 0; j <= ny-1; j++){
//     cout << phi[i][j] << " ";
//   }
//   cout << "\n";
// }
//}

//void set_phi(){
// // Assign Phi = 1 Inside A Circle Around Center Of Radiaus = 15
// // letting dx =dy =1
// // as, e.g., x_spacing =dx*diff(index), not diff(index) itself
// // phi[1][3] - phi[2][3] : spacing = dx*(2-1) = dx*(1) = 1

// double radius = 10.;
// i_c = (nx-1)/2;
// j_c = (ny-1)/2;
  
// for(i = 0; i <= nx-1; i++){
//   for(j = 0; j <= ny-1; j++){

//     if( sqrt( pow(i-i_c,2)+ pow(j-j_c,2) ) < radius ){
//	phi[i][j] =1.;}
//     else{phi[i][j] = 0.;}
//   }
// }
//}

//void save_restartfile(){
// myfileO.open("file.dat");
// for(i = 0; i <= nx-1; i++){
//   for(j = 0; j <= ny-1; j++){
//     myfileO << phi[i][j] << " ";
//   }
//   myfileO << "\n";
// }
// myfileO.close();
//}



//void read_restartfile(){
//  myfileI.open("file.dat");
//  for(i = 0; i <= nx-1; i++){
//     for(j = 0; j <= ny-1; j++){
//       myfileI >> phi[i][j];
//     }
//   }
//   myfileI.close();
//}

//void test(int nx, double phi[nx]){
//  phi[0] = phi[1] + phi[2];
//  nx = nx +1;
//}

// void test(int nx, double *p){
//   for(int i = 0; i < nx; i++){
//     cout << "p[" << i << "] = " << p[i] << "\n";
//   }
// }



void test2D(int nx, int ny, double **p){

  for(int i = 0; i < nx; i++){
    for(int j = 0; j < ny; j++){
      cout << p[i][j];
    }
  }
  cout << "\n";
}
     
int main(){
  const int nx = 10;
  const int ny = 10;
  
  double **phi;
  phi = (double **)malloc (nx * sizeof(double));

  for(int i = 0; i < nx; i++){
    phi[i] = (double *)malloc (ny * sizeof(double));
    
  }

  test2D(nx, ny, phi);






// //   const int nx = 10;
// //   double *phi;
// //   phi = (double *)malloc (nx * sizeof(double));

// //   for(int i = 0; i < nx; i++){
// //     phi[i] = double(i);
// //   }
// //   test(nx, phi);



// //   // int a =10;
// //   // cout << "a = " << a << "\n";
// //   // cout << "&a = " << &a << "\n";
// //   // &a = 0x7ffc8a36e0c4   0x <= hexadecimal


  
// //   // 
  // 7 f f c 8 a 3 6 e 0 c 4

  // int array[10];
  // for(int i=0; i<10; i++){
  //   array[i] = i;
  //   cout << "array[" << i << "] = " << array[i] << ", stored at " << &array[i] << "\n";
  // }

  // cout << "\n";

  // int *p;
  // p = &a;
  // cout << "p = "<< p << "\n";
  // cout << "*p = " << *p << "\n";






  //  const int nx = 10;
  //  double phi[nx];
  //  for(int i=0; i < nx; i++){
  //    phi[i] = double(i);
  //    cout<< "phi[" <<i<< "] = "<< phi[i] << "\n";
  //  }
  //  cout << "\n";
  //  test(nx, phi);
  //  cout << "phi[0] = " << phi[0] << "\n";
  //  cout << "\n";
  //  cout << "nx = "<< nx << "\n";
  

  
  // declare variable named "myfileI", "myfileO"
   
  //ifstream myfileI; // input file stream
  //ofstream myfileO; // output file stream   
  
  
  // MODELAR CODE <= module, modle, module
  //PORTABLE CODE
  //initialize();
  // set_phi();
  // visualize();
  // save_restartfile();
  //read_restartfile();
  //visualize();
  
 } //end main


