#include <iostream> // cout
#include <cmath> // use pow&sqrt
#include <fstream> // file I/O
using namespace std; // std = standard

// global variable
int i, j;
int i_c, j_c;
const int nx = 31;
const int ny = 31;
double phi[nx][ny];

// declare variable named "myfileI", "myfileO"
ifstream myfileI; // input file stream
ofstream myfileO; // output file stream

void initialize(){
  for(i = 0; i <= nx-1; i++){
    for(j = 0; j <= ny-1; j++){
      phi[i][j] = 0.0;
    }
  }
}

void visualize(){
  for(i = 0; i <= nx-1; i++){
    for(j = 0; j <= ny-1; j++){
      cout << phi[i][j] << " ";
    }
    cout << "\n";
  }
}
void set_phi(){
  // Assign Phi = 1 Inside A Circle Around Center Of Radiaus = 15
  // letting dx =dy =1
  // as, e.g., x_spacing =dx*diff(index), not diff(index) itself
  // phi[1][3] - phi[2][3] : spacing = dx*(2-1) = dx*(1) = 1

  double radius = 10.;
  i_c = (nx-1)/2;
  j_c = (ny-1)/2;
  
  for(i = 0; i <= nx-1; i++){
    for(j = 0; j <= ny-1; j++){

      if( sqrt( pow(i-i_c,2)+ pow(j-j_c,2) ) < radius ){
	phi[i][j] =1.;}
      else{phi[i][j] = 0.;}
    }
  }
}

void save_restartfile(){
  myfileO.open("file.dat");
  for(i = 0; i <= nx-1; i++){
    for(j = 0; j <= ny-1; j++){
      myfileO << phi[i][j] << " ";
    }
    myfileO << "\n";
  }
  myfileO.close();
}



void read_restartfile(){
  myfileI.open("file.dat");
  for(i = 0; i <= nx-1; i++){
    for(j = 0; j <= ny-1; j++){
      myfileI >> phi[i][j];
    }
  }
  myfileO.close();
}




int main(){

  initialize();
  // set_phi();
  // visualize();
  // save_restartfile();
  read_restartfile();
  visualize();
  
 } //end main


