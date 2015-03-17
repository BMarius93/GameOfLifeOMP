#include <iostream>
#include <fstream>
#include <stdlib.h> 
#include <omp.h>

using namespace std;
int main(int argc, char *argv[]){
	int nrThreads = atoi(argv[1]);
	int N = atoi(argv[2]);
	char *input = argv[3];
	char *output = argv[4];
	cout << input << endl;
	cout << output << endl;
	
	ofstream myfile;
	myfile.open (output);
	ifstream ifs;
	ifs.open(input);
	int wHarta, hHarta, w, h;
	char tip;
	int i,j,k;
	ifs >> tip >> wHarta >> hHarta >> w >> h;
	cout << tip << wHarta << hHarta << w << h << endl;
	int harta[hHarta][wHarta];
	int simulare[h][w];
	omp_set_num_threads(nrThreads);
	for(i = 0; i < hHarta; ++i){
		for(j = 0;j < wHarta; ++j){
			ifs >> harta[i][j];
		}
	}
	#pragma omp parallel for private(i,j) 
  	for(i = 0; i < h; ++i){
		for(j = 0;j < w; ++j){
			simulare[i][j] = 0;
		}
	}

	#pragma omp parallel for private(i,j)  
  	for(i = 0; i < hHarta; ++i){
		for(j = 0;j < wHarta; ++j){
			simulare[i][j] = harta[i][j];
		}
	}
	
	if(tip == 'P'){
		int copie[h + 2][w + 2];
		int copie2[h + 2][w + 2];
		#pragma omp parallel for private(i) 
		for(i = 0; i < w + 2; ++i){
			copie [0][i] = 0;
			copie [h+1][i] = 0;
			copie2 [0][i] = 0;
			copie2 [h+1][i] = 0;
		}
		#pragma omp parallel for private(i) 
		for(i = 0; i < h + 2; ++i){
			copie [i][0] = 0;
			copie [i][w+1] = 0;
			copie2 [i][0] = 0;
			copie2 [i][w+1] = 0;
		}
		//#pragma omp parallel for private(i,j) 
		for(i = 0; i < h; ++i){
			for(j = 0;j < w; ++j){
				copie [i+1][j+1] = simulare[i][j];
				copie2 [i+1][j+1] = simulare[i][j];
			}
		}
		
		

		for(k = 0; k < N; k++){
			//cout << k <<endl;
			//PARALEL
			#pragma omp parallel for private(i,j)
			for(i = 0; i < h; ++i){
				for(j = 0;j < w; ++j){
					int ok = 0;
					if(copie[i][j] == 1)
						++ok;
					if(copie[i][j+1] == 1)
						++ok;
					if(copie[i][j+2] == 1)
						++ok;
					if(copie[i+1][j] == 1)
						++ok;
					if(copie[i+1][j+2] == 1)
						++ok;
					if(copie[i+2][j] == 1)
						++ok;
					if(copie[i+2][j+1] == 1)
						++ok;
					if(copie[i+2][j+2] == 1)
						++ok;


					if(ok == 3)
						copie2[i+1][j+1] = 1;
					else if(ok == 2)
						copie2[i+1][j+1] = copie[i+1][j+1];
					else
						copie2[i+1][j+1] = 0;

					


				}
			}
			//#pragma omp parallel for private(i,j) 
			for(i = 0; i < h + 2; ++i){
				for(j = 0;j < w + 2; ++j){
					copie[i][j] = copie2[i][j];
				}
			}

		}

		int contorH = 0;
		for(i = h -1;i >= 0;i--){
			
			int ok = 0;
			for(j = 0; j < w; ++j){
				if(copie[i+1][j+1] == 1)
					ok = 1;
			}
			if(ok == 0){
				++contorH;
			}else
				break;
		}
		cout << contorH << endl;

		int contorW = 0;
		for(j = w -1;j >= 0;j--){
			
			int ok = 0;
			for(i = 0; i < h; ++i){
				if(copie[i+1][j+1] == 1)
					ok = 1;
			}
			if(ok == 0){
				++contorW;
			}else
				break;
		}
		cout << contorW << endl;
		myfile << tip << " " << w - contorW << " " <<h - contorH << " "<<w << " "<<h<<endl;
		for(i = 0; i < h - contorH; ++i){
			for(j = 0;j < w - contorW; ++j){
				myfile << copie[i+1][j+1] << " ";
			}
			myfile << endl;
		}

	}
	if(tip == 'T'){
		cout << "TOROID\n";
		int copie[h + 2][w + 2];
		int copie2[h + 2][w + 2];
		#pragma omp parallel for private(i) 
		for(i = 0; i < w; ++i){
			//cout << "for1 " <<i<<endl;

			copie [0][i+1] = simulare[h-1][i];
			copie [h+1][i+1] = simulare[0][i];
			copie2 [0][i+1] = simulare[h-1][i];
			copie2 [h+1][i+1] = simulare[0][i];
		}
		#pragma omp parallel for private(i) 
		for(i = 0; i < h; ++i){
			//cout << "for2 "<<i<<endl;

			copie [i+1][0] = simulare[i][w-1];
			copie [i+1][w+1] = simulare[i][0];
			copie2 [i+1][0] = simulare[i][w-1];
			copie2 [i+1][w+1] = simulare[i][0];
		}
		cout << "TOROID1\n";

		copie[0][0] = simulare[h-1][w-1];
		copie[0][w+1] = simulare[h-1][0];
		copie[h+1][0] = simulare[0][w-1];
		copie[h+1][w+1] = simulare[0][0];

		copie2[0][0] = simulare[h-1][w-1];
		copie2[0][w+1] = simulare[h-1][0];
		copie2[h+1][0] = simulare[0][w-1];
		copie2[h+1][w+1] = simulare[0][0];
		cout << "TOROID2\n";


		#pragma omp parallel for private(i,j) 
		for(i = 0; i < h; ++i){
			for(j = 0;j < w; ++j){
				copie [i+1][j+1] = simulare[i][j];
				copie2 [i+1][j+1] = simulare[i][j];
			}
		}


		for(k = 0; k < N; k++){
			//cout << k <<endl;
			//paralel
			#pragma omp parallel for private(i,j)

			for(i = 0; i < h; ++i){
				for(j = 0;j < w; ++j){
					int ok = 0;
					if(copie[i][j] == 1)
						++ok;
					if(copie[i][j+1] == 1)
						++ok;
					if(copie[i][j+2] == 1)
						++ok;
					if(copie[i+1][j] == 1)
						++ok;
					if(copie[i+1][j+2] == 1)
						++ok;
					if(copie[i+2][j] == 1)
						++ok;
					if(copie[i+2][j+1] == 1)
						++ok;
					if(copie[i+2][j+2] == 1)
						++ok;


					if(ok == 3)
						copie2[i+1][j+1] = 1;
					else if(ok == 2)
						copie2[i+1][j+1] = copie[i+1][j+1];
					else
						copie2[i+1][j+1] = 0;

				}
			}
			#pragma omp parallel for private(i,j) 
			for(i = 0; i < h + 2; ++i){
				for(j = 0;j < w + 2; ++j){
					copie[i][j] = copie2[i][j];
				}
			}
			#pragma omp parallel for private(i) 
			for(i = 0; i < w; ++i){
				copie [0][i+1] = copie2[h][i+1];
				copie [h+1][i+1] = copie2[1][i+1];
				
			}
			#pragma omp parallel for private(i) 
			for(i = 0; i < h; ++i){
				copie [i+1][0] = copie2[i+1][w];
				copie [i+1][w+1] = copie2[i+1][1];
				
			}

			copie[0][0] = copie2[h][w];
			copie[0][w+1] = copie2[h][1];
			copie[h+1][0] = copie2[1][w];
			copie[h+1][w+1] = copie2[1][1];

		}

		
		int contorH = 0;
		for(i = h -1;i >= 0;i--){
			
			int ok = 0;
			for(j = 0; j < w; ++j){
				if(copie[i+1][j+1] == 1)
					ok = 1;
			}
			if(ok == 0){
				++contorH;
			}else
				break;
		}
		cout << contorH << endl;

		int contorW = 0;
		for(j = w -1;j >= 0;j--){
			
			int ok = 0;
			for(i = 0; i < h; ++i){
				if(copie[i+1][j+1] == 1)
					ok = 1;
			}
			if(ok == 0){
				++contorW;
			}else
				break;
		}
		cout << contorW << endl;
		myfile << tip << " " << w - contorW << " " <<h - contorH << " "<<w << " "<<h<<endl;
		for(i = 0; i < h - contorH; ++i){
			for(j = 0;j < w - contorW; ++j){
				myfile << copie[i+1][j+1] << " ";
			}
			myfile << endl;
		}

	}
	myfile.close();
	return 0;
}
