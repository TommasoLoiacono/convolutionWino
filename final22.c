#include <stdio.h>
#include <xmmintrin.h>
#include <time.h>

#define F1 2
#define F2 2
#define NOI 1
#define CHN 3
#define NOK 1
#define IN 7


#define NUM_OF_TILES (IN-1)*(IN-1)/9
#define TD 4


float u[NOK][CHN][F1*2][F1*2]={0};
float v[NOI][CHN][NUM_OF_TILES][TD][TD]={0};
float m[NOI][NOK][NUM_OF_TILES][TD][TD]={0};
float y[NOI][NOK][NUM_OF_TILES][TD-1][TD-1]={0};

float filter[NOK][CHN][F1][F2]={0};
float inp[NOI][CHN][IN][IN]={0};
float res[NOI][NOK][IN-1][IN-1]={0};
float out[NOI][NOK][IN-1][IN-1]={0};

clock_t begin;
clock_t end;

float r1[NOK][TD][F1];
float r2[NOI][TD][TD];

//Normal convolution algorithm
void conv(float inp[NOI][CHN][IN][IN],float filter[NOK][CHN][F1][F2],float res[NOI][NOK][IN-F1+1][IN-F2+1], int numImg){
	
	float sum;
	
	for (int nk=0;nk<NOK;nk++){
		for(int ch=0;ch<CHN;ch++){
			for (int j=0;j<IN-F1+1;j++){ 
				for (int k=0;k<IN-F2+1;k++) {
					
					sum=0;

					for (int m=0;m<F1;m++){
						for(int n=0;n<F2;n++){
						
							sum=sum+inp[numImg][ch][j+m][k+n]*filter[nk][ch][m][n];
							
						}
					}	
					
					res[numImg][nk][j][k]=res[numImg][nk][j][k]+sum;
				}
			}
		}
	}
	
}

//Calculate u. Non riusciamo a tirare fuori G e G trasposta, che qui dentro fanno schifo?
void calc_u33(float filter[NOK][CHN][F1][F2],float res2[NOK][CHN][F1*2][F1*2]){

	int m=F1*2;
    float sum;


	for (int nk=0;nk<NOK;nk++){
		for(int ch=0;ch<CHN;ch++){
			for (int j=0;j<F1;j++){
                r1[nk][0][j]=filter[nk][ch][j][0];
                r1[nk][1][j]=(filter[nk][ch][j][0]+filter[nk][ch][j][1])/2;
                r1[nk][2][j]=(filter[nk][ch][j][0]-filter[nk][ch][j][1])/2;
                r1[nk][3][j]=filter[nk][ch][j][1];


			}

			for (int j=0;j<m;j++){
                res2[nk][ch][j][0]=r1[nk][j][0];
                res2[nk][ch][j][1]=(r1[nk][j][0]+r1[nk][j][1])/2;
                res2[nk][ch][j][2]=(r1[nk][j][0]-r1[nk][j][1])/2;
                res2[nk][ch][j][3]=r1[nk][j][1];
			}
		}
	}
     

}

//calculate v. Idem per v.
void calc_v33(float inp[NOI][CHN][IN][IN],float v[NOI][CHN][NUM_OF_TILES][TD][TD]){

    float sum;
    int index=0;
    
	for(int ni=0;ni<NOI;ni++){
		for(int ch=0;ch<CHN;ch++){
			for(int i=0;i<IN-1;i+=3){
				for(int j=0;j<IN-1;j+=3){
					for(int k=0;k<TD;k++){
                       r2[ni][0][k]=inp[ni][ch][i][j+k]-inp[ni][ch][i+2][j+k]; 
                        r2[ni][1][k]=inp[ni][ch][i+1][j+k]+inp[ni][ch][i+2][j+k];
                        r2[ni][2][k]=-inp[ni][ch][i+1][j+k]+inp[ni][ch][i+2][j+k]; 
                        r2[ni][3][k]=-inp[ni][ch][i+1][j+k]+inp[ni][ch][i+3][j+k]; 
					}
			
					for(int k=0;k<TD;k++){
                        v[ni][ch][index][k][0]=r2[ni][k][0]-r2[ni][k][2];
                        v[ni][ch][index][k][1]=r2[ni][k][1]+r2[ni][k][2];
                        v[ni][ch][index][k][2]=-r2[ni][k][1]+r2[ni][k][2];
                        v[ni][ch][index][k][3]=-r2[ni][k][1]+r2[ni][k][3];
						
					}
			   
					index++;
				}
			}
		
		index=0;
		
		}
    }
    
    
}

//Multiplication of u and v element wise    
void calc_Elem_wise(float u[NOK][CHN][TD][TD], float v[NOI][CHN][NUM_OF_TILES][TD][TD],float m[NOI][NOK][NUM_OF_TILES][TD][TD], int numImg){
	
    for (int nk=0;nk<NOK;nk++){
		for(int ch=0;ch<CHN;ch++){
			for(int i=0;i<NUM_OF_TILES;i++){
				for(int j=0;j<TD;j++){
					 for(int k=0;k<TD;k++){
						 
						 m[numImg][nk][i][j][k]
                             =m[numImg][nk][i][j][k]+u[nk][ch][j][k]*v[numImg][ch][i][j][k];
					
					}
				}
			}
		}
	}
}

//Calculate final result multiplying with a  
void    calc_y(float m[NOI][NOK][NUM_OF_TILES][TD][TD],float y[NOI][NOK][NUM_OF_TILES][TD-1][TD-1],int numImg){
/*
    float t[TD-2][TD];
	for(int nk=0; nk<NOK; nk++){
		for(int i=0;i<NUM_OF_TILES;i++){
			for(int j=0;j<TD;j++){
			   
				t[0][j]=m[numImg][nk][i][0][j]+m[numImg][nk][i][1][j]+m[numImg][nk][i][2][j]+m[numImg][nk][i][3][j]+m[numImg][nk][i][4][j];
				t[1][j]=m[numImg][nk][i][1][j]-m[numImg][nk][i][2][j]+2*(m[numImg][nk][i][3][j]-m[numImg][nk][i][4][j]);
				t[2][j]=m[numImg][nk][i][1][j]+m[numImg][nk][i][2][j]+4*(m[numImg][nk][i][3][j]+m[numImg][nk][i][4][j]);
				t[3][j]=m[numImg][nk][i][1][j]-m[numImg][nk][i][2][j]+8*(m[numImg][nk][i][3][j]-m[numImg][nk][i][4][j])+m[numImg][nk][i][5][j];

			}
		  
			for(int j=0;j<TD;j++){
				
				y[numImg][nk][i][j][0]=t[j][0]+t[j][1]+t[j][2]+t[j][3]+t[j][4];
				y[numImg][nk][i][j][1]=t[j][1]-t[j][2]+2*(t[j][3]-t[j][4]);
				y[numImg][nk][i][j][2]=t[j][1]+t[j][2]+4*(t[j][3]+t[j][4]);
				y[numImg][nk][i][j][3]=t[j][1]-t[j][2]+8*(t[j][3]-t[j][4])+t[j][5];
			
			}
			
		}
	}
    */
    
    float temp1,temp2,temp3,temp4;
    float t[TD-1][TD];
	for(int nk=0; nk<NOK; nk++){
		for(int i=0;i<NUM_OF_TILES;i++){
			for(int j=0;j<TD;j++){
                t[0][j]=m[numImg][nk][i][0][j]+m[numImg][nk][i][1][j]+m[numImg][nk][i][2][j];
                t[1][j]=m[numImg][nk][i][1][j]-m[numImg][nk][i][2][j];
                t[2][j]=m[numImg][nk][i][1][j]+m[numImg][nk][i][2][j]+m[numImg][nk][i][3][j];
			}		  
			for(int j=0;j<TD;j++){
                y[numImg][nk][i][j][0]=t[j][0]+t[j][1]+t[j][2];
                y[numImg][nk][i][j][1]=t[j][1]-t[j][2];
                y[numImg][nk][i][j][2]=t[j][1]+t[j][2]+t[j][3];
			
			}
			
		}
	}
    
}


void buildRes(float tiles[NOI][NOK][NUM_OF_TILES][TD-1][TD-1],int img,float out[NOI][NOK][IN-1][IN-1]){
    int ibase;
    int jbase;
    for(int nk=0; nk<NOK; nk++)
            for(int i=0;i<NUM_OF_TILES;i++){
                for(int j=0;j<TD-1;j++)
                    for(int k=0;k<TD-1;k++){
                        ibase=((i*3)/(IN-1))*3;
                        jbase=((i*3)%(IN-1));
                        out[img][nk][ibase+j][jbase+k]=tiles[img][nk][i][j][k];
                    }
		  }
    
                
            
    
    
}
//Winograd convolution
void conv33(float inp[NOI][CHN][IN][IN],float filter[NOK][CHN][F1][F2]){
	
    calc_u33(filter,u);
  
	calc_v33(inp,v);

	
	for (int i=0; i<NOI; i++){	
		calc_Elem_wise(u,v,m,i);
		calc_y(m,y,i);
        buildRes(y,i,out);
	}
    
}

int main(int argc, char *argv[]){

	//Kernel Creation
    for(int nk=0;nk<NOK;nk++){
		for(int ch=0;ch<CHN;ch++){
			for(int j=0;j<F1;j++){
				for (int k=0;k<F2;k++){
					
					filter[nk][ch][j][k]=1;//+ch+j+k+ch*nk;
				
				}
			}
		}
    }
	
	//Input Creation
	for (int ni=0;ni<NOI;ni++){
		for(int ch=0;ch<CHN;ch++){
			for(int j=0;j<IN;j++){
				for(int k=0;k<IN;k++){
				
					inp[ni][ch][j][k]=1;//j*IN+k*NOI+ni*ch;
				
				}
			}
		}
    }
	
    begin = clock();
	
	for (int i=0; i<NOI; i++){
		conv(inp,filter,res, i);
	}
	
    end = clock();
    double time_spent = (double)(end - begin)/  CLOCKS_PER_SEC;
    printf("tempo convoluzione normale = %lf\n",time_spent);

	begin = clock();
    conv33(inp,filter);
    end = clock();
	time_spent = (double)(end - begin)/  CLOCKS_PER_SEC;
    printf("tempo convoluzione wino = %lf\n",time_spent);
    
    
   for (int ni=0;ni<NOI;ni++){
		for(int ch=0;ch<NOK;ch++){
			for(int j=0;j<IN-1;j++){
				for(int k=0;k<IN-1;k++){
                    printf("%f ",out[ni][ch][j][k]);
				}
                printf("\n");
			}
            
            printf("\n\n\n");
		}
    }
    
     for (int ni=0;ni<NOI;ni++){
		for(int ch=0;ch<NOK;ch++){
			for(int j=0;j<IN-1;j++){
				for(int k=0;k<IN-1;k++){
                    printf("%f ",res[ni][ch][j][k]);	
                }
                printf("\n");
			}
            printf("\n\n\n");
		}
    }

    
    
    
     
        
                                  


    getchar();
	
    return 0;
}