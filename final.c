#include <stdio.h>
#include <xmmintrin.h>
#include <time.h>

#define F1 3
#define F2 3
#define NOI 200
#define CHN 3
#define NOK 10
#define IN 100


#define NUM_OF_TILES (IN-2)*(IN-2)/16
#define TD 6


float u[NOK][CHN][F1*2][F1*2]={0};
float v[NOI][CHN][NUM_OF_TILES][TD][TD]={0};
float m[NOI][NOK][NUM_OF_TILES][TD][TD]={0};
float y[NOI][NOK][NUM_OF_TILES][TD-2][TD-2]={0};

float filter[NOK][CHN][F1][F2]={0};
float inp[NOI][CHN][IN][IN]={0};
float res[NOI][NOK][IN-2][IN-2]={0};
float out[NOI][NOK][IN-2][IN-2];

clock_t begin;
clock_t end;

float r1[NOK][6][F1];
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

				r1[nk][0][j]=filter[nk][ch][0][j]/4;
				r1[nk][1][j]=(filter[nk][ch][0][j]+filter[nk][ch][1][j]+filter[nk][ch][2][j])/(-6);
				r1[nk][2][j]=(-filter[nk][ch][0][j]+filter[nk][ch][1][j]-filter[nk][ch][2][j])/(6);
				r1[nk][3][j]=filter[nk][ch][0][j]/24+filter[nk][ch][1][j]/12+filter[nk][ch][2][j]/6;
				r1[nk][4][j]= filter[nk][ch][0][j]/24-filter[nk][ch][1][j]/12+filter[nk][ch][2][j]/6;
				r1[nk][5][j]=filter[nk][ch][2][j];

			}

			for (int j=0;j<m;j++){

				res2[nk][ch][j][0]=r1[nk][j][0]/4;
				res2[nk][ch][j][1]=(r1[nk][j][0]+r1[nk][j][1]+r1[nk][j][2])/(-6);
				res2[nk][ch][j][2]=(-r1[nk][j][0]+r1[nk][j][1]-r1[nk][j][2])/(6);
				res2[nk][ch][j][3]=r1[nk][j][0]/24+r1[nk][j][1]/12+r1[nk][j][2]/6;
				res2[nk][ch][j][4]=r1[nk][j][0]/24-r1[nk][j][1]/12+r1[nk][j][2]/6;
				res2[nk][ch][j][5]=r1[nk][j][2];
			
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
			for(int i=0;i<IN-TD+1;i+=4){
				for(int j=0;j<IN-TD+1;j+=4){
					for(int k=0;k<TD;k++){
						
						r2[ni][0][k]=4*inp[ni][ch][i][j+k]-5*inp[ni][ch][i+2][j+k]+inp[ni][ch][i+4][j+k];
						r2[ni][1][k]=-4*(inp[ni][ch][i+1][j+k]+inp[ni][ch][i+2][j+k])+inp[ni][ch][i+3][j+k]+inp[ni][ch][i+4][j+k];
						r2[ni][2][k]=4*(inp[ni][ch][i+1][j+k]-inp[ni][ch][i+2][j+k])-inp[ni][ch][i+3][j+k]+inp[ni][ch][i+4][j+k];
						r2[ni][3][k]=2*(-inp[ni][ch][i+1][j+k]+inp[ni][ch][i+3][j+k])-inp[ni][ch][i+2][j+k]+inp[ni][ch][i+4][j+k];
						r2[ni][4][k]=2*(inp[ni][ch][i+1][j+k]-inp[ni][ch][i+3][j+k])-inp[ni][ch][i+2][j+k]+inp[ni][ch][i+4][j+k];
						r2[ni][5][k]=4*inp[ni][ch][i+1][j+k]-5*inp[ni][ch][i+3][j+k]+inp[ni][ch][i+5][j+k];
					
					}
			
					for(int k=0;k<TD;k++){
					
						v[ni][ch][index][k][0]=4*r2[ni][k][0]-5*r2[ni][k][2]+r2[ni][k][4];
						v[ni][ch][index][k][1]=-4*(r2[ni][k][1]+r2[ni][k][2])+r2[ni][k][3]+r2[ni][k][4];
						v[ni][ch][index][k][2]=4*(r2[ni][k][1]-r2[ni][k][2])-r2[ni][k][3]+r2[ni][k][4];
						v[ni][ch][index][k][3]=2*(-r2[ni][k][1]+r2[ni][k][3])-r2[ni][k][2]+r2[ni][k][4];
						v[ni][ch][index][k][4]=2*(r2[ni][k][1]-r2[ni][k][3])-r2[ni][k][2]+r2[ni][k][4];
						v[ni][ch][index][k][5]=4*r2[ni][k][1]-5*r2[ni][k][3]+r2[ni][k][5];
						
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
						 
						 m[numImg][nk][i][j][k]=m[numImg][nk][i][j][k]+u[nk][ch][j][k]*v[numImg][ch][i][j][k];
					
					}
				}
			}
		}
	}
}

//Calculate final result multiplying with a  
void    calc_y(float m[NOI][NOK][NUM_OF_TILES][TD][TD],float y[NOI][NOK][NUM_OF_TILES][TD-2][TD-2],int numImg){
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
    float t[TD-2][TD];
	for(int nk=0; nk<NOK; nk++){
		for(int i=0;i<NUM_OF_TILES;i++){
			for(int j=0;j<TD;j++){
			    temp1=m[numImg][nk][i][1][j]+m[numImg][nk][i][2][j];
                temp2=m[numImg][nk][i][1][j]-m[numImg][nk][i][2][j];
                temp3=m[numImg][nk][i][3][j]+m[numImg][nk][i][4][j];
                temp4=m[numImg][nk][i][3][j]-m[numImg][nk][i][4][j];
                t[0][j]=temp1+temp3+m[numImg][nk][i][0][j];
                t[1][j]=temp2+2*temp4;
                t[2][j]=temp1+4*temp3;
                t[3][j]=temp2+8*temp4+m[numImg][nk][i][5][j];

			}
		  
			for(int j=0;j<TD;j++){
				temp1=t[j][1]+t[j][2];
                temp2=t[j][1]-t[j][2];
                temp3=t[j][3]+t[j][4];
                temp4=t[j][3]-t[j][4];
                y[numImg][nk][i][j][0]=temp1+temp3+t[j][0];
                y[numImg][nk][i][j][1]=temp2+2*temp4;
                y[numImg][nk][i][j][2]=temp1+4*temp3;
                y[numImg][nk][i][j][3]=temp2+8*temp4+t[j][5];
			
			}
			
		}
	}
    
}


void buildRes(float tiles[NOI][NOK][NUM_OF_TILES][TD-2][TD-2],int img,float out[NOI][NOK][IN-2][IN-2]){
    int ibase;
    int jbase;
    for(int nk=0; nk<NOK; nk++)
            for(int i=0;i<NUM_OF_TILES;i++){
                for(int j=0;j<TD-2;j++)
                    for(int k=0;k<TD-2;k++){
                        ibase=((i*4)/(IN-2))*4;
                        jbase=((i*4)%(IN-2));
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
					
					filter[nk][ch][j][k]=1+ch+j+k+ch*nk;
				
				}
			}
		}
    }
	
	//Input Creation
	for (int ni=0;ni<NOI;ni++){
		for(int ch=0;ch<CHN;ch++){
			for(int j=0;j<IN;j++){
				for(int k=0;k<IN;k++){
				
					inp[ni][ch][j][k]=j*IN+k*NOI+ni*ch;
				
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
    
    /*
   for (int ni=0;ni<NOI;ni++){
		for(int ch=0;ch<NOK;ch++){
			for(int j=0;j<IN-2;j++){
				for(int k=0;k<IN-2;k++){
                    printf("%f ",out[ni][ch][j][k]);
				}
                printf("\n");
			}
            
            printf("\n\n\n");
		}
    }
    
     for (int ni=0;ni<NOI;ni++){
		for(int ch=0;ch<NOK;ch++){
			for(int j=0;j<IN-2;j++){
				for(int k=0;k<IN-2;k++){
                    printf("%f ",res[ni][ch][j][k]);				
                }
                printf("\n");
			}
            printf("\n\n\n");
		}
    }
    */                               


    getchar();
	
    return 0;
}