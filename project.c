#include <stdio.h>
#include <xmmintrin.h>
#include <time.h>

#define F1 3
#define F2 3

#define IN  10006

#define NUM_OF_TILES (IN-2)*(IN-2)/16
#define TD 6

float u[F1*2][F1*2];
float v[NUM_OF_TILES][TD][TD];
float m[NUM_OF_TILES][TD][TD];
float y[NUM_OF_TILES][TD][TD];

float filter[F1][F2];
float inp[IN][IN];
float res[IN-2][IN-2];
clock_t begin;
clock_t end;

void conv(float inp[IN][IN],float filter[F1][F2],float res[IN-2][IN-2]){
  float sum;
  for (int i = 0; i <IN-F1+1  ; i++) 
      for (int j = 0; j < IN-F2+1; j++) {
      sum=0;
        for (int k = 0; k < F1; k++)
            for(int q=0;q<F2;q++)
                sum = sum + inp[i+k][j+q]*filter[k][q];
        res[i][j] = sum; 
      }
    
    
  
}


void calc_u33(float filter[F1][F2],float res[6][6]){
    int m=6;
    float r1[6][F1];
    float sum;
    /*float g[6][3]={ {1.0/4,0,0},
            { -1.0/6,-1.0/6, -1.0/6  },
            { -1.0/6, 1.0/6, -1.0/6}, 
            {1.0/24,1.0/12,1.0/6},
            {1.0/24,-1.0/12,1.0/6},
            {0,0,1}
            };
    */
     
    for (int i = 0; i < F1 ; i++){
        r1[0][i]=filter[0][i]/4;
        r1[1][i]=(filter[0][i]+filter[1][i]+filter[2][i])/(-6);
        r1[2][i]=(-filter[0][i]+filter[1][i]-filter[2][i])/(6);
        r1[3][i]=filter[0][i]/24+filter[1][i]/12+filter[2][i]/6;
        r1[4][i]= filter[0][i]/24-filter[1][i]/12+filter[2][i]/6;
        r1[5][i]=filter[2][i];
    }
    for (int i = 0; i < m ; i++){
        res[i][0]=r1[i][0]/4;
        res[i][1]=(r1[i][0]+r1[i][1]+r1[i][2])/(-6);
        res[i][2]=(-r1[i][0]+r1[i][1]-r1[i][2])/(6);
        res[i][3]=r1[i][0]/24+r1[i][1]/12+r1[i][2]/6;
        res[i][4]=r1[i][0]/24-r1[i][1]/12+r1[i][2]/6;
        res[i][5]=r1[i][2];
    }
    
    /*
      
      for (int i = 0; i < m ; i++) 
        for (int j = 0; j < m; j++) {
        sum=0;
        for(int k = 0; k < F1; k++)
          sum += r1[i][k]*g[j][k];
        res[i][j] = sum;
        } 
        
        */ 

   
}

void calc_v33(float inp[IN][IN],float v[NUM_OF_TILES][TD][TD]){
    float bt[TD][TD]={{4,0,-5,0,1,0},
                  {0,-4,-4,1,1,0},
                  {0,4,-4,-1,1,0},
                  {0,-2,-1,2,1,0},
                  {0,2,-1,-2,1,0},
                  {0,4,0,-5,0,1}};
    float r1[TD][TD];
    float sum;
    int index=0;
    for(int i=0;i<IN-TD+1;i+=4)
        for(int j=0;j<IN-TD+1;j+=4){
            for(int k=0;k<TD;k++){
                r1[0][k]=4*inp[i][j+k]-5*inp[i+2][j+k]+inp[i+4][j+k];
                r1[1][k]=-4*(inp[i+1][j+k]+inp[i+2][j+k])+inp[i+3][j+k]+inp[i+4][j+k];
                r1[2][k]=4*(inp[i+1][j+k]-inp[i+2][j+k])-inp[i+3][j+k]+inp[i+4][j+k];
                r1[3][k]=2*(-inp[i+1][j+k]+inp[i+3][j+k])-inp[i+2][j+k]+inp[i+4][j+k];
                r1[4][k]=2*(inp[i+1][j+k]-inp[i+3][j+k])-inp[i+2][j+k]+inp[i+4][j+k];
                r1[5][k]=4*inp[i+1][j+k]-5*inp[i+3][j+k]+inp[i+5][j+k];
            }
    
            for(int k=0;k<TD;k++){
                v[index][k][0]=4*r1[k][0]-5*r1[k][2]+r1[k][4];
                v[index][k][1]=-4*(r1[k][1]+r1[k][2])+r1[k][3]+r1[k][4];
                v[index][k][2]=4*(r1[k][1]-r1[k][2])-r1[k][3]+r1[k][4];
                v[index][k][3]=2*(-r1[k][1]+r1[k][3])-r1[k][2]+r1[k][4];
                v[index][k][4]=2*(r1[k][1]-r1[k][3])-r1[k][2]+r1[k][4];
                v[index][k][5]=4*r1[k][1]-5*r1[k][3]+r1[k][5];
        }
       
            index++;
        }
    
    
}

    
void     calc_Elem_wise(float u[TD][TD], float v[NUM_OF_TILES][TD][TD],float m[NUM_OF_TILES][TD][TD]){
    for(int i=0;i<NUM_OF_TILES;i++)
        for(int j=0;j<TD;j++)
             for(int k=0;k<TD;k++)
                 m[i][j][k]=u[j][k]*v[i][j][k];
    
}
  
void    calc_y(float m[NUM_OF_TILES][TD][TD],float y[NUM_OF_TILES][TD][TD]){
    float at[TD-2][TD]={{1,1,1,1,1,0},
                       {0,1,-1,2,-2,0},
                       {0,1,1,4,4,0},
                       {0,1,-1,8,-8,1}};
    float t[TD][TD];
    for(int i=0;i<NUM_OF_TILES;i++){
        for(int j=0;j<TD;j++){
            t[0][j]=m[i][0][j]+m[i][1][j]+m[i][2][j]+m[i][3][j]+m[i][4][j];
            t[1][j]=m[i][1][j]-m[i][2][j]+2*(m[i][3][j]-m[i][4][j]);
            t[2][j]=m[i][1][j]+m[i][2][j]+4*(m[i][3][j]+m[i][4][j]);
            t[3][j]=m[i][1][j]-m[i][2][j]+8*(m[i][3][j]-m[i][4][j])+m[i][5][j];

        }
        for(int j=0;j<TD;j++){
            y[i][j][0]=t[j][0]+t[j][1]+t[j][2]+t[j][3]+t[j][4];
            y[i][j][1]=t[j][1]-t[j][2]+2*(t[j][3]-t[j][4]);
            y[i][j][2]=t[j][1]+t[j][2]+4*(t[j][3]+t[j][4]);
            y[i][j][3]=t[j][1]-t[j][2]+8*(t[j][3]-t[j][4])+t[j][5];
        }
        
    }
}
    

void conv33(float inp[IN][IN],float filter[F1][F2],float res[IN-2][IN-2]){
  

    calc_u33(filter,u);
    calc_v33(inp,v);
    calc_Elem_wise(u,v,m);
        begin = clock();

    calc_y(m,y);
    /*for (int row=0; row<4; row++)
    {
    for(int columns=0; columns<4; columns++)
         printf("%f     ", y[1][row][columns]);
    printf("\n");
    }
    */
    
}

int main(int argc, char *argv[]){

    for(int i=0;i<F1;i++){
        for(int j=0;j<F2;j++){
            filter[i][j]=1+j;
        }
    }
     for(int i=0;i<IN;i++){
        for(int j=0;j<IN;j++){
            inp[i][j]=i*IN+j;
        }
    }
    
    begin = clock();
    conv(inp,filter,res);

    end = clock();
    double time_spent = (double)(end - begin)/  CLOCKS_PER_SEC;
    printf("%lf\n",time_spent);
     
     conv33(inp,filter,res);

     end = clock();
     time_spent = (double)(end - begin)/  CLOCKS_PER_SEC;
    printf("%lf",time_spent);

    
    return 0;
}