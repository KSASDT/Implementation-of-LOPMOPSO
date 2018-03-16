#include<iostream>
#include<time.h>
#include<stdlib.h>
#include<math.h>
#include<fstream>
#include <iomanip>
#include<algorithm>
#include <random>
using namespace std;
const int tp =200;    // population size
const int nVar=3;   // problem dimesion
const int nRep=100;      //archive size
const int p = 50;    //leader particle set size
const int np=nVar+2;
double l,U,vx[tp][np],x[nVar];
double pop[tp][np],pbest[tp][np],rep[tp][np];
const int m=2;      // two objective functions
const int mu=0.1;             //Mutation Rate
const int MaxIt=60;

double z[2];
double *fun(double x[])    //functions
{
    double f1=0,f2=1;
    for(int y=0;y<nVar;y++)
	{
		f1=f1+x[y];
		f2=f2*x[y];
	}
    
	z[0]=f1;
	z[1]=f2;
  
    return z;
}


bool dominates(int x, int y)
{
  if(((pop[x][np-2])<=pop[y][np-2])&&(pop[x][np-1]<pop[x][np-1]))
  return true;
  if(((pop[x][np-2])<pop[y][np-2])&&(pop[x][np-1]<=pop[x][np-1]))
  return true;
  else
  return false;
}

int main()
{
   int i,j; 
   double *temp;
   double w, r1, r2, c1 = 1.0, c2 = 2.0;
   cout<<"enter the lower and upper bounds    ";
      cin>>l>>U;
	
	//*********Initial Popualtion generation****
	
	for(i=0;i<tp;i++)                            // initial population
                 {
                   for(j=0;j<nVar;j++)
                     {
                         x[j]=l+(U-l)*((float)rand()/RAND_MAX);
                         pop[i][j]=x[j];
                         vx[i][j] = 0.0;              //l+(U-l)*((float)rand()/RAND_MAX);
                         pbest[i][j] = pop[i][j];
                         
                     }
					 temp=fun(x); //get the 2 funct values
                     pbest[i][j]= pop[i][j] =temp[0] ;    // 1st optimization function value stored for ith population
                     pbest[i][j+1]= pop[i][j+1]=temp[1];  //value of second opt function stored   
                }
	// to set the dominated and non-dominated solutions			
	 bool pop_dominated[tp];
	 for(i=0;i<tp;i++)
	 pop_dominated[i]=false;
	  for(i=0;i<tp-1;i++)
	  {
	    for(j=i+1;j<tp;j++)
		{
		  if( dominates(i,j))
		   pop_dominated[j]=true;
		  if (dominates(j,i))
		  pop_dominated[i]=true;
		}
	  }
	  
	  // --- selecting-- for REP ----
	  int size=0;
	  
	  for(i=0;i<tp;i++)
	  {
	      if(!pop_dominated[i])
		  
    	   for(j=0;j<np;j++)
	       {
		     rep[size++][j]=pop[i][j];     // storing non dominates solutions in rep
	       }
		  
	  }
	  
	  
	  // -------Iterations------------
	 int t=0;
	 for(t=1;t<=MaxIt;t++)
	 {
	 	 double dist1[size]; 
	    double rep2[tp][np];
	     int nRep2; // for size of rep actual one may be less then nRep
   //--* * crowding_distance
		{
		   //crowding distance for each particle
		  double D[size][m];  //distance for particle wrt optm fun
		  int so_ind[size][m]; // to maintain sorted list of indices for both obj functions
		  int s=size;
		  for(int i=0;i<size;i++)
		  dist1[i]=0;                // initialize crowding distance
		  for(int k=0;k<m;k++)    // sorting rep on basis of objective function	
		  {
		     // int so_ind[size];
			 double buffer[size];
			 for(int i=0;i<size;i++)
			  { 
			     so_ind[i][k]=i;                //indices initially not sorted
				buffer[i]=rep[i][np-2+k];
			  }
			   // --- selection sort ---
			 int min, tmp ;
			 double temp2;
			 for(int i=0;i<size-1;i++)
			  {
			     min=i;
			     for(int j=i+1;j<size;j++)
				 {
				
				if(buffer[min]>buffer[j])
				      min=j;
				   
				 }
				 if(i!=min)
				 {
				    temp2 = buffer[i];
					tmp = so_ind[i][k];
					buffer[i]= buffer[min];
					so_ind[i][k]=so_ind[min][k];
					buffer[min] = temp2;
					so_ind[min][k] =tmp;
				 }
				}
				
             D[so_ind[0][k]][k]= D[so_ind[size-1][k]][k]=999999.0;	 // infinity to boundary values		
		     
			 for(int i=2;i<size-1;i++)     // distance of each objective wrt other individuals
			  D[so_ind[i][k]][k]= rep[so_ind[i+1][k]][np-2+k]-rep[so_ind[i+1][k]][np-2+k];
		      
		 }
		 for(int o;o<size;o++)
		 { 
		   dist1[o]=D[o][0]+D[o][1];   // crowding distance for all particles
		  }
		 while(size>nRep)
           { 
		     double min2=dist1[0];
			 int m_ind= 0;
			 int r;
		     for( r=0;r<size;r++)  //identify individual with minimum crowding distance
			 {
			   if(dist1[r]<min2)
			   {
			     min2= dist1[r];
				 m_ind=r;
			   }
			 }
			
			 if(size>nRep)
			 nRep2=nRep;
			 else
			 nRep2=size;
			 int index_alt,b1,b2,b;
			 
			 for(int k=0;k<m;k++)
			 {
			   for(int y=0;y<size;y++) // to find the sorted i value for r
			   {
			      if(r==so_ind[y][k])
				   {if(k==0)
				     {b1=y; b=y;}
                    else
					 {b2=y; b=y;}
					break;}
			   }
			   //updating distance of particles affected by deletion
			   
			   D[so_ind[b-1][k]][k]=rep[so_ind[b+1][k]][np-2+k]-rep[so_ind[b-2][k]][np-2+k];
			   D[so_ind[b+1][k]][k]=rep[so_ind[b+2][k]][np-2+k]-rep[so_ind[b-1][k]][np-2+k];
		       
			   // effect of removal on so_ind list
			   for(int y=b;y<size-1;y++)
			   {
			     so_ind[y][k]=so_ind[y+1][k];
			   }
			  }
			  
			  dist1[so_ind[b1-1][0]]=D[so_ind[b1-1][0]][0]+D[so_ind[b1-1][0]][1];   // crowding distance updated for neighbour partical 1
			  dist1[so_ind[b1][0]]=D[so_ind[b1][0]][0]+D[so_ind[b1][0]][1];        // crowding distance updated for neighbour particle 2
			   
			   for(int u=0;u<np;u++)
			    {rep[r][u]=rep[size-1][u];   // removal of element from rep
				 }
				dist1[r]=dist1[size-1];
			    size--;  //updating new size value of rep
            }		   
		  
		}
		
		// selection of p=7 particles
		double max = dist1[0],temp2,tmp2;
		int so_ind2[nRep2];
		int min3;
		for(int i=0;i<nRep2-1;i++)
			  {
			     min3=i;
			     for(int j=i+1;j<nRep2;j++)
				 {
				
				if(dist1[min3]>dist[j])
				      min3=j;
				   
				 }
				 if(i!=min3)
				 {
				    temp2 = dist[i];
					tmp2 = so_ind2[i];
					dist1[i]= dist1[min3];
					so_ind2[i]=so_ind2[min3];
					dist1[min3] = temp2;
					so_ind2[min3] =tmp2;
				 }
				}
				 
				 // to add 7 particles into rep2
				for(int g=0;g<7;g++)
				{
				for(int j=0;j<np;j++)
				  {
           				 // rep2[g][j] = rep[so_ind2[nRep2-g]][j];
           				 
           				 rep2[g][j] = rep[so_ind2[nRep2-g]][j];
				  }
				
				}
			
			
		for(int i= 0;i<tp;i++)
		{
		   double pg[tp][nVar];
		   //random selection of pg
		   //---
		   
		   for(int i=0;i<tp;i++)
				{
				for(int j=0;j<nVar;j++)
				  {
           				  int RandIndex = rand() % 7; // 7= P
						   pg[i][j] = rep2[RandIndex][j];
						  
				  }
				
				}
		   
		   //update velocity and position
		   
		   for(int i=0;i<tp;i++){
			   double* temp4;
		   for(int j=0;j<nVar;j++)
		   {
		     vx[i][j]=vx[i][j]+(c1*r1)*(pbest[i][j]-pop[i][j])+(c2*r2)*(pg[i][j]-pop[i][j]);
	         pop[i][j]=pop[i][j]+vx[i][j];
			 x[j] = pop[i][j];		 

			 }
			 for(int k=0;k< nVar;k++){
			 temp4 = fun(x);
			 pop[i][k] = temp4[0];
			 pop[i][k+1] = temp4[1];
			 
			 }
			 
			 
			 // update pbest through domination
			 if ((pop[i][np-2] <= pbest[i][np-2] && pop[i][np-1]< pbest[i][np-1])){
				 
				 for(int h=0; h<nVar; h++){
					 
					 pbest[i][h] = pop[i][h];
				 }
			 }
			 
			 if ((pop[i][np-2] < pbest[i][np-2] && pop[i][np-1]<=pbest[i][np-1])){
				 
				 for(int h=0; h<nVar; h++){
					 
					 pbest[i][h] = pop[i][h];
				 }
			 }else{}
		   
		}
		

	 
	 }
	 
	 
	 
	 //mutation operation 
	 double pm,dx,lb,ub;
	  pm=(1-((t-1)/(MaxIt-1)))^(1/mu);
      
	  for(int i=0;i<=tp; i++){
	  
	  // Use random_device to generate a seed for Mersenne twister engine.
      random_device rd;    

      // Use Mersenne twister engine to generate pseudo-random numbers.
      mt19937 engine(rd());
      uniform_int_distribution<int> dist(1, tp);

     int gen = dist(engine);
	  dx = pm*(U-l);
	
	  for(int j=1; j<=np;j++){
	  lb = pop[gen][i]- dx;
	  if(lb < l)
	  {
		lb = l;
      }

	  ub = pop[gen][i] + dx;
	  if(ub > U)
	  {
		ub = U;
      }
	  
	 
	 uniform_int_distribution<int> dist(lb, ub);
	 int Newgen = dist(engine);
	 pop[i][j] = pop [Newgen][j];
	} 
	}
	
	//updating REP
	int trep[tp][nVar];
	for(int i=0; i<= size-1 ; i++){
	   for(int j =0; j< nVar; j++)
	{
		trep[i][j]= rep[i][j];
		
	}
	}
	
	for(int i=size; i<= size+6 ; i++){
	   for(int j =0; j< nVar; j++)
	{
		trep[i][j]= rep2[i][j];
		
	}
	}
	
	int u = size+6;
	bool dominated[u];
	 for(i=0;i<u;i++)
	 pop_dominated[i]=false;
	  for(i=0;i<u-1;i++)
	  {
	    for(j=i+1;j<u;j++)
		{
		  if( dominates(i,j))
		   pop_dominated[j]=true;
		  if (dominates(j,i))
		  pop_dominated[i]=true;
		}
	  }
	
	  int size=0;
	  
	  for(i=0;i<u;i++)
	  {
	      if(!pop_dominated[i])
		  
    	   for(j=0;j<np;j++)
	       {
		     rep[size++][j]=trep[i][j];     // storing non dominates solutions in rep
	       }
		  
	  }
	
	
	
	
	 }
}

 





