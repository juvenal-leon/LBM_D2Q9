//
//  LBM_D2Q9.h
//  LBM_D2Q9
//
//  Created on 2/9/16.
//
//

#ifndef LBM_D2Q9_h
#define LBM_D2Q9_h


#endif /* LBM_D2Q9_h */


//=========================================================
//---------------------------------------------------------
//	----- Header file of the D2Q9 model -----
//---------------------------------------------------------
//File name: D2Q9.h
#include	<stdio.h>
#include	<stdlib.h>
#include	<math.h>
#include	<time.h>
#include	<assert.h>
#include    <string.h>




//#define	Nx	256	// number of cells in the x-direction
//#define	Ny	256	// number of cells in the y-direction
//#define	Nx	450	// number of cells in the x-direction
//#define	Ny	450	// number of cells in the y-direction

//#define Nx1	(Nx+1)
//#define Ny1	(Ny+1)
//#define L (Ny+1)	// width of the cavity
#define	Q 9		// number of discrete velocities
#define rho0	1.0	// initial density
#define ux0	0.0	// initial velocity component in x direction
#define uy0	0.0	// initial velocity component in y direction
#define uw	0.1
#define Re 400.0

int Nx, Ny, Nx1, Ny1, L;

//         0, 1, 2,  3,  4, 5,  6,  7,  8
int cx[Q]={0, 1, 0, -1,  0, 1, -1, -1,  1};
int cy[Q]={0, 0, 1,  0, -1, 1,  1, -1, -1};

//double f[Ny1][Nx1][Q]; //array of the distribution functions (DFs)
//double f_post[Ny1][Nx1][Q]; // array of the post-collision DFs
//double rho[Ny1][Nx1], ux[Ny1][Nx1], uy[Ny1][Nx1]; // arrays of fluid density and velocity
//int    boundary[Ny1][Nx1] = { };  // boolean array for boundary nodes
double ***f, ***f_post;
double **rho, **ux, **uy;
int    **boundary;

double tau; // relaxation time for BGK model
double s[Q]; // relaxation rates for MRT model
double D[Q]={9, 36, 36, 6, 12, 6, 12, 4, 4};	// D = M*M^T
double g=0.001;


double w[Q]={4.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/36,1.0/36, 1.0/36,1.0/36}; // the weights in the EDF
int rc[Q]={0,3,4,1,2,7,8,5,6}; // index of reversed velocity

void readBoundaryFromFile(const char *);

double ***new3DArray(int, int, int, size_t);
double **new2DArray(int, int, size_t);
int **new2DIntArray(int, int, size_t);

void Init_Conditions(int , char *[]);
void Init_Eq(void);	//Initialization
double feq(double RHO, double U, double V, int k);
// Equilibrium distribution function
void Coll_BGK(void);	// BGK collision
void Coll_MRT(void);	// MRT collision
double meq(double RHO, double U, double V, int k);
// Equilibrium	momenta
void Streaming(void);	// Streaming
void Den_Vel(void);	// Fluid variables
void Bounce_back(void);	// Bounce-back boundary condition
void Bounce_backV(void);


double Err(void);	// Difference in velocity field
//double u0[Ny1][Nx1],v0[Ny1][Nx1];
double **u0,**v0;
void Data_Output(void);	// Output simulation data

void mkePorous(double);
void mkeSolid(int, int, int, int);




//=========================================================


//-------------------------------------------------------------------

// Subroutine: initialization with the equilibrium method
//------------------------------------------------------------------
//
void Init_Eq()

{
    int j, i, k;
    
//    double (*rho)[Nx]=malloc(sizeof(double)*Ny*Nx);
//    double (*ux)[Nx]=malloc(sizeof(double)*Ny*Nx);
//    double (*uy)[Nx]=malloc(sizeof(double)*Ny*Nx);
    rho=new2DArray(Nx1, Ny1, sizeof(double));
    ux=new2DArray(Nx1, Ny1, sizeof(double));
    uy=new2DArray(Nx1, Ny1, sizeof(double));
    
    
    f =    new3DArray(Nx1, Ny1, (int)Q, sizeof(double));
    f_post=new3DArray(Nx1, Ny1, (int)Q, sizeof(double));
    
    
    for (int j=0; j<Ny; j++)
    {
        for (int i=0; i<Nx; i++)
        {
            rho[j][i]=0;
        }
    }
    
    for (j=0;j<=Ny;j++)
    {

        for(i=0;i<=Nx;i++)
        {
            rho[j][i]=rho0;
            
            ux[j][i]=ux0;
            
            uy[j][i]=uy0;
            
            for(k=0;k<Q;k++)
            {
                f[j][i][k]=feq(rho[j][i],ux[j][i],uy[j][i],k);
            }
            
        }
    }
}


//===========================================================================================
//
//Allocation of new arrays
//===========================================================================================

void*** newarray(int icount, int jcount, int kcount, int type_size)
{
    void*** iret = (void***)malloc(icount*sizeof(void***)+icount*jcount*sizeof(void**)+icount*jcount*kcount*type_size);
    void** jret = (void**)(iret+icount);
    char* kret = (char*)(jret+icount*jcount);
    for(int i=0;i<icount;i++)
        iret[i] = &jret[i*jcount];
    for(int i=0;i<icount;i++)
        for(int j=0;j<jcount;i++)
            jret[i*jcount+j] = &kret[i*jcount*kcount*type_size+j*kcount*type_size];
    return iret;
}

double ***new3DArray(int xCols, int yRows, int zLayers, size_t typeSize)
{
    double ***array;
    
    array=calloc(yRows, sizeof( double **));
    assert(array!=NULL);
    for (int j=0;j<yRows; j++)
    {
        array[j]=calloc(xCols, sizeof( double *));
        assert(array[j]!=NULL);
        for (int i=0; i<xCols; i++)
        {
            array[j][i]=calloc(zLayers, sizeof(double));
            assert(array[j][i]!=NULL);
        }
    }
    return array;
}

double **new2DArray(int xCols, int yRows, size_t typeSize)
{
    double **array;
    
    array=calloc(yRows, sizeof( double *));
    assert(array!=NULL);
    for (int j=0; j<yRows; j++)
    {
        array[j]=calloc(xCols, sizeof(double));
        assert(array[j]!=NULL);
    }
    return array;
}

int **new2DIntArray(int xCols, int yRows, size_t typeSize)
{
    int **array;
    
    array=calloc(yRows, sizeof( int *));
    assert(array!=NULL);
    for (int j=0; j<yRows; j++)
    {
        array[j]=calloc(xCols, sizeof(int));
        assert(array[j]!=NULL);
    }
    return array;
}
//===========================================================================================



//========================================================

void Init_Conditions(int num_Arg, char *arr_Arg[])
//=========================================================
{
    
    if (strcmp(arr_Arg[1],"0")==0)
    {
        Nx=450;
        Ny=450;
        Nx1=Nx+1;
        Ny1=Ny+1;
        boundary=new2DIntArray(Nx1, Ny1, sizeof(int));
        mkeSolid(0, 0, 0, Ny);
        mkeSolid(Nx,0,Nx,Ny);
        //mkeSolid(30, 30, 60, 60);
        mkePorous(25);
    }
    else
    {
        readBoundaryFromFile("/Users/juvenal/Documents/Tesis/LBM_D2Q9/LBM_D2Q9/mediosporosos/matrix.txt");
    }
    
    L=Ny1;
    
    u0=new2DArray(Nx1, Ny1, sizeof(double));
    v0=new2DArray(Nx1, Ny1, sizeof(double));

    


}



//-----------------------------------------------------------------

// Subroutine: calculation the equilibrium distribution
//----------------------------------------------------------------
//
double feq(double RHO, double U, double V, int k)
{
    
    double cu, U2;
    
    //V = V + (tau*g);
    
    cu=cx[k]*U+cy[k]*V; // c_k*u 
    
    U2=U*U+V*V; // u*u;
    
    return w[k]*RHO*(1.0+3.0*cu+4.5*cu*cu-1.5*U2);
}

//=========================================================

//=========================================================
//---------------------------------------------------------

// Subroutine: BGK collision
//---------------------------------------------------------
void Coll_BGK()
{
    //printf("Begin collision.... ");
    int j, i, k; 
    double FEQ;
    
    for (j=0;j<=Ny;j++) for(i=0;i<=Nx;i++) for(k=0;k<Q;k++)    
    {
        // the value of the EDF:
        FEQ=feq(rho[j][i],ux[j][i],uy[j][i],k); 
        // the post-collision distribution function: 
        f_post[j][i][k] = f[j][i][k]-(f[j][i][k]-FEQ)/tau;
        
    }
    //printf("End collision....\n");
}

//=========================================================

//=========================================================
//---------------------------------------------------------

// Subroutine: MRT collision
//---------------------------------------------------------
void Coll_MRT()
{
    
    int j, i, k; double MEQ; double m[Q];
    
    for (j=0;j<=Ny;j++) for(i=0;i<=Nx;i++)
        
    {
        
        //	Transformation from velocity space to moment space:
        m[0]=f[j][i][0]+f[j][i][1]+f[j][i][2]+f[j][i][3]+f[j][i][4]+f[j][i][5]+f[j][i][6]+f[j][i][7]+f[j][i][8];
        m[1]=-4*f[j][i][0]-f[j][i][1]-f[j][i][2]-f[j][i][3]-f[j][i][4]+2*(f[j][i][5]+f[j][i][6]+f[j][i][7]+f[j][i][8]);
        m[2]=4*f[j][i][0]-2*(f[j][i][1]+f[j][i][2]+f[j][i][3]+f[j][i][4])+f[j][i][5]+f[j][i][6]+f[j][i][7]+f[j][i][8];
        m[3]=f[j][i][1]-f[j][i][3]+f[j][i][5]-f[j][i][6]-f[j][i][7]+f[j][i][8];
        m[4]=-2*(f[j][i][1]-f[j][i][3])+f[j][i][5]-f[j][i][6]-f[j][i][7]+f[j][i][8];
        m[5]=f[j][i][2]-f[j][i][4]+f[j][i][5]+f[j][i][6]-f[j][i][7]-f[j][i][8];
        m[6]=-2*(f[j][i][2]-f[j][i][4])+f[j][i][5]+f[j][i][6]-f[j][i][7]-f[j][i][8];
        m[7]=f[j][i][1]-f[j][i][2]+f[j][i][3]-f[j][i][4]; 
        m[8]=f[j][i][5]-f[j][i][6]+f[j][i][7]-f[j][i][8];
        
        // Relaxation in moment space: for(k=0;k<Q;k++)
        { 
            MEQ = meq(rho[j][i],ux[j][i],uy[j][i],k);
            
            m[k]= m[k]-s[k]*(m[k]-MEQ);	// relaxation
            
            m[k]/=D[k];	// rescaling
        }
        
        // Transforming back to the velocity space: 
        f_post[j][i][0]= m[0]-4*(m[1]-m[2]);
        f_post[j][i][1]=m[0]-m[1]-2*(m[2]+m[4])+m[3]+m[7];
        f_post[j][i][2]=m[0]-m[1]-2*(m[2]+m[6])+m[5]-m[7]; 
        f_post[j][i][3]=m[0]-m[1]-2*(m[2]-m[4])-m[3]+m[7]; 
        f_post[j][i][4]=m[0]-m[1]-2*(m[2]-m[6])-m[5]-m[7]; 
        f_post[j][i][5]=m[0]+m[1]+m[1]+m[2]+m[3]+m[4]+m[5]+m[6]+m[8]; 
        f_post[j][i][6]=m[0]+m[1]+m[1]+m[2]-m[3]-m[4]+m[5]+m[6]-m[8]; 
        f_post[j][i][7]=m[0]+m[1]+m[1]+m[2]-m[3]-m[4]-m[5]-m[6]+m[8]; 
        f_post[j][i][8]=m[0]+m[1]+m[1]+m[2]+m[3]+m[4]-m[5]-m[6]-m[8];
        
    }
    
}
//=========================================================

//=========================================================
//---------------------------------------------------------

// Subroutine: calculation the equilibrium moment
//---------------------------------------------------------
double meq(double RHO, double U, double V, int k)
{
    double x; switch(k)
    {      
        case 0: {x=RHO; break;}
            
        case 1: {x=RHO*(-2+3*(U*U+V*V));break;} case 2: {x=RHO*(1-3*(U*U+V*V));break;} case 3: {x=RHO*U;break;}
            
        case 4: {x=-RHO*U;break;} case 5: {x=RHO*V;break;} case 6: {x=-RHO*V;break;}
            
        case 7: {x=RHO*(U*U-V*V);break;} case 8: {x=RHO*U*V;break;} default: x=0;  
    }
    return x;
}

//=========================================================

//=========================================================
//---------------------------------------------------------

// Subroutine: Streaming
//---------------------------------------------------------
void Streaming()
{
    //printf("Begin streaming.... ");
    
    int j, i, jd, id, k;
    //for (j=0;j<Ny;j++) for(i=0;i<Nx;i++) for(k=0;k<Q;k++)
    for (j=0;j<=Ny;j++) for(i=0;i<=Nx;i++) for(k=0;k<Q;k++)
    {
        jd=j-cy[k]; id=i-cx[k]; // upwind node

        if(jd>=0 && jd<=Ny && id>=0 && id<=Nx) // fluid node
        {
            if (!boundary[jd][id])
                f[j][i][k]=f_post[jd][id][k]; // streaming
            else
                // bounce-back on the boundary node:
                f[j][i][k]=f_post[jd][id][rc[k]]; //+ 6*w[k]*rho[j][i]*(cx[k]*uwx[jd][id]+cy[k]*uwy[jd][id]);
        }
        
    }
    
    //printf("End streaming\n");

}

//=========================================================

//=========================================================
//---------------------------------------------------------

// Subroutine: Bounce-back scheme
//---------------------------------------------------------
void Bounce_back()
{
    //printf("Bounce back.... ");
    int i,j;
    float u0=0.10;
    //	j=Ny: top plate
    /*
    for(i=0;i<=Nx;i++)
    {  
        f[Ny][i][4]=f_post[Ny][i][2]; 
        f[Ny][i][7]=f_post[Ny][i][5]+6*rho[Ny][i]*w[7]*cx[7]*uw; 
        f[Ny][i][8]=f_post[Ny][i][6]+6*rho[Ny][i]*w[8]*cx[8]*uw;
    }
     */
    for(i=0;i<=Nx;i++)
    {
        f[Ny][i][4]=f_post[Ny][i][2];
        f[Ny][i][7]=f_post[Ny][i][5];
        f[Ny][i][8]=f_post[Ny][i][6];
    }

    //	j=0: bottom plate
    
    for(i=0;i<=Nx;i++) 
    {  
        f[0][i][2]=f_post[0][i][4]; 
        f[0][i][5]=f_post[0][i][7]; 
        f[0][i][6]=f_post[0][i][8]; 
    }
    
/*----------------------Mohamad, flow in a channel ---------------------
    
    // i=0: left wall
    
    for(j=0;j<=Ny;j++)
    {
        double rhow=(f[j][0][0]+f[j][0][2]+f[j][0][4] + 2.0*(f[j][0][3]+f[j][0][6]+f[j][0][7]))/(1-u0);
        f[j][0][1]=f[j][0][3] + 2*rhow*u0/3;
        f[j][0][5]=f[j][0][7] + rhow*u0/6;
        f[j][0][8]=f[j][0][6] + rhow*u0/6;;
    }
    
    //	i=Nx: right wall
    for(j=0;j<=Ny;j++)
    {
        f[j][Nx][1]=2.0*f[j][Nx-1][1] - f[j][Nx-2][1];
        f[j][Nx][5]=2.0*f[j][Nx-1][5] - f[j][Nx-2][5];
        f[j][Nx][8]=2.0*f[j][Nx-1][8] - f[j][Nx-2][8];
    }
*///----------------------Mohamad, flow in a channel ---------------------


    //----------------  PBC -------------
    // i=0: left wall
    for(j=0;j<=Ny;j++)
    {
        f[j][0][1]=f[j][Nx][1];
        f[j][0][5]=f[j][Nx][5];
        f[j][0][8]=f[j][Nx][8];
    }
    
    //	i=Nx: right wall
    for(j=0;j<=Ny;j++)
    {
        f[j][Nx][3]=f[j][0][3];
        f[j][Nx][7]=f[j][0][7];
        f[j][Nx][6]=f[j][0][6];
    }
    //----------------  PBC -------------

/*
    // i=0: left wall
    for(j=0;j<=Ny;j++)
    {
        f[j][0][1]=f_post[j][0][3]; 
        f[j][0][5]=f_post[j][0][7]; 
        f[j][0][8]=f_post[j][0][6];
    }
    
    //	i=Nx: right wall
    for(j=0;j<=Ny;j++)
    {
        f[j][Nx][3]=f_post[j][Nx][1]; 
        f[j][Nx][7]=f_post[j][Nx][5]; 
        f[j][Nx][6]=f_post[j][Nx][8];
    }
*/
    //printf("End bounce back....\n");
}

void Bounce_backV()
{
    //printf("Begin bounce back V.... ");

    int i,j;
    float u0=0.10;
   
    //	j=Ny: top plate
    for(i=0;i<=Nx;i++)
    {
        f[Ny][i][4]=f_post[0][i][4];
        f[Ny][i][7]=f_post[0][i][7];
        f[Ny][i][8]=f_post[0][i][8];
    }
    
    //	j=0: bottom plate
    
    for(i=0;i<=Nx;i++)
    {
        f[0][i][2]=f_post[Ny][i][2];
        f[0][i][5]=f_post[Ny][i][5];
        f[0][i][6]=f_post[Ny][i][6];
    }
    
    
    
    
     // i=0: left wall
     for(j=0;j<=Ny;j++)
     {
         f[j][0][1]=f_post[j][0][3];
         f[j][0][5]=f_post[j][0][7];
         f[j][0][8]=f_post[j][0][6];
     }
     
     //	i=Nx: right wall
     for(j=0;j<=Ny;j++)
     {
         f[j][Nx][3]=f_post[j][Nx][1];
         f[j][Nx][7]=f_post[j][Nx][5];
         f[j][Nx][6]=f_post[j][Nx][8];
     }
    
    //printf("End bounce back V....\n");

}


//=========================================================


//=========================================================
//------------------------------------------------------------

// Subroutine: Fluid variables (density and velocity)
//------------------------------------------------------------
void Den_Vel()
{
    
    int j, i;
    
    for(j=0;j<=Ny;j++) for(i=0;i<=Nx;i++)  
    {
        if (!boundary[j][i])
        {
        
            rho[j][i]=f[j][i][0]+f[j][i][1]+f[j][i][2]+f[j][i][3]+f[j][i][4]+f[j][i][5]+f[j][i][6]+f[j][i][7]+f[j][i][8];

            ux[j][i]=(f[j][i][1]+f[j][i][5]+f[j][i][8]-f[j][i][3]-f[j][i][6]-f[j][i][7])/rho[j][i]; 
            
            uy[j][i]=(f[j][i][5]+f[j][i][6]+f[j][i][2]-f[j][i][7]-f[j][i][8]-f[j][i][4])/rho[j][i];
            
        }
            
    }
    
}

//=========================================================

double Err()	// Calculating the relative difference in velocity between two steps
//=========================================================
{
    int j, i; 
    double e1,e2;
    
    e1=e2=0.0;
    
    for(j=1;j<Ny;j++) for(i=0;i<Nx;i++)   
    {
        
        e1+=sqrt((ux[j][i]-u0[j][i])*(ux[j][i]-u0[j][i]) +(uy[j][i]-v0[j][i])*(uy[j][i]-v0[j][i]));
        
        e2+=sqrt(ux[j][i]*ux[j][i]+uy[j][i]*uy[j][i]);
        
        u0[j][i]=ux[j][i];v0[j][i]=uy[j][i];
        
    }
    
    return e1/e2;
    
}

//=========================================================

void	Data_Output() // Output data
//=========================================================


{
    
    int i,j; FILE *fp;
    
    fp=fopen("x.dat","w+");
    
    for(i=0;i<=Nx;i++) fprintf(fp,"%e \n", (i+0.5)/L); fclose(fp);
    
    fp=fopen("y.dat","w+");
    
    for(j=0;j<=Ny;j++) fprintf(fp,"%e \n", (j+0.5)/L); fclose(fp);
    
    
    fp=fopen("magnitud.txt", "w+");
    for(j=0;j<=Ny;j++)
        for(i=0;i<=Nx;i++)
        {
            //if (!boundary[j][i])
            //    fprintf(fp,"%e %e %e\n", (i+0.5)/L, (j+0.5)/L, sqrt(pow(ux[j][i],2) + pow(uy[j][i],2)));
            fprintf(fp,"%e\n", sqrt(pow(ux[j][i],2) + pow(uy[j][i],2)));
        }
    fclose(fp);

    fp=fopen("velocidad.txt", "w+");
    for(j=0;j<=Ny;j++)
        for(i=0;i<=Nx;i++)
        {
            //if (!boundary[j][i])
                fprintf(fp,"%e %e\n", ux[j][i], uy[j][i]);
        }
    fclose(fp);
    
    
    fp=fopen("ux.txt","w");
    for(j=0;j<=Ny;j++) {
        
        for (i=0; i<=Nx; i++) fprintf(fp,"%e ",ux[j][i]); fprintf(fp,"\n");
        
    }
    
    fclose(fp);
    
    fp=fopen("uy.txt","w");
    
    for(j=0;j<=Ny;j++){
        
        for (i=0; i<=Nx; i++) fprintf(fp,"%e ",uy[j][i]); fprintf(fp,"\n");
        
    }
    
    
    
    fclose(fp);
    
    fp=fopen("rho.txt","w");
    
    for(j=0;j<=Ny;j++){
        
        for (i=0; i<=Nx; i++) fprintf(fp,"%e ",rho[j][i]); fprintf(fp,"\n");
        
    }
    
    fclose(fp);
    
    
    fp=fopen("solid.txt","w");
    
    for(j=0;j<Ny1;j++){
        
        for (i=0; i<Nx1; i++) fprintf(fp,"%d ",boundary[j][i]); fprintf(fp,"\n");
        
    }
    
    fclose(fp);
        
}


//=========================================================

void	force(void)//
//=========================================================
{
    //	i=Nx: right wall
    double frce;
    double omega=1;
    double cs2=1.0/3.0;
    double uf=0.02;
    double visc=(1.0/omega-0.5)*cs2;
    double fpois=8.0*visc*uf/Ny1/Ny1;
    double rho_0=1.0;
    fpois=rho_0*fpois/6;
    frce=fpois;
    for(int j=0;j<=Ny;j++)
    {
        f[j][Nx][3]=f[j][0][3]-frce;
        f[j][Nx][7]=f[j][0][7]-frce;
        f[j][Nx][6]=f[j][0][6]-frce;
    }
    
    // i=0: left wall
    for(int j=0;j<=Ny;j++)
    {
        f[j][0][1]=f[j][Nx][1]+frce;
        f[j][0][5]=f[j][Nx][5]+frce;
        f[j][0][8]=f[j][Nx][8]+frce;
    }
}

    //=========================================================
    
    void	forceV(void)//
    //=========================================================
    {
        //	i=Nx: right wall
        double frce;
        double omega=1;
        double cs2=1.0/3.0;
        double uf=0.02;
        double visc=(1.0/omega-0.5)*cs2;
        double fpois=8.0*visc*uf/Ny1/Ny1;
        double rho_0=1.0;
        fpois=rho_0*fpois/6;
        frce=fpois;
        
        for (int i=0; i<=Nx; i++)
        {
            for (int j=0; j<=Ny; j++)
            {
                
                if (!boundary[j][i])
                {
                    f[j][i][4]+=frce;
                    f[j][i][7]+=frce;
                    f[j][i][8]+=frce;

                    f[j][i][2]-=frce;
                    f[j][i][5]-=frce;
                    f[j][i][6]-=frce;
                }
            }
        
        }
        
/*
        //	j=Ny: top plate
        for(int i=0;i<=Nx;i++)
        {
            f[Ny][i][4]=f_post[0][i][4]+frce;
            f[Ny][i][7]=f_post[0][i][7]+frce;
            f[Ny][i][8]=f_post[0][i][8]+frce;
        }

        //	j=0: bottom plate
        for(int i=0;i<=Nx;i++)
        {
            f[0][i][2]=f_post[Ny][i][2]-frce;
            f[0][i][5]=f_post[Ny][i][5]-frce;
            f[0][i][6]=f_post[Ny][i][6]-frce;
        }
 */
    
    }


void mkeSolid(int x1, int y1, int x2, int y2)
{
    for (int j=y1; j<=y2; j++)
    {
        for (int i=x1; i<=x2; i++)
        {
            boundary[j][i]=1;
        }
        
    }

}

void mkePorous(double percent)
{
    int i,j;
    srand((unsigned int)time(NULL));
    for(j=1;j<Ny;j++)
    {
        for(i=1;i<Nx;i++)
        {
            //if (!boundary[j][i])
            
            boundary[j][i] = (rand()>RAND_MAX*(percent/100))? 0 : 1;
        }
    }

}

void readBoundaryFromFile(const char * fileName)
{
    
    FILE *fp;
    int c;
    int xcols=0, yrows=0;
    int lineLength=0;
    char *line;
    int  *lineCopy;
    
    //char line[2000];
    fp=fopen(fileName,"r");
    if (fp != NULL)
    {
        do
        {
            c=fgetc(fp);
            if (c!=0x20 && c!=10)
            {
                xcols++;
            }
            
        }while (c!=10);
        printf("Número de columnas: %d\n",xcols);
        
        lineLength=xcols+xcols-1+1; //cols + blanks + LF
        line=(char *)calloc(lineLength, sizeof(char));
        
        //int res=fseek(fp, 0, SEEK_SET);
        fclose(fp);
        fp=fopen(fileName,"r");
        while (!feof(fp))
        {
            fread(line,1, lineLength, fp);
            yrows++;
        }
        yrows--;
        printf("Número de renglones: %d\n",yrows);
        
        Nx=xcols;
        Ny=yrows;
        Nx1=Nx+1;
        Ny1=Ny+1;
        
        
        lineCopy=(int *)calloc(xcols,      sizeof(int));
        int num;

        
        
        boundary=new2DIntArray(Nx1, Ny1, sizeof(int));
        int res=fseek(fp, 0, SEEK_SET);
        for (int j=yrows; j>=1; j--)
        {
            int pointer=(j-1)*lineLength;
            fseek(fp, pointer, SEEK_SET);
            fread(line, 1, lineLength, fp);
            int i=0;
            for (int n=0; n<lineLength; n++)
            {
                if(line[n]!=0x20 && line[n]!=10) //if not a space and not a LF (just 1's and 0's)
                {
                    num=line[n]-'0';
                    lineCopy[i]=num;
                    boundary[j-1][i++]=num;
                }
            }
            
            
        }
        fclose(fp);
        
        fp=fopen("matrix1.txt","w+");
        if (fp!=NULL)
        {
            for (int j=0; j<yrows; j++)
            {
                for (int i=0; i<xcols; i++)
                {
                    fprintf(fp,"%d ",boundary[j][i]);
                }
                fprintf(fp,"\n");
            }
        }
        fclose(fp);

        
        
        
        
//        int (*matrix)[xcols]=malloc(sizeof(int)*yrows*xcols);
//        
//        for (int j=0; j<yrows; j++)
//        {
//            for (int i=0; i<xcols; i++)
//            {
//                matrix[j][i]=0;
//            }
//        }
//        
//        
//        fp=fopen("matrix0.txt","w+");
//        for (int j=0; j<yrows; j++)
//        {
//            for (int i=0; i<xcols; i++)
//            {
//                fprintf(fp,"%d ",matrix[j][i]);
//            }
//            fprintf(fp,"\n");
//        }


    }
    

}



