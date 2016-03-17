//
//  main.c
//  LBM_D2Q9
//
//  Created on 2/9/16.
//  Copyright © 2016 Juvenal. All rights reserved.
//

#include <stdio.h>

/*
 * D2Q9.cpp
 *
 *  Created on: Feb 9, 2016
  */


//=========================================================

//=========================================================
#include "D2Q9.h"
int main()

{
    
    int k,M2,N2;
    double err;
    M2=Ny/2;
    N2=Nx/2;
    
    k=0;
    
    err=1.0;
   
    
    tau=3*L*uw/Re+0.5; // relaxation time for BGK
    s[7]=s[8]=1.0/tau;	s[0]=s[3]=s[5]=0.0;	s[4]=s[6]=8*(2-
    s[7])/(8-s[7]);	s[1]=1.6;	s[2]=1.8; // relaxation rates for MRT
    
    Init_Eq();
    mkeSolid(0, 0, 0, Ny);
    mkeSolid(Nx,0,Nx,Ny);
    //mkeSolid(30, 30, 60, 60);
    mkePorous(25);
    
    while(err>1.0e-6)
    //while(err>1.0e-3)
    {
        
        k++;
        Coll_BGK();	//BGK collision
        //	Coll_MRT();	//MRT collision
        Streaming();	// Streaming
        Bounce_backV(); // No-Slip boundary condition
        forceV();
        Den_Vel();	// Fluid variables
        
        if(k%1000==0)
            
        {
            
            err=Err(); // Velocity differences between two successive 1000 steps
            printf("err=%e ux_center=%e uy_center=%e k=%d\n",err,ux[M2][N2],uy[M2][N2], k);  // Display some results
        }
        
    }
    Data_Output();	// Output simulation data
    
    return 1;
}

