#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "openmx_common.h"
#include "mpi.h"
#include <omp.h>

void Calc_dDen(double *Den0,double *Den1,double *Den2,double *Den3)
{
  int MN,MN1,MN2,i,j,k,ri,ri1,ri2;
  int i1,i2,j1,j2,k1,k2,nmax,n;
  int Ng1,Ng2,Ng3;
  int DN,BN, N2D,spin,axis;
  int numprocs,myid;
  double den_min=1.0e-14;
  double ***dDen_Grid_D;
  double detA, igtv[4][4];
  int OMPID,Nthrds;
  double up_x_a,up_x_b,up_x_c;
  double up_y_a,up_y_b,up_y_c;
  double up_z_a,up_z_b,up_z_c;
  double dn_x_a,dn_x_b,dn_x_c;
  double dn_y_a,dn_y_b,dn_y_c;
  double dn_z_a,dn_z_b,dn_z_c;
  double up_a,up_b,up_c;
  double dn_a,dn_b,dn_c;

  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);
  dDen_Grid_D = (double***)malloc(sizeof(double**)*2); 
    for (k=0; k<=1; k++){
      dDen_Grid_D[k] = (double**)malloc(sizeof(double*)*3); 
      for (i=0; i<3; i++){
        dDen_Grid_D[k][i] = (double*)malloc(sizeof(double)*My_NumGridD); 
        for (j=0; j<My_NumGridD; j++) dDen_Grid_D[k][i][j] = 0.0;
      }
    }

    // dDen_Grid_B = (double***)malloc(sizeof(double**)*2); 
    // for (k=0; k<=1; k++){
    //   dDen_Grid_B[k] = (double**)malloc(sizeof(double*)*3); 
    //   for (i=0; i<3; i++){
    //     dDen_Grid_B[k][i] = (double*)malloc(sizeof(double)*My_NumGridB_AB); 
    //     for (j=0; j<My_NumGridB_AB; j++) dDen_Grid_B[k][i][j] = 0.0;
    //   }
    // }
  detA =   gtv[1][1]*gtv[2][2]*gtv[3][3]
           + gtv[1][2]*gtv[2][3]*gtv[3][1]
           + gtv[1][3]*gtv[2][1]*gtv[3][2]
           - gtv[1][3]*gtv[2][2]*gtv[3][1]
           - gtv[1][2]*gtv[2][1]*gtv[3][3]
           - gtv[1][1]*gtv[2][3]*gtv[3][2];     

  igtv[1][1] =  (gtv[2][2]*gtv[3][3] - gtv[2][3]*gtv[3][2])/detA;
  igtv[2][1] = -(gtv[2][1]*gtv[3][3] - gtv[2][3]*gtv[3][1])/detA;
  igtv[3][1] =  (gtv[2][1]*gtv[3][2] - gtv[2][2]*gtv[3][1])/detA; 

  igtv[1][2] = -(gtv[1][2]*gtv[3][3] - gtv[1][3]*gtv[3][2])/detA;
  igtv[2][2] =  (gtv[1][1]*gtv[3][3] - gtv[1][3]*gtv[3][1])/detA;
  igtv[3][2] = -(gtv[1][1]*gtv[3][2] - gtv[1][2]*gtv[3][1])/detA; 

  igtv[1][3] =  (gtv[1][2]*gtv[2][3] - gtv[1][3]*gtv[2][2])/detA;
  igtv[2][3] = -(gtv[1][1]*gtv[2][3] - gtv[1][3]*gtv[2][1])/detA;
  igtv[3][3] =  (gtv[1][1]*gtv[2][2] - gtv[1][2]*gtv[2][1])/detA; 

  #pragma omp parallel shared(My_NumGridD,Min_Grid_Index_D,Max_Grid_Index_D,igtv,dDen_Grid_D,PCCDensity_Grid_D,PCC_switch,Den0,Den1,Den2,Den3,den_min) private(OMPID,Nthrds,nmax,n,i,j,k,ri,ri1,ri2,i1,i2,j1,j2,k1,k2,MN,MN1,MN2,up_a,dn_a,up_b,dn_b,up_c,dn_c,Ng1,Ng2,Ng3)
    {

      OMPID = omp_get_thread_num();
      Nthrds = omp_get_num_threads();

      Ng1 = Max_Grid_Index_D[1] - Min_Grid_Index_D[1] + 1;
      Ng2 = Max_Grid_Index_D[2] - Min_Grid_Index_D[2] + 1;
      Ng3 = Max_Grid_Index_D[3] - Min_Grid_Index_D[3] + 1;

      for (MN=OMPID; MN<My_NumGridD; MN+=Nthrds){

        i = MN/(Ng2*Ng3);
	j = (MN-i*Ng2*Ng3)/Ng3;
	k = MN - i*Ng2*Ng3 - j*Ng3; 

      if ( i==0 || i==(Ng1-1) || j==0 || j==(Ng2-1) || k==0 || k==(Ng3-1) ){

	  dDen_Grid_D[0][0][MN] = 0.0;
	  dDen_Grid_D[0][1][MN] = 0.0;
	  dDen_Grid_D[0][2][MN] = 0.0;
	  dDen_Grid_D[1][0][MN] = 0.0;
	  dDen_Grid_D[1][1][MN] = 0.0;
	  dDen_Grid_D[1][2][MN] = 0.0;
        }

        else {

          /* set i1, i2, j1, j2, k1, and k2 */ 

          i1 = i - 1;
          i2 = i + 1;

          j1 = j - 1;
          j2 = j + 1;

          k1 = k - 1;
          k2 = k + 1;

	  /* set dDen_Grid_D */

	  if ( den_min<(Den0[MN]+Den1[MN]) ){

	    /* a-axis */

	    MN1 = i1*Ng2*Ng3 + j*Ng3 + k;
	    MN2 = i2*Ng2*Ng3 + j*Ng3 + k;

	    if (PCC_switch==0) {
	      up_a = Den0[MN2] - Den0[MN1];
	      dn_a = Den1[MN2] - Den1[MN1];
	    }
	    else if (PCC_switch==1) {
	      up_a = Den0[MN2] + PCCDensity_Grid_D[0][MN2]
	           - Den0[MN1] - PCCDensity_Grid_D[0][MN1];
	      dn_a = Den1[MN2] + PCCDensity_Grid_D[1][MN2]
	           - Den1[MN1] - PCCDensity_Grid_D[1][MN1];
	    }

	    /* b-axis */

	    MN1 = i*Ng2*Ng3 + j1*Ng3 + k; 
	    MN2 = i*Ng2*Ng3 + j2*Ng3 + k; 

	    if (PCC_switch==0) {
	      up_b = Den0[MN2] - Den0[MN1];
	      dn_b = Den1[MN2] - Den1[MN1];
	    }
	    else if (PCC_switch==1) {
	      up_b = Den0[MN2] + PCCDensity_Grid_D[0][MN2]
	           - Den0[MN1] - PCCDensity_Grid_D[0][MN1];
	      dn_b = Den1[MN2] + PCCDensity_Grid_D[1][MN2]
	           - Den1[MN1] - PCCDensity_Grid_D[1][MN1];
	    }

	    /* c-axis */

	    MN1 = i*Ng2*Ng3 + j*Ng3 + k1; 
	    MN2 = i*Ng2*Ng3 + j*Ng3 + k2; 

	    if (PCC_switch==0) {
	      up_c = Den0[MN2] - Den0[MN1];
	      dn_c = Den1[MN2] - Den1[MN1];
	    }
	    else if (PCC_switch==1) {
	      up_c = Den0[MN2] + PCCDensity_Grid_D[0][MN2]
	           - Den0[MN1] - PCCDensity_Grid_D[0][MN1];
	      dn_c = Den1[MN2] + PCCDensity_Grid_D[1][MN2]
	           - Den1[MN1] - PCCDensity_Grid_D[1][MN1];
	    }

	    /* up */

	    dDen_Grid_D[0][0][MN] = 0.5*(igtv[1][1]*up_a + igtv[1][2]*up_b + igtv[1][3]*up_c);
	    dDen_Grid_D[0][1][MN] = 0.5*(igtv[2][1]*up_a + igtv[2][2]*up_b + igtv[2][3]*up_c);
	    dDen_Grid_D[0][2][MN] = 0.5*(igtv[3][1]*up_a + igtv[3][2]*up_b + igtv[3][3]*up_c);

	    /* down */

	    dDen_Grid_D[1][0][MN] = 0.5*(igtv[1][1]*dn_a + igtv[1][2]*dn_b + igtv[1][3]*dn_c);
	    dDen_Grid_D[1][1][MN] = 0.5*(igtv[2][1]*dn_a + igtv[2][2]*dn_b + igtv[2][3]*dn_c);
	    dDen_Grid_D[1][2][MN] = 0.5*(igtv[3][1]*dn_a + igtv[3][2]*dn_b + igtv[3][3]*dn_c);

	  }

	  else{
	    dDen_Grid_D[0][0][MN] = 0.0;
	    dDen_Grid_D[0][1][MN] = 0.0;
	    dDen_Grid_D[0][2][MN] = 0.0;
	    dDen_Grid_D[1][0][MN] = 0.0;
	    dDen_Grid_D[1][1][MN] = 0.0;
	    dDen_Grid_D[1][2][MN] = 0.0;
	  }

	} /* else */
}
    #pragma omp flush(dDen_Grid_D)

    }

  Ng1 = Max_Grid_Index_D[1] - Min_Grid_Index_D[1] + 1;
  Ng2 = Max_Grid_Index_D[2] - Min_Grid_Index_D[2] + 1;
  Ng3 = Max_Grid_Index_D[3] - Min_Grid_Index_D[3] + 1;


  // printf(" dDen_Grid_D[spin][axis][DN] = %20.12f", dDen_Grid_D[0][0][100]);

  for (n=0; n<Num_Rcv_Grid_B2D[myid]; n++){
      DN = Index_Rcv_Grid_B2D[myid][n];
      BN = Index_Snd_Grid_B2D[myid][n];

      i = DN/(Ng2*Ng3);
      j = (DN-i*Ng2*Ng3)/Ng3;
      k = DN - i*Ng2*Ng3 - j*Ng3; 

      if ( !(i<=1 || (Ng1-2)<=i || j<=1 || (Ng2-2)<=j || k<=1 || (Ng3-2)<=k)){
        for (spin=0; spin<=SpinP_switch; spin++){
          for (axis = 0;axis<=2;axis++){
          dDen_Grid_B[spin][axis][BN] = dDen_Grid_D[spin][axis][DN];
          }
        }
      }
  }     

}


void Calc_g(double *Den0,double *Den1,double *Den2,double *Den3)
{
  int MN,MN1,MN2,i,j,k;
  int i1,i2,j1,j2,k1,k2;
  int Ng1,Ng2,Ng3;
  int DN,BN, N2D,spin,axis;
  int numprocs,myid;
  double den_min=1.0e-14;
  double detA, igtv[4][4];
  int OMPID,Nthrds;
  double up_x_a,up_x_b,up_x_c;
  double up_y_a,up_y_b,up_y_c;
  double up_z_a,up_z_b,up_z_c;
  double dn_x_a,dn_x_b,dn_x_c;
  double dn_y_a,dn_y_b,dn_y_c;
  double dn_z_a,dn_z_b,dn_z_c;
  double up_a,up_b,up_c;
  double dn_a,dn_b,dn_c;
  double sden[2],tden,abs_dden;
  double My_g1,My_g2;

  My_g1 = 0.0;
  My_g2 = 0.0;

  Calc_dDen(Den0,Den1,Den2,Den3);

  for (BN=0; BN<My_NumGridB_AB; BN++){
  sden[0] = Density_Grid_B[0][BN];
  sden[1] = Density_Grid_B[1][BN];
  tden = sden[0] + sden[1];
  abs_dden = 0.0;
    for (axis = 0;axis<=2;axis++){
    abs_dden += (dDen_Grid_B[0][axis][BN]+dDen_Grid_B[1][axis][BN])*(dDen_Grid_B[0][axis][BN]+dDen_Grid_B[1][axis][BN]);
    }
  
  abs_dden = sqrt(abs_dden);
  // printf("abs_dden = %20.12f",abs_dden);
  // My_g += tden;F
  My_g1 += abs_dden/tden;
  My_g2 += sqrt(abs_dden/tden);
  
  }

My_g1 *=GridVol/Cell_Volume;
My_g2 *= GridVol/Cell_Volume;
  // My_g *= GridVol/Cell_Volume;

  MPI_Allreduce(&My_g1,&g1,1,MPI_DOUBLE,MPI_SUM,mpi_comm_level1);
  MPI_Allreduce(&My_g2,&g2,1,MPI_DOUBLE,MPI_SUM,mpi_comm_level1);
}
    

void Calc_Di_Pe(int knum_i,int Knum_j,int knum_k,int SpinP_sqitch,double ***EIGEN,int ***k_op,int *T_k_op)
{
  double Unit,homo_e,lumo_e,Ave_Gap;
  int nkpoint,T_knum,i,j,k,kloop,homo_index,lumo_index;
  double VD,PF,FEFE;

  Unit = 27.2113845;
  KSGap = 0.0;

  Total_N_U = 0.0;
  T_knum = 0;

  for (i = 0; i < Kspace_grid1; i++)
  {
    for (j = 0; j < Kspace_grid2; j++)
    {
      for (k = 0; k < Kspace_grid3; k++)
      {
        if (0 < k_op[i][j][k])
        {
          T_knum++;
        }
      }
    }
  }

  for (i = 1; i <= atomnum; i++)
  {
    Total_N_U += InitN_USpin[i];
  }

  Total_N_U = round(Total_N_U);
  homo_index = (int)Total_N_U;
  lumo_index = homo_index + 1;
  homo_e = EIGEN[0][0][homo_index];
  lumo_e = EIGEN[0][0][lumo_index];

  Ave_Gap = 0.0;
  nkpoint = 0.0;

  // printf("Total_N_U = %f", Total_N_U);

  for (kloop = 0; kloop < T_knum; kloop++)
  {
    if (0 < T_k_op[kloop])
    {
      if (SpinP_switch == 0)
      {
        nkpoint += 1;
        Ave_Gap += EIGEN[0][kloop][lumo_index] - EIGEN[0][kloop][homo_index];
        if (EIGEN[0][kloop][homo_index] > homo_e)
        {
          homo_e = EIGEN[0][kloop][homo_index];
        }
        if (EIGEN[0][kloop][lumo_index] < lumo_e)
        {
          lumo_e = EIGEN[0][kloop][lumo_index];
        }
      }
      else if (SpinP_switch == 1)
      {
        // fprintf(fp_EV,"%5d  %18.14f %18.14f\n",
        //   l,EIGEN[0][kloop][l],EIGEN[1][kloop][l]);
      }
    }
  }
  if ((int)Total_N_U % 2 == 0)
  {
    if (lumo_e - homo_e > KSGap)
    {
      KSGap = (lumo_e - homo_e) * Unit;
      // printf("KSGap = %f\n", KSGap);
      // printf("homo_e = %f\n", homo_e * Unit);
      // printf("lumo_e = %f\n", lumo_e * Unit);
    }
  }
  else
  {
    // printf("KSGap = %f\n", KSGap);
  }

  Ave_Gap /= nkpoint;

  VD = Total_N_U * 2 / Cell_Volume;
  PF = sqrt(4 * M_PI * VD);
  FEFE = pow(3 * M_PI * M_PI * VD, 2 / 3) / 2;
  Di_Pe = 1 + pow(PF / Ave_Gap, 2) * (1 - Ave_Gap / 4 / FEFE);

  // printf("Penn model dielectric = %f", Di_Pe);
}