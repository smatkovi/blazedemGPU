
/* Change List
 * 1) 23/12/2012 -- adding OpenGL VBO*/

/*---------------------------------------------------------------------------*/
                    /* REQUIRED INCLUDE FILES*/
/*---------------------------------------------------------------------------*/
                       /* System C++ */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <iostream>
#include <fstream>

                       /* System CUDA */
#include <cuda_runtime.h>
#include <thrust/for_each.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/sort.h>
#include <thrust/remove.h>
#include <thrust/scan.h>
#include <thrust/count.h>
#include <thrust/sequence.h>
#include <thrust/device_vector.h>
#include <thrust/system/cuda/execution_policy.h>
/*------------------------------------------------*/

                       /* Data Structures  */
#include "../DataStructures/KSimulationObjects.h"

                       /* Public Methods  */
#include "DeviceInterface.h"

/*---------------------------------------------------------------------------*/
/*                          Device includes                                  */
/*---------------------------------------------------------------------------*/
#include "../Utilities/Device/Opertators.cuh"
#include "../Utilities/Device/Functions.cuh"

#include <device_launch_parameters.h>
#include "GPU_Memory_Transfer.cuh"

#include <string>
#include <sstream>


/* Particle Dynamics Information Global memory  */

/*---------------- SINGLE GPU---------------- */
/* Translational Parameters */
float3    *dram_position_com;
float3    *dram_velocity_com;

/* Angular Parameters */
Quaterion *dram_position_ornt;
float3    *dram_velocity_ang;

/* Particle Geometry identifiers */
uint   *dram_ObjectType;
int    *dram_P_ID;


     /* non-symmetry  PP Forces */
float3 *dram_force_PP;
float3 *dram_force_com_Wall;
float3 *dram_force_com_Lifter;

float3 *dram_force_PP_ang;
float3 *dram_force_ang_Lifter;
float3 *dram_force_ang_Wall;


 /* Symmetry PP Forces */
float *dram_force_com_X;
float *dram_force_com_Y;
float *dram_force_com_Z;

float *dram_force_ang_X;
float *dram_force_ang_Y;
float *dram_force_ang_Z;


/*----------------------------------------------------------------------*/


/* Contact point info for post processing */
float3 *dWallContactPoints;  /* Optional */
float3 *dVOContactPoints;/* Optional */
float3 *dPContactPoints;/* Optional */
uint   *dNumContacts;

/* Spatial Grid parameters */
uint  *m_dCellStart;
uint  *m_dCellEnd;
int    m_numGridCells;

uint  *dram_GridParticleHash; /* Grid hash value for each particle */
int   *dram_GridParticleIndex;/* Particle index for each particle used to sort */

/* History contact for particle Lifter/DObject */
Contact_Info *dram_Lifter_Contact_Hist;

uint  *dBroad_List;
uint  *dNumNN;


/* Temp arrays for sorting */
/* Reorder arrays based on hash */
float3       *dSortedPos;
Quaterion    *dSortedPosQ;
float3       *dSortedVel;
float3       *dSorted_velocity_ang;
uint         *dSortedPType;
int          *dSortedPID;
Contact_Info *dSortedLifter_Contact_Hist;



/*---------------- END SINGLE GPU---------------- */


/*----------------------------------------------------------------------*/
                            /* Multi-GPU */
/*----------------------------------------------------------------------*/
/* Particle Dynamics Information Global memory  */

/* Translational Parameters */
float3    *dram_position_comM  [16];
float3    *dram_velocity_comM  [16];

/* Angular Parameters */
Quaterion *dram_position_orntM [16];
float3    *dram_velocity_angM  [16];

/* Particle Geometry identifiers */
uint   *dram_ObjectTypeM       [16];
int    *dram_P_IDM             [16];


     /* non-symmetry  PP Forces */
float3 *dram_force_PPM         [16];
float3 *dram_force_com_WallM   [16];
float3 *dram_force_com_LifterM [16];

float3 *dram_force_PP_angM     [16];
float3 *dram_force_ang_LifterM [16];
float3 *dram_force_ang_WallM   [16];


 /* Symmetry PP Forces */
float *dram_force_com_XM       [16];
float *dram_force_com_YM       [16];
float *dram_force_com_ZM       [16];

float *dram_force_ang_XM       [16];
float *dram_force_ang_YM       [16];
float *dram_force_ang_ZM       [16];

/* Contact point info for post processing */
float3 *dWallContactPointsM [16]; /* Optional */
float3 *dVOContactPointsM   [16]; /* Optional */
float3 *dPContactPointsM    [16]; /* Optional */
uint   *dNumContactsM[16];

/* Spatial Grid parameters */
uint  *m_dCellStartM           [16];
uint  *m_dCellEndM             [16];

uint  *dram_GridParticleHashM  [16]; /* Grid hash value for each particle */
int   *dram_GridParticleIndexM [16];/* Particle index for each particle used to sort */

uint  *dBroad_ListM            [16];
uint  *dNumNNM                 [16];

/* History contact for particle Lifter/DObject */
Contact_Info *dram_Lifter_Contact_HistM[16];
/*----------------------------------------------------------------------*/

/* Since full c++ is not supported include here */
#include "Utilities.cuh"

#include "Physics/Collision_Response.cuh"
/*                        DEM Kernels for Polyhedra                          */
#include "Computing/Collision_Detection_BroadPhase.cuh"
/*                        DEM Kernels for Spheres                            */
#include "Computing/DEM_Kernels_Compute.cuh"
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
                         /* Local Variables */
/*---------------------------------------------------------------------------*/

ofstream OF; /* Output Log File */
ofstream        LogFileD;
ofstream        EnergyBin;
/*----------------------------------------------------------------------*/



/* Flags */

bool Device_use_multi_gpu = false;
int  run_level            = 0;
int  num_sm;

bool start_bin            = false;
bool is_cylinder_rotation = false;
bool Get_WallPoints       = false;
int  error_check          = 0;


/* GPU Information */
int  num_devices;

cudaDeviceProp DevProp;
size_t         t[16];
size_t         a[16];


int  device_num = 0;            /* Which device to use single gpu */
int  num_devices_use;           /* How many devices to use for muti-gpu */
int  num_particles_perGPU [16]; 

cudaStream_t streams [64]; /* 4 streams and events per GPU */
cudaEvent_t  events  [64];

launchParms   GlobalLaunch; /* For single GPU */
launchParms   multi[16];    /* For multi GPU */


/* Local Information */
int             NUMPARTICLES_MAX_SIM;

int             NUMWOBJECTS;
WOBJ            *m_WorldObject;
int             NUMDYOBJECTS;
SimulationInfo *SimInfo_C;


/* Data used for inital filling */
float3          *h_PosPack;           /* Initial Pos */
Quaterion       *h_OrntPack;          /* Initial Ornt */
uint            *h_Particle_TypePack; /* Initial Particle TypeID */
int             *h_Particle_IDPack;   /* Initial Particle ID */


/*---------------------------------------------*/
       /* Data for initial filling */
int    NUMPARTICLES_Current;
bool   isFilling = false;
float3 pack_pos;
int3   pack_size;

int    packNumber;
int    packNumWaves;
int    threads_perBlock;
int    fill_counter     = 0;
int    pack_clear_steps = -1;
bool   is_first         = true;
float  fill_time        = 0.0f;
/*---------------------------------------------*/

/*---------------------------------------------*/
/* Data for removing Flagged particles */
int   num_dead = 0;
float dis_mass = 0.0;
float dis_vol  = 0.0;
int   stepc    = 0;

int D_NUMPARTICLES; /* Used to kill particles */
/*---------------------------------------------*/




/*---------------------------------------------*/
/* Data for physics tallies */
int NumPhysicsTallies = 3;
/*---------------------------------------------*/


float Call_Time=0.0f;


bool NO_NN(true);

bool get_wall_forces    = false;
bool get_vobject_forces = false;
bool get_P_forces = false;



/*---------------------------------------------------------------------------*/
          /* Set the Initial distribution as a AxBxC grid */
/*---------------------------------------------------------------------------*/
void Set_Position_DefaultGrid( InitConfig *h_PosConfig, POBJ *h_POBJ )
{

	 if(error_check==1)
	 {
	    	cudaError_t errormsg=cudaGetLastError();
	    	if(errormsg>0)
	    	{
	    	cout<<"Initial Position Default Grid: Start "<<cudaGetErrorString(errormsg)<<endl;
	    	 exit(1);
	    	}
	 }

	cout<<"INFO-D: Allocating Initial positions on HOST \n";

	h_PosPack             = new float3    [NUMPARTICLES_MAX_SIM];

	 if(SimInfo_C->particle_type==1 || SimInfo_C->sphere_orient)
	 {
	   h_OrntPack            = new Quaterion [NUMPARTICLES_MAX_SIM];
	   cout<<"INFO-D: Polyhedra Setting Orient \n";
	 }

	h_Particle_TypePack   = new uint      [NUMPARTICLES_MAX_SIM];
	h_Particle_IDPack     = new int       [NUMPARTICLES_MAX_SIM];


	cout<<"INFO-D: Creating Fixed Position Grid  "<<NUMPARTICLES_MAX_SIM<<endl;

	

	cout<<" \n";

    float x,y,z;
    float x_size,y_size,z_size;


    int    NX  = h_PosConfig->num.x;
    int    NY  = h_PosConfig->num.y;
    int    NZ  = h_PosConfig->num.z;

    x_size = h_PosConfig->p_size[0].x;
    y_size = h_PosConfig->p_size[0].y;
    z_size = h_PosConfig->p_size[0].z;


    x      = h_PosConfig->start.x;
    y      = h_PosConfig->start.y;
    z      = h_PosConfig->start.z;



    cout<<"INFO-D: Start  "<<x<<", "<<y<<", "<<z<<endl;
    cout<<"  Grid  "<<NX<<", "<<NY<<", "<<NZ<<endl;

    /* Allocate positions */
    for ( int yd=0; yd< NY; yd++ )/* Height  */
      {

         for ( int zd=0; zd<NZ; zd++ )/* x B-F */
         {
      	  for ( int xd=0; xd< NX; xd++ ) /* L-R */
      	  {

      	    h_PosPack [ yd*NZ*NX + zd*NX   + xd ] = make_float3(x,y,z);

      		 if(SimInfo_C->particle_type==polyhedra || SimInfo_C->sphere_orient)
      		 {
      	       h_OrntPack[ yd*NZ*NX + zd*NX   + xd ] = make_quaterion(0,0.0,0.0,0.0);
      		 }

              x = x + x_size + h_PosConfig->space.x;
           }


      	    x = h_PosConfig->start.x;
            z = z + z_size + h_PosConfig->space.z;
          }

          z = h_PosConfig->start.z;
          y = y + y_size + h_PosConfig->space.y; /* next height level*/
       }



    cout<<"INFO-D: Allocating Initial Positions on Device: ";



	float3 *d_posPack;
	cudaMalloc( (void**) &d_posPack , sizeof(float3)*NUMPARTICLES_MAX_SIM );
	cudaDeviceSynchronize();

	cout<<" DONE! \n";

	cout<<"INFO-D: Copying  Host Positions to Device: ";
    cudaMemcpy(d_posPack, h_PosPack, sizeof(float3)*NUMPARTICLES_MAX_SIM, cudaMemcpyHostToDevice );
    cudaDeviceSynchronize();
    cout<<" DONE! \n";

	/* Set the particle type */
	int index_particle = 0;

	/* Now store particle types */
	for( int i=0; i<SimInfo_C->Num_ParticleObjects; i++ )
	{
		for(int j=0; j<SimInfo_C->Num_ParticlesPerObject[i]; j++)
		{
			h_Particle_TypePack[index_particle] = i;
			h_Particle_IDPack  [index_particle] = index_particle;
			index_particle++;
		}

	 }


	 LogFileD<<"\n";

	 /* Copy particle type to device */
	 uint *d_Particle_TypePack;
	 int  *d_Particle_IDPack;

	 LogFileD<<"INFO-D: Allocating Particle Types on Device: ";

	 cudaMalloc( (void**) &d_Particle_TypePack , sizeof(uint)*NUMPARTICLES_MAX_SIM );
	 cudaMalloc( (void**) &d_Particle_IDPack   , sizeof(int)*NUMPARTICLES_MAX_SIM );
	 cudaDeviceSynchronize();

	 if(error_check==1)
	 {
	    	cudaError_t errormsg=cudaGetLastError();
	    	if(errormsg>0)
	        {
	      	 cout<<"Initial Position Default Grid: Malloc ID "<<cudaGetErrorString(errormsg)<<endl;
	    	 exit(1);
	    	}
	 }




	 LogFileD<<" DONE! \n";
	 LogFileD<<"INFO-D: Copying  Particle Types to Device:";

	 cudaMemcpy( d_Particle_TypePack, h_Particle_TypePack, sizeof(uint)*NUMPARTICLES_MAX_SIM, cudaMemcpyHostToDevice );
	 cudaMemcpy( d_Particle_IDPack,   h_Particle_IDPack, sizeof(int)*NUMPARTICLES_MAX_SIM, cudaMemcpyHostToDevice );
	 cudaDeviceSynchronize();
	 LogFileD<<" DONE! \n";
	 LogFileD<<"\n";


	 LogFileD<<"INFO-D: Initializing arrays on Device: ";



	 Quaterion *d_orntPack;
	 if(SimInfo_C->particle_type==polyhedra || SimInfo_C->sphere_orient)
	 {
			cudaMalloc( (void**) &d_orntPack , sizeof(Quaterion)*NUMPARTICLES_MAX_SIM );
			cudaDeviceSynchronize();
		    cudaMemcpy(d_orntPack, h_OrntPack, sizeof(Quaterion)*NUMPARTICLES_MAX_SIM, cudaMemcpyHostToDevice );
		    cudaDeviceSynchronize();

	 }


	 if(error_check==1)
	 {
	    	cudaError_t errormsg=cudaGetLastError();
	    	if(errormsg>0)
	        {
	      	 cout<<"Initial Position Default Grid: Memory End "<<cudaGetErrorString(errormsg)<<endl;
	    	 exit(1);
	    	}
	 }


	 if(!Device_use_multi_gpu)
	 {
	 Set_Pos<<< GlobalLaunch.dimGrid, GlobalLaunch.dimBlock >>>( d_posPack,
			                                                     d_orntPack,
                                                                 d_Particle_TypePack,
                                                                 d_Particle_IDPack,
                                                                 dram_position_com,
			                                                     dram_position_ornt,
			                                                     dram_ObjectType,
			                                                     dram_P_ID            );

	 cudaDeviceSynchronize();
	 }
	 else
	 {
      	 Set_Pos<<< GlobalLaunch.dimGrid, GlobalLaunch.dimBlock >>>( d_posPack,
			                                                     d_orntPack,
                                                                 d_Particle_TypePack,
                                                                 d_Particle_IDPack,

                                                                 dram_position_comM  [0],
			                                                     dram_position_orntM [0],
			                                                     dram_ObjectTypeM    [0],
			                                                     dram_P_IDM          [0]  );

	 cudaDeviceSynchronize();

	 }

	 if(error_check==1)
	 {
	    cudaError_t errormsg=cudaGetLastError();
	    if(errormsg>0)
	    {
	       cout<<"Initial Position Default Grid: Set Pos Kernel "<<cudaGetErrorString(errormsg)<<endl;
	       exit(1);
	    }
	 }




	 LogFileD<<" DONE Freeing Memory ! \n";

     delete [] h_PosPack;

	 if(SimInfo_C->particle_type==polyhedra || SimInfo_C->sphere_orient)
	 {
       delete [] h_OrntPack;
	 }

     delete [] h_Particle_TypePack;
     delete [] h_Particle_IDPack;

	 cudaFree( d_posPack );
     if(SimInfo_C->particle_type==polyhedra|| SimInfo_C->sphere_orient)
	 {
	   cudaFree( d_orntPack );
	 }
	 cudaFree( d_Particle_TypePack );
	 cudaFree( d_Particle_IDPack );

	 cudaDeviceSynchronize();

	 if(error_check==1)
	 {
	    	cudaError_t errormsg=cudaGetLastError();
	    	if(errormsg>0)
	        {
	    	  cout<<"Initial Position Default Grid: Free memory "<<cudaGetErrorString(errormsg)<<endl;
	    	  exit(1);
	    	}
	 }

}
/*-----------------------------------------------------------------------------*/





/*---------------------------------------------------------------------------*/
          /* Set Grid Random Distribution of type and *Orientation  */
/*---------------------------------------------------------------------------*/
void Set_Position_Random_Grid ( InitConfig *h_PosConfig, POBJ *h_POBJ )
{


	 if(error_check==1)
	 {
	    	cudaError_t errormsg=cudaGetLastError();
	    	if(errormsg>0)
	    	{
	    	cout<<"Initial Position Random Grid: Start "<<cudaGetErrorString(errormsg)<<endl;
	    	 exit(1);
	    	}
	 }

	cout<<"INFO-D: Allocating Initial positions on HOST \n";

	h_PosPack             = new float3    [NUMPARTICLES_MAX_SIM];

	 if(SimInfo_C->particle_type==polyhedra || SimInfo_C->sphere_orient)
	 {
	   h_OrntPack            = new Quaterion [NUMPARTICLES_MAX_SIM];
	   cout<<"INFO-D: Polyhedra Setting Orient \n";
	 }

	h_Particle_TypePack   = new uint      [NUMPARTICLES_MAX_SIM];
	h_Particle_IDPack     = new int       [NUMPARTICLES_MAX_SIM];


	cout<<"INFO-D: Creating Random Position Grid  "<<NUMPARTICLES_MAX_SIM<<endl;



	cout<<" \n";

   float x,y,z;
   float x_size,y_size,z_size;


   int    NX  = h_PosConfig->num.x;
   int    NY  = h_PosConfig->num.y;
   int    NZ  = h_PosConfig->num.z;

   x_size = h_PosConfig->p_size[0].x;
   y_size = h_PosConfig->p_size[0].y;
   z_size = h_PosConfig->p_size[0].z;


   x      = h_PosConfig->start.x;
   y      = h_PosConfig->start.y;
   z      = h_PosConfig->start.z;



   cout<<"INFO-D: Start  "<<x<<", "<<y<<", "<<z<<endl;
   cout<<"  Grid  "<<NX<<", "<<NY<<", "<<NZ<<endl;

   /* Allocate positions */
   for ( int yd=0; yd< NY; yd++ )/* Height  */
     {

        for ( int zd=0; zd<NZ; zd++ )/* x B-F */
        {
     	  for ( int xd=0; xd< NX; xd++ ) /* L-R */
     	  {

     	    h_PosPack [ yd*NZ*NX + zd*NX   + xd ] = make_float3(x,y,z);

     		 if(SimInfo_C->particle_type==polyhedra || SimInfo_C->sphere_orient)
     		 {
     	       h_OrntPack[ yd*NZ*NX + zd*NX   + xd ] = make_quaterion(0,0.0,0.0,0.0);
     		 }

             x = x + x_size + h_PosConfig->space.x;
          }


     	    x = h_PosConfig->start.x;
           z = z + z_size + h_PosConfig->space.z;
         }

         z = h_PosConfig->start.z;
         y = y + y_size + h_PosConfig->space.y; /* next height level*/
      }



   cout<<"INFO-D: Allocating Initial Positions on Device: ";



	float3 *d_posPack;
	cudaMalloc( (void**) &d_posPack , sizeof(float3)*NUMPARTICLES_MAX_SIM );
	cudaDeviceSynchronize();

	cout<<" DONE! \n";

	cout<<"INFO-D: Copying  Host Positions to Device: ";
   cudaMemcpy(d_posPack, h_PosPack, sizeof(float3)*NUMPARTICLES_MAX_SIM, cudaMemcpyHostToDevice );
   cudaDeviceSynchronize();
   cout<<" DONE! \n";

	/* Set the particle type */
	int index_particle = 0;

	/* Now store particle types */
	for( int i=0; i<SimInfo_C->Num_ParticleObjects; i++ )
	{
		for(int j=0; j<SimInfo_C->Num_ParticlesPerObject[i]; j++)
		{
			h_Particle_TypePack[index_particle] = i;
			h_Particle_IDPack  [index_particle] = index_particle;
			index_particle++;
		}

	 }


    /* Now randomize grid */

		int source;
		int dest;

		int temp;

		for(int i=0; i<NUMPARTICLES_MAX_SIM;i ++)
		{
			source = rand()%(NUMPARTICLES_MAX_SIM-1);
			dest   = rand()%(NUMPARTICLES_MAX_SIM-1);

			/* Sort P_Type */
			temp                        = h_Particle_TypePack[source];
			h_Particle_TypePack[source] = h_Particle_TypePack[dest];
			h_Particle_TypePack[dest]   = temp;

			temp                      = h_Particle_IDPack[source];
			h_Particle_IDPack[source] = h_Particle_IDPack[dest];
			h_Particle_IDPack[dest]   = temp;

		}



	 LogFileD<<"\n";

	 /* Copy particle type to device */
	 uint *d_Particle_TypePack;
	 int  *d_Particle_IDPack;

	 LogFileD<<"INFO-D: Allocating Particle Types on Device: ";

	 cudaMalloc( (void**) &d_Particle_TypePack , sizeof(uint)*NUMPARTICLES_MAX_SIM );
	 cudaMalloc( (void**) &d_Particle_IDPack   , sizeof(int)*NUMPARTICLES_MAX_SIM );
	 cudaDeviceSynchronize();

	 if(error_check==1)
	 {
	    	cudaError_t errormsg=cudaGetLastError();
	    	if(errormsg>0)
	        {
	      	 cout<<"Initial Position random Grid: Malloc ID "<<cudaGetErrorString(errormsg)<<endl;
	    	 exit(1);
	    	}
	 }




	 LogFileD<<" DONE! \n";
	 LogFileD<<"INFO-D: Copying  Particle Types to Device:";

	 cudaMemcpy( d_Particle_TypePack, h_Particle_TypePack, sizeof(uint)*NUMPARTICLES_MAX_SIM, cudaMemcpyHostToDevice );
	 cudaMemcpy( d_Particle_IDPack,   h_Particle_IDPack, sizeof(int)*NUMPARTICLES_MAX_SIM, cudaMemcpyHostToDevice );
	 cudaDeviceSynchronize();
	 LogFileD<<" DONE! \n";
	 LogFileD<<"\n";


	 LogFileD<<"INFO-D: Initializing arrays on Device: ";



	 Quaterion *d_orntPack;
	 if(SimInfo_C->particle_type==polyhedra || SimInfo_C->sphere_orient)
	 {
			cudaMalloc( (void**) &d_orntPack , sizeof(Quaterion)*NUMPARTICLES_MAX_SIM );
			cudaDeviceSynchronize();
		    cudaMemcpy(d_orntPack, h_OrntPack, sizeof(Quaterion)*NUMPARTICLES_MAX_SIM, cudaMemcpyHostToDevice );
		    cudaDeviceSynchronize();

	 }


	 if(error_check==1)
	 {
	    	cudaError_t errormsg=cudaGetLastError();
	    	if(errormsg>0)
	        {
	      	 cout<<"Initial Position Default Grid: Memory End "<<cudaGetErrorString(errormsg)<<endl;
	    	 exit(1);
	    	}
	 }


	 if(!Device_use_multi_gpu)
	 {
	 Set_Pos<<< GlobalLaunch.dimGrid, GlobalLaunch.dimBlock >>>( d_posPack,
			                                                     d_orntPack,
                                                                d_Particle_TypePack,
                                                                d_Particle_IDPack,
                                                                dram_position_com,
			                                                     dram_position_ornt,
			                                                     dram_ObjectType,
			                                                     dram_P_ID            );

	 cudaDeviceSynchronize();
	 }
	 else
	 {
     	 Set_Pos<<< GlobalLaunch.dimGrid, GlobalLaunch.dimBlock >>>( d_posPack,
			                                                     d_orntPack,
                                                                d_Particle_TypePack,
                                                                d_Particle_IDPack,

                                                                 dram_position_comM  [0],
			                                                     dram_position_orntM [0],
			                                                     dram_ObjectTypeM    [0],
			                                                     dram_P_IDM          [0]  );

	 cudaDeviceSynchronize();

	 }

	 if(error_check==1)
	 {
	    cudaError_t errormsg=cudaGetLastError();
	    if(errormsg>0)
	    {
	       cout<<"Initial Position Random Grid: Set Pos Kernel "<<cudaGetErrorString(errormsg)<<endl;
	       exit(1);
	    }
	 }




	 LogFileD<<" DONE Freeing Memory ! \n";

    delete [] h_PosPack;

	 if(SimInfo_C->particle_type==polyhedra || SimInfo_C->sphere_orient)
	 {
      delete [] h_OrntPack;
	 }

    delete [] h_Particle_TypePack;
    delete [] h_Particle_IDPack;

	 cudaFree( d_posPack );
    if(SimInfo_C->particle_type==polyhedra|| SimInfo_C->sphere_orient)
	 {
	   cudaFree( d_orntPack );
	 }
	 cudaFree( d_Particle_TypePack );
	 cudaFree( d_Particle_IDPack );

	 cudaDeviceSynchronize();

	 if(error_check==1)
	 {
	    	cudaError_t errormsg=cudaGetLastError();
	    	if(errormsg>0)
	        {
	    	  cout<<"Initial Position Random Grid: Free memory "<<cudaGetErrorString(errormsg)<<endl;
	    	  exit(1);
	    	}
	 }



}
/*-----------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
          /* Fill from half-way up   */
/*---------------------------------------------------------------------------*/
void Set_Position_Fill_Rectangle ( POBJ *h_POBJ,WOBJ *h_WOBJ, float vel, 
							       float3 fill_planeStart, float3 fill_planeEnd )
{

	cout<<"INFO-D: Allocating Fill positions on HOST \n";

	h_PosPack             = new float3    [NUMPARTICLES_MAX_SIM];

	if(SimInfo_C->particle_type==polyhedra || SimInfo_C->sphere_orient)
	{
	  h_OrntPack            = new Quaterion [NUMPARTICLES_MAX_SIM];
	}

	h_Particle_TypePack   = new uint      [NUMPARTICLES_MAX_SIM];
	h_Particle_IDPack     = new int       [NUMPARTICLES_MAX_SIM];

	 if(error_check==1)
	 {
	    	cudaError_t errormsg=cudaGetLastError();
	    	if(errormsg>0)
	    	{
	    	  cout<<"Initial Position Fill: Start "<<cudaGetErrorString(errormsg)<<endl;
	    	 exit(1);
	    	}
	 }

	cout<<"INFO-D: Creating Fill Config \n";

	cout<<" \n";

    float x,y,z;
    float x_size,y_size,z_size;

    float3 space;
    float3 start;
    space = make_float3(0.00100f,0.0010f,0.0000f);

    /* Biggest particle is first */


      x_size = 2.0*h_POBJ[0].radius;
      y_size = 2.0*h_POBJ[0].radius;
      z_size = 2.0*h_POBJ[0].radius;


    /* For mill box width is diameter */
    int    NX,NY,NZ;

    float radius;
    if(NUMDYOBJECTS>0)
    {
     radius = SimInfo_C->VolumeObject_CylinderBoundR;
    }
    else
    {
    	radius = h_WOBJ[0].cyl_geo.radius;
    }


    if(SimInfo_C->Simulation_Type==ballmill && !SimInfo_C->Rotating_WObject)
    {
      printf("Mill detected\n");

      if(NUMDYOBJECTS>0)
      {
        x      =  x_size + (h_WOBJ[0].cyl_geo.radius - SimInfo_C->VolumeObject_CylinderBoundR);
      }
      else
      {
    	  x    =  x_size + space.x;
      }

      y      = h_WOBJ[0].cyl_geo.radius;

      z      =  h_WOBJ[0].cyl_geo.center_bot_cap.z
    		    + z_size/2.0f + space.z;

      NX  = floor( ( radius*2 - x_size )/(x_size + space.x) );

      NZ  = floor((h_WOBJ[0].cyl_geo.height)/(z_size + space.z));

      NY  = ceil( NUMPARTICLES_MAX_SIM/(float)(NX*NZ) );

    }
    else
    {

      printf("Planes Detected \n");
      x      = fill_planeStart.x
    		    + x_size;

      y      = fill_planeStart.y;

      z      = fill_planeStart.z
    		     + z_size;

      NX  = floor((fill_planeEnd.x - fill_planeStart.x)/(x_size + space.x));
      NZ  = floor((fill_planeEnd.z - fill_planeStart.z) /(z_size + space.z));

      NY  = ceilf( NUMPARTICLES_MAX_SIM/(float)(NX*NZ) );

    }

    start = make_float3(x,y,z);

    pack_pos = start;


    printf("Start %f %f %f  \n",x,y,z);

    printf("Grid Size %d %d %d Tot %d   SIMMAX %d \n",NX,NY,NZ, NX*NY*NZ,NUMPARTICLES_MAX_SIM);

    pack_size = make_int3(NX,NY,NZ);

    packNumber = NX*NZ;

    packNumWaves= NY;



	pack_clear_steps = fabsl((int)ceil((2.0f*y_size/(SimInfo_C->InitalDelta_t*vel))));
    printf("clear steps %d \n",pack_clear_steps);

    int selectcount[6];

    /* Get prob for each type */
 	for(int i=0; i<SimInfo_C->Num_ParticleObjects; i++)
 	{
 	  selectcount[i] = 0;
 	}


 	int selectT=-1;
 	int pty=0;



    for ( int yd=0; yd< NY; yd++ )/* Height  */
    {

       for ( int zd=0; zd<NZ; zd++ )/* x B-F */
       {
    	  for ( int xd=0; xd< NX; xd++ ) /* L-R */
    	  {
    		if(yd*NZ*NX + zd*NX   + xd <NUMPARTICLES_MAX_SIM)
    		{
    	      h_PosPack[ yd*NZ*NX + zd*NX   + xd ]= make_float3(x,y,z);


//    	 	 if(SimInfo_C->particle_type==polyhedra || SimInfo_C->sphere_orient)
//    	 	 {	      /* Random Orientation */
//      	      float rand_w = ((float)rand()/RAND_MAX)*30.0;
//      	      float rand_x = (float)rand()/RAND_MAX;
//      	      float rand_y = (float)rand()/RAND_MAX;
//      	      float rand_z = (float)rand()/RAND_MAX;
//      	      h_OrntPack[ yd*NZ*NX + zd*NX   + xd ]= make_quaterion(rand_w,rand_x,rand_y,rand_z);//(90.0,0.0,1.0,0.0);//
//    	     }

    	       if(SimInfo_C->particle_type==polyhedra || SimInfo_C->sphere_orient)
    	     {
    	      h_OrntPack[ yd*NZ*NX + zd*NX   + xd ]= make_quaterion(0.0,0.0,0.0,0.0);
    	     }
      	    selectT++;

    	      /* Reset after each alternating */
    	      if(selectT >= SimInfo_C->Num_ParticleObjects )
    	      {
    	    	  selectT=0;
   	          }

    	      pty=selectT;


              /* Make sure we dont exceed count */
    	      if( selectcount[pty]>=SimInfo_C->Num_ParticlesPerObject[pty] )
    	      {
    	       /* Find someone who is free */
    	    	 	for(int i=0; i<SimInfo_C->Num_ParticleObjects; i++)
    	    	 	{
    	      	      if( selectcount[i]<=(SimInfo_C->Num_ParticlesPerObject[i]-1) )
    	      	      {
    	      	    	pty=i;
    	      	    	break;
    	      	      }
    	    	 	}
    	      }

    	    selectcount[pty]++;


    	    h_Particle_TypePack[yd*NZ*NX + zd*NX   + xd] = pty;
    	    h_Particle_IDPack[yd*NZ*NX + zd*NX   + xd] = yd*NZ*NX + zd*NX   + xd;

            x = x + x_size + space.x;
          }
    	  else
    	  {
    	     break;
    	   }
       }


    	  x = start.x;
          z = z + z_size + space.z;
        }

        z = start.z;

       //y = y + y_size + space.y; /* next height level*/
     }


    uint  *d_Particle_TypePack;
    int   *d_Particle_IDPack;
     /* Now randomize grid */

		int source;
		int dest;

		int temp;

		for(int i=0; i<NUMPARTICLES_MAX_SIM;i ++)
		{
			source = rand()%(NUMPARTICLES_MAX_SIM-1);
			dest   = rand()%(NUMPARTICLES_MAX_SIM-1);

			temp                        = h_Particle_TypePack[source];
			h_Particle_TypePack[source] = h_Particle_TypePack[dest];
			h_Particle_TypePack[dest]   = temp;

			temp                        = h_Particle_IDPack[source];
			h_Particle_IDPack[source]   = h_Particle_IDPack[dest];
			h_Particle_IDPack[dest]     = temp;
		}



    cout<<"INFO-D: Allocating Initial Positions on Device: ";

    float3 *d_posPack;
	cudaMalloc( (void**) &d_posPack , sizeof(float3)*NUMPARTICLES_MAX_SIM );
	cudaDeviceSynchronize();


	cout<<" DONE! \n";


	cout<<"INFO-D: Copying  Positions to Device: ";
    cudaMemcpy(d_posPack, h_PosPack, sizeof(float3)*NUMPARTICLES_MAX_SIM, cudaMemcpyHostToDevice );
    cudaDeviceSynchronize();
    cout<<" DONE! \n";


	LogFileD<<"\n";
	 /* Copy particle type to device */


	 LogFileD<<"INFO-D: Allocating Particle Types on Device: ";

	 cudaMalloc( (void**) &d_Particle_TypePack , sizeof(uint)*NUMPARTICLES_MAX_SIM );
	 cudaMalloc( (void**) &d_Particle_IDPack , sizeof(int)*NUMPARTICLES_MAX_SIM );
	 cudaDeviceSynchronize();
	 LogFileD<<" DONE! \n";
	 LogFileD<<"INFO-D: Copying  Particle Types to Device:";

	 cudaMemcpy( d_Particle_TypePack, h_Particle_TypePack, sizeof(uint)*NUMPARTICLES_MAX_SIM, cudaMemcpyHostToDevice );
	 cudaMemcpy( d_Particle_IDPack,   h_Particle_IDPack, sizeof(int)*NUMPARTICLES_MAX_SIM, cudaMemcpyHostToDevice );
	 cudaDeviceSynchronize();

	 Quaterion *d_orntPack;
	 if(SimInfo_C->particle_type==polyhedra|| SimInfo_C->sphere_orient)
	 {

			cudaMalloc( (void**) &d_orntPack , sizeof(Quaterion)*NUMPARTICLES_MAX_SIM );
			cudaDeviceSynchronize();
		    cudaMemcpy(d_orntPack, h_OrntPack, sizeof(Quaterion)*NUMPARTICLES_MAX_SIM, cudaMemcpyHostToDevice );
		    cudaDeviceSynchronize();

	 }

	 LogFileD<<" DONE! \n";
	 LogFileD<<"\n";


	 LogFileD<<"INFO-D: Initializing arrays on Device: ";
	 cudaDeviceSynchronize();


	 if(error_check==1)
	 {
	    	cudaError_t errormsg=cudaGetLastError();
	    	if(errormsg>0)
	    	{
	    	cout<<"Initial Position Fill: Mem End "<<cudaGetErrorString(errormsg)<<endl;
	    	 exit(1);
	    	}
	 }


	 if(!Device_use_multi_gpu)
	 {
	 Set_Pos<<< GlobalLaunch.dimGrid, GlobalLaunch.dimBlock >>>( d_posPack,
			                                                     d_orntPack,
                                                                 d_Particle_TypePack,
                                                                 d_Particle_IDPack,
                                                                 dram_position_com,
			                                                     dram_position_ornt,
			                                                     dram_ObjectType,
			                                                     dram_P_ID            );

	 cudaDeviceSynchronize();
	 }
	 else
	 {
      	 Set_Pos<<< GlobalLaunch.dimGrid, GlobalLaunch.dimBlock >>>( d_posPack,
			                                                     d_orntPack,
                                                                 d_Particle_TypePack,
                                                                 d_Particle_IDPack,
																 dram_position_comM[0],
			                                                     dram_position_orntM[0],
			                                                     dram_ObjectTypeM[0],
			                                                     dram_P_IDM[0]            );

	 cudaDeviceSynchronize();

	 }



	 if(error_check==1)
	 {
	    	cudaError_t errormsg=cudaGetLastError();
	    	if(errormsg>0)
	    	{
	    	cout<<"Initial Position Fill: Set_Pos Kernel "<<cudaGetErrorString(errormsg)<<endl;
	    	 exit(1);
	    	}
	 }

     cout<<" DONE Freeing Memory ! \n";

     NUMPARTICLES_Current = 0;

     delete [] h_PosPack;

	 if(SimInfo_C->particle_type==1 || SimInfo_C->sphere_orient)
	 {
       delete [] h_OrntPack;
	 }
     delete [] h_Particle_TypePack;
     delete [] h_Particle_IDPack;

	 cudaFree( d_posPack );
	 if(SimInfo_C->particle_type==1 || SimInfo_C->sphere_orient)
	 {
	  cudaFree( d_orntPack );
	 }

	 cudaFree( d_Particle_TypePack );
	 cudaFree( d_Particle_IDPack );

	 if(error_check==1)
	 {
	    	cudaError_t errormsg=cudaGetLastError();
	    	if(errormsg>0)
	    	{
	    	cout<<"Initial Position Fill: Mem Free "<<cudaGetErrorString(errormsg)<<endl;
	    	 exit(1);
	    	}
	 }
}
/*-----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------*/
    /* Adds particles to the simulation until NUMPARTICLES_Max is reached */
/*-----------------------------------------------------------------------------*/
void Fill_State()
{
	int num_rem      = NUMPARTICLES_MAX_SIM - NUMPARTICLES_Current;
	int num_particle = min(pack_size.x*pack_size.z,num_rem);


	 if(error_check==1)
	 {
	    	cudaError_t errormsg=cudaGetLastError();
	    	if(errormsg>0)
	    	{
	    	 cout<<" Fill State: Start "<<cudaGetErrorString(errormsg)<<endl;
	    	 exit(1);
	    	}
	 }



	   NUMPARTICLES_Current += num_particle;


	   int num_particles_SM = (int)ceil(NUMPARTICLES_Current/(float)num_sm);

	   /* Now get the number of blocks per SM */
	   int num_blocks_SM =  (int)ceil(num_particles_SM/(float)threads_perBlock);

	   int num_threads_block=threads_perBlock;

	   /* block size is too big */
	   if(num_particles_SM < threads_perBlock)
	   {
		   num_threads_block = num_particles_SM; /* Single block is sufficient */
	   }

	   GlobalLaunch.dimBlock = make_uint3( num_threads_block,1,1);
	   GlobalLaunch.dimGrid  = make_uint3( num_blocks_SM*num_sm,1,1 );



	    SimInfo_C->Num_Particles = NUMPARTICLES_Current;
	    
		if(!Device_use_multi_gpu)
	    {
		  cudaMemcpyToSymbol( SimParms,       SimInfo_C,     sizeof( SimulationInfo ) );
	      cudaDeviceSynchronize();
		}
		else
		{
          for ( int i=0; i < num_devices_use; i++ )
		  {
			cudaSetDevice(i);
            cudaMemcpyToSymbol( SimParms,       SimInfo_C,     sizeof( SimulationInfo ) );
	        cudaDeviceSynchronize();
		  }
		  cudaSetDevice(0);
		}



	/* Only modify subset */
	    if(!Device_use_multi_gpu)
	    {
	       Fill_Plane<<< GlobalLaunch.dimGrid, GlobalLaunch.dimBlock >>>(pack_pos.y,dram_position_com,dram_velocity_com);
	       cudaDeviceSynchronize();
	    }
	    else
	    {
	    	   int num_perGPU = (int)floor(NUMPARTICLES_Current/(float)num_devices_use);

			   int num_odd = NUMPARTICLES_Current - num_perGPU*num_devices_use;

				for ( int i=0; i < num_devices_use; i++ )
				{
					num_particles_perGPU[i] = num_perGPU;

					if(i==0)
					{
							num_particles_perGPU[0] = num_perGPU + num_odd;

					}


				num_particles_SM = (int)ceil(num_particles_perGPU[i]/(float)num_sm);

				/* Now get the number of blocks per SM */
				num_blocks_SM      =  (int)ceil(num_particles_SM/(float)threads_perBlock);
				num_threads_block = threads_perBlock ;

				/* block size is too big */
				if( num_particles_SM < threads_perBlock )
				{
				 num_threads_block = num_particles_SM; /* Single block is sufficient */
				}

				multi[i].dimBlock = make_uint3( num_threads_block,1,1);
				multi[i].dimGrid  = make_uint3( num_blocks_SM*num_sm,1,1 );

				//printf(" GPU %d Particles %d per GPU threads = %d \n",i, num_particles_perGPU[i],multi[i].dimGrid.x*multi[i].dimBlock.x);

				}/* end loop over device */

	    	    Fill_Plane<<< GlobalLaunch.dimGrid, GlobalLaunch.dimBlock >>>(pack_pos.y,dram_position_comM[0],dram_velocity_comM[0]);
	    		cudaDeviceSynchronize();
	    }

	 if(error_check==1)
	 {
	    	cudaError_t errormsg=cudaGetLastError();
	    	if(errormsg>0)
	    	{
	    	 cout<<"Fill State: Fill Kernel "<<cudaGetErrorString(errormsg)<<endl;
	    	 exit(1);
	    	}
	 }

}
/*-----------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------*/
            /* Set the positions only on device */
/*-----------------------------------------------------------------------------*/
void Device_Set_Pos( float3 *h_pos,  Quaterion *h_Ornt,
		             uint *h_ptype,   int       *h_pid /*TODO add number to set */  )
{

    LogFileD<<"INFO-D: Allocating Only File Positions on Device: ";

	float3    *d_pos;
	Quaterion *d_Ornt;
	uint      *d_ptype;
	int       *d_pid;


	cudaMalloc( (void**) &d_pos , sizeof(float3)*NUMPARTICLES_MAX_SIM );
	cudaMalloc( (void**) &d_ptype, sizeof(uint) *NUMPARTICLES_MAX_SIM  );
	cudaMalloc( (void**) &d_pid,   sizeof(int)*NUMPARTICLES_MAX_SIM    );
	cudaDeviceSynchronize();


	LogFileD<<" DONE! \n";
	LogFileD<<"INFO-D: Copying  Positions to Device: ";

	cudaMemcpy(d_pos, h_pos, sizeof(float3)*NUMPARTICLES_MAX_SIM, cudaMemcpyHostToDevice );

	if(SimInfo_C->particle_type==polyhedra|| SimInfo_C->sphere_orient)
	{
		cudaMalloc( (void**) &d_Ornt , sizeof(Quaterion)*NUMPARTICLES_MAX_SIM );
		cudaDeviceSynchronize();
		cudaMemcpy(d_Ornt, h_Ornt, sizeof(Quaterion)*NUMPARTICLES_MAX_SIM, cudaMemcpyHostToDevice );
	}

    cudaMemcpy(d_ptype, h_ptype, sizeof(uint)*NUMPARTICLES_MAX_SIM,  cudaMemcpyHostToDevice );
    cudaMemcpy(d_pid,   h_pid,  sizeof(int)*NUMPARTICLES_MAX_SIM,  cudaMemcpyHostToDevice );
    cudaDeviceSynchronize();

	 if(error_check==1)
	 {
	    	cudaError_t errormsg=cudaGetLastError();
	    	if(errormsg>0)
	    	{
	    	cout<<"Set Position: Mem  End "<<cudaGetErrorString(errormsg)<<endl;
	    	 exit(1);
	    	}
	 }



    LogFileD<<" DONE! \n";

    if(!Device_use_multi_gpu)
    {
      Set_Pos<<< GlobalLaunch.dimGrid, GlobalLaunch.dimBlock >>>( d_pos,d_Ornt, d_ptype, d_pid,
    		                                                    dram_position_com, dram_position_ornt,
    		                                                     dram_ObjectType, dram_P_ID   );
	  cudaDeviceSynchronize();
    }
    else
    {
    	Set_Pos<<< GlobalLaunch.dimGrid, GlobalLaunch.dimBlock >>>( d_pos,d_Ornt, d_ptype, d_pid,
    	    		                                                    dram_position_comM[0], dram_position_orntM[0],
    	    		                                                     dram_ObjectTypeM[0], dram_P_IDM[0]   );
        cudaDeviceSynchronize();
    }

	 if(error_check==1)
	 {
	    	cudaError_t errormsg=cudaGetLastError();
	    	if(errormsg>0)
	    	{
	    	cout<<"Set Position: Set_Pos Kernel "<<cudaGetErrorString(errormsg)<<endl;
	    	 exit(1);
	    	}
	 }


	cudaFree(d_pos);
	if(SimInfo_C->particle_type==polyhedra || SimInfo_C->sphere_orient)
	{
	  cudaFree(d_Ornt);
	}

	cudaFree(d_ptype);
	cudaFree(d_pid);

	 if(error_check==1)
	 {
	    	cudaError_t errormsg=cudaGetLastError();
	    	if(errormsg>0)
	    	{
	    	cout<<"Set Position: Mem Free "<<cudaGetErrorString(errormsg)<<endl;
	    	 exit(1);
	    	}
	 }


	LogFileD<<" DONE Freeing Memory ! \n";

}
/*-----------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
                 /* Set SPEC TYPE:2  Pos and Velocity */
/*---------------------------------------------------------------------------*/
void Device_Set_Pos_Velocity( float3 *h_pos, Quaterion *h_Ornt, float3 *h_vel,
		                      float3 *h_Avel, uint *h_ptype,  int     *h_pid)
{

    LogFileD<<"INFO-D: Allocating File Pos Vel on Device: ";


	if(error_check==1)
	 {
	    	cudaError_t errormsg=cudaGetLastError();
	    	if(errormsg>0)
	    	{
	    	cout<<"Set Position Velocity: Start "<<cudaGetErrorString(errormsg)<<endl;
	    	 exit(1);
	    	}
	 }


	float3    *d_pos;
	float3    *d_vel;
	float3    *d_Avel;
	Quaterion *d_Ornt;
	uint      *d_ptype;
	int       *d_pid;


	cudaMalloc( (void**) &d_pos , sizeof(float3)*NUMPARTICLES_MAX_SIM );

	if(SimInfo_C->particle_type==polyhedra|| SimInfo_C->sphere_orient)
	{
		cudaMalloc( (void**) &d_Ornt , sizeof(Quaterion)*NUMPARTICLES_MAX_SIM );
	}

	cudaMalloc( (void**) &d_vel , sizeof(float3)*NUMPARTICLES_MAX_SIM );
	cudaMalloc( (void**) &d_Avel, sizeof(float3)*NUMPARTICLES_MAX_SIM );
	cudaMalloc( (void**) &d_ptype, sizeof(uint) *NUMPARTICLES_MAX_SIM  );
	cudaMalloc( (void**) &d_pid,   sizeof(int)*NUMPARTICLES_MAX_SIM    );
	cudaDeviceSynchronize();


	LogFileD<<" DONE! \n";
	LogFileD<<"INFO-D: Copying  Positions to Device: ";

	cudaMemcpy(d_pos, h_pos, sizeof(float3)*NUMPARTICLES_MAX_SIM, cudaMemcpyHostToDevice );

	if(SimInfo_C->particle_type==polyhedra)
	{
		cudaMemcpy(d_Ornt, h_Ornt, sizeof(Quaterion)*NUMPARTICLES_MAX_SIM, cudaMemcpyHostToDevice );
	}

    cudaMemcpy(d_vel,  h_vel,   sizeof(float3)*NUMPARTICLES_MAX_SIM, cudaMemcpyHostToDevice );
    cudaMemcpy(d_Avel, h_Avel,  sizeof(float3)*NUMPARTICLES_MAX_SIM, cudaMemcpyHostToDevice );
    cudaMemcpy(d_ptype, h_ptype, sizeof(uint)*NUMPARTICLES_MAX_SIM,  cudaMemcpyHostToDevice );
    cudaMemcpy(d_pid,   h_pid,   sizeof(int)*NUMPARTICLES_MAX_SIM,  cudaMemcpyHostToDevice );
    cudaDeviceSynchronize();

    LogFileD<<" DONE! \n";


    if(!Device_use_multi_gpu)
    {
      Set_Pos_Vel<<< GlobalLaunch.dimGrid, GlobalLaunch.dimBlock >>>( d_pos,d_Ornt,d_vel,d_Avel,
    		                                                        d_ptype,d_pid,
    		                                                        dram_position_com, dram_position_ornt,
    		                                                         dram_velocity_com, dram_velocity_ang, dram_ObjectType, dram_P_ID   );
	  cudaDeviceSynchronize();
    }
    else
    {
      Set_Pos_Vel<<< GlobalLaunch.dimGrid, GlobalLaunch.dimBlock >>>( d_pos,d_Ornt,d_vel,d_Avel,
      		                                                        d_ptype,d_pid,
      		                                                        dram_position_comM[0], dram_position_orntM[0],
      		                                                         dram_velocity_comM[0], dram_velocity_angM[0], 
																	 dram_ObjectTypeM[0], dram_P_IDM[0]   );
  	  cudaDeviceSynchronize();

    }


	if(error_check==1)
	 {
	    	cudaError_t errormsg=cudaGetLastError();
	    	if(errormsg>0)
	    	{
	    	cout<<"Set Position Velocity: Set_Pos_Vel Kernel "<<cudaGetErrorString(errormsg)<<endl;
	    	 exit(1);
	    	}
	 }



	cudaFree(d_pos);
	cudaFree(d_vel);
	cudaFree(d_Avel);

	if(SimInfo_C->particle_type==polyhedra || SimInfo_C->sphere_orient)
	{
	  cudaFree(d_Ornt);
	}

	cudaFree(d_ptype);
	cudaFree(d_pid);

	LogFileD<<" DONE Freeing Memory ! \n";

}
/*-----------------------------------------------------------------------------*/




/*---------------------------------------------------------------------------*/
          /* Spec type 3 all particle properties */
/*---------------------------------------------------------------------------*/
void Device_Set_System_State( float3    *h_pos,
		                      Quaterion *h_Ornt,
		                      float3 *h_vel,
		                      float3 *h_Avel,
		                      float3 *h_acc,
		                      uint  *h_ptype,
		                      int      *h_pid)
{

    cout<<"INFO-D: Allocating File State on Device: ";

	float3    *d_pos;
	float3    *d_vel;
	Quaterion *d_Ornt;
	float3    *d_Avel;
	uint      *d_ptype;
	int       *d_pid;

	cudaMalloc( (void**) &d_pos  , sizeof(float3)   *NUMPARTICLES_MAX_SIM );

	if(SimInfo_C->particle_type==polyhedra || SimInfo_C->sphere_orient)
    {
      cudaMalloc( (void**) &d_Ornt , sizeof(Quaterion)*NUMPARTICLES_MAX_SIM );
    }

	cudaMalloc( (void**) &d_vel  , sizeof(float3)   *NUMPARTICLES_MAX_SIM );
	cudaMalloc( (void**) &d_Avel , sizeof(float3)   *NUMPARTICLES_MAX_SIM );
	cudaMalloc( (void**) &d_ptype , sizeof(uint)     *NUMPARTICLES_MAX_SIM );
	cudaMalloc( (void**) &d_pid   , sizeof(int)     *NUMPARTICLES_MAX_SIM );
	cudaDeviceSynchronize();

	 if(error_check==1)
	 {
	    	cudaError_t errormsg=cudaGetLastError();
	    	if(errormsg>0)
	    	{
	    	cout<<"Set System State: Start "<<cudaGetErrorString(errormsg)<<endl;
	    	 exit(1);
	    	}
	 }


	cout<<" DONE! \n";

	cout<<"INFO-D: Copying  State to Device: ";

	cudaMemcpy( d_pos,  h_pos,  sizeof(float3)*NUMPARTICLES_MAX_SIM, cudaMemcpyHostToDevice );

	if(SimInfo_C->particle_type==polyhedra || SimInfo_C->sphere_orient)
    {
      cudaMemcpy(d_Ornt, h_Ornt, sizeof(Quaterion)*NUMPARTICLES_MAX_SIM, cudaMemcpyHostToDevice );
    }
    cudaMemcpy( d_vel,  h_vel,  sizeof(float3)*NUMPARTICLES_MAX_SIM, cudaMemcpyHostToDevice );
    cudaMemcpy( d_Avel, h_Avel, sizeof(float3)*NUMPARTICLES_MAX_SIM, cudaMemcpyHostToDevice );
    cudaMemcpy( d_ptype, h_ptype, sizeof(uint)*NUMPARTICLES_MAX_SIM,   cudaMemcpyHostToDevice );
    cudaMemcpy( d_pid, h_pid, sizeof(int)*NUMPARTICLES_MAX_SIM,   cudaMemcpyHostToDevice );
    cudaDeviceSynchronize();
    cout<<" DONE! \n";

	if ( !Device_use_multi_gpu)
	{
	Set_PSystem_State<<< GlobalLaunch.dimGrid, GlobalLaunch.dimBlock >>>( d_pos,d_Ornt,
			                                                             d_vel,
			                                                            d_Avel,
			                                                            d_ptype,d_pid,
			                                                            dram_position_com,
			                                                            dram_position_ornt,
			                                                            dram_velocity_com,
			                                                            dram_velocity_ang,
					                        		                	dram_ObjectType,
					                        		                	dram_P_ID);
	cudaDeviceSynchronize();
	}
	else
	{
     	Set_PSystem_State<<< GlobalLaunch.dimGrid, GlobalLaunch.dimBlock >>>( d_pos,d_Ornt,
			                                                             d_vel,
			                                                            d_Avel,
			                                                            d_ptype,d_pid,
			                                                            dram_position_comM[0],
			                                                            dram_position_orntM[0],
			                                                            dram_velocity_comM[0],
			                                                            dram_velocity_angM[0],
					                        		                	dram_ObjectTypeM[0],
					                        		                	dram_P_IDM[0]);

     	cudaDeviceSynchronize();
	}



	 if(error_check==1)
	 {
	    	cudaError_t errormsg=cudaGetLastError();
	    	if(errormsg>0)
	    	{
	    	cout<<"Set System State: Set_PSystem_State Kernel "<<cudaGetErrorString(errormsg)<<endl;
	    	 exit(1);
	    	}
	 }



	 cudaFree(d_pos);

	 if(SimInfo_C->particle_type==polyhedra|| SimInfo_C->sphere_orient)
	 {
	   cudaFree(d_Ornt);
	 }
	 cudaFree(d_vel);

	 cudaFree(d_Avel);
	 cudaFree(d_ptype);
	 cudaFree(d_pid);

	 cout<<" DONE Freeing Memory ! \n";

	 if(error_check==1)
	 {
	    	cudaError_t errormsg=cudaGetLastError();
	    	if(errormsg>0)
	    	{
	    	cout<<"Set System State: End MemFree "<<cudaGetErrorString(errormsg)<<endl;
	    	 exit(1);
	    	}
	 }

}
/*-----------------------------------------------------------------------------*/

void Device_Set_Tallys(float *h_Tallys)
{

    if(!Device_use_multi_gpu)
    {
	  float *d_Tallys;

      cudaMalloc( (void**) &d_Tallys, sizeof(float)*6 );
      cudaDeviceSynchronize();
      cudaMemcpy( d_Tallys, h_Tallys, sizeof(float)*6, cudaMemcpyHostToDevice );
      cudaDeviceSynchronize();


	  Set_Tallys<<< 1,1 >>>(d_Tallys);
	  cudaDeviceSynchronize();

	  if(error_check==1)
	  {
	    	cudaError_t errormsg=cudaGetLastError();
	    	if(errormsg>0)
	    	{
	    	cout<<"Set Tallies: Set_Tallys Kernel "<<cudaGetErrorString(errormsg)<<endl;
	    	 exit(1);
	    	}
	  }

	 cudaFree(d_Tallys);
    }
    else
    {
    	for (int i=0;i<num_devices_use;i++)
    	{
    		cudaSetDevice(i);

    		  float *d_Tallys;

    	      cudaMalloc( (void**) &d_Tallys, sizeof(float)*6 );
    	      cudaDeviceSynchronize();
    	      cudaMemcpy( d_Tallys, h_Tallys, sizeof(float)*6, cudaMemcpyHostToDevice );
    	      cudaDeviceSynchronize();


    		  Set_Tallys<<< 1,1 >>>(d_Tallys);
    		  cudaDeviceSynchronize();

    		  if(error_check==1)
    		  {
    		    	cudaError_t errormsg=cudaGetLastError();
    		    	if(errormsg>0)
    		    	{
    		    	cout<<"Set Tallies: Set_Tallys Kernel "<<cudaGetErrorString(errormsg)<<endl;
    		    	 exit(1);
    		    	}
    		  }

    		 cudaFree(d_Tallys);
    	}
    	cudaSetDevice(0);
    }

}





/*---------------------------------------------------------------------------*/
               /* 1. Get the properties of all CUDA DEVICES */
/*---------------------------------------------------------------------------*/
void Get_DeviceInfo()
{


    cudaGetDeviceCount(&num_devices);

    cout<<"INFO-D: Total Devices: "<<num_devices<<endl;

    for( int i=0;i<num_devices; i++ )
	{
    	cudaSetDevice(i);
		cudaDeviceReset();
		cudaDeviceSynchronize();
		cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);

	    cudaGetDeviceProperties(&DevProp,i);
    	cout<<"INFO-D: Device: "<<i<<" with "<<DevProp.multiProcessorCount<<" SM "<<DevProp.name<<endl;
	    cudaMemGetInfo(&a[i],&t[i]);
	    cout<<" Clock (MHZ): "<<DevProp.clockRate/1000 << " Total Memory (MB): "
			  <<t[i]*9.53674e-7<<" Free (MB): "<<a[i]*9.53674e-7<< " Clock(MHZ): "<<
			   DevProp.memoryClockRate/1000<< " BUS (Bit) "<<
			   DevProp.memoryBusWidth<<" Copy Engines "<<DevProp.asyncEngineCount<<" MultiBoard "<<DevProp.isMultiGpuBoard<<endl;

		if(i==device_num)
	    {
	     num_sm = DevProp.multiProcessorCount;
	    }

		if( Device_use_multi_gpu  && (i < num_devices_use) )
        {
	        cudaSetDevice(i);
	        for (int j = 0; j < 4; j++)
		    {
	          cudaStreamCreate(&streams[i*4+j]);
	          cudaEventCreate(&events[i*4+j]);
		    }
         }
		cudaDeviceSynchronize();
      }

	   cout<<" \n";
	   cout<<" \n";

       if(num_devices_use==0)
       {
		  cudaSetDevice(device_num);


          cudaGetDeviceProperties(&DevProp,device_num);
			  
	      cout<<"INFO-D: Using Device: SM "<<num_sm  <<" , "<<DevProp.name<<endl;
	     
		  for (int j = 0; j < 4; j++)
		  {
	        cudaStreamCreate(&streams[j]);
	        cudaEventCreate(&events[j]);
		  }
        }

       cout<<" \n";

	   cudaSetDevice(0);


}
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
   /* 1) ENTRY POINT Copies all data to the device sets kernel parameters*/
/*---------------------------------------------------------------------------*/
void Device_Set_SimulationData(  WOBJ*           h_WorldOBJ,
								 POBJ*           h_ParticleOBJ,
								 DOBJ*           h_DynamicOBJ,
								 SimulationInfo* h_SimInfo,
								 int             h_Num_WorldOBJ,
								 int             h_Num_ParticleOBJ,
								 int             h_Num_DynamicOBJ,
								 InitConfig      *h_PosConfig             )

{

	m_WorldObject = h_WorldOBJ;

	error_check = 1;
	LogFileD.open( "../SimulationDevice.Log" );

	 cout<<" \n";
	 cout<<" \n";
	cout<<"INFO-D: Starting Device Configuration \n";
	
	 cout<<" \n";
	/* Create a local copy of the SimOBJ */
	SimInfo_C            = h_SimInfo;
	Device_use_multi_gpu = h_PosConfig->multi_gpu;
	device_num           = h_PosConfig->use_device;
    num_devices_use      = h_PosConfig->num_gpus;

	Get_DeviceInfo();
	
	cudaDeviceSynchronize();
	
	cudaSetDevice(0);


	NUMPARTICLES_MAX_SIM = h_SimInfo->Num_Particles;
	NUMDYOBJECTS         = h_Num_DynamicOBJ;
	NUMWOBJECTS          = h_Num_WorldOBJ;
	D_NUMPARTICLES       = NUMPARTICLES_MAX_SIM;/* Used for killing Silo */
	NUMPARTICLES_Current = NUMPARTICLES_MAX_SIM;/* Used for filling */


	cout<<"INFO-D: Total Particles = "<<NUMPARTICLES_MAX_SIM<<endl;
	cout<<"INFO-D: Size Types = "<<h_SimInfo->Num_ParticleObjects<<endl;

	for(int i=0; i< h_SimInfo->Num_ParticleObjects;i++)
	{
		cout<<"INFO-D: Num per type = "<<h_SimInfo->Num_ParticlesPerObject[i]<<endl;
	}


	/* Set thread info on the GPU */
    printf( "RUNLEVEL 0: Max SM occupy \n");

    int num_particles_SM = (int)ceil(NUMPARTICLES_MAX_SIM/(float)num_sm);

    /* Now get the number of blocks per SM */
    int num_blocks_SM =  (int)ceil(num_particles_SM/(float)h_PosConfig->threads_perBlock);
    printf(" blocks per SM %d \n",num_blocks_SM);

    threads_perBlock = h_PosConfig->threads_perBlock;
    int num_threads_block=h_PosConfig->threads_perBlock;

    /* block size is too big */
    if(num_particles_SM < h_PosConfig->threads_perBlock)
    {
	   num_threads_block = num_particles_SM; /* Single block is sufficient */
    }

    GlobalLaunch.dimBlock = make_uint3( num_threads_block,1,1);
    GlobalLaunch.dimGrid  = make_uint3( num_blocks_SM*num_sm,1,1 );

	printf("INFO-D: Launching %d Blocks with %d Threads: Total threads = %d \n",GlobalLaunch.dimGrid.x,
								 GlobalLaunch.dimBlock.x,GlobalLaunch.dimGrid.x*GlobalLaunch.dimBlock.x);


    /* Dynamically allocate Dynamics arrays */
	cudaSetDevice(device_num);
	cudaDeviceSynchronize();

	cout<<"INFO-D: Allocating Dynamic Arrays : \n";

    m_numGridCells = SimInfo_C->num_NNCells.x*SimInfo_C->num_NNCells.y
    		                                 *SimInfo_C->num_NNCells.z;


    cudaMalloc( (void**) &dram_position_com,  sizeof(float3)*NUMPARTICLES_MAX_SIM );
    cudaMalloc( (void**) &dram_velocity_com,  sizeof(float3)*NUMPARTICLES_MAX_SIM );

    if (SimInfo_C->particle_type==1 || SimInfo_C->sphere_orient)
    {
      cudaMalloc( (void**) &dram_position_ornt, sizeof(Quaterion)*NUMPARTICLES_MAX_SIM );
    }


    cudaMalloc( (void**) &dram_velocity_ang,  sizeof(float3)*NUMPARTICLES_MAX_SIM );


    cudaMalloc( (void**) &dram_ObjectType  ,  sizeof(uint)*NUMPARTICLES_MAX_SIM );
    cudaMalloc( (void**) &dram_P_ID        ,  sizeof(int)*NUMPARTICLES_MAX_SIM );



    cudaMalloc( (void**) &dram_Lifter_Contact_Hist ,      sizeof(Contact_Info)*NUMPARTICLES_MAX_SIM );


    if(SimInfo_C->use_symmetry)
    {
      cudaMalloc( (void**) &dram_force_com_X ,  sizeof(float)*NUMPARTICLES_MAX_SIM );
      cudaMalloc( (void**) &dram_force_com_Y ,  sizeof(float)*NUMPARTICLES_MAX_SIM );
      cudaMalloc( (void**) &dram_force_com_Z ,  sizeof(float)*NUMPARTICLES_MAX_SIM );

      cudaMalloc( (void**) &dram_force_ang_X ,  sizeof(float)*NUMPARTICLES_MAX_SIM );
      cudaMalloc( (void**) &dram_force_ang_Y ,  sizeof(float)*NUMPARTICLES_MAX_SIM );
      cudaMalloc( (void**) &dram_force_ang_Z ,  sizeof(float)*NUMPARTICLES_MAX_SIM );
    }
    else
    {
      cudaMalloc( (void**) &dram_force_PP    , sizeof(float3)   *NUMPARTICLES_MAX_SIM );
      cudaMalloc( (void**) &dram_force_PP_ang, sizeof(float3)   *NUMPARTICLES_MAX_SIM );
    }


    /* NN arrays */
    /* Allocation of cell data */
    cudaMalloc( (void**) &m_dCellStart,         sizeof(uint)*m_numGridCells);
    cudaMalloc( (void**) &m_dCellEnd,           sizeof(uint)*m_numGridCells);

    cudaMalloc( (void**) &dram_GridParticleHash , sizeof(uint)*NUMPARTICLES_MAX_SIM );
    cudaMalloc( (void**) &dram_GridParticleIndex, sizeof(int)*NUMPARTICLES_MAX_SIM );


    cudaMalloc( (void**) &dram_force_com_Wall       , sizeof(float3)*NUMPARTICLES_MAX_SIM );
    cudaMalloc( (void**) &dram_force_com_Lifter       , sizeof(float3)*NUMPARTICLES_MAX_SIM );

    cudaMalloc( (void**) &dram_force_ang_Wall   , sizeof(float3)*NUMPARTICLES_MAX_SIM );
    cudaMalloc( (void**) &dram_force_ang_Lifter , sizeof(float3)*NUMPARTICLES_MAX_SIM );

	get_wall_forces    = h_PosConfig->get_wall_points;
	get_vobject_forces = h_PosConfig->get_vobject_points;
	get_P_forces = h_PosConfig->get_particle_points;

	if(get_wall_forces)
	{
	  printf(" devvice allco W\n");
      cudaMalloc( (void**) &dWallContactPoints, sizeof(float3)*NUMPARTICLES_MAX_SIM );
	}

	if(get_vobject_forces)
	{
	  printf(" devvice allco V\n");
      cudaMalloc( (void**) &dVOContactPoints, sizeof(float3)*NUMPARTICLES_MAX_SIM );
	}

	if(get_P_forces)
	{
	 printf(" devvice allco P\n");
     cudaMalloc( (void**) &dPContactPoints, sizeof(float3)*NUMPARTICLES_MAX_SIM*64 );
     cudaMalloc( (void**) &dNumContacts, sizeof(uint)*NUMPARTICLES_MAX_SIM);
	}

	                 /* Allocation of PDA data */
    cudaMalloc( (void**) &dSortedPos  , sizeof(float3)   *NUMPARTICLES_MAX_SIM );     /* 120MB */
    cudaMalloc( (void**) &dSortedVel  , sizeof(float3)   *NUMPARTICLES_MAX_SIM );     /* 120MB */

    if(SimInfo_C->particle_type==polyhedra || SimInfo_C->sphere_orient)
	{
      cudaMalloc( (void**) &dSortedPosQ , sizeof(Quaterion)*NUMPARTICLES_MAX_SIM ); /* 160MB */
	}

    cudaMalloc( (void**) &dSorted_velocity_ang , sizeof(float3)   *NUMPARTICLES_MAX_SIM );    /* 120MB */
    cudaMalloc( (void**) &dSortedPType, sizeof(uint)     *NUMPARTICLES_MAX_SIM );        /* 12MB */
    cudaMalloc( (void**) &dSortedPID  , sizeof(int)      *NUMPARTICLES_MAX_SIM );        /* 12MB */

    cudaMalloc( (void**) &dSortedLifter_Contact_Hist     , sizeof(Contact_Info)*NUMPARTICLES_MAX_SIM );


	//cudaMallocHost((void**)&host_pinned_pos,  sizeof(float3)   *NUMPARTICLES_MAX_SIM);
	//cudaMallocHost((void**)&host_pinned_vel,  sizeof(float3)   *NUMPARTICLES_MAX_SIM);
	//cudaMallocHost((void**)&host_pinned_Type, sizeof(uint)     *NUMPARTICLES_MAX_SIM );  
	//cudaMallocHost((void**)&host_pinned_ID  , sizeof(int)      *NUMPARTICLES_MAX_SIM );


	 if(!NO_NN)
	 {
     cudaMalloc( (void**) &dBroad_List, sizeof(uint)*NUMPARTICLES_MAX_SIM*32 );
	 cudaMalloc( (void**) &dNumNN     , sizeof(uint)*NUMPARTICLES_MAX_SIM    );
	 }
     cudaDeviceSynchronize();


	float Tallys [6];

	for(int i=0;i<6;i++)
	{
	  Tallys[i]=0.0f;
	}

	float *d_Tally;
	cudaMalloc((void **)&d_Tally, sizeof(float)*6);
	cudaDeviceSynchronize();

	cudaMemcpy(d_Tally,Tallys, sizeof(float)*6,          cudaMemcpyHostToDevice );
	cudaDeviceSynchronize();

	Set_Tallys<<< 1,1 >>>(d_Tally);
	cudaDeviceSynchronize();
	cudaFree(d_Tally);



    cudaDeviceSynchronize();

	 if(error_check==1)
	 {
		cudaError_t errormsg=cudaGetLastError();
		if(errormsg>0)
		{
		 cout<<"Set Simulation Data: Dram memory alloc "<<cudaGetErrorString(errormsg)<<endl;
		 exit(1);
		}
	 }


    cout<<"INFO-D: Copying Constant Data Device: \n";

	/* 1) Copy Objects to Device */
	cudaMemcpyToSymbol( SimParms,       h_SimInfo,     sizeof( SimulationInfo ) );
	cudaMemcpyToSymbol( WorldObject,    h_WorldOBJ,    sizeof(WOBJ)*h_Num_WorldOBJ );
	cudaMemcpyToSymbol( ParticleObject, h_ParticleOBJ, sizeof(POBJ)*h_SimInfo->Num_ParticleObjects);


	if(h_Num_DynamicOBJ>0)
	{
		cudaMemcpyToSymbol( VolumeObject,  h_DynamicOBJ,  sizeof(DOBJ)*h_Num_DynamicOBJ);
	}

	if(!h_PosConfig->use_file)
	{
	if(h_PosConfig->grid_type<2)
	{
	 cudaMemcpyToSymbol( InitVel,        h_PosConfig->velocity,     sizeof(float3)*2);
	}
	else
	{
	  h_PosConfig->velocity[0]= make_float3(h_PosConfig->launch_Vel.x,h_PosConfig->launch_Vel.y,h_PosConfig->launch_Vel.z);
	  cudaMemcpyToSymbol( InitVel,        h_PosConfig->velocity,     sizeof(float3)*2);
	}
	cudaDeviceSynchronize();
	}
	else
	{
		  h_PosConfig->velocity[0]= make_float3(0.0,0.0,0.0f);
		  cudaMemcpyToSymbol( InitVel,        h_PosConfig->velocity,     sizeof(float3)*2);
	}

	 if(error_check==1)
	 {
		cudaError_t errormsg=cudaGetLastError();
		if(errormsg>0)
		{
		 cout<<"Set Simulation Data: Constant memory alloc "<<cudaGetErrorString(errormsg)<<endl;
		 exit(1);
		}
	 }


	cout<<" GPU Malloc DONE!\n";




   if(SimInfo_C->use_symmetry)
   {
	/* Launch Kernel to set Initial parameters on device */
    Set_Init_Config_Device<<< GlobalLaunch.dimGrid, GlobalLaunch.dimBlock >>>
	   		                  ( SimInfo_C->Num_ParticleObjects,
	   		                	SimInfo_C->Num_WorldObjects,
	   		                	NUMDYOBJECTS,

	   		                	dram_force_com_X, dram_force_com_Y, dram_force_com_Z,
	   		                	dram_force_com_Wall,dram_force_com_Lifter,
	   		                	dram_force_ang_X, dram_force_ang_Y, dram_force_ang_Z,
	   		                	dram_force_ang_Wall,dram_force_ang_Lifter,

	   		                	dram_velocity_com,
			                    dram_velocity_ang,
			                    dram_Lifter_Contact_Hist,SimInfo_C->Num_Particles,0);
   }
   else
   {
		/* Launch Kernel to set Initial parameters on device */
	    Set_Init_Config_Device<<< GlobalLaunch.dimGrid, GlobalLaunch.dimBlock >>>
		   		                  ( SimInfo_C->Num_ParticleObjects,
		   		                	SimInfo_C->Num_WorldObjects,
		   		                	NUMDYOBJECTS,

		   		                	dram_force_PP,
		   		                	dram_force_com_Wall,dram_force_com_Lifter,
		   		                	dram_force_PP_ang,
		   		                	dram_force_ang_Wall,dram_force_ang_Lifter,

		   		                	dram_velocity_com,
				                    dram_velocity_ang,
									
				                    dram_Lifter_Contact_Hist,SimInfo_C->Num_Particles,0 );
   }
    cudaDeviceSynchronize();

	 if(error_check==1)
	 {
		cudaError_t errormsg=cudaGetLastError();
		if(errormsg>0)
		{
		 cout<<"Set Simulation Data: Set_Init_Config_Device Kernel "<<cudaGetErrorString(errormsg)<<endl;
		 exit(1);
		}
	 }

    /* Select Initial Position */
	if(!h_PosConfig->use_file)
	{
	  if( h_PosConfig->grid_type==0 )
	  {
	    Set_Position_DefaultGrid(h_PosConfig,h_ParticleOBJ);
	  }
	  else if( h_PosConfig->grid_type==1 )
	  {
		Set_Position_Random_Grid(h_PosConfig,h_ParticleOBJ);
	  }
	  else if( h_PosConfig->grid_type==2 )
	  {
		  Set_Position_Fill_Rectangle(h_ParticleOBJ,h_WorldOBJ,h_PosConfig->launch_Vel.y,h_PosConfig->fill_plane_start,h_PosConfig->fill_plane_end);
		  isFilling = true;
	  }

	}
	cudaDeviceSynchronize();

	LogFileD<<"INFO-D: All data copied to Device "<<std::endl;


    cudaMemGetInfo(&a[device_num],&t[device_num]);
    cout<<" Clock (MHZ): "<<" Total Memory (MB): "
		  <<t[device_num]*9.53674e-7<<" Free (MB): "<<a[device_num]*9.53674e-7 <<endl;

    LogFileD.close();

	 if(error_check==1)
	 {
		cudaError_t errormsg=cudaGetLastError();
		if(errormsg>0)
		{
		 cout<<"Set Simulation Data: Error at end: Debug further "<<cudaGetErrorString(errormsg)<<endl;
		 exit(1);
		}
	 }

}
/*-----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------*/
void Filling_Logic_Simulation()
{
	/* Dynamically add particles to simulation by updating number of particles
	 * and changing the downstream kernel launches */
	if(isFilling)
	{
		if((fill_counter>pack_clear_steps || NUMPARTICLES_Current==0) )
		{
		  fill_counter=0;

		  if(NUMPARTICLES_Current<NUMPARTICLES_MAX_SIM)
		  {
			Fill_State();
		  }
		  if(NUMPARTICLES_Current>=NUMPARTICLES_MAX_SIM)
		  {
			  printf("fill ended\n");
			  if(SimInfo_C->Simulation_Type==silo)
			  {
				if(!Device_use_multi_gpu)
				{
			      Packing_Vel<<< GlobalLaunch.dimGrid, GlobalLaunch.dimBlock >>>(dram_velocity_com);
				}
				else
				{
                  Packing_Vel<<< GlobalLaunch.dimGrid, GlobalLaunch.dimBlock >>>(dram_velocity_comM[0]);
				}
				cudaDeviceSynchronize();
					 if(error_check==1)
	 {
	    	cudaError_t errormsg=cudaGetLastError();
	    	if(errormsg>0)
	    	{
	    	cout<<"DEM_UpdateSim: Filling "<<cudaGetErrorString(errormsg)<<endl;
	    	 exit(1);
	    	}
	 }

			  }
			  isFilling = false;
		  }
		}
		else
		{
			fill_counter++;
		}
	}
	else
	{
		NUMPARTICLES_Current = NUMPARTICLES_MAX_SIM;
	}

}
/*-----------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------*/
 /* Flags particles to be removed from the sim and gives them a position such 
    that the hash will be large */
/*-----------------------------------------------------------------------------*/
void Flag_DeadParticles()
{

       stepc++;
       /* Dont check every step to save computations */
       if( stepc >= 250 )
       {

		   if(SimInfo_C->use_valid_zone)
    	  {
        	 Remove_PL_Volume_Particles<<< GlobalLaunch.dimGrid, GlobalLaunch.dimBlock >>>(dram_position_com,dram_P_ID,dram_ObjectType,NUMPARTICLES_Current);
        	 cudaDeviceSynchronize();
    	  }
    	  else /* Silo */
    	  {
    	    Remove_Silo_Discharged_Particles<<< GlobalLaunch.dimGrid, GlobalLaunch.dimBlock >>>(dram_position_com,dram_P_ID,NUMPARTICLES_Current);
    	    cudaDeviceSynchronize();
    	  }

    	  /* Copy number killed to Host (Check if can do this without memcpy */
    	  int *d_kill;
    	  int  h_kill[2];


    	  float *d_stats;
    	  float  h_stats[2];

    	  cudaMalloc( (void**) &d_kill        , sizeof(int)*2   );
    	  cudaMalloc( (void**) &d_stats        , sizeof(float)*2   );

    	  Kill_copy<<<1,1>>>(d_kill,d_stats);
    	  cudaDeviceSynchronize();

    	  cudaMemcpy( h_kill,        d_kill        ,sizeof(int)*2,
    	    		                                     cudaMemcpyDeviceToHost );

    	  cudaMemcpy( h_stats,        d_stats        ,sizeof(float)*2,
    	    		                                     cudaMemcpyDeviceToHost );
    	  cudaDeviceSynchronize();
    	  

    	  num_dead=h_kill[0];
    	  cudaFree(d_kill);

          dis_mass = h_stats[0];
          dis_vol = h_stats[1];
    	  cudaFree(d_stats);
    	  stepc = 0;
       }
     

	 if(error_check==1)
	 {
	    	cudaError_t errormsg=cudaGetLastError();
	    	if(errormsg>0)
	    	{
	    	 cout<<"DEM_UpdateSim: KillParticles "<<cudaGetErrorString(errormsg)<<endl;
	    	 exit(1);
	    	}
	 }

	 
}

/* Adjusts the block size so that threads do not access dead particles */
void Remove_Flagged_Particles()
{
	//printf(" Start %d \n",NUMPARTICLES_Current);
		NUMPARTICLES_Current -= num_dead;
		//printf(" Removing %d Particles current %d \n",num_dead,NUMPARTICLES_Current);
		//printf("sim type %d \n",SimInfo_C->Simulation_Type);
		if( SimInfo_C->Simulation_Type==Pulp_Lift)
		{
		  printf(" %f %d %f %f \n", Call_Time, num_dead, dis_mass, dis_vol);
		}
		num_dead=0;

	    int num_particles_SM = (int)ceil(NUMPARTICLES_Current/(float)num_sm);

	    /* Now get the number of blocks per SM */
	   	int num_blocks_SM =  (int)ceil(num_particles_SM/(float)threads_perBlock);

	   	   int num_threads_block=threads_perBlock;

	   	   /* block size is too big */
	   	   if(num_particles_SM < threads_perBlock)
	   	   {
	   		   num_threads_block = num_particles_SM; /* Single block is sufficient */
	   	   }


	   	   GlobalLaunch.dimBlock = make_uint3( num_threads_block,1,1);
	   	   GlobalLaunch.dimGrid  = make_uint3( num_blocks_SM*num_sm,1,1 );

	   	   SimInfo_C->Num_Particles = NUMPARTICLES_Current;
	   	   cudaMemcpyToSymbol( SimParms,       SimInfo_C,     sizeof( SimulationInfo ) );
	   	   cudaDeviceSynchronize();
	
		   
}


/*---------------------------------------------------------------------------*/
               /* Performs the simulation using non-symmetry */
/*---------------------------------------------------------------------------*/
void Device_DEM_UpdateSim()
{

	 if(error_check==1)
	 {
	    	cudaError_t errormsg=cudaGetLastError();
	    	if(errormsg>0)
	    	{
	    	cout<<"DEM_UpdateSim: Start "<<cudaGetErrorString(errormsg)<<endl;
	    	 exit(1);
	    	}
	 }

	if(isFilling)
	{
	 Filling_Logic_Simulation();
	}
	                /* Reorder All arrays based on hash */
    /*-----------------------------------------------------------------------*/
           /* 1. Calculate the GridIndex of Each Particle (Grid Hash) */
    /*-----------------------------------------------------------------------*/

    /* Calculate the Hash for All particles based on its Spatial Position */
    Kernel_SpatialDecomp_CalcPHash<<< GlobalLaunch.dimGrid, GlobalLaunch.dimBlock >>>
    		                          ( dram_GridParticleHash, dram_GridParticleIndex,
    		                            dram_position_com,NUMPARTICLES_Current );
    /*-----------------------------------------------------------------------*/

	 if(error_check==1)
	 {
	    	cudaError_t errormsg=cudaGetLastError();
	    	if(errormsg>0)
	    	{
	    	cout<<"DEM_UpdateSim: Hash "<<cudaGetErrorString(errormsg)<<endl;
	    	 exit(1);
	    	}
	 }
    /*-----------------------------------------------------------------------*/
           /* 2. Sort Particles according to (Grid Hash) */
    /*-----------------------------------------------------------------------*/

    /* In index 0 is lowest and last is highest */
    SpatialDecomp_ThrustSort_ByHash(dram_GridParticleHash, dram_GridParticleIndex,NUMPARTICLES_Current);

    /*-----------------------------------------------------------------------*/

    /* Reorder arrays based on hash */

    /* m_dGridParticleIndex contains the order which the dynamics arrays must be reordered */

    /* Sort arrays using temp copies */
    Kernel_SortArrays
      <<< GlobalLaunch.dimGrid, GlobalLaunch.dimBlock >>>
      ( dram_GridParticleIndex,
    	dSortedPos,
    	dSortedPosQ,
    	dSortedVel,
    	dSorted_velocity_ang,

    	dSortedPType,
    	dSortedPID,

        dSortedLifter_Contact_Hist,

        dram_position_com,
        dram_position_ornt,
        dram_velocity_com,
        dram_velocity_ang,
        dram_ObjectType,
        dram_P_ID,

        dram_Lifter_Contact_Hist,

        NUMPARTICLES_Current );


	 if(error_check==1)
	 {
	    	cudaError_t errormsg=cudaGetLastError();
	    	if(errormsg>0)
	    	{
	    	cout<<"DEM_UpdateSim: Sort "<<cudaGetErrorString(errormsg)<<endl;
	    	 exit(1);
	    	}
	 }

	 cudaMemsetAsync(m_dCellStart, 0xffffffff, m_numGridCells*sizeof(uint));
     cudaMemsetAsync(m_dCellEnd,   0xffffffff, m_numGridCells*sizeof(uint));


    /*-----------------------------------------------------------------------*/
                /* 4. Update the pos and vel based on sorting in 3. */
      /*-----------------------------------------------------------------------*/
      Kernel_SpatialDecomp_ReorderPDA
                  <<< GlobalLaunch.dimGrid, GlobalLaunch.dimBlock,0, streams[1] >>>
      		    ( dSortedPos,dSortedPosQ, dSortedVel,dSorted_velocity_ang, dSortedPType,dSortedPID,
      		      dSortedLifter_Contact_Hist,

      		      dram_position_com,
      		      dram_position_ornt,
      		      dram_velocity_com,
                  dram_velocity_ang,
                  dram_ObjectType,
                  dram_P_ID,
                  dram_Lifter_Contact_Hist,

                  NUMPARTICLES_Current );


      /*-----------------------------------------------------------------------*/
                 /* 3. Bin data into Cells based on the Hash and index */
      /*-----------------------------------------------------------------------*/


          /* Bins data into cells - Creates from scratch at each step */
          Kernel_SpatialDecomp_BinData
            <<< GlobalLaunch.dimGrid, GlobalLaunch.dimBlock,0, streams[2] >>>
            ( dram_GridParticleHash,m_dCellStart, m_dCellEnd,NUMPARTICLES_Current);
          /*-----------------------------------------------------------------------*/

      /*-----------------------------------------------------------------------*/

	 if(error_check==1)
	 {
	    	cudaError_t errormsg=cudaGetLastError();
	    	if(errormsg>0)
	    	{
	    	cout<<"DEM_UpdateSim: Reorder "<<cudaGetErrorString(errormsg)<<endl;
	    	 exit(1);
	    	}
	 }


    cudaDeviceSynchronize();

    /* Remove flagged particles from the sim by reducing num threads */
	if(num_dead>0)
	{
	  Remove_Flagged_Particles();
	}


	 if(error_check==1)
	 {
	    	cudaError_t errormsg=cudaGetLastError();
	    	if(errormsg>0)
	    	{
	    	cout<<"DEM_UpdateSim: BinData "<<cudaGetErrorString(errormsg)<<endl;
	    	 exit(1);
	    	}
	 }




 if(!SimInfo_C->use_symmetry)
 {
    /*-----------------------------------------------------------------------*/
                       /* 5. Scan ALL NN PAIRS */
    /*-----------------------------------------------------------------------*/
	 cudaDeviceSynchronize();



	 if(NO_NN)	  
     { 



	    Kernel_Spheres_PP
                         <<< GlobalLaunch.dimGrid, GlobalLaunch.dimBlock,0, streams[1] >>>
    		             (get_P_forces,
				  dPContactPoints, m_dCellStart,m_dCellEnd, dram_position_com,
				  dram_velocity_com,
				  dram_velocity_ang,
				  dram_ObjectType,
				  dram_P_ID, 
				  dram_force_PP, dram_force_PP_ang,
				  NUMPARTICLES_Current,dNumContacts);
		
		if(error_check==1)
 	    {
 	    	cudaError_t errormsg=cudaGetLastError();
 	    	if(errormsg>0)
 	    	{
 	    	cout<<"DEM_UpdateSim: SpheresPP "<<cudaGetErrorString(errormsg)<<endl;
 	    	 exit(1);
 	    	}
 	     }
	 }
	 else
	 {
      
        Kernel_BroadCollisionDetection_NonSymmetry
                         <<< GlobalLaunch.dimGrid, GlobalLaunch.dimBlock>>>
    		             ( m_dCellStart,m_dCellEnd,dBroad_List,dNumNN,dram_position_com,dram_ObjectType,dram_P_ID,NUMPARTICLES_Current);
 	    
		if(error_check==1)
 	    {
 	    	cudaError_t errormsg=cudaGetLastError();
 	    	if(errormsg>0)
 	    	{
 	    	cout<<"DEM_UpdateSim: NN Search NonSymm "<<cudaGetErrorString(errormsg)<<endl;
 	    	 exit(1);
 	    	}
 	     }
		cudaDeviceSynchronize();
    /*-----------------------------------------------------------------------*/
                       /* 6. Apply Particle Particle Forces */
    /*-----------------------------------------------------------------------*/


         if(SimInfo_C->particle_type==0)
	     {
			Kernel_ParticleInteraction_Spheres_NonSymmetry<<< GlobalLaunch.dimGrid, GlobalLaunch.dimBlock,0, streams[1] >>>
				( get_P_forces,
				  dPContactPoints,
				  dNumNN,dBroad_List,
				  dram_position_com,
				  dram_velocity_com,
				  dram_velocity_ang,
				  dram_ObjectType,
				  dram_P_ID,dram_force_PP,dram_force_PP_ang,NUMPARTICLES_Current   );


		 	 if(error_check==1)
		 	 {
		 	    	cudaError_t errormsg=cudaGetLastError();
		 	    	if(errormsg>0)
		 	    	{
		 	    	cout<<"DEM_UpdateSim: Force spheres NonSymm "<<cudaGetErrorString(errormsg)<<endl;
		 	    	 exit(1);
		 	    	}
		 	 }
		  }
      else
      {
			Kernel_ParticleInteraction_Polyhedra_NonSymmetry<<< GlobalLaunch.dimGrid, GlobalLaunch.dimBlock,0, streams[1]  >>>
				( dNumNN,dBroad_List,
				  dram_position_com,
				  dram_position_ornt,
				  dram_velocity_com,
				  dram_velocity_ang,
				  dram_ObjectType,
				  dram_P_ID,dram_force_PP,dram_force_PP_ang,
				  NUMPARTICLES_Current);

      }
	}

  }/* symmetry not verified */
  else
  {

	    /*-----------------------------------------------------------------------*/
	                       /* 5. Scan 1/2 NN PAIRS */
	    /*-----------------------------------------------------------------------*/

	      Kernel_BroadCollisionDetection_Symmetry
	                         <<< GlobalLaunch.dimGrid, GlobalLaunch.dimBlock>>>
	    		             ( m_dCellStart,m_dCellEnd,dBroad_List,dNumNN,dram_position_com,dram_ObjectType,dram_P_ID,NUMPARTICLES_Current);
	      cudaDeviceSynchronize();

	      //PrintNN<<< 1,1>>>(dNumNN,dBroad_List,dram_P_ID);

		 	 if(error_check==1)
		 	 {
		 	    	cudaError_t errormsg=cudaGetLastError();
		 	    	if(errormsg>0)
		 	    	{
		 	    	 cout<<"DEM_UpdateSim: NN search Symm "<<cudaGetErrorString(errormsg)<<endl;
		 	    	 exit(1);
		 	    	}
		 	 }

   } /* End symmetry check */

    if(NUMDYOBJECTS>0)
    {

      if(SimInfo_C->particle_type==spheres)
      {
	    VolumeObject_InteractionSpheres<<< GlobalLaunch.dimGrid, GlobalLaunch.dimBlock,0, streams[2] >>>(
			     get_vobject_forces,
				 dVOContactPoints,
			     
				 dram_force_com_Lifter,
			     dram_force_ang_Lifter,
			     dram_Lifter_Contact_Hist,
			     dram_position_com,
			     dram_velocity_com,
			     dram_velocity_ang,
	             dram_ObjectType, dram_P_ID,NUMPARTICLES_Current);
	   // cudaDeviceSynchronize();
      }
      else
      {
  	    VolumeObject_InteractionPolyhedra<<< GlobalLaunch.dimGrid, GlobalLaunch.dimBlock,0, streams[2] >>>(
  			     dram_force_com_Lifter,
  			     dram_force_ang_Lifter,
  			     dram_Lifter_Contact_Hist,
  			     dram_position_com,
  			     dram_position_ornt,
  			     dram_velocity_com,
  			     dram_velocity_ang,
  	             dram_ObjectType, dram_P_ID,NUMPARTICLES_Current);
  	    //cudaDeviceSynchronize();

      }

	 	 if(error_check==1)
	 	 {
	 	    	cudaError_t errormsg=cudaGetLastError();
	 	    	if(errormsg>0)
	 	    	{
	 	    	 cout<<"DEM_UpdateSim: Lifter Inter "<<cudaGetErrorString(errormsg)<<endl;
	 	    	 exit(1);
	 	    	}
	 	 }
    }


    if(SimInfo_C->particle_type==spheres)
    {
		WorldObject_InteractionSpheres_Planar<<< GlobalLaunch.dimGrid, GlobalLaunch.dimBlock,0, streams[3] >>>( get_wall_forces,
			                                                                        dram_force_com_Wall,
			                                                                        dram_force_ang_Wall,
			                                                                        dWallContactPoints,
			                                                                        dram_position_com,
			                                                                        dram_position_ornt,
			                                                                        dram_velocity_com,
			                                                                        dram_velocity_ang,
                                                                                    dram_ObjectType,
                                                                                    dram_P_ID,NUMPARTICLES_Current,0                   );
    }
    else
    {
    	WorldObject_InteractionPolyhedra_Planar<<< GlobalLaunch.dimGrid, GlobalLaunch.dimBlock,0, streams[3]  >>>( get_wall_forces, is_cylinder_rotation,
    			                                                                        dram_force_com_Wall,
    			                                                                        dram_force_ang_Wall,
    			                                                                        dWallContactPoints,
    			                                                                        dram_position_com,
    			                                                                        dram_position_ornt,
    			                                                                        dram_velocity_com,
    			                                                                        dram_velocity_ang,
                                                                                        dram_ObjectType,
                                                                                        dram_P_ID,NUMPARTICLES_Current                   );
    }


	 if(error_check==1)
	 {
	    	cudaError_t errormsg=cudaGetLastError();
	    	if(errormsg>0)
	    	{
	    	 cout<<"DEM_UpdateSim: Surface Interaction "<<cudaGetErrorString(errormsg)<<endl;
	    	 exit(1);
	    	}
	 }

	 cudaDeviceSynchronize();
    
	if(SimInfo_C->use_symmetry)
    {
      Integrate_Euler_Symmetry<<< GlobalLaunch.dimGrid, GlobalLaunch.dimBlock >>>
            ( dram_force_com_X, dram_force_com_Y,dram_force_com_Z,
              dram_force_ang_X, dram_force_ang_Y, dram_force_ang_Z,
              dram_force_com_Wall,   dram_force_ang_Wall,
              dram_force_com_Lifter, dram_force_ang_Lifter,

              dram_position_com, dram_position_ornt,
              dram_velocity_com, dram_velocity_ang,
              dram_ObjectType, dram_P_ID,NUMPARTICLES_Current );

	 	 if(error_check==1)
	 	 {
	 	    	cudaError_t errormsg=cudaGetLastError();
	 	    	if(errormsg>0)
	 	    	{
	 	    	 cout<<"DEM_UpdateSim: Integration Euler Symm "<<cudaGetErrorString(errormsg)<<endl;
	 	    	 exit(1);
	 	    	}
	 	 }

  }
	   else
	   {
         if(SimInfo_C->Integration_Type==0)
         {
            Integrate_Euler_NonSymmetry<<< GlobalLaunch.dimGrid, GlobalLaunch.dimBlock >>>
              ( dram_force_PP,          dram_force_PP_ang,
              dram_force_com_Wall,   dram_force_ang_Wall,
              dram_force_com_Lifter, dram_force_ang_Lifter,

              dram_position_com, dram_position_ornt,
              dram_velocity_com, dram_velocity_ang,
              dram_ObjectType, dram_P_ID,NUMPARTICLES_Current,Get_WallPoints );


	 	  if(error_check==1)
	 	  {
	 	    	cudaError_t errormsg=cudaGetLastError();
	 	    	if(errormsg>0)
	 	    	{
	 	    	 cout<<"DEM_UpdateSim: Integration Euler Non Symm "<<cudaGetErrorString(errormsg)<<endl;
	 	    	 exit(1);
	 	    	}
	 	   }
         }
	   }



	 if(error_check==1)
	 {
	    	cudaError_t errormsg=cudaGetLastError();
	    	if(errormsg>0)
	    	{
	    	 cout<<"DEM_UpdateSim: END "<<cudaGetErrorString(errormsg)<<endl;
	    	 exit(1);
	    	}
	 }

	if( SimInfo_C->Kill_Particles)
    {
	  Flag_DeadParticles();
	}
	 Call_Time+=SimInfo_C->InitalDelta_t;



	 if (isFilling)
	 {
		 fill_time+=SimInfo_C->InitalDelta_t;
		 if(fill_time>0.25f)
		 {
			 printf(" Filled Particles %d \n",NUMPARTICLES_Current);
			 fill_time = 0.0f;
		 }
	 }


}
/*---------------------------------------------------------------------------*/











void Sync_Killed()
{

    /* Calculate the Hash for All particles based on its Spatial Position */
    Kernel_SpatialDecomp_CalcPHash<<< GlobalLaunch.dimGrid, GlobalLaunch.dimBlock >>>
    		                          ( dram_GridParticleHash, dram_GridParticleIndex,
    		                            dram_position_com,NUMPARTICLES_Current );
    /*-----------------------------------------------------------------------*/
    cudaDeviceSynchronize();

	 if(error_check==1)
	 {
	    	cudaError_t errormsg=cudaGetLastError();
	    	if(errormsg>0)
	    	{
	    	cout<<"DEM_UpdateSim: Hash "<<cudaGetErrorString(errormsg)<<endl;
	    	 exit(1);
	    	}
	 }
    /*-----------------------------------------------------------------------*/
           /* 2. Sort Particles according to (Grid Hash) */
    /*-----------------------------------------------------------------------*/

    /* In index 0 is lowest and last is highest */
    SpatialDecomp_ThrustSort_ByHash(dram_GridParticleHash, dram_GridParticleIndex,NUMPARTICLES_Current);
    cudaDeviceSynchronize();
    /*-----------------------------------------------------------------------*/


    /* m_dGridParticleIndex contains the order which the dynamics arrays must be reordered */

    /* Sort arrays using temp copies */
    Kernel_SortArrays
      <<< GlobalLaunch.dimGrid, GlobalLaunch.dimBlock >>>
      ( dram_GridParticleIndex,
    	dSortedPos,
    	dSortedPosQ,
    	dSortedVel,
    	dSorted_velocity_ang,

    	dSortedPType,
    	dSortedPID,

        dSortedLifter_Contact_Hist,

        dram_position_com,
        dram_position_ornt,
        dram_velocity_com,
        dram_velocity_ang,
        dram_ObjectType,
        dram_P_ID,

        dram_Lifter_Contact_Hist,

        NUMPARTICLES_Current );

    cudaDeviceSynchronize();

	 if(error_check==1)
	 {
	    	cudaError_t errormsg=cudaGetLastError();
	    	if(errormsg>0)
	    	{
	    	cout<<"DEM_UpdateSim: Sort "<<cudaGetErrorString(errormsg)<<endl;
	    	 exit(1);
	    	}
	 }
    /*-----------------------------------------------------------------------*/
                /* 4. Update the pos and vel based on sorting in 3. */
      /*-----------------------------------------------------------------------*/
      Kernel_SpatialDecomp_ReorderPDA
                  <<< GlobalLaunch.dimGrid, GlobalLaunch.dimBlock >>>
      		    ( dSortedPos,dSortedPosQ, dSortedVel,dSorted_velocity_ang, dSortedPType,dSortedPID,
      		      dSortedLifter_Contact_Hist,

      		      dram_position_com,
      		      dram_position_ornt,
      		      dram_velocity_com,
                  dram_velocity_ang,
                  dram_ObjectType,
                  dram_P_ID,
                  dram_Lifter_Contact_Hist,

                  NUMPARTICLES_Current );
     cudaDeviceSynchronize();
      /*-----------------------------------------------------------------------*/

	 if(error_check==1)
	 {
	    	cudaError_t errormsg=cudaGetLastError();
	    	if(errormsg>0)
	    	{
	    	cout<<"DEM_UpdateSim: Reorder "<<cudaGetErrorString(errormsg)<<endl;
	    	 exit(1);
	    	}
	 }



	  Remove_Flagged_Particles();
	  cudaDeviceSynchronize();

}



/*---------------------------------------------------------------------------*/
/*    Returns the new positions to the host to render per time-step/frame    */
/*---------------------------------------------------------------------------*/
void Device_Get_P_Positions( float3 *Pos, Quaterion *Quart_ornt,uint *PType,
		                     int *hP_ID, int *Number_Particles )
{
	cudaDeviceSynchronize();/* Make Sure Nothing running of Device */
	
	if(num_dead>0)
	{
		if(!Device_use_multi_gpu)
		{
		  Sync_Killed();
		  cudaDeviceSynchronize();
		}
		else
		{
			Sync_KilledMulti();
		    cudaDeviceSynchronize();
		}
	}

	Number_Particles[0] =  NUMPARTICLES_Current;

	if( NUMPARTICLES_Current > 0 )
	{

	  if(!Device_use_multi_gpu)
	  {
		cudaMemcpyAsync( Pos,        dram_position_com        ,sizeof(float3)*NUMPARTICLES_Current,
			cudaMemcpyDeviceToHost );


		cudaMemcpyAsync( PType,      dram_ObjectType      ,sizeof(uint)*NUMPARTICLES_Current,
    		                                     cudaMemcpyDeviceToHost );

		cudaMemcpyAsync( hP_ID,       dram_P_ID         ,sizeof(int)*NUMPARTICLES_Current,
    		                                     cudaMemcpyDeviceToHost );

		if(SimInfo_C->particle_type==polyhedra || SimInfo_C->sphere_orient)
		{
			cudaMemcpyAsync( Quart_ornt, dram_position_ornt ,sizeof(Quaterion)*NUMPARTICLES_Current,
    												 cudaMemcpyDeviceToHost );
		}
	  }
	  else
	  {
       		cudaMemcpyAsync( Pos,        dram_position_comM[0]        ,sizeof(float3)*NUMPARTICLES_Current,
			cudaMemcpyDeviceToHost );

		cudaMemcpyAsync( PType,      dram_ObjectTypeM[0]      ,sizeof(uint)*NUMPARTICLES_Current,
    		                                     cudaMemcpyDeviceToHost );

		cudaMemcpyAsync( hP_ID,       dram_P_IDM[0]         ,sizeof(int)*NUMPARTICLES_Current,
    		                                     cudaMemcpyDeviceToHost );

		if(SimInfo_C->particle_type==polyhedra || SimInfo_C->sphere_orient)
		{
			cudaMemcpyAsync( Quart_ornt, dram_position_orntM[0] ,sizeof(Quaterion)*NUMPARTICLES_Current,
    												 cudaMemcpyDeviceToHost );
		}
	  }

	  cudaDeviceSynchronize();
	  //memcpy(Pos,host_pinned_pos,sizeof(float3)*NUMPARTICLES_Current);
	  //memcpy(PType,host_pinned_Type,sizeof(uint)*NUMPARTICLES_Current);
	  //memcpy(hP_ID,host_pinned_ID,sizeof(int)*NUMPARTICLES_Current);

	 if(error_check==1)
	 {
		cudaError_t errormsg=cudaGetLastError();
		if(errormsg>0)
		{
		 cout<<"Get Positions : Error End "<<cudaGetErrorString(errormsg)<<endl;
		 exit(1);
		}
	 }

	}
		else /* Start up during filling */
	{
	 hP_ID[0] =-2;
	}

}
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*    Returns the new positions to the host to render per time-step/frame    */
/*---------------------------------------------------------------------------*/
void Device_Get_PID( uint *PType,
		                     int *hP_ID, int *Number_Particles )
{
	cudaDeviceSynchronize();/* Make Sure Nothing running of Device */
	
	if(num_dead>0)
	{
		if(!Device_use_multi_gpu)
		{
		  Sync_Killed();
		  cudaDeviceSynchronize();
		}
		else
		{
			Sync_KilledMulti();
		    cudaDeviceSynchronize();
		}
	}

	Number_Particles[0] =  NUMPARTICLES_Current;

	if( NUMPARTICLES_Current > 0 )
	{

	  if(!Device_use_multi_gpu)
	  {
		
		cudaMemcpyAsync( PType,      dram_ObjectType      ,sizeof(uint)*NUMPARTICLES_Current,
    		                                     cudaMemcpyDeviceToHost );

		cudaMemcpyAsync( hP_ID,       dram_P_ID         ,sizeof(int)*NUMPARTICLES_Current,
    		                                     cudaMemcpyDeviceToHost );
	  }
	  else
	  {
		cudaMemcpyAsync( PType,      dram_ObjectTypeM[0]      ,sizeof(uint)*NUMPARTICLES_Current,
    		                                     cudaMemcpyDeviceToHost );

		cudaMemcpyAsync( hP_ID,       dram_P_IDM[0]         ,sizeof(int)*NUMPARTICLES_Current,
    		                                     cudaMemcpyDeviceToHost );

		
	  }

	  cudaDeviceSynchronize();
	  
	 if(error_check==1)
	 {
		cudaError_t errormsg=cudaGetLastError();
		if(errormsg>0)
		{
		 cout<<"Get ID : Error End "<<cudaGetErrorString(errormsg)<<endl;
		 exit(1);
		}
	 }

	}
		else /* Start up during filling */
	{
	 hP_ID[0] =-2;
	}

}
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/*    Returns the new positions to the host to render per time-step/frame    */
/*---------------------------------------------------------------------------*/
void Device_Get_System_State( float3 *h_pos, Quaterion *h_Ornt, float3 *h_vel,
		                      float3 *h_Avel, uint *h_type,
		                      int *h_pid, int *Number_Particles )
{

	cudaDeviceSynchronize();/* Make Sure Nothing running of Device */
	
	if(num_dead>0)
	{
		if(!Device_use_multi_gpu)
		{
		  Sync_Killed();
		  cudaDeviceSynchronize();
		}
		else
		{
			Sync_KilledMulti();
		    cudaDeviceSynchronize();
		}
	}


	Number_Particles[0] =  NUMPARTICLES_Current;

	if( NUMPARTICLES_Current > 0 )
	{

      if(!Device_use_multi_gpu)
	  {
		cudaMemcpy( h_pos,        dram_position_com        ,sizeof(float3)*NUMPARTICLES_Current,
			cudaMemcpyDeviceToHost);

		cudaMemcpy( h_type,      dram_ObjectType      ,sizeof(uint)*NUMPARTICLES_Current,
    		                                     cudaMemcpyDeviceToHost );

		cudaMemcpy( h_pid,       dram_P_ID         ,sizeof(int)*NUMPARTICLES_Current,
    		                                     cudaMemcpyDeviceToHost);

		cudaMemcpy( h_vel,  dram_velocity_com  ,sizeof(float3)*NUMPARTICLES_Current,
		cudaMemcpyDeviceToHost  );


		cudaMemcpy( h_Avel, dram_velocity_ang, sizeof(float3)*NUMPARTICLES_Current,
			                                          cudaMemcpyDeviceToHost  );

        if(SimInfo_C->particle_type==polyhedra || SimInfo_C->sphere_orient)
	    {
		  cudaMemcpy( h_Ornt, dram_position_ornt ,sizeof(Quaterion)*NUMPARTICLES_Current,
    		                                     cudaMemcpyDeviceToHost );
	    }
	  }
	  else
	  {
        cudaMemcpy( h_pos,        dram_position_comM[0]        ,sizeof(float3)*NUMPARTICLES_Current,
		cudaMemcpyDeviceToHost);

		cudaMemcpy( h_type,      dram_ObjectTypeM[0]      ,sizeof(uint)*NUMPARTICLES_Current,
    		                                     cudaMemcpyDeviceToHost );

		cudaMemcpy( h_pid,       dram_P_IDM[0]         ,sizeof(int)*NUMPARTICLES_Current,
    		                                     cudaMemcpyDeviceToHost);

		cudaMemcpy( h_vel,       dram_velocity_comM[0]  ,sizeof(float3)*NUMPARTICLES_Current,
		cudaMemcpyDeviceToHost  );

		cudaMemcpy( h_Avel,      dram_velocity_angM[0] ,sizeof(float3)*NUMPARTICLES_Current,
			                                          cudaMemcpyDeviceToHost  );

        if(SimInfo_C->particle_type==polyhedra || SimInfo_C->sphere_orient)
	    {
		  cudaMemcpy( h_Ornt, dram_position_orntM[0] ,sizeof(Quaterion)*NUMPARTICLES_Current,
    		                                     cudaMemcpyDeviceToHost );
	    }
	  
	  }

	  cudaDeviceSynchronize();
	    if(error_check==1)
	    {
		   cudaError_t errormsg=cudaGetLastError();
		   if(errormsg>0)
		   {
		     cout<<"Get system end : Error End "<<cudaGetErrorString(errormsg)<<endl;
		     exit(1);
		   }
	    }

	 }
	 else /* Start up during filling */
	 {
	   h_pid[0] =-2;
	 }

}
/*---------------------------------------------------------------------------*/

void Device_Get_P_Velocity( float3 *H_Vel, float3 *H_Rvel, int Flag,
                            uint *h_ptype, int *h_pid, int *Number_Particles )
{

	cudaDeviceSynchronize();/* Make Sure Nothing running of Device */
	if( num_dead > 0 )
	{
		if(!Device_use_multi_gpu)
		{
		  Sync_Killed();
		  cudaDeviceSynchronize();
		}
		else
		{
			Sync_KilledMulti();
		    cudaDeviceSynchronize();
		}
	}

	Number_Particles[0] =  NUMPARTICLES_Current;

	 if(error_check==1)
	 {
		cudaError_t errormsg=cudaGetLastError();
		if(errormsg>0)
		{
		 cout<<"Get Velocity : Error Start "<<cudaGetErrorString(errormsg)<<endl;
		 exit(1);
		}
	 }
	 
	 if(!Device_use_multi_gpu)
	 {
	    cudaMemcpyAsync( H_Vel,  dram_velocity_com  ,sizeof(float3)*NUMPARTICLES_Current,
		cudaMemcpyDeviceToHost  );

	    if(Flag==1)
	    {
			cudaMemcpyAsync( H_Rvel, dram_velocity_ang ,sizeof(float3)*NUMPARTICLES_Current,
			                                          cudaMemcpyDeviceToHost  );
	    }

	    cudaMemcpyAsync( h_ptype,      dram_ObjectType      ,sizeof(uint)*NUMPARTICLES_Current,
    		                                     cudaMemcpyDeviceToHost  );

	     cudaMemcpyAsync( h_pid,       dram_P_ID        ,sizeof(int)*NUMPARTICLES_Current,
    		                                     cudaMemcpyDeviceToHost  );
	 }
	 else
	 {
        cudaMemcpyAsync( H_Vel,  dram_velocity_comM[0]  ,sizeof(float3)*NUMPARTICLES_Current,
		cudaMemcpyDeviceToHost  );

	    if(Flag==1)
	    {
			cudaMemcpyAsync( H_Rvel, dram_velocity_angM[0] ,sizeof(float3)*NUMPARTICLES_Current,
			                                          cudaMemcpyDeviceToHost );
	    }

	    cudaMemcpyAsync( h_ptype,      dram_ObjectTypeM[0]      ,sizeof(uint)*NUMPARTICLES_Current,
    		                                     cudaMemcpyDeviceToHost  );

	    cudaMemcpyAsync( h_pid,       dram_P_IDM[0]        ,sizeof(int)*NUMPARTICLES_Current,
    		                                     cudaMemcpyDeviceToHost  );
	 }

	 cudaDeviceSynchronize();

	 if(error_check==1)
	 {
		cudaError_t errormsg=cudaGetLastError();
		if(errormsg>0)
		{
		 cout<<"Get Velocity : Error End "<<cudaGetErrorString(errormsg)<<endl;
		 exit(1);
		}
	 }
}





void Device_Get_W_Forces(float3 *WallForces,float3 *WallContact_Point)
{
		  if (Device_use_multi_gpu)
      {
        int offset=0;
	    for(int i=0;i<num_devices_use;i++)
	    {
		  cudaSetDevice(i);
		  cudaDeviceSynchronize();

		  cudaMemcpyAsync( WallContact_Point+offset,   dWallContactPointsM[i]        ,sizeof(float3)*num_particles_perGPU[i],
		              		                                     cudaMemcpyDeviceToHost,streams[i*4+1] );


		  cudaMemcpyAsync( WallForces+offset,   dram_force_com_WallM[i]        ,sizeof(float3)*num_particles_perGPU[i],
			  cudaMemcpyDeviceToHost,streams[i*4]);


		  offset += num_particles_perGPU[i];
	     }

		for(int i=0;i<num_devices_use;i++)
	    {
		  cudaSetDevice(i);
		  cudaDeviceSynchronize();
		}
	  }
	  else
	  {
	 cudaDeviceSynchronize();
	 cudaMemcpy( WallForces,   dram_force_com_Wall        ,sizeof(float3)*NUMPARTICLES_Current,
		  cudaMemcpyDeviceToHost);

      cudaMemcpy( WallContact_Point,   dWallContactPoints        ,sizeof(float3)*NUMPARTICLES_Current,
            		                                     cudaMemcpyDeviceToHost );
	  
	  cudaDeviceSynchronize();

     }
	      if(error_check==1)
    {
    	cudaError_t errormsg=cudaGetLastError();
    	 if(errormsg>0)
    	 {
    	   cout<<"Get Wforces: "<<cudaGetErrorString(errormsg)<<endl;
    	   exit(1);
    	 }
    }

}



void Device_Get_P_Forces(uint *PnumContacts,float3 *PContact_Point)
{

	  if (Device_use_multi_gpu)
      {
        int offset=0;
	    for(int i=0;i<num_devices_use;i++)
	    {
		  cudaSetDevice(i);
		  cudaDeviceSynchronize();

		  cudaMemcpyAsync( PnumContacts+offset,   dNumContactsM[i]        ,sizeof(uint)*num_particles_perGPU[i],
			  cudaMemcpyDeviceToHost,streams[i*4]);

		  cudaMemcpyAsync( PContact_Point+offset,   dPContactPointsM[i]        ,sizeof(float3)*num_particles_perGPU[i]*64,
            		                                     cudaMemcpyDeviceToHost,streams[i*4+1] );
	  
		  offset += num_particles_perGPU[i];
	     }

		for(int i=0;i<num_devices_use;i++)
	    {
		  cudaSetDevice(i);
		  cudaDeviceSynchronize();
		}
	  }
	  else
	  {
	 cudaDeviceSynchronize();

	 cudaMemcpy( PnumContacts,   dNumContacts        ,sizeof(uint)*NUMPARTICLES_Current,
		  cudaMemcpyDeviceToHost);

	 cudaMemcpy( PContact_Point,   dPContactPoints        ,sizeof(float3)*NUMPARTICLES_Current*64,
            		                                     cudaMemcpyDeviceToHost );
	  
	  cudaDeviceSynchronize();
	  }
	      if(error_check==1)
    {
    	cudaError_t errormsg=cudaGetLastError();
    	 if(errormsg>0)
    	 {
    	   cout<<"Get Pforces: "<<cudaGetErrorString(errormsg)<<endl;
    	   exit(1);
    	 }
    }


	 
}


void Device_Get_V_Forces(float3 *VForces,float3 *VContact_Point)
{

	  if (Device_use_multi_gpu)
      {
        int offset=0;
	    for(int i=0;i<num_devices_use;i++)
	    {
		  cudaSetDevice(i);
		  cudaDeviceSynchronize();
		  cudaMemcpyAsync( VForces+offset,   dram_force_com_LifterM[i]        ,sizeof(float3)*num_particles_perGPU[i],
			  cudaMemcpyDeviceToHost,streams[i*4]);

		  cudaMemcpyAsync( VContact_Point+offset,   dVOContactPointsM[i]        ,sizeof(float3)*num_particles_perGPU[i],
            		                                     cudaMemcpyDeviceToHost,streams[i*4+1] );
	  
		  offset += num_particles_perGPU[i];
	     }

		for(int i=0;i<num_devices_use;i++)
	    {
		  cudaSetDevice(i);
		  cudaDeviceSynchronize();
		}
	  }
	  else
	  {
	    cudaDeviceSynchronize();
	    cudaMemcpy( VForces,   dram_force_com_Lifter        ,sizeof(float3)*NUMPARTICLES_Current,
		  cudaMemcpyDeviceToHost);

         cudaMemcpy( VContact_Point,   dVOContactPoints        ,sizeof(float3)*NUMPARTICLES_Current,
            		                                     cudaMemcpyDeviceToHost );
	  
	     cudaDeviceSynchronize();
	   }

    if(error_check==1)
    {
    	cudaError_t errormsg=cudaGetLastError();
    	 if(errormsg>0)
    	 {
    	   cout<<"Get Vobject: "<<cudaGetErrorString(errormsg)<<endl;
    	   exit(1);
    	 }
    }
	 
}


void Start_Bin()
{
	start_bin=true;
}


void Device_Reset_Energy(int Num)
{

	ResetEnergy<<< 1,1 >>>(Num);
	cudaDeviceSynchronize();
}


void Device_Get_Energy(float *Energy)
{
	float *d_energy;

	cudaMalloc( (void**) &d_energy  , sizeof(float)*NumPhysicsTallies);
	cudaDeviceSynchronize();

	GetEnergy<<< 1,1 >>>(d_energy,NumPhysicsTallies);
	cudaDeviceSynchronize();

	cudaMemcpy( Energy,     d_energy     ,sizeof(float)*NumPhysicsTallies,
  		                                     cudaMemcpyDeviceToHost );
    cudaDeviceSynchronize();
    cudaFree( d_energy );

    if(error_check==1)
    {
    	cudaError_t errormsg=cudaGetLastError();
    	 if(errormsg>0)
    	 {
    	   cout<<"Get Energy: "<<cudaGetErrorString(errormsg)<<endl;
    	   exit(1);
    	 }
    }
}



/*---------------------------------------------------------------------------*/
/*                    Passes updated Dynamic object List                     */
/*---------------------------------------------------------------------------*/
void Device_Update_DObjectPostions( DOBJ* h_DynamicOBJ)
{
	if(!Device_use_multi_gpu)
	{
	  cudaMemcpyToSymbol( VolumeObject, h_DynamicOBJ,    sizeof(DOBJ)*NUMDYOBJECTS );
	  cudaDeviceSynchronize();
	}
	else
	{
	   for(int i=0;i<num_devices_use;i++)
	   {
		 cudaSetDevice(i);
		 cudaMemcpyToSymbol( VolumeObject, h_DynamicOBJ,    sizeof(DOBJ)*NUMDYOBJECTS );
	     cudaDeviceSynchronize();
	   }
	   cudaSetDevice(0);
	}
}
/*-----------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
/*                    Passes updated Dynamic object List                     */
/*---------------------------------------------------------------------------*/
void Device_Update_WorldObjectPostions( WOBJ *h_WorldOBJ)
{
	if(!Device_use_multi_gpu)
	{
	 cudaMemcpyToSymbol( WorldObject, h_WorldOBJ,    sizeof(WOBJ)*NUMWOBJECTS );
	 cudaDeviceSynchronize();
	}
	else
	{
	   for(int i=0;i<num_devices_use;i++)
	   {
		 cudaSetDevice(i);
		 cudaMemcpyToSymbol( WorldObject, h_WorldOBJ,    sizeof(WOBJ)*NUMWOBJECTS );
	     cudaDeviceSynchronize();
	   }
	   cudaSetDevice(0);
	}
}
/*-----------------------------------------------------------------------------*/



void drum_rotation(bool cylinder_rotation)
{
	is_cylinder_rotation = cylinder_rotation;
}




/*-----------------------------------------------------------------------------*/
            /* Passes a new object removing the last object */
/*-----------------------------------------------------------------------------*/
void Device_OpenHatch( WOBJ* h_WorldOBJ,int h_Num_WorldOBJ)
{

	if(!Device_use_multi_gpu)
	{
	cudaMemcpyToSymbol( WorldObject,    h_WorldOBJ,    sizeof(WOBJ)*h_Num_WorldOBJ );
	cudaDeviceSynchronize();
    Hatch<<< 1,1 >>>();
    cudaDeviceSynchronize();
	}
	else
	{
	   for(int i=0;i<num_devices_use;i++)
	   {
		  cudaSetDevice(i);
		  cudaMemcpyToSymbol( WorldObject,    h_WorldOBJ,    sizeof(WOBJ)*h_Num_WorldOBJ );
		  cudaDeviceSynchronize();
		  Hatch<<< 1,1 >>>();
		  cudaDeviceSynchronize();
		}
	   cudaSetDevice(0);
	}

	printf("Device Hatch Opened\n");

}
/*-----------------------------------------------------------------------------*/

void Device_AddForce( )
{

    AddForce<<< 1,1 >>>();
    cudaDeviceSynchronize();
}


void Device_SubForce( )
{

    SubForce<<< 1,1 >>>();
    cudaDeviceSynchronize();
}


void Device_ZeroForce( )
{

    ZeroForce<<< 1,1 >>>();
    cudaDeviceSynchronize();
}





/*-----------------------------------------------------------------------------*/
                      /* Frees allocated arrays*/
/*-----------------------------------------------------------------------------*/
void Device_Clean()
{

  if (!Device_use_multi_gpu)
  {
	if(!SimInfo_C->use_symmetry)
	{
	  cudaFree(dram_force_com_X);
	  cudaFree(dram_force_com_Y);
	  cudaFree(dram_force_com_Z);

	  cudaFree(dram_force_ang_X);
	  cudaFree(dram_force_ang_Y);
	  cudaFree(dram_force_ang_Z);
	}
	else
	{
	    cudaFree(dram_force_PP);
	    cudaFree(dram_force_PP_ang);
	}


	cudaFree(dram_force_com_Wall);
	cudaFree(dram_force_com_Lifter);

	cudaFree(dram_force_ang_Wall);
	cudaFree(dram_force_ang_Lifter);

	if(Get_WallPoints)
	{
	  cudaFree(dWallContactPoints);
	}

	if(get_vobject_forces)
	{
	  cudaFree(dVOContactPoints);
	}

	if(get_P_forces)
	{
	  cudaFree(dPContactPoints);
	  cudaFree(dNumContacts);
	}

	cudaFree(m_dCellStart);
	cudaFree(m_dCellEnd);

	cudaFree(dram_GridParticleHash);
	cudaFree(dram_GridParticleIndex);

	cudaFree(dram_Lifter_Contact_Hist);
	cudaFree(dram_ObjectType);
	cudaFree(dram_P_ID);
	cudaFree(dram_position_com);

	if(SimInfo_C->particle_type==polyhedra || SimInfo_C->sphere_orient)
    {
	 cudaFree(dram_position_ornt);
    }
	cudaFree(dram_velocity_ang);
	cudaFree(dram_velocity_com);

	if(SimInfo_C->particle_type==polyhedra)
	{
	cudaFree(dBroad_List);
    cudaFree(dNumNN);
	}

  }

	  /* Free reamaining temp arrays*/
      /* Free All but NNList  */
    cudaFree(dSortedPos);

	if(SimInfo_C->particle_type==polyhedra || SimInfo_C->sphere_orient)
	{
      cudaFree(dSortedPosQ);
	}

    cudaFree(dSortedVel);
    cudaFree(dSorted_velocity_ang);
    cudaFree(dSortedPType);
    cudaFree(dSortedPID);
    cudaFree(dSortedLifter_Contact_Hist);



  if (Device_use_multi_gpu)
  {
	for(int i=0;i<num_devices_use;i++)
	{
		cudaSetDevice(i);
       for (int j = 0; j < 4; j++)
	   {
		   cudaStreamDestroy(streams[i*4+j]);
		   cudaEventDestroy(events[i*4+j]);
	   }

    if(!SimInfo_C->use_symmetry)
	{
	  cudaFree(dram_force_com_XM[i]);
	  cudaFree(dram_force_com_YM[i]);
	  cudaFree(dram_force_com_ZM[i]);

	  cudaFree(dram_force_ang_XM[i]);
	  cudaFree(dram_force_ang_YM[i]);
	  cudaFree(dram_force_ang_ZM[i]);
	}
	else
	{
	    cudaFree(dram_force_PPM[i]);
	    cudaFree(dram_force_PP_angM[i]);
	}


	cudaFree(dram_force_com_WallM[i]);
	cudaFree(dram_force_com_LifterM[i]);

	cudaFree(dram_force_ang_WallM[i]);
	cudaFree(dram_force_ang_LifterM[i]);

	if(Get_WallPoints)
	{
	  cudaFree(dWallContactPointsM[i]);
	}

	if(get_vobject_forces)
	{
	  cudaFree(dVOContactPointsM[i]);
	}

	if(get_P_forces)
	{
	  cudaFree(dPContactPointsM[i]);
	  cudaFree(dNumContactsM[i]);
	}

	cudaFree(m_dCellStartM[i]);
	cudaFree(m_dCellEndM[i]);

	cudaFree(dram_GridParticleHashM[i]);
	cudaFree(dram_GridParticleIndexM[i]);

	cudaFree(dram_Lifter_Contact_HistM[i]);
	cudaFree(dram_ObjectTypeM[i]);
	cudaFree(dram_P_IDM[i]);
	cudaFree(dram_position_comM[i]);

	if(SimInfo_C->particle_type==polyhedra || SimInfo_C->sphere_orient)
    {
	 cudaFree(dram_position_orntM[i]);
    }
	cudaFree(dram_velocity_angM[i]);
	cudaFree(dram_velocity_comM[i]);


  cudaFree(dBroad_ListM[i]);
  cudaFree(dNumNNM[i]);
  }

	cudaDeviceSynchronize();
	cudaDeviceReset();

}

  	cudaDeviceSynchronize();
	cudaDeviceReset();
	printf(" Freed GPU Memory!\n");
}

//Macro for checking cuda errors following a cuda launch or api call
#define cudaCheckError() {                                          \
        cudaError_t e=cudaGetLastError();                                 \
        if(e!=cudaSuccess) {                                              \
            printf("Cuda failure %s:%d: '%s'\n",__FILE__,__LINE__,cudaGetErrorString(e));           \
            exit(EXIT_FAILURE);                                           \
        }                                                                 \
    }
/* Start og MultiGPU */

void outputBandwidthMatrix(int numGPUs)
{
    int numElems=10000000;
    int repeat=5;
    vector<int *> buffers(numGPUs);
    vector<cudaEvent_t> start(numGPUs);
    vector<cudaEvent_t> stop(numGPUs);

    for (int d=0; d<numGPUs; d++)
    {
        cudaSetDevice(d);
        cudaMalloc(&buffers[d],numElems*sizeof(int));
        cudaCheckError();
        cudaEventCreate(&start[d]);
        cudaCheckError();
        cudaEventCreate(&stop[d]);
        cudaCheckError();
    }

    vector<double> bandwidthMatrix(numGPUs*numGPUs);

    for (int i=0; i<numGPUs; i++)
    {
        cudaSetDevice(i);

        for (int j=0; j<numGPUs; j++)
        {

            cudaDeviceSynchronize();
            cudaCheckError();
            cudaEventRecord(start[i]);

            for (int r=0; r<repeat; r++)
            {
                cudaMemcpyPeerAsync(buffers[i],i,buffers[j],j,sizeof(int)*numElems);
            }

            cudaEventRecord(stop[i]);
            cudaDeviceSynchronize();
            cudaCheckError();

            float time_ms;
            cudaEventElapsedTime(&time_ms,start[i],stop[i]);
            double time_s=time_ms/1e3;

            double gb=numElems*sizeof(int)*repeat/(double)1e9;
            bandwidthMatrix[i*numGPUs+j]=gb/time_s;
        }
    }

    printf("   D\\D");

    for (int j=0; j<numGPUs; j++)
    {
        printf("%6d ", j);
    }

    printf("\n");

    for (int i=0; i<numGPUs; i++)
    {
        printf("%6d ",i);

        for (int j=0; j<numGPUs; j++)
        {
            printf("%6.02f ", bandwidthMatrix[i*numGPUs+j]);
        }

        printf("\n");
    }

    for (int d=0; d<numGPUs; d++)
    {
        cudaSetDevice(d);
        cudaFree(buffers[d]);
        cudaCheckError();
        cudaEventDestroy(start[d]);
        cudaCheckError();
        cudaEventDestroy(stop[d]);
        cudaCheckError();
    }
}

void outputBidirectionalBandwidthMatrix(int numGPUs)
{
    int numElems=10000000;
    int repeat=5;
    vector<int *> buffers(numGPUs);
    vector<cudaEvent_t> start(numGPUs);
    vector<cudaEvent_t> stop(numGPUs);
    vector<cudaStream_t> stream0(numGPUs);
    vector<cudaStream_t> stream1(numGPUs);

    for (int d=0; d<numGPUs; d++)
    {
        cudaSetDevice(d);
        cudaMalloc(&buffers[d],numElems*sizeof(int));
        cudaCheckError();
        cudaEventCreate(&start[d]);
        cudaCheckError();
        cudaEventCreate(&stop[d]);
        cudaCheckError();
        cudaStreamCreate(&stream0[d]);
        cudaCheckError();
        cudaStreamCreate(&stream1[d]);
        cudaCheckError();
    }

    vector<double> bandwidthMatrix(numGPUs*numGPUs);

    for (int i=0; i<numGPUs; i++)
    {
        cudaSetDevice(i);

        for (int j=0; j<numGPUs; j++)
        {

            cudaDeviceSynchronize();
            cudaCheckError();
            cudaEventRecord(start[i]);

            for (int r=0; r<repeat; r++)
            {
                cudaMemcpyPeerAsync(buffers[i],i,buffers[j],j,sizeof(int)*numElems,stream0[i]);
                cudaMemcpyPeerAsync(buffers[j],j,buffers[i],i,sizeof(int)*numElems,stream1[i]);
            }

            cudaEventRecord(stop[i]);
            cudaDeviceSynchronize();
            cudaCheckError();

            float time_ms;
            cudaEventElapsedTime(&time_ms,start[i],stop[i]);
            double time_s=time_ms/1e3;

            double gb=2.0*numElems*sizeof(int)*repeat/(double)1e9;
            bandwidthMatrix[i*numGPUs+j]=gb/time_s;
        }
    }

    printf("   D\\D");

    for (int j=0; j<numGPUs; j++)
    {
        printf("%6d ", j);
    }

    printf("\n");

    for (int i=0; i<numGPUs; i++)
    {
        printf("%6d ",i);

        for (int j=0; j<numGPUs; j++)
        {
            printf("%6.02f ", bandwidthMatrix[i*numGPUs+j]);
        }

        printf("\n");
    }

    for (int d=0; d<numGPUs; d++)
    {
        cudaSetDevice(d);
        cudaFree(buffers[d]);
        cudaCheckError();
        cudaEventDestroy(start[d]);
        cudaCheckError();
        cudaEventDestroy(stop[d]);
        cudaCheckError();
        cudaStreamDestroy(stream0[d]);
        cudaCheckError();
        cudaStreamDestroy(stream1[d]);
        cudaCheckError();
    }
}

void outputLatencyMatrix(int numGPUs)
{
    int repeat=10000;
    vector<int *> buffers(numGPUs);
    vector<cudaEvent_t> start(numGPUs);
    vector<cudaEvent_t> stop(numGPUs);

    for (int d=0; d<numGPUs; d++)
    {
        cudaSetDevice(d);
        cudaMalloc(&buffers[d],1);
        cudaCheckError();
        cudaEventCreate(&start[d]);
        cudaCheckError();
        cudaEventCreate(&stop[d]);
        cudaCheckError();
    }

    vector<double> latencyMatrix(numGPUs*numGPUs);

    for (int i=0; i<numGPUs; i++)
    {
        cudaSetDevice(i);

        for (int j=0; j<numGPUs; j++)
        {

            cudaDeviceSynchronize();
            cudaCheckError();
            cudaEventRecord(start[i]);

            for (int r=0; r<repeat; r++)
            {
                cudaMemcpyPeerAsync(buffers[i],i,buffers[j],j,1);
            }

            cudaEventRecord(stop[i]);
            cudaDeviceSynchronize();
            cudaCheckError();

            float time_ms;
            cudaEventElapsedTime(&time_ms,start[i],stop[i]);

            latencyMatrix[i*numGPUs+j]=time_ms*1e3/repeat;
        }
    }

    printf("   D\\D");

    for (int j=0; j<numGPUs; j++)
    {
        printf("%6d ", j);
    }

    printf("\n");

    for (int i=0; i<numGPUs; i++)
    {
        printf("%6d ",i);

        for (int j=0; j<numGPUs; j++)
        {
            printf("%6.02f ", latencyMatrix[i*numGPUs+j]);
        }

        printf("\n");
    }

    for (int d=0; d<numGPUs; d++)
    {
        cudaSetDevice(d);
        cudaFree(buffers[d]);
        cudaCheckError();
        cudaEventDestroy(start[d]);
        cudaCheckError();
        cudaEventDestroy(stop[d]);
        cudaCheckError();
    }
}



/*---------------------------------------------------------------------------*/
   /* 1) ENTRY POINT Copies all data to the device sets kernel parameters*/
/*---------------------------------------------------------------------------*/
void Device_Set_SimulationData_MultiGPU( WOBJ*           h_WorldOBJ,
								         POBJ*           h_ParticleOBJ,
								         DOBJ*           h_DynamicOBJ,
								         SimulationInfo* h_SimInfo,
								         int             h_Num_WorldOBJ,
								         int             h_Num_ParticleOBJ,
								         int             h_Num_DynamicOBJ,
								         InitConfig      *h_PosConfig             )

{

	m_WorldObject = h_WorldOBJ;

	error_check = 1;
	LogFileD.open( "../SimulationDevice.Log" );

	cout<<"INFO-D: Starting Device Configuration \n";
	 cout<<" \n";
	 cout<<" \n";

	device_num = h_PosConfig->use_device;
	/* Create a local copy of the SimOBJ */
	SimInfo_C = h_SimInfo;
	Device_use_multi_gpu = h_PosConfig->multi_gpu;

	num_devices_use = h_PosConfig->num_gpus;
	Get_DeviceInfo();

	cudaSetDevice(0);
	cudaDeviceSynchronize();

	printf( " num devices %d \n",num_devices_use);

	cout<<" Using P2P copy "<<endl;

	    for (int i=0; i<num_devices_use; i++)
	    {
	        cudaSetDevice(i);

	        for (int j=0; j<num_devices_use; j++)
	        {
	            int access;
	            if (i!=j)
	            {
	                cudaDeviceCanAccessPeer(&access,i,j);
	                printf("Device=%d %s Access Peer Device=%d\n", i, access ? "CAN" : "CANNOT", j);
	            }
	        }
	    }


	    for (int i=0; i<num_devices_use; i++)
	    {
	        cudaSetDevice(i);

	        for (int j=0; j<num_devices_use; j++)
	        {
	            int access;
	            cudaDeviceCanAccessPeer(&access,i,j);

	            if (access)
	            {
	                cudaDeviceEnablePeerAccess(j,0);
	       		 if(error_check==1)
	       		 {
	       			cudaError_t errormsg=cudaGetLastError();
	       			if(errormsg>0)
	       			{
	       			 cout<<"Set Simulation Data: PeerAccess "<<cudaGetErrorString(errormsg)<<endl;
	       			 exit(1);
	       			}
	       		 }
	            }
	        }
	    }

    int numGPUs;
    cudaGetDeviceCount(&numGPUs);

    //compute cliques
    vector<vector<int> > cliques;

    vector<bool> added(numGPUs,false);

    for (int i=0; i<numGPUs; i++)
    {
        if (added[i]==true)
            continue;         //already processed

        //create new clique with i
        vector<int> clique;
        added[i]=true;
        clique.push_back(i);

        for (int j=i+1; j<numGPUs; j++)
        {
            int access;
            cudaDeviceCanAccessPeer(&access,i,j);

            if (access)
            {
                clique.push_back(j);
                added[j]=true;
            }
        }

        cliques.push_back(clique);
    }

    printf("P2P Cliques: \n");

    for (int c=0; c<(int)cliques.size(); c++)
    {
        printf("[");

        for (int j=0; j<(int)cliques[c].size()-1; j++)
        {
            printf("%d ",cliques[c][j]);
        }

        printf("%d]\n",cliques[c][cliques[c].size()-1]);
    }

    printf("Unidirectional P2P=Enabled Bandwidth Matrix (GB/s)\n");
    outputBandwidthMatrix(numGPUs);

	printf("Bidirectional P2P=Enabled Bandwidth Matrix (GB/s)\n");
    outputBidirectionalBandwidthMatrix(numGPUs);


    printf("P2P=Enabled Latency Matrix (us)\n");
    outputLatencyMatrix(numGPUs);



    cudaSetDevice(0);

	NUMPARTICLES_MAX_SIM = h_SimInfo->Num_Particles;
	NUMDYOBJECTS         = h_Num_DynamicOBJ;
	NUMWOBJECTS          = h_Num_WorldOBJ;
	
	D_NUMPARTICLES       = NUMPARTICLES_MAX_SIM;/* Used for killing Silo */
	NUMPARTICLES_Current = NUMPARTICLES_MAX_SIM;/* Used for filling */

	cout<<" \n";
	 cout<<" \n";

	cout<<"INFO-D: Total Particles = "<<NUMPARTICLES_MAX_SIM<<endl;
	cout<<"INFO-D: Size Types = "<<h_SimInfo->Num_ParticleObjects<<endl;

	for( int i=0; i< h_SimInfo->Num_ParticleObjects;i++ )
	{
		cout<<"INFO-D: Num particles per type = "<<h_SimInfo->Num_ParticlesPerObject[i]<<endl;
	}


	/* Set thread info on the GPU */
    printf( "RUNLEVEL 0: Max SM occupy \n");

    int num_particles_SM = (int)ceil(NUMPARTICLES_MAX_SIM/(float)num_sm);

    /* Now get the number of blocks per SM */
    int num_blocks_SM =  (int)ceil(num_particles_SM/(float)h_PosConfig->threads_perBlock);
    printf(" blocks per SM %d \n",num_blocks_SM);

    threads_perBlock = h_PosConfig->threads_perBlock;
    int num_threads_block=h_PosConfig->threads_perBlock;

    /* block size is too big */
    if(num_particles_SM < h_PosConfig->threads_perBlock)
    {
	   num_threads_block = num_particles_SM; /* Single block is sufficient */
    }

    GlobalLaunch.dimBlock = make_uint3( num_threads_block,1,1);
    GlobalLaunch.dimGrid  = make_uint3( num_blocks_SM*num_sm,1,1 );

	printf("INFO-D: Total %d Blocks with %d Threads: Total threads = %d \n",GlobalLaunch.dimGrid.x,
								 GlobalLaunch.dimBlock.x,GlobalLaunch.dimGrid.x*GlobalLaunch.dimBlock.x);
   

	

	/* Start Multi GPU allocations */

		printf(" Multi-GPU Performance Mode Splitting over %d GPUS \n",num_devices_use);

		int num_perGPU = (int)floor(NUMPARTICLES_MAX_SIM/(float)num_devices_use);
	    
		int num_rem = NUMPARTICLES_MAX_SIM - num_perGPU*num_devices_use;

	
		for ( int i=0; i < num_devices_use; i++ )
		{
			num_particles_perGPU[i] = num_perGPU;
	
			if(i==0)
			{
					num_particles_perGPU[0] = num_perGPU + num_rem;

			}
		
		    num_particles_SM = (int)ceil(num_particles_perGPU[i]/(float)num_sm);

	        /* Now get the number of blocks per SM */
	        num_blocks_SM      =  (int)ceil(num_particles_SM/(float)h_PosConfig->threads_perBlock);

	        threads_perBlock      = h_PosConfig->threads_perBlock;
	        num_threads_block = h_PosConfig->threads_perBlock;

	        /* block size is too big */
	        if( num_particles_SM < h_PosConfig->threads_perBlock )
	        {
		      num_threads_block = num_particles_SM; /* Single block is sufficient */
	        }

	        multi[i].dimBlock = make_uint3( num_threads_block,1,1);
	        multi[i].dimGrid  = make_uint3( num_blocks_SM*num_sm,1,1 );

	        printf(" Particles %d per GPU threads = %d \n",num_particles_perGPU[i],multi[i].dimGrid.x*multi[i].dimBlock.x);

		}
		
		 cout<<" Starting GPU Malloc \n";

		int mem_loc=0;
		/* Sets information on each GPU */

	    m_numGridCells = SimInfo_C->num_NNCells.x*SimInfo_C->num_NNCells.y
											 *SimInfo_C->num_NNCells.z;

	    printf(" num grid cells %d \n",m_numGridCells);


        for (int i=0;i<num_devices_use;i++)
        {
        	cout<<" \n";

	        cudaSetDevice(i);
            cudaDeviceSynchronize();

		     if(error_check==1)
			 {
				cudaError_t errormsg=cudaGetLastError();
				if(errormsg>0)
				{
				 cout<<"Set Simulation Data: Start Dram memory alloc "<<cudaGetErrorString(errormsg)<<" Device "<<i<<endl;
				 exit(1);
				}
			 }


	        cout<<"INFO-D: Allocating Dynamic Arrays GPU "<<i<<endl;



			/* Primary arrays for particle information */
		    cudaMalloc( (void**) &dram_ObjectTypeM[i]  ,  sizeof(uint)*NUMPARTICLES_MAX_SIM );
			cudaMalloc( (void**) &dram_P_IDM      [i]  ,  sizeof(int)*NUMPARTICLES_MAX_SIM );

			cudaMalloc( (void**) &dram_position_comM[i],  sizeof(float3)*NUMPARTICLES_MAX_SIM );
			cudaMalloc( (void**) &dram_velocity_comM[i],  sizeof(float3)*NUMPARTICLES_MAX_SIM );
			cudaMalloc( (void**) &dram_velocity_angM[i],  sizeof(float3)*NUMPARTICLES_MAX_SIM );

			if (SimInfo_C->particle_type==polyhedra || SimInfo_C->sphere_orient)
			{
			  cudaMalloc( (void**) &dram_position_orntM[i], sizeof(Quaterion)*NUMPARTICLES_MAX_SIM );
			}

			cudaMalloc( (void**) &dram_Lifter_Contact_HistM[i], sizeof(Contact_Info)*NUMPARTICLES_MAX_SIM ); 
			
			if(SimInfo_C->particle_type==polyhedra)
			{
			cudaMalloc( (void**) &dBroad_ListM[i], sizeof(uint)*num_particles_perGPU[i]*32 );
	        cudaMalloc( (void**) &dNumNNM[i]     , sizeof(uint)*num_particles_perGPU[i]    );


			printf(" Size of Primary Arrays %f MB \n", ( NUMPARTICLES_MAX_SIM*4*2 + NUMPARTICLES_MAX_SIM*12*3 +
				NUMPARTICLES_MAX_SIM*16 + num_particles_perGPU[i]*32*4 + num_particles_perGPU[i]*4 )*1E-6 ); 
			}
			else
			{
              			printf(" Size of Primary Arrays %f MB \n", ( NUMPARTICLES_MAX_SIM*4*2 + NUMPARTICLES_MAX_SIM*12*3 +
				NUMPARTICLES_MAX_SIM*16  )*1E-6 ); 
			}
		     if(error_check==1)
			 {
				cudaError_t errormsg=cudaGetLastError();
				if(errormsg>0)
				{
				 cout<<"Set Simulation Data: Prime Dram memory alloc "<<cudaGetErrorString(errormsg)<<" GPU "<<i<<endl;
				 exit(1);
				}
			 }

			/* Secondary arrays allocation based on particles per GPU */

						/* Allocation of cell data, changes with particle distn */
			cudaMalloc( (void**) &m_dCellStartM[i],           sizeof(uint)*m_numGridCells);
			cudaMalloc( (void**) &m_dCellEndM[i],             sizeof(uint)*m_numGridCells);
			cudaMalloc( (void**) &dram_GridParticleHashM[i],  sizeof(uint)*NUMPARTICLES_MAX_SIM         );
			cudaMalloc( (void**) &dram_GridParticleIndexM[i], sizeof(int)*NUMPARTICLES_MAX_SIM          );



			printf(" Size of Cell %f MB \n", ( m_numGridCells*4*2.0 + NUMPARTICLES_MAX_SIM*4*2.0 )*1E-6);

			/*Arrays for Force info */
			if(SimInfo_C->use_symmetry)
			{
			  cudaMalloc( (void**) &dram_force_com_XM[i] ,  sizeof(float)*num_particles_perGPU[i] );
			  cudaMalloc( (void**) &dram_force_com_YM[i] ,  sizeof(float)*num_particles_perGPU[i] );
			  cudaMalloc( (void**) &dram_force_com_ZM[i] ,  sizeof(float)*num_particles_perGPU[i] );

			  cudaMalloc( (void**) &dram_force_ang_XM[i] ,  sizeof(float)*num_particles_perGPU[i] );
			  cudaMalloc( (void**) &dram_force_ang_YM[i] ,  sizeof(float)*num_particles_perGPU[i] );
			  cudaMalloc( (void**) &dram_force_ang_ZM[i] ,  sizeof(float)*num_particles_perGPU[i] );
			}
			else
			{
			  cudaMalloc( (void**) &dram_force_PPM[i]       , sizeof(float3)*num_particles_perGPU[i] );
			  cudaMalloc( (void**) &dram_force_PP_angM[i]   , sizeof(float3)*num_particles_perGPU[i] );
			}

			cudaMalloc( (void**) &dram_force_com_WallM[i]   , sizeof(float3)*num_particles_perGPU[i] );
			cudaMalloc( (void**) &dram_force_com_LifterM[i] , sizeof(float3)*num_particles_perGPU[i] );

			cudaMalloc( (void**) &dram_force_ang_WallM[i]   , sizeof(float3)*num_particles_perGPU[i] );
			cudaMalloc( (void**) &dram_force_ang_LifterM[i] , sizeof(float3)*num_particles_perGPU[i] );

			printf(" Size of Forces %f MB \n", ( num_particles_perGPU[i]*4*6.0 )*1E-6);


	
	get_wall_forces    = h_PosConfig->get_wall_points;
	get_vobject_forces = h_PosConfig->get_vobject_points;
	get_P_forces = h_PosConfig->get_particle_points;



	if(get_wall_forces)
	{
		printf(" devvice allco W\n");
     cudaMalloc( (void**) &dWallContactPointsM[i], sizeof(float3)*num_particles_perGPU[i] );
	}

	if(get_vobject_forces)
	{
		printf(" devvice allco V\n");
     cudaMalloc( (void**) &dVOContactPointsM[i], sizeof(float3)*num_particles_perGPU[i] );
	}

	if(get_P_forces)
	{
	 printf(" devvice allco P %f \n",12*sizeof(float3)*num_particles_perGPU[i]*32.0*1E-6);
     cudaMalloc( (void**) &dPContactPointsM[i], sizeof(float3)*num_particles_perGPU[i]*64 );
     cudaMalloc( (void**) &dNumContactsM[i], sizeof(uint)*num_particles_perGPU[i]);
	}

			/* End of secondary array allocation */
      

			float Tallys [6];

			for(int j=0;j<6;j++)
			{
			  Tallys[j]=0.0f;
			}

			float *d_Tally;
			cudaMalloc((void **)&d_Tally, sizeof(float)*6);
			cudaDeviceSynchronize();

			cudaMemcpy(d_Tally,Tallys, sizeof(float)*6,          cudaMemcpyHostToDevice );
			cudaDeviceSynchronize();

			Set_Tallys<<< 1,1 >>>(d_Tally);
			cudaDeviceSynchronize();
			cudaFree(d_Tally);

			 if(error_check==1)
			 {
				cudaError_t errormsg=cudaGetLastError();
				if(errormsg>0)
				{
				 cout<<"Set Simulation Data: Dram memory alloc "<<cudaGetErrorString(errormsg)<<" Device "<<i<<endl;
				 exit(1);
				}
			 }


			cout<<"INFO-D: Copying Constant Data Device: "<<i<<endl;

			/* 1) Copy Objects to Device */
			cudaMemcpyToSymbol( SimParms,       h_SimInfo,     sizeof( SimulationInfo ) );
			cudaMemcpyToSymbol( WorldObject,    h_WorldOBJ,    sizeof(WOBJ)*h_Num_WorldOBJ );
			cudaMemcpyToSymbol( ParticleObject, h_ParticleOBJ, sizeof(POBJ)*h_SimInfo->Num_ParticleObjects);


			if( h_Num_DynamicOBJ > 0 )
			{
				cudaMemcpyToSymbol( VolumeObject,  h_DynamicOBJ,  sizeof(DOBJ)*h_Num_DynamicOBJ);
			}

			if( !h_PosConfig->use_file )
			{
			  if( h_PosConfig->grid_type < 2 )
			  {
			    cudaMemcpyToSymbol( InitVel,        h_PosConfig->velocity,     sizeof(float3)*2);
			  }
			  else
			  {
			    h_PosConfig->velocity[0]= make_float3(h_PosConfig->launch_Vel.x,h_PosConfig->launch_Vel.y,h_PosConfig->launch_Vel.z);
			    cudaMemcpyToSymbol( InitVel,        h_PosConfig->velocity,     sizeof(float3)*2);
			  }
			  cudaDeviceSynchronize();
			}
			else
			{
				  h_PosConfig->velocity[0]= make_float3(0.0,0.0,0.0f);
				  cudaMemcpyToSymbol( InitVel,        h_PosConfig->velocity,     sizeof(float3)*2);
			}


			cout<<" GPU Malloc DONE! GPU "<<i<<endl;
			cudaDeviceSynchronize();


			cout<<" Launching Init Kernel GPU "<<i<<endl;
			/* Launch Kernel to set Initial parameters on  device */
			 Set_Init_Config_Device<<< multi[i].dimGrid, multi[i].dimBlock >>>
		   								  ( SimInfo_C->Num_ParticleObjects,
		   		                			SimInfo_C->Num_WorldObjects,
		   		                			NUMDYOBJECTS,

		   		                			dram_force_PPM[i],
		   		                			dram_force_com_WallM[i],dram_force_com_LifterM[i],
		   		                			dram_force_PP_angM[i],
		   		                			dram_force_ang_WallM[i], dram_force_ang_LifterM[i],

		   		                			dram_velocity_comM[i],
											dram_velocity_angM[i],
											dram_Lifter_Contact_HistM[i],
											num_particles_perGPU[i],
											mem_loc );
			 if(error_check==1)
			 {
				cudaError_t errormsg=cudaGetLastError();
				if(errormsg>0)
				{
				 cout<<"Set Simulation Data: Init_Config Kernel "<<cudaGetErrorString(errormsg)<<" Device "<<i<<endl;
				 exit(1);
				}
			 }

			cudaDeviceSynchronize();

			mem_loc += num_particles_perGPU[i];

			cout<<" GPU Init Done "<<i<<endl;

			/* Master GPU */
			if (i==0)
			{
				cout<<" GPU 0 Only Tasks ! "<<endl;
				                 /* Allocation of PDA data */
                cudaMalloc( (void**) &dSortedPos  , sizeof(float3)   *NUMPARTICLES_MAX_SIM );     /* 12bytes*N */
                cudaMalloc( (void**) &dSortedVel  , sizeof(float3)   *NUMPARTICLES_MAX_SIM );     /* 12bytes*N */

                if(SimInfo_C->particle_type==polyhedra || SimInfo_C->sphere_orient)
	            {
                  cudaMalloc( (void**) &dSortedPosQ , sizeof(Quaterion)*NUMPARTICLES_MAX_SIM ); /* 16bytes*N */
	            }

                 cudaMalloc( (void**) &dSorted_velocity_ang, sizeof(float3)   *NUMPARTICLES_MAX_SIM );  /* 12bytes*N */
   
				 cudaMalloc( (void**) &dSortedPType,         sizeof(uint)     *NUMPARTICLES_MAX_SIM );  /* 12bytes*N */
    
				 cudaMalloc( (void**) &dSortedPID,           sizeof(int)      *NUMPARTICLES_MAX_SIM );  /* 12bytes*N */
    
				 cudaMalloc( (void**) &dSortedLifter_Contact_Hist     , sizeof(Contact_Info)*NUMPARTICLES_MAX_SIM ); /* 12bytes*N */

				printf(" GPU 0 Sorted Arrays %f MB \n", ( NUMPARTICLES_MAX_SIM*4*2 + NUMPARTICLES_MAX_SIM*12*3 +
				NUMPARTICLES_MAX_SIM*16 + NUMPARTICLES_MAX_SIM*16 )*1E-6 ); 

				cout<<" Setting positions on GPU 0 only  "<<endl;
			     /* Select Initial Position */
			     if(!h_PosConfig->use_file)
			     {
			       if( h_PosConfig->grid_type==0 )
			       {
				     Set_Position_DefaultGrid(h_PosConfig,h_ParticleOBJ);
			       }
			       else if( h_PosConfig->grid_type==1 )
			       {
				     Set_Position_Random_Grid(h_PosConfig,h_ParticleOBJ);
			        }
			        else if( h_PosConfig->grid_type==2 )
			        {
				      Set_Position_Fill_Rectangle(h_ParticleOBJ,h_WorldOBJ,h_PosConfig->launch_Vel.y,
						  h_PosConfig->fill_plane_start,h_PosConfig->fill_plane_end );
				      isFilling = true;
				    }
			     }
			    
			}

			cudaDeviceSynchronize();

			cudaMemGetInfo(&a[i],&t[i]);
			cout<<" Clock (MHZ): "<<" Total Memory (MB): "
				  <<t[i]*9.53674e-7<<" Free (MB): "<<a[i]*9.53674e-7 <<endl;

			cout<<"INFO-D: All data copied to Device "<<i<<endl;
			cout<<"\n"<<endl;
			cout<<"\n"<<endl;
      
     } /* End Loop over devices */
        
	    cudaSetDevice(0);

        mem_loc=0;

    	for (int i=1;i<num_devices_use;i++)
    	{
          mem_loc += num_particles_perGPU[i-1];
    	  cudaMemcpyPeer( dram_velocity_comM[0]       + mem_loc,0, dram_velocity_comM[i]       + mem_loc, i,sizeof(float3)*num_particles_perGPU[i]);
    	  cudaMemcpyPeer( dram_velocity_angM[0]       + mem_loc,0, dram_velocity_angM[i]       + mem_loc, i,sizeof(float3)*num_particles_perGPU[i]);
		  cudaMemcpyPeer( dram_Lifter_Contact_HistM[0]+ mem_loc,0, dram_Lifter_Contact_HistM[i]+ mem_loc, i,sizeof(float3)*num_particles_perGPU[i]);
	      cudaDeviceSynchronize();
    	}
    	cudaDeviceSynchronize();


		/* MULTI_1: Overlap memory transactions */
	for (int i=0;i<num_devices_use;i++)
	{
	 cudaSetDevice(i);

	 cudaMemset(m_dCellStartM[i], 0xffffffff, m_numGridCells*sizeof(uint));
     cudaMemset(m_dCellEndM[i],   0xffffffff, m_numGridCells*sizeof(uint));

     if(error_check==1)
    	 {
    	    	cudaError_t errormsg=cudaGetLastError();
    	    	if(errormsg>0)
    	    	{
    	    	cout<<"Start: MEMSET "<<cudaGetErrorString(errormsg)<<Call_Time<<" device :"<<i<<endl;
    	    	 exit(1);
    	    	}
    	 }
    	/*-----------------------------------------------------------------------*/
	}

	cudaDeviceSynchronize();


      LogFileD.close();

	 if(error_check==1)
	 {
		cudaError_t errormsg=cudaGetLastError();
		if(errormsg>0)
		{
		 cout<<"Set Simulation Data: Error at end: Debug further "<<cudaGetErrorString(errormsg)<<endl;
		 exit(1);
		}
	 }

}
/*-----------------------------------------------------------------------------*/




/* Flags particles to be removed from the sim and gives them a position such that the hash will be large */
void Flag_DeadParticlesMulti()
{
	cudaSetDevice(0);
	 stepc++;
	       /* Dont check every step to save computations */
	       if( stepc >= 250 )
	       {
			   if(SimInfo_C->use_valid_zone)
	    	  {
	        	 Remove_PL_Volume_Particles<<< GlobalLaunch.dimGrid, GlobalLaunch.dimBlock >>>
	        			 (dram_position_comM[0],dram_P_IDM[0],dram_ObjectTypeM[0],NUMPARTICLES_Current);
	        	 cudaDeviceSynchronize();
	    	  }
	    	  else /* Silo */
	    	  {
	    	    Remove_Silo_Discharged_Particles<<< GlobalLaunch.dimGrid, GlobalLaunch.dimBlock >>>
	    	    		(dram_position_comM[0],dram_P_IDM[0],NUMPARTICLES_Current);
	    	    cudaDeviceSynchronize();
	    	  }

	    	  /* Copy number killed to Host (Check if can do this without memcpy */
	    	  int *d_kill;
	    	  int  h_kill[2];


	    	  float *d_stats;
	    	  float  h_stats[2];

	    	  cudaMalloc( (void**) &d_kill        , sizeof(int)*2   );
	    	  cudaMalloc( (void**) &d_stats        , sizeof(float)*2   );

	    	  Kill_copy<<<1,1>>>(d_kill,d_stats);
	    	  cudaDeviceSynchronize();

	    	  cudaMemcpy( h_kill,        d_kill        ,sizeof(int)*2,
	    	    		                                     cudaMemcpyDeviceToHost );

	    	  cudaMemcpy( h_stats,        d_stats        ,sizeof(float)*2,
	    	    		                                     cudaMemcpyDeviceToHost );
	    	  cudaDeviceSynchronize();


	    	  num_dead=h_kill[0];
	    	  cudaFree(d_kill);

	          dis_mass = h_stats[0];
	          dis_vol = h_stats[1];
	    	  cudaFree(d_stats);
	    	  stepc = 0;
	       }


		 if(error_check==1)
		 {
		    	cudaError_t errormsg=cudaGetLastError();
		    	if(errormsg>0)
		    	{
		    	 cout<<"DEM_UpdateSim: KillParticles "<<cudaGetErrorString(errormsg)<<endl;
		    	 exit(1);
		    	}
		 }


}


/* Adjusts the block size so that threads do not access dead particles */
void Remove_Flagged_ParticlesMulti()
{
	//printf(" NUM BEFORE %d \n",NUMPARTICLES_Current);
	NUMPARTICLES_Current -= num_dead;
	//printf(" Removing %d Particles current %d \n",num_dead,NUMPARTICLES_Current);

	//printf("sim type %d \n",SimInfo_C->Simulation_Type);
	if( SimInfo_C->Simulation_Type==Pulp_Lift)
	{
	  printf(" %f %d %f %f \n", Call_Time, num_dead, dis_mass, dis_vol);
	}
	num_dead=0;

   int num_particles_SM = (int)ceil(NUMPARTICLES_Current/(float)num_sm);

   /* Now get the number of blocks per SM */
  	int num_blocks_SM =  (int)ceil(num_particles_SM/(float)threads_perBlock);

  	   int num_threads_block=threads_perBlock;

  	   /* block size is too big */
  	   if(num_particles_SM < threads_perBlock)
  	   {
  		   num_threads_block = num_particles_SM; /* Single block is sufficient */
  	   }


  	   GlobalLaunch.dimBlock = make_uint3( num_threads_block,1,1);
  	   GlobalLaunch.dimGrid  = make_uint3( num_blocks_SM*num_sm,1,1 );

  	   SimInfo_C->Num_Particles = NUMPARTICLES_Current;

  	   for (int i=0;i<num_devices_use;i++)
   	   {
  		   cudaSetDevice(i);
  	     cudaMemcpyToSymbol( SimParms,       SimInfo_C,     sizeof( SimulationInfo ) );
  	     cudaDeviceSynchronize();
   	   }
  	 cudaSetDevice(0);

  	   int num_perGPU = (int)floor(NUMPARTICLES_Current/(float)num_devices_use);

	   int num_rem = NUMPARTICLES_Current - num_perGPU*num_devices_use;


		for ( int i=0; i < num_devices_use; i++ )
		{
			num_particles_perGPU[i] = num_perGPU;

			if(i==0)
			{
					num_particles_perGPU[0] = num_perGPU + num_rem;

			}


		num_particles_SM = (int)ceil(num_particles_perGPU[i]/(float)num_sm);

	    /* Now get the number of blocks per SM */
	    num_blocks_SM      =  (int)ceil(num_particles_SM/(float)threads_perBlock);
	    num_threads_block = threads_perBlock ;

	    /* block size is too big */
	    if( num_particles_SM < threads_perBlock )
	    {
		 num_threads_block = num_particles_SM; /* Single block is sufficient */
	    }

	    multi[i].dimBlock = make_uint3( num_threads_block,1,1);
	    multi[i].dimGrid  = make_uint3( num_blocks_SM*num_sm,1,1 );

	    //printf(" GPU %d Particles %d per GPU threads = %d \n",i, num_particles_perGPU[i],multi[i].dimGrid.x*multi[i].dimBlock.x);

		}

}



void Sync_KilledMulti()
{
	cudaDeviceSynchronize();
	cudaSetDevice(0);

   /* Calculate the Hash for All particles based on its Spatial Position */
   Kernel_SpatialDecomp_CalcPHash<<< GlobalLaunch.dimGrid, GlobalLaunch.dimBlock >>>
   		                          ( dram_GridParticleHashM[0], dram_GridParticleIndexM[0],
   		                            dram_position_comM[0], NUMPARTICLES_Current );
   /*-----------------------------------------------------------------------*/
   cudaDeviceSynchronize();

	 if(error_check==1)
	 {
	    	cudaError_t errormsg=cudaGetLastError();
	    	if(errormsg>0)
	    	{
	    	cout<<"DEM_UpdateSim: Hash "<<cudaGetErrorString(errormsg)<<endl;
	    	 exit(1);
	    	}
	 }
   /*-----------------------------------------------------------------------*/
          /* 2. Sort Particles according to (Grid Hash) */
   /*-----------------------------------------------------------------------*/

   /* In index 0 is lowest and last is highest */
   SpatialDecomp_ThrustSort_ByHash(dram_GridParticleHashM[0], dram_GridParticleIndexM[0], NUMPARTICLES_Current);
   cudaDeviceSynchronize();
   /*-----------------------------------------------------------------------*/

   /* m_dGridParticleIndex contains the order which the dynamics arrays must be reordered */

   /* Sort arrays using temp copies */
   Kernel_SortArrays
     <<< GlobalLaunch.dimGrid, GlobalLaunch.dimBlock >>>
     ( dram_GridParticleIndex,
   	dSortedPos,
   	dSortedPosQ,
   	dSortedVel,
   	dSorted_velocity_ang,

   	dSortedPType,
   	dSortedPID,

       dSortedLifter_Contact_Hist,

       dram_position_comM[0],
       dram_position_orntM[0],
       dram_velocity_comM[0],
       dram_velocity_angM[0],
       dram_ObjectTypeM[0],
       dram_P_IDM[0],

       dram_Lifter_Contact_HistM[0],

       NUMPARTICLES_Current );

   cudaDeviceSynchronize();

	 if(error_check==1)
	 {
	    	cudaError_t errormsg=cudaGetLastError();
	    	if(errormsg>0)
	    	{
	    	cout<<"DEM_UpdateSim: Sort "<<cudaGetErrorString(errormsg)<<endl;
	    	 exit(1);
	    	}
	 }
   /*-----------------------------------------------------------------------*/
               /* 4. Update the pos and vel based on sorting in 3. */
     /*-----------------------------------------------------------------------*/
     Kernel_SpatialDecomp_ReorderPDA
                 <<< GlobalLaunch.dimGrid, GlobalLaunch.dimBlock >>>
     		    ( dSortedPos,dSortedPosQ, dSortedVel,dSorted_velocity_ang, dSortedPType,dSortedPID,
     		      dSortedLifter_Contact_Hist,

     		      dram_position_comM[0],
     		      dram_position_orntM[0],
     		      dram_velocity_comM[0],
                 dram_velocity_angM[0],
                 dram_ObjectTypeM[0],
                 dram_P_IDM[0],
                 dram_Lifter_Contact_HistM[0],

                 NUMPARTICLES_Current );
    cudaDeviceSynchronize();
     /*-----------------------------------------------------------------------*/

	 if(error_check==1)
	 {
	    	cudaError_t errormsg=cudaGetLastError();
	    	if(errormsg>0)
	    	{
	    	cout<<"DEM_UpdateSim: Reorder "<<cudaGetErrorString(errormsg)<<endl;
	    	 exit(1);
	    	}
	 }

	  Remove_Flagged_ParticlesMulti();
	  cudaDeviceSynchronize();

}


/*---------------------------------------------------------------------------*/
                      /* Brute Force Multi-GPU  */
/*---------------------------------------------------------------------------*/
void Device_DEM_UpdateSimMulti()
{
	int mem_loc = 0;

	cudaSetDevice(0);

	 if(error_check==1)
	 {
	    	cudaError_t errormsg=cudaGetLastError();
	    	if(errormsg>0)
	    	{
	    	 cout<<"DEM_UpdateSim: Start  "<<cudaGetErrorString(errormsg)<<Call_Time<<endl;
	    	 exit(1);
	    	}
	 }

	if( isFilling )
	{
	 Filling_Logic_Simulation();
	}

	                /* Reorder All arrays based on hash */
    /*-----------------------------------------------------------------------*/
           /* 1. Calculate the GridIndex of Each Particle (Grid Hash) */
    /*-----------------------------------------------------------------------*/

    

    /* Calculate the Hash for All particles based on its Spatial Position */
    Kernel_SpatialDecomp_CalcPHash<<< GlobalLaunch.dimGrid, GlobalLaunch.dimBlock,0,streams[0] >>>
    		                          ( dram_GridParticleHashM[0], dram_GridParticleIndexM[0],
    		                            dram_position_comM[0],NUMPARTICLES_Current );

	 if(error_check==1)
	 {
	    	cudaError_t errormsg=cudaGetLastError();
	    	if(errormsg>0)
	    	{
	    	cout<<"DEM_UpdateSim: Hash "<<cudaGetErrorString(errormsg)<<Call_Time<<endl;
	    	 exit(1);
	    	}
	 }
	/*-----------------------------------------------------------------------*/
	
	 
	 /*-----------------------------------------------------------------------*/
           /* 2. Sort Particles according to (Grid Hash) */
    /*-----------------------------------------------------------------------*/
    
	cudaSetDevice(0);
    /* In index 0 is lowest and last is highest */
    SpatialDecomp_ThrustSort_ByHashMulti(dram_GridParticleHashM[0], dram_GridParticleIndexM[0],NUMPARTICLES_Current,streams[0]);
			   if(error_check==1)
	       {
	    	 cudaError_t errormsg=cudaGetLastError();
	    	 if(errormsg>0)
	    	 {
	    	  cout<<"DEM_UpdateSim: SortHash  "<<cudaGetErrorString(errormsg)<<Call_Time<< "Device 0 "<<endl;
	    	  exit(1);
	    	 }
	       }
	
	cudaDeviceSynchronize();
	



	for (int i=1;i<num_devices_use;i++)
	{
	/* Copy Hash List to all GPUS */
	 cudaMemcpyPeerAsync(dram_GridParticleHashM [i], i,  dram_GridParticleHashM[0],0,sizeof(uint)*NUMPARTICLES_Current,streams[1]);
	 cudaMemcpyPeerAsync(dram_GridParticleIndexM[i], i,  dram_GridParticleIndexM[0],0,sizeof(uint)*NUMPARTICLES_Current,streams[1]);
	}
	cudaEventRecord(events[0],streams[1]);


	
	/*-----------------------------------------------------------------------*/


    /* m_dGridParticleIndex contains the order which the dynamics arrays must be reordered */

    /* Sort arrays using temp copies */
    Kernel_SortArrays
      <<< GlobalLaunch.dimGrid, GlobalLaunch.dimBlock,0,streams[0]  >>>
	  ( dram_GridParticleIndexM[0],
    	dSortedPos,
    	dSortedPosQ,
    	dSortedVel,
    	dSorted_velocity_ang,

    	dSortedPType,
    	dSortedPID,

        dSortedLifter_Contact_Hist,

        dram_position_comM[0],
        dram_position_orntM[0],
        dram_velocity_comM[0],
        dram_velocity_angM[0],
        dram_ObjectTypeM[0],
        dram_P_IDM[0],

        dram_Lifter_Contact_HistM[0],
        NUMPARTICLES_Current);

		   if(error_check==1)
	       {
	    	 cudaError_t errormsg=cudaGetLastError();
	    	 if(errormsg>0)
	    	 {
	    	  cout<<"DEM_UpdateSim: Integrate  "<<cudaGetErrorString(errormsg)<<Call_Time<< "Device Sort "<<endl;
	    	  exit(1);
	    	 }
	       }

	
    /*-----------------------------------------------------------------------*/
           /* 3. Bin data into Cells based on the Hash 
		         Each GPU has its own copy of the full cell data */
    /*-----------------------------------------------------------------------*/

    cudaEventSynchronize(events[0]);/* Check if Copy is Done use hash to bin into cells */
    for (int i=0;i<num_devices_use;i++)
	{
     
	   cudaSetDevice(i); 
       /* Bins data into cells - Creates from scratch at each step */
       Kernel_SpatialDecomp_BinData
       <<< GlobalLaunch.dimGrid, GlobalLaunch.dimBlock,0,streams[i*4+1] >>>
	   ( dram_GridParticleHashM[i],m_dCellStartM[i], m_dCellEndM[i],NUMPARTICLES_Current);
	 
		   if(error_check==1)
	       {
	    	 cudaError_t errormsg=cudaGetLastError();
	    	 if(errormsg>0)
	    	 {
	    	  cout<<"DEM_UpdateSim: Bin  "<<cudaGetErrorString(errormsg)<<Call_Time<< "Device "<<i<<endl;
	    	  exit(1);
	    	 }
	       }

	}/* End binning on all GPUS */
    /*-----------------------------------------------------------------------*/
   
    cudaSetDevice(0); 

	 if( error_check == 1 )
	 {
	    	cudaError_t errormsg=cudaGetLastError();
	    	if(errormsg>0)
	    	{
	    	cout<<"DEM_UpdateSim: Sort "<<cudaGetErrorString(errormsg)<<endl;
	    	 exit(1);
	    	}
	 }

    /*-----------------------------------------------------------------------*/
                /* 4. Update the pos and vel based on sorting in 3. */
      /*-----------------------------------------------------------------------*/
      Kernel_SpatialDecomp_ReorderPDA
                  <<< GlobalLaunch.dimGrid, GlobalLaunch.dimBlock,0,streams[0] >>>
      		    ( dSortedPos,dSortedPosQ, dSortedVel,dSorted_velocity_ang, dSortedPType,dSortedPID,
      		      dSortedLifter_Contact_Hist,

      		      dram_position_comM[0],
      		      dram_position_orntM[0],
      		      dram_velocity_comM[0],
                  dram_velocity_angM[0],
                  dram_ObjectTypeM[0],
                  dram_P_IDM[0],
                  dram_Lifter_Contact_HistM[0],
                  NUMPARTICLES_Current );

 		   if(error_check==1)
	       {
	    	 cudaError_t errormsg=cudaGetLastError();
	    	 if(errormsg>0)
	    	 {
	    	  cout<<"DEM_UpdateSim: Reorder  "<<cudaGetErrorString(errormsg)<<Call_Time<< "Device 0 "<<endl;
	    	  exit(1);
	    	 }
	       }
      /*-----------------------------------------------------------------------*/

	 if(error_check==1)
	 {
	    	cudaError_t errormsg=cudaGetLastError();
	    	if(errormsg>0)
	    	{
	    	cout<<"DEM_UpdateSim: Reorder "<<cudaGetErrorString(errormsg)<<endl;
	    	 exit(1);
	    	}
	 }

	/* Remove flagged particles from the sim by reducing num threads */
	if(num_dead>0)
	{
		Remove_Flagged_ParticlesMulti();
	}


	/* Copy all dynamics arrays to other GPUS only after reorder is done */
	cudaDeviceSynchronize();
	for (int i=1;i<num_devices_use;i++)
	{
	/* Copy some  primary info to All GPUS */ 
	 cudaMemcpyPeerAsync(dram_position_comM[i],i, dram_position_comM[0],0, sizeof(float3)*NUMPARTICLES_Current,streams[0]);
	 cudaMemcpyPeerAsync(dram_P_IDM        [i],i, dram_P_IDM        [0],0, sizeof(int)*NUMPARTICLES_Current,   streams[0]);
	 cudaMemcpyPeerAsync(dram_ObjectTypeM  [i],i, dram_ObjectTypeM  [0],0, sizeof(uint)*NUMPARTICLES_Current,  streams[0]);
	 
	 /* Need all data if we not keeping a list */
	 if(NO_NN)
	 {
	    /* Copy remaining data to other GPUS */
		cudaMemcpyPeerAsync(dram_velocity_comM[i],i,dram_velocity_comM[0],0,sizeof(float3)*NUMPARTICLES_Current,streams[1]);
	    cudaMemcpyPeerAsync(dram_velocity_angM[i],i,dram_velocity_angM[0],0,sizeof(float3)*NUMPARTICLES_Current,streams[1]);

	   if(SimInfo_C->particle_type==polyhedra || SimInfo_C->sphere_orient)
	   {
		 cudaMemcpyPeerAsync(dram_position_orntM[i],i,dram_position_orntM[0],0,sizeof(Quaterion)*NUMPARTICLES_Current,streams[1]);
	   }

	  }

	}

	for (int i=0;i<num_devices_use;i++)
	{
	  cudaSetDevice(i);
	  cudaDeviceSynchronize();
	}
	cudaSetDevice(0);
    /*-----------------------------------------------------------------------*/
                       /* 5. Scan ALL NN PAIRS */
    /*-----------------------------------------------------------------------*/
	
	if( SimInfo_C->particle_type==spheres)
	{
	   mem_loc=0;
	   for (int i=0;i<num_devices_use;i++)
	   {
 		
	     cudaSetDevice(i);
	     Kernel_Spheres_PP
                         <<< multi[i].dimGrid, multi[i].dimBlock,0,streams[i*4]>>>
    		             (get_P_forces,
				  dPContactPointsM[i], 
				  m_dCellStartM[i],m_dCellEndM[i], dram_position_comM[i],
				  dram_velocity_comM[i],
				  dram_velocity_angM[i],
				  dram_ObjectTypeM[i],
				  dram_P_IDM[i], dram_force_PPM[i],dram_force_PP_angM[i],
				  num_particles_perGPU[i],dNumContactsM[i],mem_loc);
	      /* Register events */	   
	      mem_loc += num_particles_perGPU[i];
	  	  if(error_check==1)
	       {
	    	 cudaError_t errormsg=cudaGetLastError();
	    	 if(errormsg>0)
	    	 {
	    	  cout<<"DEM_UpdateSim: Spheres NN  "<<cudaGetErrorString(errormsg)<<Call_Time<< "Device "<<i<<endl;
	    	  exit(1);
	    	 }
	       }
	   }
	}
	else
	{ 
	   mem_loc=0;
	   for (int i=0;i<num_devices_use;i++)
	   {
 		
	     cudaSetDevice(i);
         Kernel_BroadCollisionDetection_NonSymmetry
		    <<< multi[i].dimGrid, multi[i].dimBlock,0,streams[i*4] >>>
    		             ( m_dCellStartM[i],
						   m_dCellEndM[i],
						   dBroad_ListM[i],
						   dNumNNM[i],
						   dram_position_comM[i],
						   dram_ObjectTypeM[i],
						   dram_P_IDM[i],
						   num_particles_perGPU[i], 
						   mem_loc);

         /* Register events */	   
	     mem_loc += num_particles_perGPU[i];
	  	 if(error_check==1)
	       {
	    	 cudaError_t errormsg=cudaGetLastError();
	    	 if(errormsg>0)
	    	 {
	    	  cout<<"DEM_UpdateSim: BroadPhase  "<<cudaGetErrorString(errormsg)<<Call_Time<< "Device "<<i<<endl;
	    	  exit(1);
	    	 }
	       }
	 }
    cudaSetDevice(0);
	for (int i=1;i<num_devices_use;i++)
	{
	 /* Copy remaining data to other GPUS */

		cudaMemcpyPeerAsync(dram_velocity_comM[i],i,dram_velocity_comM[0],0,sizeof(float3)*NUMPARTICLES_Current,streams[1]);

	   cudaMemcpyPeerAsync(dram_velocity_angM[i],i,dram_velocity_angM[0],0,sizeof(float3)*NUMPARTICLES_Current,streams[1]);


	 
	   if(SimInfo_C->particle_type==polyhedra || SimInfo_C->sphere_orient)
	   {
		 cudaMemcpyPeerAsync(dram_position_orntM[i],i,dram_position_orntM[0],0,sizeof(Quaterion)*NUMPARTICLES_Current,streams[1]);
	   }

	}

  }
	 


	/* Wait for All GPUS to finish */
	for (int i=0;i<num_devices_use;i++)
    {
      cudaSetDevice(i);
	  cudaDeviceSynchronize();
    }

    /*-----------------------------------------------------------------------*/
                       /* 6. Apply Particle Particle Forces */
    /*-----------------------------------------------------------------------*/


      if(SimInfo_C->particle_type==polyhedra)
	  {

         mem_loc = 0;
		 for (int i=0;i<num_devices_use;i++)
	     {
	        cudaSetDevice(i);
	        Kernel_ParticleInteraction_Polyhedra_NonSymmetry<<< multi[i].dimGrid, multi[i].dimBlock,0,streams[i*4] >>>
				( dNumNNM[i],dBroad_ListM[i],
				  dram_position_comM[i],
				  dram_position_orntM[i],
				  dram_velocity_comM[i],
				  dram_velocity_angM[i],
				  dram_ObjectTypeM[i],
				  dram_P_IDM[i],
				  dram_force_PPM[i],
				  dram_force_PP_angM[i],
				  num_particles_perGPU[i], mem_loc   );

			mem_loc += num_particles_perGPU[i];

			  if(error_check==1)
		       {
		    	 cudaError_t errormsg=cudaGetLastError();
		    	 if(errormsg>0)
		    	 {
		    	  cout<<"DEM_UpdateSim: Particle Forces  "<<cudaGetErrorString(errormsg)<<Call_Time<< "Device "<<i<<endl;
		    	  exit(1);
		    	 }
		       }
		 }
      }
	
	
		 
  

    if(NUMDYOBJECTS>0)
    {
	  mem_loc=0;
      if(SimInfo_C->particle_type==spheres)
      {
        
       	for (int i=0;i<num_devices_use;i++)
	    {
	      cudaSetDevice(i);
	      VolumeObject_InteractionSpheres<<< multi[i].dimGrid, multi[i].dimBlock, 0,streams[i*4 + 1 ]  >>>(
			     get_vobject_forces,
				 dVOContactPointsM[i],
			     dram_force_com_LifterM[i],
			     dram_force_ang_LifterM[i],
			     dram_Lifter_Contact_HistM[i],
			     dram_position_comM[i],
			     dram_velocity_comM[i],
			     dram_velocity_angM[i],
	             dram_ObjectTypeM[i], 
				 dram_P_IDM[i], 
				 num_particles_perGPU[i],
				 mem_loc);
		  mem_loc+=num_particles_perGPU[i];

		  if(error_check==1)
	       {
	    	 cudaError_t errormsg=cudaGetLastError();
	    	 if(errormsg>0)
	    	 {
	    	  cout<<"DEM_UpdateSim: Volume  "<<cudaGetErrorString(errormsg)<<Call_Time<< "Device "<<i<<endl;
	    	  exit(1);
	    	 }
	       }
		}
      }
      else
      {
        for (int i=0;i<num_devices_use;i++)
	    {
	       cudaSetDevice(i);
  	       VolumeObject_InteractionPolyhedra<<< multi[i].dimGrid, multi[i].dimBlock, 0,streams[i*4 + 1 ] >>>(
			     dram_force_com_LifterM[i],
			     dram_force_ang_LifterM[i],
			     dram_Lifter_Contact_HistM[i],
			     dram_position_comM[i],
                 dram_position_orntM[i],
				 dram_velocity_comM[i],
			     dram_velocity_angM[i],
	             dram_ObjectTypeM[i], 
				 dram_P_IDM[i], 
				 num_particles_perGPU[i],
				 mem_loc);

		   mem_loc += num_particles_perGPU[i];
		   if(error_check==1)
		 	       {
		 	    	 cudaError_t errormsg=cudaGetLastError();
		 	    	 if(errormsg>0)
		 	    	 {
		 	    	  cout<<"DEM_UpdateSim: Volume  "<<cudaGetErrorString(errormsg)<<Call_Time<< "Device "<<i<<endl;
		 	    	  exit(1);
		 	    	 }
		 	       }
		}

       }

    }/* End volume objects check */


	mem_loc = 0;
    if(SimInfo_C->particle_type==spheres)
    {
      for (int i=0;i<num_devices_use;i++)
	  {
	    cudaSetDevice(i);

	    /* We do one world object at time so less divergence make sure its run one at a time */

	    	WorldObject_InteractionSpheres_Planar<<< multi[i].dimGrid, multi[i].dimBlock,0,streams[i*4 +2 ]  >>>( get_wall_forces,
			                                                                        dram_force_com_WallM[i],
			                                                                        dram_force_ang_WallM[i],
			                                                                        dWallContactPointsM[i],
			                                                                        dram_position_comM[i],
			                                                                        dram_position_orntM[i],
			                                                                        dram_velocity_comM[i],
			                                                                        dram_velocity_angM[i],
                                                                                    dram_ObjectTypeM[i],
                                                                                    dram_P_IDM[i],
																					num_particles_perGPU[i],
																					mem_loc);



		    	WorldObject_InteractionSpheres_Macro<<< multi[i].dimGrid, multi[i].dimBlock,0,streams[i*4 +2 ]  >>>( get_wall_forces,
		    			                                                                is_cylinder_rotation,
				                                                                        dram_force_com_WallM[i],
				                                                                        dram_force_ang_WallM[i],
				                                                                        dWallContactPointsM[i],
				                                                                        dram_position_comM[i],
				                                                                        dram_position_orntM[i],
				                                                                        dram_velocity_comM[i],
				                                                                        dram_velocity_angM[i],
	                                                                                    dram_ObjectTypeM[i],
	                                                                                    dram_P_IDM[i],
																						num_particles_perGPU[i],

																						mem_loc);



		     mem_loc += num_particles_perGPU[i];

		   if(error_check==1)
	       {
	    	 cudaError_t errormsg=cudaGetLastError();
	    	 if(errormsg>0)
	    	 {
	    	  cout<<"DEM_UpdateSim: World  "<<cudaGetErrorString(errormsg)<<Call_Time<< "Device "<<i<<endl;
	    	  exit(1);
	    	 }
	       }
	  }

      /* Check Macros */

    }
    else
    {
      for (int i=0;i<num_devices_use;i++)
	  {
	    cudaSetDevice(i); 


    	WorldObject_InteractionPolyhedra_Planar<<< multi[i].dimGrid, multi[i].dimBlock,0,streams[i*4 +2 ]   >>>( get_wall_forces, is_cylinder_rotation,
			                                                                        dram_force_com_WallM[i],
			                                                                        dram_force_ang_WallM[i],
			                                                                        dWallContactPointsM[i],
			                                                                        dram_position_comM[i],
			                                                                        dram_position_orntM[i],
			                                                                        dram_velocity_comM[i],
			                                                                        dram_velocity_angM[i],
                                                                                    dram_ObjectTypeM[i],
																					dram_P_IDM[i],
																					num_particles_perGPU[i],
																					mem_loc);

	    	WorldObject_InteractionPolyhedra_Macro<<< multi[i].dimGrid, multi[i].dimBlock,0,streams[i*4 +2 ]   >>>( get_wall_forces, is_cylinder_rotation,
	    				                                                                        dram_force_com_WallM[i],
	    				                                                                        dram_force_ang_WallM[i],
	    				                                                                        dWallContactPointsM[i],
	    				                                                                        dram_position_comM[i],
	    				                                                                        dram_position_orntM[i],
	    				                                                                        dram_velocity_comM[i],
	    				                                                                        dram_velocity_angM[i],
	    	                                                                                    dram_ObjectTypeM[i],
	    																						dram_P_IDM[i],
	    																						num_particles_perGPU[i],
	    																						mem_loc);



		 mem_loc += num_particles_perGPU[i];

		   if(error_check==1)
 {
	 cudaError_t errormsg=cudaGetLastError();
	 if(errormsg>0)
	 {
	  cout<<"DEM_UpdateSim: World  "<<cudaGetErrorString(errormsg)<<Call_Time<< "Device "<<i<<endl;
	  exit(1);
	 }
 }

		}
    }

	/* Copy lifter history  */
	 
	 if( NUMDYOBJECTS > 0 )
	 {

		mem_loc = 0;
		for (int i=1;i<num_devices_use;i++)
		{
		  cudaSetDevice(i);
		  cudaDeviceSynchronize();
	      mem_loc += num_particles_perGPU[i-1];

	 	  /* Copy lifter hist to GPU 0 */
	 	  cudaMemcpyPeerAsync( dram_Lifter_Contact_HistM[0] + mem_loc,0,
			                   dram_Lifter_Contact_HistM[i] + mem_loc,i, 
							   sizeof(Contact_Info)*num_particles_perGPU[i],
							   streams[i*4 + 1]                             );
		 }
	    cudaSetDevice(0);
	    cudaDeviceSynchronize();

	 }
	 else
	 {
       	for (int i=0;i<num_devices_use;i++)
		{
		  cudaSetDevice(i);
		  cudaDeviceSynchronize();
		}
	 }

	 
	  mem_loc = 0;
	  for ( int i=0; i < num_devices_use; i++ )
	  {
	    cudaSetDevice(i);



            Integrate_Euler_NonSymmetry<<< multi[i].dimGrid, multi[i].dimBlock >>>
              ( dram_force_PPM[i],          dram_force_PP_angM[i],
              dram_force_com_WallM[i],   dram_force_ang_WallM[i],
              dram_force_com_LifterM[i], dram_force_ang_LifterM[i],

              dram_position_comM[i], dram_position_orntM[i],
              dram_velocity_comM[i], dram_velocity_angM[i],
              dram_ObjectTypeM[i], dram_P_IDM[i],
              num_particles_perGPU[i], get_wall_forces, mem_loc );

		      mem_loc += num_particles_perGPU[i];

		   if(error_check==1)
	       {
	    	 cudaError_t errormsg=cudaGetLastError();
	    	 if(errormsg>0)
	    	 {
	    	  cout<<"DEM_UpdateSim: Integrate  "<<cudaGetErrorString(errormsg)<<Call_Time<< "Device "<<i<<endl;
	    	  exit(1);
	    	 }
	       }
	  }
   

	/* Copy data from all gpus back to 0 */



	  for (int i=0;i<num_devices_use;i++)
	  {
	    cudaSetDevice(i);
	    cudaDeviceSynchronize();
	  }


	    for (int i=0;i<num_devices_use;i++)
		{

		   cudaSetDevice(i);
	       /* Bins data into cells - Creates from scratch at each step */
		   Kernel_SpatialDecomp_BinDataR
	       <<< GlobalLaunch.dimGrid, GlobalLaunch.dimBlock,0,streams[i*4+1] >>>
		   ( dram_GridParticleHashM[i],m_dCellStartM[i], m_dCellEndM[i],NUMPARTICLES_Current);

			   if(error_check==1)
		       {
		    	 cudaError_t errormsg=cudaGetLastError();
		    	 if(errormsg>0)
		    	 {
		    	  cout<<"DEM_UpdateSim: BinReset  "<<cudaGetErrorString(errormsg)<<Call_Time<< "Device "<<i<<endl;
		    	  exit(1);
		    	 }
		       }

		}/* End binning on all GPUS */
	    /*-----------------------------------------------------------------------*/


	  cudaSetDevice(0);

	 mem_loc = 0;
	
	for ( int i=1; i < num_devices_use; i++ )
	{
      mem_loc += num_particles_perGPU[i-1];
	  /* Copy All primary info to All GPUS */ 
	  cudaMemcpyPeerAsync(dram_position_comM[0] + mem_loc,0, dram_position_comM[i]+ mem_loc, i,sizeof(float3)*num_particles_perGPU[i]);
	  cudaMemcpyPeerAsync(dram_velocity_comM[0] + mem_loc,0, dram_velocity_comM[i]+ mem_loc, i,sizeof(float3)*num_particles_perGPU[i]);

      cudaMemcpyPeerAsync(dram_velocity_angM[0] + mem_loc,0, dram_velocity_angM[i] + mem_loc, i,sizeof(float3)*num_particles_perGPU[i]);
	
	  if(SimInfo_C->particle_type==polyhedra)
      {
		 cudaMemcpyPeerAsync(dram_position_orntM[0]+mem_loc,0,dram_position_orntM[i]+ mem_loc, i,sizeof(Quaterion)*num_particles_perGPU[i]);
	  }
	}

	  /* Wait for all copies to be done */
	  for (int i=0;i<num_devices_use;i++)
	  {
	    cudaSetDevice(i);
	    cudaDeviceSynchronize();
	  }
	  cudaSetDevice(0);


     if( SimInfo_C->Kill_Particles)
     {
		Flag_DeadParticlesMulti();
	 }
	 Call_Time+=SimInfo_C->InitalDelta_t;

	 	 if (isFilling)
	 {
		 fill_time+=SimInfo_C->InitalDelta_t;
		 if(fill_time>0.25f)
		 {
			 printf(" Filled Particles %d \n",NUMPARTICLES_Current);
			 fill_time = 0.0f;
		 }
	 }
     
}
/*---------------------------------------------------------------------------*/
