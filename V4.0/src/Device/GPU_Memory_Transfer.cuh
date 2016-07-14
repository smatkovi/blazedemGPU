

/* Contains the Kernels and methods that facilitate
 * the transfer between Device and Host  */


/* Structure Definitions */
/* Classification of contact Type */
enum COL_TYPE   { UNDEF       = -1,
	              SP_PLANE    =  0,
	              Vertex_Face =  1,
	              Edge_Face   =  2,
	              Face_Face   =  3,
	              Edge_Edge   =  4,
	              MEdge_MEdge =  5,
	              EdgeF_EdgeF =  6 };


/* Return type for polyhedra contact */
struct PointInPoly
{
  bool  is_point_inside;
  int   Close_Face_num [3];
  float distance   [3];
};





/* Particle Object for rotation */
struct POBJ_ROT
{
  float3 Vertex [32];
  Plane  face   [32];
};


/*------------------------------------------------------*/
             /* Contact Info */
/*------------------------------------------------------*/

/* Generic Contact Info */
struct CollData
{
  bool      is_collision;
  float     dis;
  float3    contact_point;
  float3    Normal;

  COL_TYPE  collision_class;
};
/*------------------------------------------------------*/


struct CollData_Sphere
{
  bool      is_collision;
  float     dis;
  float3    contact_point;
  float3    Normal;
};


struct NN_Selected_Cell
{
	int NN_Index[32];
};





/*---------------------------------------------------------------------------*/
                      /* Simulation Data */
/*---------------------------------------------------------------------------*/

__constant__  SimulationInfo     SimParms;

__constant__ WOBJ    WorldObject    [8];
__constant__ POBJ    ParticleObject [8];
__constant__ DOBJ    VolumeObject  [16];


__constant__ float3  InitVel[2];


__device__ int Num_WorldObjects;    /* Max 10 World Objects    */
__device__ int Num_ParticleObjects; /* Max 10 Particle Objects */
__device__ int Num_DynamicObjects;  /* Max 10 Particle Objects */

__device__ int Total_NN;

__device__ float fc;

/*---------------------------------------------------------------------------*/

__device__ float Tally_EnergyDis_PP;
__device__ float Tally_EnergyDisLift;
__device__ float Tally_EnergyDisSurf;
__device__ float Tally_EnergyDis_Norm;
__device__ float Tally_EnergyDis_Shear;


__device__ float Tally_EnergyTot;
/*---------------------------------------------------------------------------*/



__device__ int Kill_count=0;
__device__ float Mass_count=0.0f;
__device__ float Vol_count=0.0f;

__device__ float Device_Time=0.0f;

__device__ int TotalNN=0;


/*---------------------------------------------------------------------------*/
        /* Start SIM: Sets the initial configuration values for particles */
/*---------------------------------------------------------------------------*/
__global__ void Set_Init_Config_Device(  int    num_Particle_Types,
		                                 int    num_world_Objects,
		                                 int    num_dynamic_Objects,

		                                 float  *force_com_PP_X,
		                                 float  *force_com_PP_Y,
		                                 float  *force_com_PP_Z,
		                                 float3 *force_com_PW,
		                                 float3 *force_com_PL,

		                                 float  *force_ang_PP_X,
		                                 float  *force_ang_PP_Y,
		                                 float  *force_ang_PP_Z,
		                                 float3  *force_ang_PW,
		                                 float3  *force_ang_PL,

		                                 float3  *velocity_com,
		                                 float3  *velocity_ang,
		  		                         Contact_Info    *Lifter_Contact_Hist,
										 int     Num_Particles,
		  		                         uint gpu_offset = 0          )
{
	uint index = blockIdx.x*blockDim.x  + threadIdx.x + gpu_offset;



	if( index < Num_Particles)
   {

	 /* First thread store number of objects */
	 if( index==0 )
	 {

	   Num_WorldObjects    = num_world_Objects;
	   Num_ParticleObjects = num_Particle_Types;
	   Num_DynamicObjects  = num_dynamic_Objects;
	   /*TODO COUNTERS HERE*/
	   Total_NN = 0;
	 }
	 __syncthreads();

	 velocity_com [index]   = InitVel[0];
     velocity_ang [index]   = make_float3 (0.000000f,0.000000f,0.000000f); /*rad/sec*/

	 /* Particle Particle force split into components due to atomic adds */
     force_com_PP_X [index] = 0.000000f;
     force_com_PP_Y [index] = 0.000000f;
     force_com_PP_Z [index] = 0.000000f;

	 force_ang_PP_X [index] = 0.000000f;
     force_ang_PP_Y [index] = 0.000000f;
     force_ang_PP_Z [index] = 0.000000f;



     force_com_PW [index]   = make_float3 (0.000000f,0.000000f,0.000000f);
     force_ang_PW [index]   = make_float3 (0.000000f,0.000000f,0.000000f);

     /* If we have lifters */
     if(Num_DynamicObjects>0)
     {
        force_com_PL [index]   = make_float3 (0.000000f,0.000000f,0.000000f);
        force_ang_PL [index]   = make_float3 (0.000000f,0.000000f,0.000000f);

        Lifter_Contact_Hist[index].num_contact = 0;
        Lifter_Contact_Hist[index].cont_type   = -1;
     }

   }

}
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
        /* Start SIM: Sets the initial configuration values for particles */
/*---------------------------------------------------------------------------*/
__global__ void Set_Init_Config_Device(  int    num_Particle_Types,
		                                 int    num_world_Objects,
		                                 int    num_dynamic_Objects,
		                                 float3  *force_com_PP,
		                                 float3  *force_com_PW,
		                                 float3  *force_com_PL,
		                                 float3  *force_ang_PP,
		                                 float3  *force_ang_PW,
		                                 float3  *force_ang_PL,
		                                 float3  *velocity_com,
		                                 float3  *velocity_ang,
		  		                         Contact_Info    *Lifter_Contact_Hist,
										 int     Num_Particles,
		  		                         uint gpu_offset = 0)
{
	uint index = blockIdx.x*blockDim.x  + threadIdx.x ;



   if( index < Num_Particles )
   {

	 /* First thread store number of objects */
	 if( index==0 )
	 {

	   Num_WorldObjects    = num_world_Objects;
	   Num_ParticleObjects = num_Particle_Types;
	   Num_DynamicObjects  = num_dynamic_Objects;
	   /*TODO COUNTERS HERE*/
	   Total_NN = 0;

	   printf(" hello  %f %f %f \n",InitVel[0].x,InitVel[0].y,InitVel[0].z);
	 }




	 
	 velocity_com [index+ gpu_offset]   = InitVel[0];

//	 if(gpu_offset/Num_Particles==1)
//	 {
//	   velocity_com [index+ gpu_offset] = make_float3 (0.000000f,300.000000f,0.000000f);
//
//	 }
//	 if(gpu_offset/Num_Particles==2)
//	 {
//	   velocity_com [index+ gpu_offset] = make_float3 (0.000000f,400.000000f,0.000000f);
//	 }

	// printf(" GPU %d index %d Val %f \n",gpu_offset/Num_Particles ,index+ gpu_offset,velocity_com[index+gpu_offset ].y);

	 velocity_ang [index + gpu_offset]   = make_float3 (0.000000f,0.000000f,0.000000f); /*rad/sec*/

//	 if(gpu_offset/Num_Particles==1)
//	 {
//		 velocity_ang [index+ gpu_offset] = make_float3 (100.000000f,200.000000f,300.000000f);
//
//	 }

     force_com_PP [index]   = make_float3 (0.000000f,0.000000f,0.000000f);
     force_com_PW [index]   = make_float3 (0.000000f,0.000000f,0.000000f);

     force_ang_PP [index]   = make_float3 (0.000000f,0.000000f,0.000000f);
     force_ang_PW [index]   = make_float3 (0.000000f,0.000000f,0.000000f);

     /* If we have lifters */
     if(Num_DynamicObjects>0)
     {
        force_com_PL [index]   = make_float3 (0.000000f,0.000000f,0.000000f);
        force_ang_PL [index]   = make_float3 (0.000000f,0.000000f,0.000000f);
        Lifter_Contact_Hist[index+ gpu_offset].num_contact = 0;
        Lifter_Contact_Hist[index+ gpu_offset].cont_type = -1;
     }


   }

}
/*---------------------------------------------------------------------------*/







__global__ void Set_Tallys( float *Tallys )
{
	   Tally_EnergyDis_PP    = Tallys[0];
	   Tally_EnergyDisLift   = Tallys[1];
	   Tally_EnergyDisSurf   = Tallys[2];
}



/*---------------------------------------------------------------------------*/
        /* Only updates subset of particles */
/*---------------------------------------------------------------------------*/
__global__ void Fill_Plane(float pack_pos ,float3 *position_com, float3 *velocity_com)
{
   uint index = blockIdx.x*blockDim.x  + threadIdx.x;

   if( index <SimParms.Num_Particles )
   {

     /* To Fill From reverse direction change here */
     if( position_com[index].y > pack_pos )
     {
    	 position_com[index].y = pack_pos;
		 velocity_com[index] = InitVel[0];
		 

     }
	 	//	 if(position_com[index].x>230.0)
		 //{
			//
   //        velocity_com[index].x = 150.0f;
		 //}

   }


}

/*---------------------------------------------------------------------------*/
               /* Update velocity after PP collision*/
/*---------------------------------------------------------------------------*/
__global__ void Packing_Vel( float3 *velocity_com)
{
  uint index = blockIdx.x*blockDim.x  + threadIdx.x ;

  if( index <SimParms.Num_Particles )
  {
	  
     velocity_com[index]*=0.25f;
	  
  }

}






/*---------------------------------------------------------------------------*/
        /* Start SIM: Sets the initial configuration values for particles */
/*---------------------------------------------------------------------------*/
__global__ void Set_Pos( float3 *Host_Pos,     Quaterion *Host_ORNT,
		                 uint   *Host_Ptype,   int   *Host_P_ID,
		                 float3 *position_com, Quaterion *position_ornt,
		                 uint   *P_ObjectType, int    *P_ID)
{
	uint index = blockIdx.x*blockDim.x  + threadIdx.x;

   if( index <SimParms.Num_Particles )
   {
	  if(SimParms.particle_type==polyhedra)
	  {
	    position_com  [index] = Host_Pos[index];// +ParticleObject[P_ObjectType[index]].COM;
	    position_ornt       [index] = Host_ORNT[index];
	  }
	  else
	  {
		position_com  [index] = Host_Pos[index];
	  }

	  if(SimParms.sphere_orient)
	  {
	   position_ornt       [index] = make_IquaterionD();
	  }


	   P_ID          [index] = Host_P_ID [index];
	   P_ObjectType  [index] = Host_Ptype[index];

   }

}


__global__ void Init_OldPos(float3 *position_com,float3 *position_comOld)
{
	   uint index = blockIdx.x*blockDim.x  + threadIdx.x;

	   if( index <SimParms.Num_Particles )
	   {
		   position_comOld[index] = position_com[index];
	   }

}



/*---------------------------------------------------------------------------*/
        /* Sets Pos and Velocity */
/*---------------------------------------------------------------------------*/
__global__ void Set_Pos_Vel( float3 *Host_Pos,     Quaterion *Host_ORNT,
		                     float3 *Host_Vel,     float3    *Host_AVel,
		                     uint   *Host_Ptype,   int       *Host_PID,
		                     float3 *position_com, Quaterion *position_ornt,
		                     float3 *velocity_com, float3    *velocity_ang,
		                     uint   *P_ObjectType,
		                     int    *P_ID)

{
   uint index = blockIdx.x*blockDim.x  + threadIdx.x;

   if( index <SimParms.Num_Particles )
   {

	  if(SimParms.particle_type==polyhedra)
	  {
		position_com  [index] = Host_Pos[index];//+ParticleObject[P_ObjectType[index]].COM;
		position_ornt [index] = Host_ORNT[index];
	  }
	  else
	  {
		position_com  [index] = Host_Pos[index];
	  }

	  velocity_com  [index] = Host_Vel  [index];
	  velocity_ang  [index] = Host_AVel [index];
	  P_ObjectType  [index] = Host_Ptype[index];
	  P_ID          [index] = Host_PID  [index];

   }

}

/*---------------------------------------------------------------------------*/
        /* Start SIM: Sets the initial configuration values for particles */
/*---------------------------------------------------------------------------*/
__global__ void Set_PSystem_State( float3 *Host_Pos,     Quaterion *Host_ORNT,
		                        float3 *Host_Vel,     float3    *Host_AVel,
		                        uint   *Host_Ptype,   int       *Host_PID,
		                        float3 *position_com, Quaterion *position_ornt,
		                        float3 *velocity_com, float3    *velocity_ang,
 		                        uint   *P_ObjectType,
 		                        int    *P_ID)

{
   uint index = blockIdx.x*blockDim.x  + threadIdx.x;

   if( index <SimParms.Num_Particles )
   {

	  if(SimParms.particle_type==polyhedra)
	  {
		position_com  [index] = Host_Pos[index] ;//+ ParticleObject[P_ObjectType[index]].COM;
		position_ornt       [index] = Host_ORNT[index];
	  }
	  else
	  {
		  position_com  [index] = Host_Pos[index];
	  }

	  if(SimParms.sphere_orient)
	  {
	   position_ornt       [index] = make_IquaterionD();
	  }



	 velocity_com  [index] = Host_Vel   [index];
	 velocity_ang  [index] = Host_AVel  [index];
	 P_ObjectType  [index] = Host_Ptype [index];
	 P_ID          [index] = Host_PID[index];

   }

}


/*---------------------------------------------------------------------------*/
           /* Returns the COM Position and Orientation to Host */
/*---------------------------------------------------------------------------*/
__global__ void Get_Particle_PosArray( float3 *Host_Pos,     Quaterion *Host_Quart,
		                               uint   *Host_PType,   int *Host_PID,
		                               float3 *position_com, Quaterion *position_ornt,
		  		                       uint   *P_ObjectType, int *P_ID, int Num_Particles                 )
{
	uint index = blockIdx.x*blockDim.x  + threadIdx.x;

	if( index <Num_Particles )
	{
	   if(SimParms.particle_type==spheres)
	   {
	     Host_Pos   [index] = position_com  [index];
	   }
	   else
	   {
		 Host_Pos   [index] = (position_com  [index]);
		 Host_Quart [index] = position_ornt [index];
	   }

		  if(SimParms.sphere_orient)
		  {
			  Host_Quart [index] = position_ornt [index];
		  }


	   Host_PType [index] = P_ObjectType  [index];
	   Host_PID   [index] = P_ID [index];
	}

}
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
           /* Returns the COM Position and Orientation to Host */
/*---------------------------------------------------------------------------*/
__global__ void Get_SystemState( float3   *Host_Pos,  Quaterion *Host_Quart,
		                         float3  *Host_Vel,   float3  *Host_AVel,
		                         uint    *Host_PType, int       *Host_PID,

		                         float3    *position_com, Quaterion *position_ornt,
		                         float3    *velocity_com, float3    *velocity_ang,
  		                         uint      *P_ObjectType, int *P_ID, int Num_Particles              )
{
	uint index = blockIdx.x*blockDim.x  + threadIdx.x;

	if( index < Num_Particles )
	{
	   if(SimParms.particle_type==spheres)
	   {
	     Host_Pos   [index] = position_com  [index];
	   }
	   else
	   {
		   Host_Pos   [index] = (position_com  [index]);
		   Host_Quart [index] = position_ornt [index];
	   }

	   Host_Vel      [index] = (velocity_com  [index]);
	   Host_AVel     [index] = (velocity_ang  [index]);
	  // Host_Acc      [index] = (accelrat_com  [index]);

	   Host_PType [index] = P_ObjectType      [index];
	   Host_PID   [index] = P_ID [index];


	}

}
/*---------------------------------------------------------------------------*/




/*---------------------------------------------------------------------------*/
           /* Returns the COM Position and Orientation to Host */
/*---------------------------------------------------------------------------*/
__global__ void Get_Particle_VelArray( float3 *Host_Vel, float3 *Host_VelR,
		                               uint    *Host_PType, int       *Host_PID,
		                               int Flag,
		                               float3 *velocity_com, float3 *velocity_ang,
		                               uint      *P_ObjectType, int *P_ID, int Num_Particles        )
{
	uint index = blockIdx.x*blockDim.x  + threadIdx.x;

	if( index < Num_Particles )
	{
	   Host_Vel      [index] = (velocity_com  [index]);

	   if (Flag==1)
	   {
	     Host_VelR     [index] = (velocity_ang  [index]);
	   }

	   Host_PType [index] = P_ObjectType [index];
	   Host_PID   [index] = P_ID         [index];
	}

}
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
           /* Returns the COM Position and Orientation to Host */
/*---------------------------------------------------------------------------*/
__global__ void Get_Wall_Forces( float3 *dWall_Forces,float3 *Wall_Forces,
		                         float3 *Contact_Points )
{
	uint index = blockIdx.x*blockDim.x  + threadIdx.x;

	if( index <SimParms.Num_Particles )
	{
         Wall_Forces[index] = -1.0f*dWall_Forces[index];
	}

}
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
               /* Update velocity after PP collision*/
/*---------------------------------------------------------------------------*/
__global__ void Update_VelArray(float3 *velocity_comT,float3 *velocity_com)
{
  uint index = blockIdx.x*blockDim.x  + threadIdx.x ;

  if( index < SimParms.Num_Particles )
  {
    velocity_com[index] = velocity_comT[index];
  }

}




__global__ void ResetEnergy(int num)
{

	if(num<=3)
	{
	  Tally_EnergyDis_PP  = 0.0f;
	  Tally_EnergyDisLift = 0.0f;
	  Tally_EnergyDisSurf = 0.0f;
	}

	if(num>3)
	{
	  Tally_EnergyDis_Norm  = 0.0f;
	  Tally_EnergyDis_Shear = 0.0f;
	}

}

__global__ void GetEnergy(float *H_Energy, int num)
{
	if(num<=3)
    {
      /*0.5 since we add energy for each collision twice */
	  H_Energy[0]=0.5f*Tally_EnergyDis_PP*1E-4;
	  H_Energy[1]=Tally_EnergyDisLift*1E-4;
	  H_Energy[2]=Tally_EnergyDisSurf*1E-4;
    }

	if(num>3)
	{
	 H_Energy[3]=Tally_EnergyDis_Norm*1E-4;
	 H_Energy[4]=Tally_EnergyDis_Shear*1E-4;
	}

}



/* Flags particles which meet the criteria below */
__global__ void Silo_Kill(float3 *position_com, int *P_ID)
{
    uint index = blockIdx.x*blockDim.x  + threadIdx.x ;

	if( index < SimParms.Num_Particles  && P_ID[index]>-1)
	{
		if( position_com[index].y < SimParms.Silo_Kill_Height )
		{
			 position_com[index]=SimParms.max_size;
			 P_ID[index]=-1;
			 atomicAdd(&Kill_count,1);
		}

	}
}


__global__ void Kill_copy(int *killed, float *stats)
{

    killed[0]=Kill_count;
    stats[0] = Mass_count;
    stats[1] = Vol_count;
}




/* Flags particles which meet the criteria below */
__global__ void Remove_Silo_Discharged_Particles(float3  *position_com, int *P_ID,int Num_Particles )
{
    uint index = blockIdx.x*blockDim.x  + threadIdx.x ;

	if( (index < Num_Particles) && P_ID[index]!=-1 )
	{
		if( position_com[index].y < SimParms.Silo_Kill_Height )
		{
			 position_com[index]=SimParms.max_size;
			 P_ID[index]=-1;
			 atomicAdd(&Kill_count,1);
		}

	}
}





/* Flags particles which meet the criteria below */
__global__ void Remove_PL_Volume_Particles(float3 *position_com, int *P_ID,uint *P_Type,int Num_Particles )
{
    uint index = blockIdx.x*blockDim.x  + threadIdx.x ;

	if( (index < Num_Particles) && P_ID[index]!=-1 )
	{
	 if( position_com[index].z <= SimParms.valid_zone_Z_Start || position_com[index].z >= SimParms.valid_zone_Z_End 
			|| position_com[index].y <= SimParms.valid_zone_Y_Start || position_com[index].y >= SimParms.valid_zone_Y_End ||
			position_com[index].x <= SimParms.valid_zone_X_Start || position_com[index].x >= SimParms.valid_zone_X_End )
		//if( position_com[index].z <= 41.0f ||position_com[index].z >= 300  || position_com[index].y <= 50.0f || position_com[index].y >= 1500.0f|| position_com[index].x <= 50.0f || position_com[index].x >= 1500.0f)

		{
			 position_com[index] = SimParms.max_size;
			 P_ID[index]=-1;
			 
			 atomicAdd(&Kill_count,1);
			 
			 if(SimParms.Simulation_Type==Pulp_Lift)
			 {
			  atomicAdd(&Mass_count,ParticleObject[P_Type[index]].mass);
			  atomicAdd(&Vol_count, ParticleObject[P_Type[index]].volume);
			 }
		}

	}
}





__global__ void AddForce()
{
	fc+=10.0f;
}

__global__ void SubForce()
{
	fc-=10.0f;
}

__global__ void ZeroForce()
{
	fc=0.0f;
}

__global__ void Hatch()
{
	Num_WorldObjects = Num_WorldObjects -1;
}

