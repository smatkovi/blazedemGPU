

/*---------------------------------------------------------------------------*/
/* (5) */
/*---------------------------------------------------------------------------*/
__global__ void Integrate_Euler_NonSymmetry(
		                       float3    *force_com_PP,
		                       float3    *force_ang_PP,

		                       float3    *force_com_wall,
		                       float3    *force_ang_wall,

		                       float3    *force_com_lifter,
		                       float3    *force_ang_lifter,

		                       float3    *position_com,
		                       Quaterion *position_ornt,

		                       float3    *velocity_com,
		                       float3    *velocity_ang,
		                       uint      *P_ObjectType,
		                       int       *P_ID,
		                       uint       Num_Particles,
		                       bool       get_surfacePoints,
							   int        gpu_offset=0 )
{

    uint index = blockIdx.x*blockDim.x  + threadIdx.x ;

	if( (index < Num_Particles) && P_ID[index + gpu_offset]>-1 )
	{

	float delta = SimParms.InitalDelta_t;



	/* acc is in cm so F/m= (N/kg)*/
    float Mass = ParticleObject[P_ObjectType[index + gpu_offset]].mass;


    /* 1. Net acceleration acting on particle for this step */
    float3 accelrat = (force_com_wall[index] + force_com_PP[index]) +SimParms.Force_Field*Mass;
//
//	if(P_ID[index+gpu_offset]==1932)
//			    	  {
//		float3 val = (accelrat/Mass)*SimParms.InitalDelta_t ;
//		               printf("Vel change PP %f %f %f \n", val.x, val.y, val.z);
//			    	  }
    

    if(Num_DynamicObjects>0)
    {
    	accelrat += force_com_lifter[index] ;
    }

    accelrat*=(1.0f/Mass);



	/* 2. Integrate the Translational Velocity "cm/s" a(t-1) + a(t) */
	if(gpu_offset==0)
	{
	 //PrintVectorD(position_com [index]);

      velocity_com [index] += (accelrat)*SimParms.InitalDelta_t;
      position_com [index] += velocity_com[index]*delta;


	}
	else
	{
      /* Start from the first index in the array */
	  velocity_com [index + gpu_offset]  +=  (accelrat)*SimParms.InitalDelta_t;
      position_com [index + gpu_offset]  +=  velocity_com[index+ gpu_offset]*delta;
	}

    
    if(SimParms.Mode2D)
    {
      velocity_com[ index ].z = 0.0f;
    }



    if( SimParms.Rotation )
    {

        float3 forceL = force_ang_PP[index] + force_ang_wall[index];

        if( Num_DynamicObjects > 0 )
        {
        	forceL += force_ang_lifter[index];
        }


      /* Compute acceleration and update angular velocity */
      if( dot(forceL,forceL) > 0.0f )
      {

 	     /* Collision Response Rotation */
          float3 ang_acc;

          if(SimParms.particle_type==spheres)
          {
             ang_acc   = ParticleObject[P_ObjectType[index + gpu_offset]].InertiaT[0]*forceL;
          }
          else
          {
            ang_acc = Get_Angular_Acc(index + gpu_offset, P_ObjectType[index + gpu_offset],forceL,position_ornt[index + gpu_offset]);
          }


		  if(gpu_offset==0)
	      {
             velocity_ang [index]  += ((ang_acc*SimParms.InitalDelta_t));

             //if(SimParms.particle_type==spheres)
             {
               /* Rolling resistance */
               velocity_ang[index]   -= velocity_ang[index]*SimParms.Roll_Res;
             }
		  }
		  else
		  {
			 velocity_ang [index+ gpu_offset]  += ((ang_acc*SimParms.InitalDelta_t));

             //if(SimParms.particle_type==spheres)
             {
               /* Rolling resistance */
               velocity_ang[index+ gpu_offset]   -= velocity_ang[index+ gpu_offset]*SimParms.Roll_Res;
             }
		  }


      }


      /* Angular integration */
      if( SimParms.particle_type==polyhedra || SimParms.sphere_orient )
      {
	    float3   omega    = velocity_ang[index+ gpu_offset ];
        float    magOmega = dot(omega,omega);

        Quaterion dQ = make_IquaterionD();

        if( magOmega > 0.000000f )
        {


      	    magOmega     = sqrt(magOmega);
      	    float3 dir   = omega/magOmega;
            dQ.w         = cos(0.5000f*magOmega*delta);
            double sinv  = sin(0.5000f*magOmega*delta);
            dQ.x         = sinv*dir.x;
            dQ.y         = sinv*dir.y;
            dQ.z         = sinv*dir.z;
        }
         /* load current orientation */
	     Quaterion C_Ornt       = position_ornt[index + gpu_offset];
	     Quaterion Q            = normaliseD((dQ*C_Ornt));/* Integrate */

	     /* Check for singularity at 360% */
	     if( Q.w <= -1.00000f)
	     {
	    	 position_ornt[ index] = dQ;
	     }
	     else
	     {
	      /* 2. Integrate the Angular Position */
	      position_ornt[ index] = Q;
	     }

	   	if(SimParms.Mode2D)
	   	{
	   	   position_ornt[ index].x = 0.0f;
	   	   position_ornt[ index].y = 0.0f;
	   	}

      }/* End  */

  	  if(SimParms.Mode2D)
  	  {
  	    velocity_ang [index].x   = 0.0f;
  	    velocity_ang [index].y   = 0.0f;
  	  }

    } /* End if rotation */


	 force_com_wall    [index] = make_float3(0.0f);

	 if(SimParms.Rotation)
	 {
	   force_ang_wall[index] = make_float3(0.0f);
	 }



 

   }/* End particles */

}
/*---------------------------------------------------------------------------*/



/*___________________________________________________________________________*/

__global__ void Integrate_Euler_Symmetry( float     *force_com_X,
		                                  float     *force_com_Y,
		                                  float     *force_com_Z,
		                                  float     *force_ang_X,
		                                  float     *force_ang_Y,
		                                  float     *force_ang_Z,

		                                  float3    *force_com_wall,
										  float3    *force_ang_wall,
										  float3    *force_com_lifter,
										  float3    *force_ang_lifter,

										  float3    *position_com,
										  Quaterion *position_ornt,
										  float3    *velocity_com,
										  float3    *velocity_ang,
										  uint      *P_ObjectType,
										  int       *P_ID,
										  uint       Num_Particles,
										  int  gpu_offset=0)
{

    
	    uint index = blockIdx.x*blockDim.x  + threadIdx.x + gpu_offset ;

	if( (index < Num_Particles) && P_ID[index]>-1 )
	{

	float delta = SimParms.InitalDelta_t;



	/* acc is in cm so F/m= (N/kg)*/
    float Mass = ParticleObject[P_ObjectType[index]].mass;


	    float3 force_com_PP = make_float3( force_com_X[index- gpu_offset], force_com_Y[index- gpu_offset], force_com_Z[index- gpu_offset]);

	    float3 force_ang_PP = make_float3( force_ang_X[index- gpu_offset], force_ang_Y[index- gpu_offset], force_ang_Z[index- gpu_offset]);


    /* 1. Net acceleration acting on particle for this step */
    float3 accelrat = (force_com_wall[index- gpu_offset] + force_com_PP) +SimParms.Force_Field*Mass;


    if(Num_DynamicObjects>0)
    {
    	accelrat += force_com_lifter[index- gpu_offset] ;
    }

    accelrat*=(1.0f/Mass);




	/* 2. Integrate the Translational Velocity "cm/s" a(t-1) + a(t) */
	if(gpu_offset==0)
	{
      velocity_com [index]  += (accelrat)*SimParms.InitalDelta_t;
      position_com [index] += velocity_com[index]*delta;
	}
	else
	{
      /* Start from the first index in the array */
	  velocity_com [index-gpu_offset]  = velocity_com [index] + (accelrat)*SimParms.InitalDelta_t;
      position_com [index-gpu_offset]  = position_com [index] + velocity_com[index-gpu_offset]*delta;
	}

    
    if(SimParms.Mode2D)
    {
      velocity_com[ index - gpu_offset ].z = 0.0f;
    }



    if( SimParms.Rotation )
    {

        float3 forceL = force_ang_PP + force_ang_wall[index];

        if( Num_DynamicObjects > 0 )
        {
        	forceL += force_ang_lifter[index];
        }


      /* Compute acceleration and update angular velocity */
      if( dot(forceL,forceL) > 0.0f )
      {

 	     /* Collision Response Rotation */
          float3 ang_acc;

          if(SimParms.particle_type==spheres)
          {
             ang_acc   = ParticleObject[P_ObjectType[index]].InertiaT[0]*forceL;
          }
          else
          {
            ang_acc = Get_Angular_Acc(index,P_ObjectType[index],forceL,position_ornt[index]);
          }

		  if(gpu_offset==0)
	      {
             velocity_ang [index]  += ((ang_acc*SimParms.InitalDelta_t));

             if(SimParms.particle_type==spheres)
             {
               /* Rolling resistance */
               velocity_ang[index]   -= velocity_ang[index]*SimParms.Roll_Res;
             }
		  }
		  else
		  {
			  velocity_ang [index-gpu_offset]  = velocity_ang [index] + ((ang_acc*SimParms.InitalDelta_t));

             if(SimParms.particle_type==spheres)
             {
               /* Rolling resistance */
               velocity_ang[index-gpu_offset]   -= velocity_ang[index-gpu_offset]*SimParms.Roll_Res;
             }
		  }


      }


      /* Angular integration */
      if( SimParms.particle_type==polyhedra || SimParms.sphere_orient )
      {
	    float3   omega    = velocity_ang[index];
        float    magOmega = dot(omega,omega);

        Quaterion dQ = make_IquaterionD();

        if( magOmega > 0.000000f )
        {


      	    magOmega     = sqrt(magOmega);
      	    float3 dir   = omega/magOmega;
            dQ.w         = cos(0.5000f*magOmega*delta);
            double sinv  = sin(0.5000f*magOmega*delta);
            dQ.x         = sinv*dir.x;
            dQ.y         = sinv*dir.y;
            dQ.z         = sinv*dir.z;
        }
         /* load current orientation */
	     Quaterion C_Ornt       = position_ornt[index];
	     Quaterion Q            = normaliseD((dQ*C_Ornt));/* Integrate */

	     /* Check for singularity at 360% */
	     if( Q.w <= -1.00000f)
	     {
	    	 position_ornt[ index - gpu_offset] = dQ;
	     }
	     else
	     {
	      /* 2. Integrate the Angular Position */
	      position_ornt[ index - gpu_offset] = Q;
	     }

	   	if(SimParms.Mode2D)
	   	{
	   	   position_ornt[ index - gpu_offset ].x = 0.0f;
	   	   position_ornt[ index  - gpu_offset].y = 0.0f;
	   	}

      }/* End  */

  	  if(SimParms.Mode2D)
  	  {
  	    velocity_ang [index- gpu_offset].x   = 0.0f;
  	    velocity_ang [index- gpu_offset].y   = 0.0f;
  	  }

    } /* End if rotation */


		/* Set PP forces to zero */
		 force_com_X[index]=0.0f;
		 force_com_Y[index]=0.0f;
		 force_com_Z[index]=0.0f;

		 force_ang_X[index]=0.0f;
		 force_ang_Y[index]=0.0f;
		 force_ang_Z[index]=0.0f;

		 
    force_com_wall  [index- gpu_offset] = make_float3(0.0f);
    force_ang_wall  [index- gpu_offset] = make_float3(0.0f);

    if( Num_DynamicObjects > 0 )
    {
      force_com_lifter    [index- gpu_offset] = make_float3(0.0f);
      force_ang_lifter    [index- gpu_offset] = make_float3(0.0f);
    }

   }/* End particles */


}
