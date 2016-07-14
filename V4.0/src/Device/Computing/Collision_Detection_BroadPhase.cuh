/* Contains the Kernels and methods that perform
 * Spatial sub-division of the problem  */


/*-----------------------------------------------------------------------------*/
                     /* DEVICE COMPUTING METHODS */
/*-----------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------*/
                      /* Given a Point: p
                       * Returns it position in the
                       * computational grid */

/* !!Cost per Call!!
 *    Storage       : 6 Bytes
 *    Memory        : 6 Reads (Const)
 *    Computational : 3 Sub,3 Div, 1 Floor */
/*-----------------------------------------------------------------------------*/
__device__ int3 SpatialDecomp_CalcGridPos(float3 p)
{
    int3 gridPos;

    gridPos.x = floor((p.x - SimParms.worldOrigin.x) / SimParms.cellSize.x);
    gridPos.y = floor((p.y - SimParms.worldOrigin.y) / SimParms.cellSize.y);
    gridPos.z = floor((p.z - SimParms.worldOrigin.z) / SimParms.cellSize.z);

    return gridPos;
}
/*-----------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------*/
                     /* Given a grid position
                      * Returns a Hash based on the spatial location */
/* !!Cost per Call!!
 *    Storage       : 0 Bytes
 *    Memory        : 3 Reads (Const)
 *    Computational : 3 Add, 3 Mul */

/*-----------------------------------------------------------------------------*/
__device__ uint SpatialDecomp_CalcGridHash(int3 gridPos)
{

    return ( (gridPos.x + gridPos.z*SimParms.num_NNCells.z) +
			 (gridPos.y*SimParms.num_NNCells.x)*SimParms.num_NNCells.z);


}
/*-----------------------------------------------------------------------------*/





/*-----------------------------------------------------------------------------*/
                       /* GLOBAL Kernels */
/*-----------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------*/
         /* Calculate grid hash value for all particles */
         /* 1) Input   : numberParticles
          * 2) Returns : Hash of each particle : gridParticleHash_R
          * 3) Returns : Index of each particle: gridParticleIndex_R
          * This linked list allows sorting to get better hit rate */


/* !!Cost per thread!!
 *    Storage       : 9 + 6 Bytes
 *    Memory        : Reads:# 1 (global-device) 9 (constant), Writes:# 2 (global-malloc)
 *    Computational : ADD:# 1+3-3, MUL:# 1+3, DIV:# +1,FLOOR:# +1 */
/* Constant Memory IO */
/*-----------------------------------------------------------------------------*/

__global__ void Kernel_SpatialDecomp_CalcPHash( uint   *gridParticleHash,
		                                        int    *gridParticleIndex,
		                                        float3 *position_com,
		                                        uint   Num_Particles,
												int    gpu_offset=0)

{
	uint index = blockIdx.x*blockDim.x  + threadIdx.x;



    if (index < Num_Particles )/* stop at max_index-1 */
    {

		if(index==0)
		{
          Kill_count = 0;
		}

    	//printf("%d %f \n",index,position_com[index].x );

    	//int3 GridPos = SpatialDecomp_CalcGridPos( position_com[index+ gpu_offset]);

    /* store grid hash and particle index */
      gridParticleHash [index + gpu_offset]   = SpatialDecomp_CalcGridHash( SpatialDecomp_CalcGridPos( position_com[index + gpu_offset]));
      gridParticleIndex[index + gpu_offset]   = index;



    }

}
/*-----------------------------------------------------------------------------*/


__global__ void Kernel_SpatialDecomp_CalcPHash( uint   *gridParticleHash,
		                                        int    *gridParticleIndex,
		                                        float3 *position_com,
		                                        uint   Num_Particles,
		                                        int    *P_ID,
												int    gpu_offset=0)

{
	uint index = blockIdx.x*blockDim.x  + threadIdx.x;



    if (index < Num_Particles )/* stop at max_index-1 */
    {

		if(index==0)
		{
          Kill_count = 0;
		}

    	//printf("%d %f \n",index,position_com[index].x );

    	int3 GridPos = SpatialDecomp_CalcGridPos( position_com[index+ gpu_offset]);

    /* store grid hash and particle index */
      gridParticleHash [index + gpu_offset]   = SpatialDecomp_CalcGridHash(GridPos);
      gridParticleIndex[index + gpu_offset]   = index;

      /* Pre-selection for certain planar geometries */
      if( GridPos.x<2 || GridPos.y<2  || GridPos.z<2   )
      {
    	  P_ID[index + gpu_offset]+= Num_Particles;
      }


    }

}









/*-----------------------------------------------------------------------------*/
             /*   Uses Thrust Sort Method(uses existing device arrays)
                * 1) Ascending based on: dGridParticleHash
                * 2) Linked Sorted     : dGridParticleIndex  */
/*-----------------------------------------------------------------------------*/
void SpatialDecomp_ThrustSort_ByHash( uint *dGridParticleHash, int *dGridParticleIndex,
		                                                       uint numParticles     )
{
    thrust::stable_sort_by_key(thrust:: device_ptr<uint>(dGridParticleHash),
                                thrust:: device_ptr<uint>(dGridParticleHash + numParticles),
                                thrust:: device_ptr<int>(dGridParticleIndex));

}
/*-----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------*/
             /*   Uses Thrust Sort Method(uses existing device arrays)
                * 1) Ascending based on: dGridParticleHash
                * 2) Linked Sorted     : dGridParticleIndex  */
/*-----------------------------------------------------------------------------*/
void SpatialDecomp_ThrustSort_ByHashMulti( uint *dGridParticleHash, int *dGridParticleIndex,
		                                                       uint numParticles, cudaStream_t stream )
{
    thrust::stable_sort_by_key(thrust::cuda::par.on(stream), thrust:: device_ptr<uint>(dGridParticleHash),
                                thrust:: device_ptr<uint>(dGridParticleHash + numParticles),
                                thrust:: device_ptr<int>(dGridParticleIndex));

}
/*-----------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------*/
/* Note we keep history between particles and lifter for smooth edge contact */
__global__ void Kernel_SortArrays   ( int      *dGridParticleIndex,
					                  float3    *SortedPos_R,
					                  Quaterion *SortedPos_RQ,
					                  float3    *SortedVel_R,
					                  float3    *SortedVel_Rang,
					                  uint      *SortedPType_R,
					                  int       *SortedPID_R,

					                  Contact_Info      *SortedLifter_Contact_Hist_R,
									  float3    *position_com,
									  Quaterion *position_ornt,
									  float3    *velocity_com,
									  float3    *velocity_ang,
									  uint      *P_ObjectType,
									  int       *P_ID,
									  Contact_Info     *Lifter_Contact_Hist,
									  uint   Num_Particles)
{
    uint index = blockIdx.x*blockDim.x + threadIdx.x;

    if (index < Num_Particles )/* stop at max_index-1 */
    {

      /* dGridParticleIndex contains the address of the element in
          * the Particle Dynamics Array (PDA) */

       /* Get the sorted particle address from PDA */
       uint   sortedIndex = dGridParticleIndex[index];

       /* Create a sorted PDA in temp storage as we cant do a
        * direct assignment (similar to temp when sorting)  */
       SortedPos_R    [index]  = position_com  [sortedIndex];
       SortedVel_R    [index]  = velocity_com  [sortedIndex];
       SortedVel_Rang [index]  = velocity_ang  [sortedIndex];
       SortedPType_R  [index]  = P_ObjectType  [sortedIndex];
       SortedPID_R    [index]  = P_ID          [sortedIndex];

       if(SimParms.particle_type==polyhedra || SimParms.sphere_orient)
       {
         SortedPos_RQ  [index]  = position_ornt [sortedIndex];
       }

       SortedLifter_Contact_Hist_R    [index] = Lifter_Contact_Hist   [sortedIndex];

    }

}
/*-----------------------------------------------------------------------------*/




/*---------------------------------------------------------------------------*/
               /* Updates PDA will sorted list to
                * make memory accesses more efficient
                * 1) Input: All previously allocated*/
/* !!Cost per thread!!
 *    Storage       : 2 Bytes
 *    Memory        : Reads:# 4 (global-malloc), Writes:# 4 (global-device)
 *    Computational : ADD:# 6-1, MUL:# 1,
 *    */
/*---------------------------------------------------------------------------*/
__global__ void Kernel_SpatialDecomp_ReorderPDA( float3 *SortedPos,    Quaterion *SortedPosQ,
		                                         float3 *SortedVel,    float3 *SortedVelQ,
		                                         uint   *SortedPType,  int *SortedPID,
		                                         Contact_Info   *SortedLifter_Contact_Hist,

		                                         float3 *position_com, Quaterion *position_ornt,
		                                         float3 *velocity_com, float3 *velocity_ang,
		                                         uint   *P_ObjectType, int *P_ID,
		                                         Contact_Info   *Lifter_Contact_Hist,
		                                         uint   Num_Particles)
{
   uint index = blockIdx.x*blockDim.x  + threadIdx.x ;

   if( index < Num_Particles )
   {

	 /* update */
     position_com  [index] = SortedPos   [index];
     velocity_com  [index] = SortedVel   [index];

     if(SimParms.particle_type==1 || SimParms.sphere_orient )
     {
       position_ornt [index] = SortedPosQ  [index];
     }



     velocity_ang [index] = SortedVelQ  [index];

     P_ObjectType [index] = SortedPType [index];
     P_ID         [index] = SortedPID   [index];

     Lifter_Contact_Hist    [index] = SortedLifter_Contact_Hist[index];


   }

}

/*---------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------*/
               /* Bins particles into Grid Cells based on the Hash
                * and reorders based on sorting
                * 1) Input  : dGridParticleHash, dGridParticleIndex
                * 2) Output : Cell bin address: cellStart_R, cellEnd_R
                * 3) Output : Sorted Dynamics Info */
/*-----------------------------------------------------------------------------*/
__global__ void Kernel_SpatialDecomp_BinData( uint *dGridParticleHash,
		                                      uint *cellStart_R, uint *cellEnd_R,
		                                      uint   Num_Particles)
{
    uint mem_index = blockIdx.x*blockDim.x + threadIdx.x;


    if (mem_index < Num_Particles )/* stop at max_index-1 */
    {


      /* index gives the index of a particle in the array */
      uint hash;
      uint nexthash;

      hash    = dGridParticleHash[mem_index];

      if ( (mem_index < Num_Particles-1) )/* stop at max_index-1 */
      {
         /* load next hash in list */
         nexthash = dGridParticleHash[mem_index+1];

         /* First entry so cell start */
         if( mem_index==0 )
         {
           cellStart_R[hash] = mem_index;
         }

         /* if the hashes are different its a new cell*/
         if ( hash != nexthash )
         {
    	   /* End of will be at address (index+1) in sorted List */
    	   cellEnd_R[hash]       = mem_index+1;
    	   /* Start will be at address (index+1) in sorted List */
    	   cellStart_R[nexthash] = mem_index+1;
         }

      }
        /* Last particle must be the cell end */
      if ( mem_index == SimParms.Num_Particles - 1 )
      {
         cellEnd_R[hash] = mem_index + 1;
      }

    }

}
/*-----------------------------------------------------------------------------*/



__global__ void Kernel_SpatialDecomp_BinDataR( uint *dGridParticleHash,
		                                      uint *cellStart_R, uint *cellEnd_R,
		                                      uint   Num_Particles)
{
    uint mem_index = blockIdx.x*blockDim.x + threadIdx.x;


    if (mem_index < Num_Particles )/* stop at max_index-1 */
    {


      /* index gives the index of a particle in the array */
      uint hash;
      uint nexthash;

      hash    = dGridParticleHash[mem_index];


      if ( (mem_index < Num_Particles-1) )/* stop at max_index-1 */
      {
         /* load next hash in list */
         nexthash = dGridParticleHash[mem_index+1];

         /* First entry so cell start */
         if( mem_index==0 )
         {
           cellStart_R[hash] = 0;
         }

         /* if the hashes are different its a new cell*/
         if ( hash != nexthash )
         {
    	   /* End of will be at address (index+1) in sorted List */
    	   cellEnd_R[hash]       = 0;
    	   /* Start will be at address (index+1) in sorted List */
    	   cellStart_R[nexthash] = 0;
         }

      }
        /* Last particle must be the cell end */
      if ( mem_index == SimParms.Num_Particles - 1 )
      {
         cellEnd_R[hash] = 0;
      }

    }

}




/*-----------------------------------------------------------------------------*/
               /* 1) Check if two spheres are colliding  */
/*-----------------------------------------------------------------------------*/
__device__ NN_Selected_Cell CollisionDetection_Filter_SphereNNCell(
		                                                      int3   NNgridPos,
		                                                      uint   index,
                                                              uint   *cellStart,
                                                              uint   *cellEnd,
                                                              float3 *position_com,
                                                              uint   *P_ObjectType )
{
    int gridHash = SpatialDecomp_CalcGridHash(NNgridPos);

    /* get start of bucket for this cell */
    int startIndex = cellStart[gridHash];

     NN_Selected_Cell Selected_mem_index;

    int num_sel=0;

    /* check cell is not empty */
    if (startIndex != 0xffffffff)
    {


        /* Get end of the cell */
        int endIndex = cellEnd[gridHash];

        /* Get the COM Position of Current Particle */

        float3 PosA = position_com[index];

        /* Get the Bound Radius Current Particle Type */
        float  R_A = ParticleObject[P_ObjectType[index]].radius;



        /* Loop over all entries  in the cell and compute NN */
        for (int j=startIndex; j<endIndex; j++)
        {

           if (j != index) /* check not colliding with self */
           {
               /* Check to see if there is an overlap */
			   if( (dot(PosA-position_com[j],PosA-position_com[j]) -
				   (( R_A + ParticleObject[ P_ObjectType[j]].radius)*( R_A + ParticleObject[ P_ObjectType[j]].radius)) ) < -0.0000100f )
               {
                  Selected_mem_index.NN_Index[num_sel+1] = j; /* mem index of neighbour */
                  num_sel++;
               }

            }

        }/* End loop over Current cell */

    }/* End checking cell */


    /* First value is how many selected */
    Selected_mem_index.NN_Index[0] = num_sel;

    return Selected_mem_index;
}
/*-----------------------------------------------------------------------------*/




/*-----------------------------------------------------------------------------*/
               /* 1) Check if two spheres are colliding  */
/*-----------------------------------------------------------------------------*/
__device__ Forces CollisionDetection_Filter_SphereNNCell_Force(
	                                                      
		                                                      int3   NNgridPos,
		                                                      uint   index,
                                                              uint   *cellStart,
                                                              uint   *cellEnd,
                                                              float3 Pos_A,
                                                              float3 Vel_A,
                                                              float3 VelAng_A,
                                                              int    P_typeA,
                                                              float  R_A,
                                                              float3 *position_com,
													          float3    *velocity_com,
							                                  float3    *velocity_ang,
							                                  uint      *P_ObjectType,
															  bool    get_surfacePoints,
                                                              float3  *dContactPoints,
                                                              int *ncon)
{



    int gridHash = SpatialDecomp_CalcGridHash(NNgridPos);

    /* get start of bucket for this cell */
    int startIndex = cellStart[gridHash];

	Forces Force;

	Force.trans  = make_float3(0.0f);
	Force.torque = make_float3(0.0f);

	Forces ForceC;

	ForceC.trans  = make_float3(0.0f);
	ForceC.torque = make_float3(0.0f);

	int num_contacts=0;

	/* check cell is not empty */
    if (startIndex != 0xffffffff)
    {

        /* Get end of the cell */
        int endIndex = cellEnd[gridHash];

        /* Get the COM Position of Current Particle */



		

        /* Loop over all entries  in the cell and compute NN */
        for (int j=startIndex; j<endIndex; j++) /* mem-locations of particles in the cell */
        {
        
           if (j != index) /* check not colliding with self */
           {


			   int    P_typeB  = P_ObjectType[j];
			   float  R_B      = ParticleObject[P_typeB].radius; 

		        float3 relPos     = Pos_A - position_com[j]  ;

               /* Check to see if there is an overlap */
			   if( (dot(relPos,relPos) -
				   (( R_A + R_B)*( R_A + R_B))) < -0.0000100f )
               {

	                 {

				 		 ForceC = Collision_Response_Particle_Particle(   Vel_A-velocity_com[j],
                                  relPos/length(relPos),
                                   (abs( (length(relPos) -
                      					 (R_A + R_B ) ) ) + 0.000010f),
                                  (-1.0f*(relPos/length(relPos))*R_A),
                                  ((relPos/length(relPos))*R_B),
					                  VelAng_A,
					              velocity_ang[j],
			  					      P_typeA,
			  				          P_typeB,
			  					      Vel_A,
			  					      velocity_com[j],
			  					      ParticleObject[P_typeA].mass,
			  					      ParticleObject[P_typeB].mass,
			  					relPos);



				 		Force += ForceC;

				 		 dContactPoints[num_contacts] = Pos_A + ( (relPos/length(relPos))*R_A);

				 		 dContactPoints[num_contacts+8] = ForceC.trans;

		             }
				 	num_contacts++;


               }

            }

        }/* End loop over Current cell */

    }/* End checking cell */

    ncon[0]=num_contacts;

	return Force;
}
/*-----------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------*/
       /* 1) Check if two spheres are colliding and return force only */
/*-----------------------------------------------------------------------------*/
__device__ Forces CollisionDetection_Filter_SphereNNCell_Force(

		                                                      int3   NNgridPos,
		                                                      uint   index,
                                                              uint   *cellStart,
                                                              uint   *cellEnd,
                                                              float3 Pos_A,
                                                              float3 Vel_A,
                                                              float3 VelAng_A,
                                                              int    P_typeA,
                                                              float  R_A,
                                                              float3 *position_com,
													          float3    *velocity_com,
							                                  float3    *velocity_ang,
							                                  uint      *P_ObjectType)
{
    int gridHash = SpatialDecomp_CalcGridHash(NNgridPos);

    /* get start of bucket for this cell */
    int startIndex = cellStart[gridHash];

	Forces Force;

	Force.trans  = make_float3(0.0f);
	Force.torque = make_float3(0.0f);

	/* check cell is not empty */
    if (startIndex != 0xffffffff)
    {

        /* Get end of the cell */
        int endIndex = cellEnd[gridHash];

        /* Loop over all entries  in the cell and compute NN */
        for (int j=startIndex; j<endIndex; j++) /* mem-locations of particles in the cell */
        {

           if (j != index) /* check not colliding with self */
           {


			   int    P_typeB  = P_ObjectType[j];
			   float  R_B      = ParticleObject[P_typeB].radius;

		       float3 relPos     = Pos_A - position_com[j]  ;

               /* Check to see if there is an overlap */
			   if( (dot(relPos,relPos) -
				   (( R_A + R_B)*( R_A + R_B))) < -0.0000100f )
               {


				 	 {
				         /* Compute force exerted on Particle A */
						 Force += Collision_Response_Particle_Particle(   Vel_A-velocity_com[j],
				          		                                      relPos/length(relPos),
				          		                                    (fabs( (length(relPos) -
				          		                                    					 (R_A + R_B ) ) ) + 0.000010f),
				          		                                      (-1.0f*(relPos/length(relPos))*R_A),
				          		                                      ((relPos/length(relPos))*R_B),
				          		  					                  VelAng_A,
				          		  					              velocity_ang[j],
				          		  			  					      P_typeA,
				          		  			  				          P_typeB,
				          		  			  					      Vel_A,
				          		  			  					      velocity_com[j],
				          		  			  					      ParticleObject[P_typeA].mass,
				          		  			  					      ParticleObject[P_typeB].mass,
				          		  			  					relPos);


				 	 }

               }

            }

        }/* End loop over Current cell */

    }/* End checking cell */

	return Force;
}
/*-----------------------------------------------------------------------------*/



/*-----------------------------------------------------------------------------*/
       /* Broad Phase Collision Test1: Bounding Sphere Distance
        * Stores List of NN for each Particle, checks in parallel */
/*-----------------------------------------------------------------------------*/
__global__ void Kernel_BroadCollisionDetection_NonSymmetry(
                          uint   *cellStart,
                          uint   *cellEnd ,
                          uint   *Broad_List,
                          uint   *NN_NUM,
                          float3 *position_com,
                          uint   *p_ObjectID,
                          int    *P_ID,
                          uint   Num_Particles,
						  uint gpu_offset=0)
{

	/* get the unique memory location index  */
    uint mem_index = blockIdx.x*blockDim.x + threadIdx.x ;

   if( ( mem_index < Num_Particles) && P_ID[mem_index+gpu_offset] > -1 )
   {

	   NN_NUM[mem_index] = 0;
	   Broad_List[(mem_index)*32] = 0;

    /* Read particle pos from sorted arrays */
    float3 pos = position_com[mem_index + gpu_offset];

    /* Get hashed address in grid */
    int3 gridPos = SpatialDecomp_CalcGridPos(pos);

    /* examine neighboring cells */
    /* check NN -1 : +1 */

    int zstart = -1;
    int zend   =  1;
    int ystart = -1;
    int yend   =  1;
    int xstart = -1;
    int xend   =  1;

    if(gridPos.x==0)
    {
      xstart = 0;
    }
    if(gridPos.y==0)
    {
      ystart = 0;
    }
    if(gridPos.z==0)
    {
      zstart = 0;
    }

    if(gridPos.x==(SimParms.num_NNCells.x-1))
    {
      xend = 0;
    }
    if(gridPos.y==(SimParms.num_NNCells.y-1))
    {
      yend = 0;
    }
    if(gridPos.z==(SimParms.num_NNCells.z-1))
    {
      zend = 0;
    }

    NN_Selected_Cell Selected_ParticleMemIndex;

    int counter=0;

    /* Search through each of the NN Cells */
    for (int y = ystart; y <= yend; y++)
    {
        for (int z = zstart; z <= zend ; z++)
        {
            for (int x=  xstart; x <= xend; x++)
            {

                Selected_ParticleMemIndex =
                		CollisionDetection_Filter_SphereNNCell( gridPos + make_int3(x, y, z), mem_index+ gpu_offset,cellStart, cellEnd, position_com, p_ObjectID);


                 int3 cellc =  gridPos + make_int3(x, y, z);


                for(int i=0; i <Selected_ParticleMemIndex.NN_Index[0];i++)
                {


                  Broad_List[(mem_index)*32 + (counter+i)]= Selected_ParticleMemIndex.NN_Index[i+1];

                }

                counter = counter + Selected_ParticleMemIndex.NN_Index[0]; /*location 0 is num being returned*/

            }

        }
    }

    NN_NUM[(mem_index)] = counter;


    if( counter > 32 )
    {
      printf("ERROR: MAX num neighbors(32) exceeded \n");
      /* Throw an Exception when supported by Device */
    }

    //printf("Pos %f %f %f  Grid Pos %d %d %d  : Particle %d has %d NN \n ", position_com[mem_index].x, position_com[mem_index].y, position_com[mem_index].z,
  		//  gridPos.x,gridPos.y,gridPos.z,  P_ID[mem_index],NN_NUM[mem_index]);

  }/* thread check */

}



/*-----------------------------------------------------------------------------*/
       /* Broad Phase Collision Test1: Bounding Sphere Distance
        * Stores List of NN for each Particle, checks in parallel */
/*-----------------------------------------------------------------------------*/
__global__ void Kernel_Spheres_PP(
		                  bool      get_surfacePoints,
                          float3    *dWallContactPoints,
                          uint   *cellStart,
                          uint   *cellEnd ,
						  float3    *position_com,
					      float3    *velocity_com,
						  float3    *velocity_ang,
						  uint      *P_ObjectType,
						  int       *P_ID  ,
						  float3    *force_com,
						  float3    *force_ang,
                          uint   Num_Particles,
                          uint   *Num_contacts,
						  uint gpu_offset=0)
{

	/* get the unique memory location index  */
    uint mem_index = blockIdx.x*blockDim.x + threadIdx.x ;

   if( ( mem_index < Num_Particles) && P_ID[mem_index+gpu_offset] > -1 )
   {

    /* Read particle pos from sorted arrays */
    float3 Pos_A    = position_com[mem_index + gpu_offset];
	float3 Vel_A    = velocity_com[mem_index + gpu_offset];
	float3 VelAng_A = velocity_ang[mem_index + gpu_offset];

    /* Get the Bound Radius Current Particle Type */

	int    P_typeA = P_ObjectType[mem_index + gpu_offset];
	float  R_A = ParticleObject[P_typeA].radius;


    /* Get hashed address in grid */
    int3 gridPos = SpatialDecomp_CalcGridPos(Pos_A);

    /* examine neighboring cells */
    /* check NN -1 : +1 */

    int zstart = -1;
    int zend   =  1;
    int ystart = -1;
    int yend   =  1;
    int xstart = -1;
    int xend   =  1;

    if(gridPos.x==0)
    {
      xstart = 0;
    }
    if(gridPos.y==0)
    {
      ystart = 0;
    }
    if(gridPos.z==0)
    {
      zstart = 0;
    }

    if(gridPos.x==(SimParms.num_NNCells.x-1))
    {
      xend = 0;
    }
    if(gridPos.y==(SimParms.num_NNCells.y-1))
    {
      yend = 0;
    }
    if(gridPos.z==(SimParms.num_NNCells.z-1))
    {
      zend = 0;
    }


	Forces Force_P;
	Force_P.trans  = make_float3(0.0f);
	Force_P.torque = make_float3(0.0f);



	float3 contact_point[16];

	int n_con[1];
	uint n_contacts=0;



    /* Search through each of the NN Cells */
    for (int y = ystart; y <= yend; y++)
    {
        for (int z = zstart; z <= zend ; z++)
        {
            for (int x=  xstart; x <= xend; x++)
            {

            	if(get_surfacePoints)
            	{
				 Force_P += CollisionDetection_Filter_SphereNNCell_Force( gridPos + make_int3(x, y, z),
					 mem_index + gpu_offset,
					           cellStart,
					           cellEnd,
					           Pos_A,
					           Vel_A,
					           VelAng_A,
					           P_typeA,
					           R_A,
							   position_com,
							   velocity_com,
							   velocity_ang,
							   P_ObjectType,
							   get_surfacePoints,
							   contact_point,
							   n_con);

				if(n_con[0]>0 )
				{
					for(int i=0;i<n_con[0];i++)
					{
						dWallContactPoints[mem_index*64 + n_contacts ]      = contact_point[i]; /* Point */
						dWallContactPoints[mem_index*64 + n_contacts + 32   ] = contact_point[i+8]; /* Force */
						n_contacts++;
					}
				 }
            	}
            	else
            	{


      				 Force_P += CollisionDetection_Filter_SphereNNCell_Force( gridPos + make_int3(x, y, z),
      					 mem_index + gpu_offset,
      					           cellStart,
      					           cellEnd,
   					           Pos_A,
   					           Vel_A,
   					           VelAng_A,
   					           P_typeA,
   					           R_A,
      							   position_com,
      							   velocity_com,
      							   velocity_ang,
      							   P_ObjectType );

            	}

            }

        }
    }


	force_com[mem_index] = Force_P.trans;

    if( SimParms.Rotation )
	{
	    force_ang[mem_index] = Force_P.torque;
	}
	else
	{
	    force_ang[mem_index] = make_float3(0.0f);
	}

    if(get_surfacePoints)
    {
    	Num_contacts[mem_index] = n_contacts;
		 if(n_contacts>27)
		 {
			 printf("too many links %d \n", n_contacts);

		 }
    }


	
  }/* thread check */

}









                          /* Need to verify this */

/*-----------------------------------------------------------------------------*/
               /* 1) Check if two spheres are colliding  */
/*-----------------------------------------------------------------------------*/
__device__ NN_Selected_Cell CollisionDetection_Filter_2SphereSameCell( int3 OwngridPos,    uint index,
									 uint *cellStart, uint *cellEnd,float3 *position_com,uint *P_ObjectType, int *P_ID)
{
    uint gridHash = SpatialDecomp_CalcGridHash(OwngridPos);

    /* get start of bucket for this cell */
    uint startIndex = cellStart[gridHash];
    /* Get end of the cell */
    uint endIndex = cellEnd[gridHash];


     NN_Selected_Cell Selected_ID;

    int num_sel=0;
   // float delta = SimParms.InitalDelta_t;


   /* Only if we are the first particle in the cell */
   int smallest= SimParms.Num_Particles+1;

   for(uint j=startIndex; j<endIndex; j++)
   {
      if (P_ID[j]<smallest)
      {
         smallest =P_ID[j];
       }
    }

    if(P_ID[index]==smallest)
    {
        /* Get the COM Position of Current Particle */

        float3 PosA = position_com[index];

        /* Get the Bound Radius Current Particle Type */
        float  R_A = ParticleObject[P_ObjectType[index]].radius;


        /* Loop over all entries  in the cell and compute NN */
        for (uint j=startIndex; j<endIndex; j++)
        {

           if (j != index) /* check not colliding with self */
            {

               /* Check to see if there is an overlap */
			   if( (dot(PosA-position_com[j],PosA-position_com[j]) -
				   (( R_A + ParticleObject[ P_ObjectType[j]].radius)*( R_A + ParticleObject[ P_ObjectType[j]].radius))) < -0.0000100f )
               {
				   //printf(" Contact %d : with %d  Val %f  dis %f rad %f \n", P_ObjectType[index],P_ObjectType[j],val,dis_Ps,dis_rads );
                  Selected_ID.NN_Index[num_sel+1] = j; /* mem index of neighbour */
                  num_sel++;
               }

            }

        }/* End loop over Current cell */
    	}
    /* First value is how many selected */
    Selected_ID.NN_Index[0] = num_sel;

    return Selected_ID;
}
/*-----------------------------------------------------------------------------*/



__device__ NN_Selected_Cell CollisionDetection_Filter_SphereSameCellMany( int3 OwngridPos,    uint index,
									 uint *cellStart, uint *cellEnd,float3 *position_com,uint *P_ObjectType, int *P_ID)
{

NN_Selected_Cell  Selected_ID;
int   num_sel=0;
//float delta = SimParms.InitalDelta_t;

uint gridHash = SpatialDecomp_CalcGridHash(OwngridPos);

    /* Get end of the cell */
uint endIndex = cellEnd[gridHash];


/* Get the COM Position of Current Particle */
float3 PosA = position_com[index];// + delta*velocity_com[index] ;
/* Get the Bound Radius Current Particle Type */
float R_A = ParticleObject[P_ObjectType[index]].radius;

for(int j=index+1; j<endIndex;j++)
{

               /* Check to see if there is an overlap */
			   if( (dot(PosA-position_com[j],PosA-position_com[j]) -
				   (( R_A + ParticleObject[ P_ObjectType[j]].radius)*( R_A + ParticleObject[ P_ObjectType[j]].radius))) < -0.0000100f )
               {
				   //printf(" Contact %d : with %d  Val %f  dis %f rad %f \n", P_ObjectType[index],P_ObjectType[j],val,dis_Ps,dis_rads );
                  Selected_ID.NN_Index[num_sel+1] = j; /* mem index of neighbour */
                  num_sel++;
               }

}/* End loop thru cell*/

    /* First value is how many selected */
Selected_ID.NN_Index[0] = num_sel;

    return Selected_ID;
}
/*-----------------------------------------------------------------------------*/



/*-----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------*/
__global__ void Kernel_BroadCollisionDetection_Symmetry(
									  uint   *cellStart,
									  uint   *cellEnd ,
									  uint   *Broad_List,
									  uint   *NN_NUM,
									  float3 *position_com,
									  uint   *P_ObjectID,
									  int    *P_ID,
									  int     Num_Particles  )
{

/* get the unique index */
uint index = blockIdx.x*blockDim.x + threadIdx.x;

if( index < Num_Particles && P_ID[index]>-1)
{

		   NN_NUM[index] = 0;
	   Broad_List[index*32] = 0;

/* read particle data from sorted arrays */
float3 pos_A = position_com[index];


/* get address in grid */
int3 gridPos = SpatialDecomp_CalcGridPos(pos_A);


/* examine neighboring cells */
/* in a forward direction */

int zstart =  0;
int zend   =  1;

int ystart =  0;
int yend   =  1;

int xstart =  0;
int xend   =  1;

/* Snap to corners */
if(gridPos.x==(SimParms.num_NNCells.x-1))
{
xend = 0;
}
if(gridPos.y==(SimParms.num_NNCells.y-1))
{
yend = 0;
}
if(gridPos.z==(SimParms.num_NNCells.z-1))
{
zend = 0;
}

NN_Selected_Cell Selected_Index;


int counter=0;


for (int y = ystart; y <= yend; y++)
{
for (int z = zstart; z <= zend ; z++)
{
    for (int x=  xstart; x <= xend; x++)
    {

      if((x+y+z)>0)/* dont check own cell */
      {
        int3 neighbourPos = gridPos + make_int3(x, y, z);


        Selected_Index = CollisionDetection_Filter_SphereNNCell( neighbourPos, index,cellStart, cellEnd, position_com, P_ObjectID);

        for(int i=0; i <Selected_Index.NN_Index[0];i++)
        {
        	Broad_List[index*32 + (counter+i) ] = Selected_Index.NN_Index[i+1];

        }
        counter = counter + Selected_Index.NN_Index[0];
      }




    }
}
}


/* Need to include the diagonal neighbors */

Selected_Index.NN_Index[0]=0;


if( gridPos.x!=0 && gridPos.z!=(SimParms.num_NNCells.z-1))
{
int3 neighbourPos = gridPos + make_int3(-1, 0, 1);

Selected_Index = CollisionDetection_Filter_SphereNNCell( neighbourPos, index,cellStart, cellEnd,position_com, P_ObjectID );


for(int i=0; i <Selected_Index.NN_Index[0];i++)
{
   Broad_List[index*32 + (counter+i)] = Selected_Index.NN_Index[i+1];

}
counter = counter + Selected_Index.NN_Index[0];

}




Selected_Index.NN_Index[0]=0;
if( gridPos.x!=0 && gridPos.z!=(SimParms.num_NNCells.z-1) && gridPos.y!=(SimParms.num_NNCells.y-1) )
{
int3 neighbourPos = gridPos + make_int3(-1, 1, 1);

Selected_Index = CollisionDetection_Filter_SphereNNCell( neighbourPos, index,cellStart, cellEnd, position_com, P_ObjectID);


for(int i=0; i <Selected_Index.NN_Index[0];i++)
{
   Broad_List[index*32 + (counter+i)] = Selected_Index.NN_Index[i+1];
}
counter = counter + Selected_Index.NN_Index[0];

}


Selected_Index.NN_Index[0]=0;
if( gridPos.x!=(SimParms.num_NNCells.x-1) && gridPos.y!=(SimParms.num_NNCells.y-1) && gridPos.z!=0 )
{
int3 neighbourPos = gridPos + make_int3(1, 1, -1);

Selected_Index = CollisionDetection_Filter_SphereNNCell( neighbourPos, index,cellStart, cellEnd, position_com, P_ObjectID);


for(int i=0; i <Selected_Index.NN_Index[0];i++)
{
   Broad_List[index*32 + (counter+i)] = Selected_Index.NN_Index[i+1];

}
counter = counter + Selected_Index.NN_Index[0];

}



Selected_Index.NN_Index[0]=0;
if( gridPos.x!=0 && gridPos.y!=(SimParms.num_NNCells.y-1) && gridPos.z!=0 )
{
int3 neighbourPos = gridPos + make_int3(-1, 1, -1);

Selected_Index = CollisionDetection_Filter_SphereNNCell( neighbourPos, index,cellStart, cellEnd, position_com, P_ObjectID );


 for(int i=0; i < Selected_Index.NN_Index[0];i++)
{
   Broad_List[index*32 + (counter+i)] = Selected_Index.NN_Index[i+1];
}
counter = counter + Selected_Index.NN_Index[0];

}


Selected_Index.NN_Index[0]=0;
if( gridPos.x!=(SimParms.num_NNCells.x-1) && gridPos.y!=0 )
{
int3 neighbourPos = gridPos + make_int3(1, -1, 0);

Selected_Index = CollisionDetection_Filter_SphereNNCell( neighbourPos, index,cellStart, cellEnd, position_com, P_ObjectID );


 for(int i=0; i <Selected_Index.NN_Index[0];i++)
{
   Broad_List[index*32 + (counter+i)] = Selected_Index.NN_Index[i+1];

}
counter = counter + Selected_Index.NN_Index[0];

}



Selected_Index.NN_Index[0]=0;
if( gridPos.y!=0 && gridPos.z!=(SimParms.num_NNCells.z-1)  )
{
int3 neighbourPos = gridPos + make_int3(0, -1, 1);

Selected_Index = CollisionDetection_Filter_SphereNNCell( neighbourPos, index,cellStart, cellEnd, position_com , P_ObjectID);


 for(int i=0; i < Selected_Index.NN_Index[0];i++)
{
   Broad_List[index*32 + (counter+i)] = Selected_Index.NN_Index[i+1];

}
counter = counter + Selected_Index.NN_Index[0];

}


/*  Own cell Logic */
uint gridHash = SpatialDecomp_CalcGridHash(gridPos);

/* get start of bucket for this cell */
uint startIndex = cellStart[gridHash];
/* Get end of the cell */
uint endIndex = cellEnd[gridHash];


int num_particlesCell = endIndex-startIndex;

/* 2 Particles in cell */
if (num_particlesCell>1 && num_particlesCell<=2)
{
Selected_Index.NN_Index[0]=0;
Selected_Index = CollisionDetection_Filter_2SphereSameCell( gridPos, index,cellStart, cellEnd, position_com, P_ObjectID, P_ID);


for(int i=0; i <Selected_Index.NN_Index[0];i++)
{
	Broad_List[index*32 + (counter+i) ] = Selected_Index.NN_Index[i+1];//make_uint2(index,Selected_Index[i+1]);

}
counter = counter + Selected_Index.NN_Index[0];
}/* get a list of pairs*/
else if (num_particlesCell>2)
{
Selected_Index.NN_Index[0]=0;
        Selected_Index = CollisionDetection_Filter_SphereSameCellMany( gridPos, index,cellStart, cellEnd, position_com, P_ObjectID, P_ID);


        for(int i=0; i <Selected_Index.NN_Index[0];i++)
        {
        	Broad_List[index*32 + (counter+i) ] = Selected_Index.NN_Index[i+1];//make_uint2(index,Selected_Index[i+1]);

        }
        counter = counter + Selected_Index.NN_Index[0];

}



/* Store total number of neighbours */
NN_NUM[index] = counter;
//printf("stored %d  %d \n",(index),counter);


    //printf("Pos %f %f %f  Grid Pos %d %d %d  : Particle %d has %d NN \n ", position_com[index].x, position_com[index].y, position_com[index].z,
  		//  gridPos.x,gridPos.y,gridPos.z,  P_ID[index],NN_NUM[index]);
} /* thread check */

/*multi gpu */
/* get the unique index */
//uint index = blockIdx.x*blockDim.x + threadIdx.x + gpu_offset;

}
/*-----------------------------------------------------------------------------*/


__global__ void PrintNN (uint *NumNN, uint *NNList, int *P_ID)
{
//  for ( int i=0; i<SimParms.Num_Particles; i++ )
//  {
//	  //printf( "particle %d has %d NN \n",P_ID[i], NumNN[i]);
//
//	  for ( int j=0; j<NumNN[i]; j++ )
//	  {
//		printf(" %d %d \n",i,NNList[i*32+j]);
//	  }
//  }

printf("\n");

}



