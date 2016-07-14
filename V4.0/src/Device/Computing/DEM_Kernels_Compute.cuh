#include "Collision_Detection_General.cuh"
#include "Collision_Detection_Spheres.cuh"
#include "Collision_Detection_Polyhedra.cuh"

#include "Integration.cuh"



/*---------------------------------------------------------------------------*/
/* (1) Performs collision detection for objects that have there positions
                         Modified in constant memory
        Contact History is used for each volume object so that we dont
        jump contacts an the edge: It is assumed that each particle
        can be in contact with only one volume object at a time  */
/*---------------------------------------------------------------------------*/
__global__ void VolumeObject_InteractionPolyhedra(
		                                    float3       *dLifterForces,
		                                    float3       *dLifterForces_ang,
		                                    Contact_Info *VObject_Contact_Hist,
		                                    float3       *position_com,
		                                    Quaterion    *position_ornt,
		                                    float3       *velocity_com,
		                                    float3       *velocity_ang,
		                                    uint         *P_ObjectType,
		                                    int          *P_ID ,
		                                    int           Num_Particles,
											int           gpu_offset=0         )
{

	uint index = blockIdx.x*blockDim.x  + threadIdx.x ;


    if( (index < Num_Particles) && P_ID[index + gpu_offset]>-1)
	{
      uint d_index=0;

      dLifterForces     [index] = make_float3(0.0f);
      dLifterForces_ang [index] = make_float3(0.0f);

	  CollData  Result_D;

	  Result_D.is_collision  =  false;

	  uint   P_type = P_ObjectType[index+ gpu_offset];
	  float3 Pos_P  = position_com[index+ gpu_offset];
	  float3 Vel_P  = velocity_com[index+ gpu_offset];


     bool has_rotated=false;

     Edge Edges_P[32];

     Contact_Info  local;

      /* For a ball only check if the particle is close to the mill shell */
      if( ( VolumeObject[d_index].is_translating || !SimParms.isMillSim || SimParms.Rotating_WObject) ||
    	 /* ball mill case */
       ( length(make_float3(Pos_P.x*WorldObject[0].cyl_geo.AXIS.x,Pos_P.y*WorldObject[0].cyl_geo.AXIS.y,Pos_P.z*WorldObject[0].cyl_geo.AXIS.z) - SimParms.cent) > (SimParms.VolumeObject_CylinderBoundR-2.0f*ParticleObject[P_type].radius ))   )
      {


    	  POBJ_ROT POBJ_RT;

	     /* Loop over all Dynamic Objects */
	     for ( d_index=0; d_index<Num_DynamicObjects; d_index++ )
	     {
		    float3 Pos_D = VolumeObject[d_index].COM;
		    float3 BV    = (Pos_P - Pos_D)*VolumeObject[d_index].Axis;

		    /* Check against bounding cylinder of a VObject */
		    if( length(BV) < ( VolumeObject[d_index].boundR +
		    		ParticleObject[P_type].radius) )
		    {

 				/* Update Lifter Edges  */
 			    Edge EdgesD[20];
 			    for( uint i=0; i<VolumeObject[ d_index ].num_edges; i++ )
 			    {

 			      EdgesD[i].Point = VolumeObject[d_index].vertex
 			    		 [VolumeObject[d_index].edge[i].point_vindex[0]]
 			                                                      ;

 			      EdgesD[i].dir = (VolumeObject[ d_index].vertex
 			    	  [VolumeObject[ d_index ].edge[i].point_vindex[1]]
 			                                  ) - EdgesD[i].Point;

 			     }



                /* Rotate particle only if bounds overlap*/
				 if(!has_rotated)
				 {
				    /* 1. Update Current Orientation of Particle in global cords */
				   POBJ_RT = RotatePolyhedra(P_type, Pos_P,SimParms.Rotation,position_ornt[index + gpu_offset] );
				   has_rotated =true;

				   /* 1.1 Make Particle Edges */

					for( uint i=0; i<ParticleObject[ P_type ].num_edges; i++ )
					{

					 Edges_P[i].Point = POBJ_RT.Vertex[ParticleObject[P_type].edge[i].point_vindex[0]]
														 ;
					 Edges_P[i].dir   = (POBJ_RT.Vertex[ParticleObject[P_type].edge[i].point_vindex[1]]
													) - Edges_P[i].Point;
					}


				 }



					 /* 1.Quick check for Face SP  */

					 /* Have to use a different method call since Vobject faces are surfaces not planes */
					 bool SP_Check = is_SP_FacesVO_A_VertexB( VolumeObject[d_index].num_faces,
													     VolumeObject[d_index].faces,
													     Pos_P,
													     ParticleObject[P_type].num_vertex,
													     POBJ_RT.Vertex                      ) ;

					 if(!SP_Check)
					 {
						 /* Check for SP between Faces of Particle and Volume object */
						 SP_Check = is_SP_FacesA_VertexB(  ParticleObject[P_type].num_faces,
													  POBJ_RT.face,
													  Pos_D,
													  VolumeObject[d_index].num_vertex,
													  VolumeObject[d_index].vertex               ) ;
					 }



					 /* So either there contact or an Edge SP */
					 if( !SP_Check )
					 {
                            /* Can only be in contact with one face of a VObject *change for concave */
						    for ( uint i=0; i < VolumeObject[d_index].num_faces; i++ )
				            {
										Result_D = Collision_Detection_Polyhedra_Surface
												 ( POBJ_RT.Vertex, P_type, Pos_P, Vel_P,
												   VolumeObject[d_index].faces[i],
												   VolumeObject[d_index].vertex                );


										if(Result_D.is_collision)
										{
											local.snum = i;
											local.cont_type=1;
											break;
										}

						     }

					      /* Check Edge */
						  if ( !Result_D.is_collision )
						  {
								  Pos_D.z = Pos_P.z;

								  float3 RelPos = (Pos_P - Pos_D)/length( Pos_P-Pos_D);

								  Result_D = Collision_Detection_Polyhedra_Polyhedra_EdgesN( EdgesD,   Edges_P,
										  VolumeObject[d_index].num_edges,   ParticleObject[P_type].num_edges,
								          RelPos,  Pos_D, Pos_P,Vel_P  );
								  if(Result_D.is_collision)
								  {
									  local.cont_type=2;
								  }

						  }



					}



			   if(Result_D.is_collision)
			   {

				   float3 RelPos = Pos_P - Pos_D*VolumeObject[d_index].Axis;


					/* Indicates a contact error */
					  if( Result_D.dis > 0.1500f*ParticleObject[P_type].radius )
					  {
						  //printf("%f !!Dwarning large dis  Coll type %d Coll claf %d : P %d %d  vel %f  \n ",
						//									  					  Result_D.dis,Result_D.colltype, Result_D.collision_class, P_ID[index + gpu_offset], P_ID[index], length((velocity_com[index + gpu_offset])/100.0f));

						   //Result_D.Normal = RelPos/length(RelPos);
						   Result_D.dis = 0.025f*ParticleObject[P_type].radius;
						   Result_D.contact_point = Pos_P;
					  }

					  /* Limits velocity in case of large penetrations */
					  else if( Result_D.dis > 0.0500f*ParticleObject[P_type].radius)
					  {

						 Result_D.dis = 0.020f*ParticleObject[P_type].radius;
                         Result_D.contact_point = Pos_P;
					  }



				   Forces Force;
				   float3 Vel_Lifter = make_float3(0.0f);

				   if(VolumeObject[d_index].is_attached_DWorld)
				   {
				     Vel_Lifter = cross( make_float3(0.0,0.0,-SimParms.MillRotVel),
   		                                        Result_D.contact_point -
   		                                        make_float3(WorldObject[0].cyl_geo.radius,
   		            		                    WorldObject[0].cyl_geo.radius,0.0f)       );
				   }
				   else if(VolumeObject[d_index].is_translating)
				   {

					   Vel_Lifter = VolumeObject[d_index].translate_vel ;
				   }


					  Force = Collision_Response_Particle_VolumeObject( Result_D.Normal,
                                                                 Result_D.dis,
                                                                (Result_D.contact_point-Pos_P),
                                                                 make_float3(0.0f),
                                                                 velocity_ang[index+ gpu_offset],
                                                                 make_float3(0.0f),
		  					                                     P_type,
		  					                                     Vel_P,
		  						                                 Vel_Lifter,
		  					                                     ParticleObject[P_type].mass    );


				  dLifterForces[index] += Force.trans;

				  /* Rotational Velocity */
				  if( SimParms.Rotation )
				  {
					  dLifterForces_ang[index] += Force.torque;

				  }

					 /* Assume only one lifter at a time */
					 break; /* *EARLY EXIT MIGHT not Help */
			   }/* End If collision */

		 } /* End Bound Collision */

	   }/* End Loop over Dynamic Objects */

     }/* End loop over lifters bound */

     /* We use contact history so we dont jump from edge to face contact at
      * the corner */
      local.obj_id    = d_index;
      VObject_Contact_Hist[index+ gpu_offset]  = local;

   }


}

/*---------------------------------------------------------------------------*/
/* (1) Performs collision detection for objects that have there positions
                         Modified in constant memory                         */
/*---------------------------------------------------------------------------*/
__global__ void VolumeObject_InteractionSpheres(
	                                        bool      get_surfacePoints,
                                            float3    *dWallContactPoints,
		                                    float3       *dLifterForces,
		                                    float3       *dLifterForces_ang,
		                                    Contact_Info *VObject_Contact_Hist,
		                                    float3       *position_com,
		                                    float3       *velocity_com,
		                                    float3       *velocity_ang,
		                                    uint         *P_ObjectType,
		                                    int          *P_ID ,
		                                    int           Num_Particles,
											int           gpu_offset=0)
{

	uint index = blockIdx.x*blockDim.x  + threadIdx.x;


    if( (index < Num_Particles) && P_ID[index + gpu_offset]>-1)
	{
      uint d_index=0;

      dLifterForces     [index] = make_float3(0.0f);
      dLifterForces_ang [index] = make_float3(0.0f);

	  CollData  Result_D;

	  Result_D.is_collision  =  false;

	  uint   P_type = P_ObjectType [index+ gpu_offset];
	  float3 Pos_P  = position_com [index+ gpu_offset];
	  float3 Vel_P  = velocity_com [index+ gpu_offset];

	  if(get_surfacePoints)
      {
		dWallContactPoints[index] = make_float3(0.0f);
	  }

	  			   
	  Contact_Info  local;
	  local = VObject_Contact_Hist[index+ gpu_offset];



      /* Check only outer particles for mill */
      /* For a ball only check if the particle is close to the mill shell */
      if( ( VolumeObject[d_index].is_translating || !SimParms.isMillSim || SimParms.Rotating_WObject) ||
    	 /* ball mill case */
       ( length(make_float3(Pos_P.x*WorldObject[0].cyl_geo.AXIS.x,Pos_P.y*WorldObject[0].cyl_geo.AXIS.y,Pos_P.z*WorldObject[0].cyl_geo.AXIS.z) - SimParms.cent) > (SimParms.VolumeObject_CylinderBoundR-2.0f*ParticleObject[P_type].radius ))   )
      {
	     /* Loop over all Dynamic Objects */
	     for ( d_index=0; d_index<Num_DynamicObjects; d_index++ )
	     {
		    float3 Pos_D = VolumeObject[d_index].COM;
		    float3 BV    = (Pos_P - Pos_D)*VolumeObject[d_index].Axis;

		    /* Check against bounding cylinder of a DObject */
		    if( length(BV) < ( VolumeObject[d_index].boundR +
		    		ParticleObject[P_type].radius) )
		    {

				/* Make local copy Lifter Edges */
			    Edge EdgesD[16];
			    for( uint i=0; i<VolumeObject[ d_index ].num_edges; i++ )
			    {

			      EdgesD[i].Point = VolumeObject[d_index].vertex
			    		 [VolumeObject[d_index].edge[i].point_vindex[0]]
			                                                      ;

			      EdgesD[i].dir = (VolumeObject[ d_index].vertex
			    	  [VolumeObject[ d_index ].edge[i].point_vindex[1]]
			                                  ) - EdgesD[i].Point;

			     }

			     /* Sphere particle */
                 if (SimParms.particle_type==0)
                 {

			        Result_D = Collision_Detection_Sphere_DObject
							( index, d_index ,Pos_P,
							  Vel_P, EdgesD,
							  VolumeObject[d_index].num_edges,ParticleObject[P_type].radius,local);

                 }


			   if(Result_D.is_collision)
			   {
                   /* No rotation for sphere edge contact */
				   if(SimParms.particle_type==0 && Result_D.collision_class==Edge_Face)
				   {
					   Result_D.contact_point = Pos_P;
				   }

				    if(SimParms.unit_test==1)
				    {
					    printf("PV %d pen dis %f ",P_ID[index],Result_D.dis);
				        PrintVectorD(Vel_P);
				    }

				   Forces Force;
				   float3 Vel_Lifter = make_float3(0.0f);

				   if(VolumeObject[d_index].is_attached_DWorld)
				   {
				     Vel_Lifter = cross( make_float3(0.0,0.0,-SimParms.MillRotVel),
   		                                        Result_D.contact_point -
   		                                        make_float3(WorldObject[0].cyl_geo.radius,
   		            		                    WorldObject[0].cyl_geo.radius,0.0f)       );
				   }
				   else if(VolumeObject[d_index].is_translating)
				   {

					   Vel_Lifter = VolumeObject[d_index].translate_vel ;
				   }



					  Force = Collision_Response_Particle_VolumeObject( Result_D.Normal,
                                                                 Result_D.dis,
                                                                (Result_D.contact_point-Pos_P),
                                                                 make_float3(0.0f),
                                                                 velocity_ang[index],
                                                                 make_float3(0.0f),
		  					                                     P_type,
		  					                                     Vel_P,
		  						                                 Vel_Lifter,
		  					                                     ParticleObject[P_type].mass    );


				  dLifterForces[index] += Force.trans;

				  /* Rotational Velocity */
				  if( SimParms.Rotation )
				  {
					  dLifterForces_ang[index] += Force.torque;

				  }

				   if(get_surfacePoints)
				   {
				     dWallContactPoints[index] = Result_D.contact_point;
				   }

					 /* Assume only one lifter at a time */
					 break;
			   }/* End If collision */

		 } /* End Bound Collision */

	   }/* End Loop over Dynamic Objects */

     }/* End loop over lifters bound */

     /* We use contact history so we dont jump from edge to face contact at
      * the corner */


      local.obj_id    = d_index;
      VObject_Contact_Hist[index]  = local;

   }


}



/*---------------------------------------------------------------------------*/
/*        (2)        World interaction for Polyhedral particles               */
/*---------------------------------------------------------------------------*/
__global__ void WorldObject_InteractionPolyhedra_Planar( bool      get_surfacePoints,
		                                     bool      is_cylinder_rotation,
		                                     float3    *dWall_Forces_com,
		                                     float3    *dWall_Forces_ang,
		                                     float3    *dWallContactPoints,
		                                     float3    *position_com,
		                                     Quaterion *position_ornt,
		                                     float3    *velocity_com,
		                                     float3    *velocity_ang,
                                             uint      *P_ObjectType,
                                             int       *P_ID,
                                             int        Num_Particles,
											int           gpu_offset=0     )
{

	uint index = blockIdx.x*blockDim.x  + threadIdx.x;

   if( (index < Num_Particles) && P_ID[index + gpu_offset]>-1)
   {


	 CollData Result;

	 uint     P_type = P_ObjectType[index+ gpu_offset];

	 float3   Pos_P  = position_com[index+ gpu_offset];
	 float3   Vel_P  = velocity_com[index+ gpu_offset];
	 POBJ_ROT POBJ_RT;

	 Result.is_collision    = false;
	 Result.collision_class = UNDEF;


	 if(get_surfacePoints)
	 {
		 dWallContactPoints[index] = make_float3(0.0f);
	 }


	      /* 1. Update Current Orientation of Particle in global cords */
	     POBJ_RT = RotatePolyhedra( P_type, Pos_P,true,position_ornt[index+ gpu_offset] );



     /* Loop over all World Objects */
	 for ( uint k=0; k<Num_WorldObjects; k++ )
	 {

		 /* If the particle has a very high velocity it can
		  * collied with a surface and its new velocity take
		  * it past another surface */

		if( WorldObject[k].surface_type==plane)
		{
		    /* Check for possible collisions 1 surface at a time */
		    for ( uint i=0; i < WorldObject[k].num_surfaces; i++ )
            {

		     /* Treat as an infinite plane */
		     if(WorldObject[k].surfaces[i].area==-1.0f)
		     {

		    		Result = Collision_Detection_Polyhedra_Plane
		    									 ( POBJ_RT.Vertex, P_type, Pos_P, Vel_P,
		    									   WorldObject[k].surfaces[i],
		    									   WorldObject[k].vertex                );

		     }
			 else
		     {

						Result = Collision_Detection_Polyhedra_Surface
								 ( POBJ_RT.Vertex, P_type, Pos_P, Vel_P,
								   WorldObject[k].surfaces[i],
								   WorldObject[k].vertex                );

		     }

			   if(Result.is_collision)
			   {

			     Forces Force;

				   if(SimParms.Rotating_WObject && is_cylinder_rotation)
				   {
				     float3 Vel_Lifter = cross( make_float3(0.0,0.0,-SimParms.MillRotVel),
 		                                        Result.contact_point -
 		                                        make_float3(WorldObject[0].cyl_geo.radius,
 		                                        		WorldObject[0].cyl_geo.radius,0.0f)       );


					  Force = Collision_Response_Particle_MovingSurface( Result.Normal,
                                                               Result.dis,
                                                              (Result.contact_point-Pos_P),
                                                               velocity_ang[index+ gpu_offset],
		  					                                     P_type,
		  					                                     Vel_P,
		  						                                 Vel_Lifter   );
				   }
				   else
				   {
				    Force = Collision_Response_Particle_StaticSurface(
											Vel_P,
											Result.Normal,
											Result.dis,
										   (Result.contact_point -Pos_P),
											velocity_ang[index+ gpu_offset],
											P_type  );
				   }

				   /* Update output array */
				   if(get_surfacePoints)
				   {
				     dWallContactPoints[index] = Result.contact_point;
				   }

					dWall_Forces_com[index] += Force.trans;

					/* Rotational Velocity */
					 if( SimParms.Rotation )
					 {
						 if(Result.collision_class!=Face_Face)
						 {
						   dWall_Forces_ang[index] += Force.torque;
						 }/* When we have face face particle should stop rotating */
                     }
			   }
            }
		}



	 }/* End Loop over World Objects */




   } /* thread check */

}
/*---------------------------------------------------------------------------*/





/*---------------------------------------------------------------------------*/
/*        (2)        World interaction for Polyhedral particles               */
/*---------------------------------------------------------------------------*/
__global__ void WorldObject_InteractionPolyhedra_Macro( bool      get_surfacePoints,
		                                     bool      is_cylinder_rotation,
		                                     float3    *dWall_Forces_com,
		                                     float3    *dWall_Forces_ang,
		                                     float3    *dWallContactPoints,
		                                     float3    *position_com,
		                                     Quaterion *position_ornt,
		                                     float3    *velocity_com,
		                                     float3    *velocity_ang,
                                             uint      *P_ObjectType,
                                             int       *P_ID,
                                             int        Num_Particles,
											int           gpu_offset=0     )
{

	uint index = blockIdx.x*blockDim.x  + threadIdx.x;

   if( (index < Num_Particles) && P_ID[index + gpu_offset]>-1)
   {


	 CollData Result;

	 uint     P_type = P_ObjectType[index+ gpu_offset];

	 float3   Pos_P  = position_com[index+ gpu_offset];
	 float3   Vel_P  = velocity_com[index+ gpu_offset];
	 POBJ_ROT POBJ_RT;

	 Result.is_collision    = false;
	 Result.collision_class = UNDEF;



	      /* 1. Update Current Orientation of Particle in global cords */
	     POBJ_RT = RotatePolyhedra( P_type, Pos_P,SimParms.Rotation,position_ornt[index+ gpu_offset] );



     /* Loop over all World Objects */
	 for ( uint k=0; k<Num_WorldObjects; k++ )
	 {

	      if ( WorldObject[k].surface_type == cylinder )
	      {
	         Macro_Cylinder Cyl =  WorldObject[k].cyl_geo;



      	     float3 Velnorm = velocity_com[index+ gpu_offset];
      	            Velnorm = Velnorm/length(Velnorm);

      	     /* Check particle Cylinder Surface */


        	   Result = Collision_Detection_Cylinder_PolyF
        	    		   ( index,P_type, Pos_P, Vel_P, Cyl,POBJ_RT,k);


			   if(Result.is_collision)
			   {
			     Forces Force;

                   /* Limit Max Pen Distance */
			       if(Result.dis>ParticleObject[P_type].radius*0.05f)
			       {
			          Result.dis=ParticleObject[P_type].radius*0.015f;
			       }


			      /* Add Movement of drum */

			      if (is_cylinder_rotation)
			      {


				     float3 Vel_Drum = cross( make_float3(0.0,0.0,-SimParms.MillRotVel),
		                                        Result.contact_point -
		                                        make_float3(WorldObject[0].cyl_geo.radius,
		            		                    WorldObject[0].cyl_geo.radius,0.0f)       );


					  Force = Collision_Response_Particle_MovingSurface( Result.Normal,
                                                              Result.dis,
                                                             (Result.contact_point-Pos_P),
                                                              velocity_ang[index + gpu_offset],
		  					                                     P_type,
		  					                                     Vel_P,
		  						                                 Vel_Drum   );
			      }
			      else
			      {

			      Force = Collision_Response_Particle_StaticSurface(
			    		                            Vel_P,
                                                    Result.Normal,
                                                    Result.dis,
                                                    (Result.contact_point -Pos_P),
                                                    velocity_ang[index + gpu_offset ],
 										            P_type  );

			      }


				   /* Update output array */
				   if(get_surfacePoints)
				   {
				     dWallContactPoints[index] = Result.contact_point;
				   }


				   dWall_Forces_com[index] += Force.trans;

					/* Rotational Velocity if its polyhedra with face flat on cylinder dont induce a moment*/
		           if( SimParms.Rotation )
				   {
		        	 if(!(Result.collision_class==Face_Face ))
		        	 {
	                   dWall_Forces_ang[index] += Force.torque;
		        	 }
			       }
			   }


      	     if(Cyl.has_top_cap)
      	     {


    			   Result = Collision_Detection_Vertcies_Plane( ParticleObject[P_type].num_vertex, POBJ_RT.Vertex,
    						                        Cyl.normal_top_cap,
    												Cyl.center_top_cap       );


 			   if(Result.is_collision)
 			   {
 				  //printf("top cap contact %d \n",P_ID[index]);
 			     Forces Force;

			      Force = Collision_Response_Particle_StaticSurface(
			    		                            Vel_P,
                                                   Result.Normal,
                                                   Result.dis,
                                                   (Result.contact_point -Pos_P),
                                                   velocity_ang[index + gpu_offset ],
										            P_type  );

				   /* Update output array */
				   if(get_surfacePoints)
				   {
				     dWallContactPoints[index] = Result.contact_point;
				   }

				   dWall_Forces_com[index]+= (Force.trans  );


 					/* Rotational Velocity */
 					 if( SimParms.Rotation )
 					 {
 						 dWall_Forces_ang[index] += Force.torque;
 					 }
 			   }
      	     }

      	     if(Cyl.has_bot_cap)
      	     {
      	    	//printf("bot cap contact %d \n",P_ID[index]);

    			   Result = Collision_Detection_Vertcies_Plane( ParticleObject[P_type].num_vertex, POBJ_RT.Vertex,
    						                        Cyl.normal_bot_cap,
    												Cyl.center_bot_cap       );


 			   if(Result.is_collision)
 			   {
 			     Forces Force;



			      Force = Collision_Response_Particle_StaticSurface(
			    		                            Vel_P,
                                                   Result.Normal,
                                                   Result.dis,
                                                   (Result.contact_point -Pos_P),
                                                   velocity_ang[index + gpu_offset ],
										            P_type  );




 					   /* Update output array */
 					   if(get_surfacePoints)
 					   {

 					     dWallContactPoints[index] =Result.contact_point;
 					   }

 					  dWall_Forces_com[index] += (Force.trans  );


 					/* Rotational Velocity */
 					 if( SimParms.Rotation )
 					 {
 						 dWall_Forces_ang[index] += Force.torque;

 					  }/* End if rotation */

 			   }/* End if cap collision */

      	     }   /* End bot cap */

	      }/* End cylinder check */



	 }/* End Loop over World Objects */


   } /* thread check */

}
/*---------------------------------------------------------------------------*/




/*---------------------------------------------------------------------------*/
/*  (2)  World interaction for Spherical particles per world object          */
/*---------------------------------------------------------------------------*/
__global__ void WorldObject_InteractionSpheres_Planar( bool      get_surfacePoints,
		                                     float3    *dWall_Forces,
		                                     float3    *dWall_Forces_ang,
		                                     float3    *dWallContactPoints,
		                                     float3    *position_com,
		                                     Quaterion *position_ornt,
		                                     float3    *velocity_com,
		                                     float3    *velocity_ang,
                                             uint      *P_ObjectType,
                                             int       *P_ID,
                                             int        Num_Particles,
											 int        gpu_offset=0     )
{
    bool is_cylinder_rotation=false;
	uint index = blockIdx.x*blockDim.x  + threadIdx.x;

   if( (index < Num_Particles) && P_ID[index+ gpu_offset]>-1)
   {

	  // P_ID[index+ gpu_offset]-= Num_Particles;

	 CollData_Sphere Result;

	 float3   Pos_P  = position_com[index+ gpu_offset];
	 float3   Vel_P  = velocity_com[index+ gpu_offset];
	 uint     P_type = P_ObjectType[index+ gpu_offset];

	 Result.is_collision    = false;

	 Forces Force;
	 Force.trans  = make_float3(0.0f);
	 Force.torque = make_float3(0.0f);

	 if(get_surfacePoints)
	 {
		 dWallContactPoints[index] = make_float3(0.0f);
	 }

     /* Loop over all World Objects */
	 for ( uint k=0; k<Num_WorldObjects; k++ )
	 {

		 /* If the particle has a very high velocity it can
		  * collied with a surface and its new velocity take
		  * it past another surface */

		//if( P_ID[index+ gpu_offset] > Num_Particles )
         if ( WorldObject[k].surface_type == plane )

         {
		    /* Check for possible collisions 1 surface at a time */
		    for ( uint i=0; i < WorldObject[k].num_surfaces; i++ )
            {

		     /* Treat as an infinite plane */
		     if(WorldObject[k].surfaces[i].area<0.0f)
		     {

				 if(SimParms.particle_type==0)
				 {

					 Result = Collision_Detection_Sphere_Plane
								(P_type, Pos_P,
								  WorldObject[k].surfaces[i].normal,
								  WorldObject[k].surfaces[i].centroid);

				 }
				 else
				 {
					 printf( "warning no collision method found \n");
                    /*TODO*/
				 }
		     }
			 else
		     {
				 if(SimParms.particle_type==0)
				 {

					 Result = Collision_Detection_Sphere_Surface
								( P_type, Pos_P,
								  WorldObject[k].surfaces[i],
								  WorldObject[k].vertex     );

				 }

		     }
			   if(Result.is_collision)
			   {

                   /* For a world surface that is rotating */
				   if(SimParms.Rotating_WObject && is_cylinder_rotation)
				   {

					  Force = Collision_Response_Particle_MovingSurface( Result.Normal,
                                                               Result.dis,
                                                              (Result.contact_point-Pos_P),
                                                               velocity_ang[index + gpu_offset],
		  					                                     P_type,
		  					                                     Vel_P,
		  					                                   cross( make_float3(0.0,0.0,-SimParms.MillRotVel),
		  					                                    		                                        Result.contact_point -
		  					                                    		                                        make_float3(WorldObject[0].cyl_geo.radius,
		  					                                    		                                        		WorldObject[0].cyl_geo.radius,0.0f)       )   );
				   }
				   else
				   {
				    Force = Collision_Response_Particle_StaticSurface(
											Vel_P,
											Result.Normal,
											Result.dis,
										   (Result.contact_point -Pos_P),
											velocity_ang[index + gpu_offset],
											P_type  );
				   }

				   /* Update output array */
				   if(get_surfacePoints)
				   {
				     dWallContactPoints[index] = Result.contact_point;
				   }

					dWall_Forces[index] += Force.trans;

					/* Rotational Velocity */
					 if( SimParms.Rotation )
					 {
						   dWall_Forces_ang[index] += Force.torque;
                     }
			   }
            }
		}


	 }/* End Loop over World Objects */


   } /* thread check */

}
/*---------------------------------------------------------------------------*/





/*---------------------------------------------------------------------------*/
/*        (2)        World interaction for Spherical particles               */
/*---------------------------------------------------------------------------*/
__global__ void WorldObject_InteractionSpheres_Macro( bool      get_surfacePoints,
		                                     bool      is_cylinder_rotation,
		                                     float3    *dWall_Forces,
		                                     float3    *dWall_Forces_ang,
		                                     float3    *dWallContactPoints,
		                                     float3    *position_com,
		                                     Quaterion *position_ornt,
		                                     float3    *velocity_com,
		                                     float3    *velocity_ang,
                                             uint      *P_ObjectType,
                                             int       *P_ID,
                                             int        Num_Particles,

											int           gpu_offset=0     )
{

	uint index = blockIdx.x*blockDim.x  + threadIdx.x;

   if( (index < Num_Particles) && P_ID[index+ gpu_offset]>-1)
   {


	 CollData_Sphere Result;

	 float3   Pos_P  = position_com[index+ gpu_offset];
	 float3   Vel_P  = velocity_com[index+ gpu_offset];
	 uint     P_type = P_ObjectType[index+ gpu_offset];


	 Result.is_collision    = false;




	 Forces Force;
	 Force.trans  = make_float3(0.0f);
	 Force.torque = make_float3(0.0f);



     /* Loop over all World Objects */
	 for ( uint k=0; k<Num_WorldObjects; k++ )
	 {

		 /* If the particle has a very high velocity it can
		  * collied with a surface and its new velocity take
		  * it past another surface */

          if ( WorldObject[k].surface_type == cylinder || WorldObject[k].surface_type == cone_trunc )
          {
	         Macro_Cylinder Cyl =  WorldObject[k].cyl_geo;


      	     /* Check particle Cylinder Surface */
			 if( SimParms.particle_type==0 && WorldObject[k].surface_type==cylinder )
      	     {
				// printf("check cylinder\n");
      	       Result = Collision_Detection_Sphere_Cylinder
      	    		   ( index + gpu_offset,P_type, Pos_P, Vel_P, Cyl,k);
   	         }
			 else if( SimParms.particle_type==0 && WorldObject[k].surface_type==cone_trunc )
			{
				Result = Collision_Detection_Sphere_ConeLL
      	    		   ( index + gpu_offset,P_type, Pos_P, Vel_P, Cyl);
			}



			   if(Result.is_collision)
			   {

			      /* Add Movement of drum */

			      if (is_cylinder_rotation)
			      {



					  Force = Collision_Response_Particle_MovingSurface( Result.Normal,
                                                              Result.dis,
                                                             (Result.contact_point-Pos_P),
                                                              velocity_ang[index + gpu_offset],
		  					                                     P_type,
		  					                                     Vel_P,
		  					                                   cross( make_float3(0.0,0.0,-SimParms.MillRotVel),
		  					                                   		                                        Result.contact_point -
		  					                                   		                                        make_float3(WorldObject[0].cyl_geo.radius,
		  					                                   		            		                    WorldObject[0].cyl_geo.radius,0.0f)       )   );
			      }
			      else
			      {

			      Force = Collision_Response_Particle_StaticSurface(
			    		                            Vel_P,
                                                    Result.Normal,
                                                    Result.dis,
                                                    (Result.contact_point -Pos_P),
                                                    velocity_ang[index + gpu_offset],
 										            P_type  );

			      }


				   /* Update output array */
				   if(get_surfacePoints)
				   {
				     dWallContactPoints[index] = Result.contact_point;
				   }


				   dWall_Forces[index] += Force.trans;

					/* Rotational Velocity if its polyhedra with face flat on cylinder dont induce a moment*/
		           if( SimParms.Rotation )
				   {
	                   dWall_Forces_ang[index] += Force.torque;

			       }
			   }


      	     if(Cyl.has_top_cap)
			 {

 			       Result = Collision_Detection_Sphere_Plane(P_type, Pos_P,
 			    		                            Cyl.normal_top_cap,
 			    		                            Cyl.center_top_cap  );



 			   if(Result.is_collision)
			   {
 			     Forces Force;
				 if(SimParms.Simulation_Type==ballmill)
				 {
 			       Force = Collision_Response_EndPlate_Particle(
 			    		                            Vel_P,
                                                     Result.Normal,
                                                     Result.dis,
                                                    (Result.contact_point -Pos_P),
                                                     velocity_ang[index + gpu_offset],
                                                     P_type                       );
				 }
				 else
				 {
					 Force = Collision_Response_Particle_StaticSurface(
 			    		        Vel_P,
                                    Result.Normal,
                                    Result.dis,
                                (Result.contact_point -Pos_P),
                                    velocity_ang[index + gpu_offset],
                                    P_type                       );
				 }
				   /* Update output array */
				   if(get_surfacePoints)
				   {
				     dWallContactPoints[index] = Result.contact_point;
				   }

				   dWall_Forces[index] += Force.trans ;


 					/* Rotational Velocity */
 					 if( SimParms.Rotation )
 					 {
 						 dWall_Forces_ang[index] += Force.torque;
 					 }
 			   }
      	     }

      	     if(Cyl.has_bot_cap)
      	     {

 			       Result = Collision_Detection_Sphere_Plane(P_type, Pos_P,
 			    		                            Cyl.normal_bot_cap,
 			    		                            Cyl.center_bot_cap  );


 			   if(Result.is_collision)
 			   {
 			     Forces Force;

				 if(SimParms.Simulation_Type==ballmill)
				 {
 			       Force = Collision_Response_EndPlate_Particle(
 			    		                            Vel_P,
                                                     Result.Normal,
                                                     Result.dis,
                                                    (Result.contact_point -Pos_P),
                                                     velocity_ang[index + gpu_offset],
                                                     P_type                       );
				 }
				 else
				 {
					 Force = Collision_Response_Particle_StaticSurface(
 			    		        Vel_P,
                                    Result.Normal,
                                    Result.dis,
                                (Result.contact_point -Pos_P),
                                    velocity_ang[index + gpu_offset],
                                    P_type                       );
				 }

 					   /* Update output array */
 					   if(get_surfacePoints)
 					   {

 					     dWallContactPoints[index] = Result.contact_point;
 					   }

 					  dWall_Forces[index]+= Force.trans;


 					/* Rotational Velocity */
 					 if( SimParms.Rotation )
 					 {
 						 dWall_Forces_ang[index] += Force.torque;

 					  }/* End if rotation */

 			   }/* End if cap collision */

      	     }   /* End bot cap */

	      }/* End cylinder check */



	 }/* End Loop over World Objects */


   } /* thread check */

}
/*---------------------------------------------------------------------------*/








/*---------------------------------------------------------------------------*/
/*            (3) Particle Particle Collision Response Spheres:
 *                     each particle has its own thread  */
/*---------------------------------------------------------------------------*/
__global__ void Kernel_ParticleInteraction_Spheres_NonSymmetry
							 ( 
							 
	                           bool      get_surfacePoints,
                               float3    *dWallContactPoints,
							   uint      *NumNN,
							   uint      *Broad_List,
							   float3    *position_com,
							   float3    *velocity_com,
							   float3    *velocity_ang,
							   uint      *P_ObjectType,
							   int       *P_ID  ,
							   float3    *force_com,
							   float3    *force_ang,
							   int       Num_Particles,
						       uint gpu_offset=0)
{
	uint mem_index_A = blockIdx.x*blockDim.x + threadIdx.x;

	/* Check that the mem location is valid and the particle is alive */
	if( (mem_index_A < Num_Particles) && P_ID[mem_index_A + gpu_offset]>-1 )
	{

	//printf("Contact ID %d NumNN %d NN ID %d \n",P_ObjectType[mem_index_A],NumNN[mem_index_A],P_ObjectType[Broad_List[mem_index_A*32]] );

	  /*-----------------------------------------------------------------*/
	                      /* Local data for kernel */
	  CollData Result;
	  float3   relPos;
	  float3   RVel;
	  Forces   Force_PairA;
	 /*-----------------------------------------------------------------*/


	  /* Load information for particle in location A */
	  uint   P_typeA  = P_ObjectType [mem_index_A + gpu_offset];
	  float3 Vel_A    = velocity_com [mem_index_A + gpu_offset];
	  float3 Pos_A    = position_com [mem_index_A + gpu_offset];
	  float3 VelAng_A = velocity_ang [mem_index_A + gpu_offset];

	  /* Translational Force on Particle A */
	  float3 force  = make_float3(0.0f);
	  /* Rotational Force on particle A */
	  float3 forceL = make_float3(0.0f);

	 if(get_surfacePoints)
	 {
		 dWallContactPoints[mem_index_A]=make_float3(0.0f);
	 }


	 float3 avg_point = make_float3(0.0f);

      /* Loop over all the NN memlocations  of particle A */
	  for ( uint j=0; j< NumNN[mem_index_A]; j++)
	  {
		/* load the memory location of particle who is in a NN */
		int    mem_index_B  = Broad_List[(mem_index_A)*32 + j];

		/* Load information for particle in location B */
		uint   P_typeB  = P_ObjectType [mem_index_B];
		float3 Pos_B    = position_com [mem_index_B];
        float3 Vel_B    = velocity_com [mem_index_B];
		float3 VelAng_B = velocity_ang [mem_index_B];

		/* Compute local kernel info */
		relPos     = Pos_A - Pos_B  ;
		RVel       = Vel_A - Vel_B ;

		/* Compute Contact information */
		Result.Normal        = relPos/length(relPos); /* Move A away from B */

		Result.dis           = fabs( (length(relPos) -
				(ParticleObject[P_typeA].radius + ParticleObject[P_typeB].radius ) ) ) + 0.000010f;


		/* Compute force exerted on Particle A */
		Force_PairA = Collision_Response_Particle_Particle(   RVel,
		          		                                      Result.Normal,
		          		                                      Result.dis,
		          		                                      (-1.0f*Result.Normal*ParticleObject[P_typeA].radius),
		          		                                      (Result.Normal*ParticleObject[P_typeB].radius),
		          		  					                  VelAng_A,
		          		  					                  VelAng_B,
		          		  			  					      P_typeA,
		          		  			  				          P_typeB,
		          		  			  					      Vel_A,
		          		  			  						  Vel_B,
		          		  			  					      ParticleObject[P_typeA].mass,
		          		  			  					      ParticleObject[P_typeB].mass,
		          		  			  					      relPos);

	    if(get_surfacePoints)
	    {
			avg_point += Pos_A + (Result.Normal*ParticleObject[P_typeA].radius);
		}


	    if(SimParms.unit_test==1)
	    {
		  if(P_ID[mem_index_A + gpu_offset]==0)
		  {
	        printf("PP pen dis %f ",Result.dis);
	        PrintVectorND(Vel_A);
	        PrintVectorD(Vel_B);
		  }
	    }

         /* Add to force A normal + friction */
         force  += Force_PairA.trans;
         forceL += Force_PairA.torque;

         	      //printf("pen dis %f,  Object_ID %d NumNN %d Object_ID %d Force %f %f %f ang %f %f %f Pos %f %f %f Vel %f %f %f  Mass %f \n",Result.dis, P_ObjectType[mem_index_A],NumNN[mem_index_A],P_ObjectType[Broad_List[mem_index_A*32]],
         	     	//    		  force.x,force.y,force.z,forceL.x,forceL.y,forceL.z, position_com[mem_index_A].x, position_com[mem_index_A].y, position_com[mem_index_A].z,
         	     	//    		  velocity_com[mem_index_A].x,velocity_com[mem_index_A].y,velocity_com[mem_index_A].z, ParticleObject[P_ObjectType[mem_index_A]].mass);
	  }

        /* Update net PP force on A + Gravity */
	      force_com[mem_index_A] = force;

	      /* Update net PP ang force on A */
	      if( SimParms.Rotation )
	      {
	        force_ang[mem_index_A] = forceL;
	      }
	      else
	      {
	    	force_ang[mem_index_A] = make_float3(0.0f);
	      }

         
	 if(get_surfacePoints)
	 {
		 dWallContactPoints[mem_index_A]=avg_point;
	 }

	}/* End */


}
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
/*  (4) Polyhedra Polyhedra Collision Response using point in poly test */
/*---------------------------------------------------------------------------*/
__global__ void Kernel_ParticleInteraction_Polyhedra_NonSymmetry
                                             ( uint      *NumNN,
		                                       uint      *Broad_List,
		                                       float3    *position_com,
		                                       Quaterion *position_ornt,
		                                       float3    *velocity_com,
		                                       float3    *velocity_ang,
		                                       uint      *P_ObjectType,
		                                       int       *P_ID  ,
                                               float3    *force_com,
                                               float3    *force_ang,
                                               int        Num_Particles,
						                       uint gpu_offset=0         )
{
	uint mem_index_A = blockIdx.x*blockDim.x + threadIdx.x ;



	if( (mem_index_A < Num_Particles) && P_ID[mem_index_A + gpu_offset]>-1)
	{


		  Forces Force_PairA;

	      uint  P_typeA   = P_ObjectType [mem_index_A + gpu_offset];
	      float3 Vel_A    = velocity_com [mem_index_A + gpu_offset];
	      float3 Pos_A    = position_com [mem_index_A + gpu_offset];
	      float3 VelAng_A = velocity_ang [mem_index_A + gpu_offset];

		  /* Translational Force on Particle A */
			  float3 force  = make_float3(0.0f);
		   /* Rotational Force */
			  float3 forceL = make_float3(0.0f);

			  bool heru=false;

		  /* Rotate Particle A */
		  POBJ_ROT POBJA;
		  Edge     EdgesA[32];

		  if ( NumNN[mem_index_A] > 0 )
		  {
			 POBJA = RotatePolyhedra( P_typeA, Pos_A,SimParms.Rotation ,position_ornt[mem_index_A + gpu_offset]);

			 /* Make Edges A */
			 for( uint i=0; i < ParticleObject[ P_typeA ].num_edges; i++ )
			 {
					EdgesA[i].Point = POBJA.Vertex[ParticleObject[P_typeA].
													  edge[i].point_vindex[0] ];

					EdgesA[i].dir = (POBJA.Vertex[ParticleObject[P_typeA].
									edge[i].point_vindex[1]] ) - EdgesA[i].Point;


			  }
		   }

		   CollData Result_P;

		   float3 RVel;
		   uint   num_coll = 0;

		   /* Info for NN particle B */
		  int    mem_index_B;
		  uint   P_typeB;

		  float3 Pos_B;
		  float3 Vel_B;
		  float3 VelAng_B;
		  POBJ_ROT POBJB;

		  /* Loop over all the NN of particle A */
		  for ( uint j=0; j< NumNN[mem_index_A]; j++)
		  {

			  Result_P.is_collision = false;

			  mem_index_B  = Broad_List[(mem_index_A)*32 + j];
			  P_typeB      = P_ObjectType[mem_index_B];
			  Vel_B        = velocity_com[mem_index_B];
			  VelAng_B     = velocity_ang [mem_index_B];
			  Pos_B        = position_com[mem_index_B];
			  RVel         = Vel_A - Vel_B ;

			  /* Rotate Particle B */
		      POBJB = RotatePolyhedra( P_typeB, Pos_B,SimParms.Rotation,position_ornt[mem_index_B]);

			  Edge EdgesB[32];

			  for( uint i=0; i<ParticleObject[ P_typeB ].num_edges; i++ )
			  {
					EdgesB[i].Point = POBJB.Vertex[ParticleObject[P_typeB].edge[i].
															   point_vindex[0] ];

					EdgesB[i].dir = (POBJB.Vertex[ParticleObject[P_typeB].edge[i].
											 point_vindex[1]]) - EdgesB[i].Point;
			  }

			  float3 RelPos = (Pos_A-Pos_B)/length( Pos_A-Pos_B);


			   /* Check for vertex face SP  AB */
			   bool SP_Check = is_SP_FacesA_VertexB( ParticleObject[P_typeA].num_faces,
												  POBJA.face, Pos_B,
												  ParticleObject[P_typeB].num_vertex,
												  POBJB.Vertex                        );

				 /* Check for vertex face SP  BA */
				 if(!SP_Check)
				 {
					 SP_Check = is_SP_FacesA_VertexB( ParticleObject[P_typeB].num_faces,
												 POBJB.face, Pos_A,
												 ParticleObject[P_typeA].num_vertex,
												 POBJA.Vertex                        );
				 }

				 /* Check for vertex face contact or Edge SP */
				 if(!SP_Check)
				 {
					  /* Check Edges for contact or SP */

					  {

						  Result_P = Collision_Detection_Polyhedra_Polyhedra_EdgesN(EdgesA, EdgesB,
							     			 ParticleObject[P_typeA].num_edges,
							     			 ParticleObject[P_typeB].num_edges,
							                  RelPos,Pos_A,Pos_B,Vel_A);

					      if(Result_P.is_collision)
						  {
	  						 Result_P.collision_class = Edge_Edge;
							 Result_P.is_collision=true;
						   }
						   else
						   {
	                          Result_P.is_collision    = false;
							  Result_P.collision_class = SP_PLANE;
						   }
					  }

					  if (!Result_P.is_collision)
					  {
                         /* 1. Check Vertex A Face B */
					     Result_P = Collision_Detection_Polyhedara_Polyhedra_VertexFaces
							      ( Pos_A, ParticleObject[P_typeA].num_vertex,
							        POBJA.Vertex,
							        Pos_B ,ParticleObject[P_typeB].num_faces,
							        POBJB.face, ParticleObject[P_typeB].num_vertex,
							        POBJB.Vertex, ParticleObject[P_typeB],
							        P_typeB, (mem_index_A + gpu_offset),mem_index_B);

					  }
					  /* 2. Check Vertex B Face A */
					  if (!Result_P.is_collision)
					  {

						  Result_P = Collision_Detection_Polyhedara_Polyhedra_VertexFaces
								  ( Pos_B, ParticleObject[P_typeB].num_vertex,
								   POBJB.Vertex,
                                   Pos_A ,ParticleObject[P_typeA].num_faces,
                                   POBJA.face, ParticleObject[P_typeA].num_vertex,
                                   POBJA.Vertex, ParticleObject[P_typeA],
                                   P_typeA,  mem_index_B, (mem_index_A + gpu_offset)  );
						  /* B-A so reverse normal */
						  Result_P.Normal *= -1.0f;
					  }
                   }





				  if(Result_P.is_collision)
				  {
                      /* Since we sometimes make an error this corrects very large penetrations */
					  heru = false;
					  /* Indicates a contact error */
					  if( Result_P.dis > 0.1500f*ParticleObject[P_typeA].radius )
					  {

						   Result_P.Normal = RelPos/length(RelPos);
						   Result_P.dis = 0.025f*ParticleObject[P_typeA].radius;
						   Result_P.contact_point = Pos_A;
						   heru = true;

					  }
					  else if( Result_P.dis > 0.0500f*ParticleObject[P_typeA].radius)
					  {
						   Result_P.Normal = RelPos/length(RelPos);

						   Result_P.dis    = 0.020f*ParticleObject[P_typeA].radius;

						   Result_P.contact_point = Pos_A;

						   heru = true;
					  }

					  if(!heru)
					  {

				      Force_PairA = Collision_Response_Particle_Particle(     RVel,
						          		                                      Result_P.Normal,
						          		                                      Result_P.dis,
						          		                                      Result_P.contact_point -Pos_A,
						          		                                      Result_P.contact_point -Pos_B,
						          		  					                  VelAng_A,
						          		  					                  VelAng_B,
						          		  			  					      P_typeA,
						          		  			  				          P_typeB,
						          		  			  					      Vel_A,
						          		  			  						  Vel_B,
						          		  			  					      ParticleObject[P_typeA].mass,
					            		  			  					      ParticleObject[P_typeB].mass,RelPos );
					  }
					  else/* No rotation if we had to apply a distance herustic */
					  {
					      Force_PairA = Collision_Response_Particle_Particle(     RVel,
							          		                                      Result_P.Normal,
							          		                                      Result_P.dis,
							          		                                      make_float3(0.0f),
							          		                                      make_float3(0.0f),
							          		  					                  VelAng_A,
							          		  					                  VelAng_B,
							          		  			  					      P_typeA,
							          		  			  				          P_typeB,
							          		  			  					      Vel_A,
							          		  			  						  Vel_B,
							          		  			  					      ParticleObject[P_typeA].mass,
						            		  			  					      ParticleObject[P_typeB].mass,RelPos );
					  }


			     /* We can sum contributions in local memory */
			      force  += Force_PairA.trans;/* normal + tangential */
			      forceL += Force_PairA.torque;/* Rotational */

			      num_coll++; /* Increment number of collisions */

				}/* End if collision */

			  }/* End checking neighbours of A */


		  /* Write total forces to global memory */
	      force_com[mem_index_A] = force;

	      if( SimParms.Rotation )
	      {
	        force_ang[mem_index_A] = forceL;
	      }
	      else
	      {
	    	force_ang[mem_index_A] = make_float3(0.0f);
	      }

	}/* End */


}
/*---------------------------------------------------------------------------*/



