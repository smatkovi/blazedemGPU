/* govender.nicolin@gmail.com */


/*---------------------------------------------------------------------------*/
/*     (1) Checks if Polyhedra is colliding with an infinite Plane           */
/*---------------------------------------------------------------------------*/
__device__ CollData Collision_Detection_Polyhedra_Plane( float3 *P_Vertex,
		                                                 uint     P_type,
		                                                 float3   Pos_P,
		                                                 float3   Vel_P,
		                                                 Surface  surface,
		                                                 float3  *surface_vertex  )
{

	CollData Result;

	Result.is_collision   = false;


	 /* 1) Apply Equation 1 to check position of COM */
	 if( dot( surface.normal, ( Pos_P - surface.centroid ) ) > 0.000000f )
	 {
         							   ;
          /* 3) Check distance between bounding sphere */
          if(  dot( surface.normal,( Pos_P - surface.centroid)) <= ParticleObject[P_type].radius )
          {
               /* 1. Check if a vertex is colliding */

        	    Result = Collision_Detection_Vertcies_Plane( ParticleObject[P_type].num_vertex,
        	    		                                       P_Vertex,
        	    		                                       surface.normal,
        	    		                                       surface.centroid       );

        	    /* No moment for face face contact */
        	    if( Result.collision_class==Face_Face )
        	    {
        	    	Result.contact_point = Pos_P;
        	    }

        	    return Result;

           } /* End Bound Collision */

       }/* End which side we on  */

   return Result;
}
/*---------------------------------------------------------------------------*/





              /* START Object Collision detection */

/*---------------------------------------------------------------------------*/
/*     (5) Checks if Polyhedra is colliding with a surface (area)            */
/*---------------------------------------------------------------------------*/
__device__ CollData Collision_Detection_Polyhedra_Surface( float3 *P_Vertex,
                                                           uint     P_type,
                                                           float3   Pos_P,
                                                           float3   Vel_P,
                                                           Surface  surface,
                                                           float3  *surface_vertex )
{

	CollData Result;

	Result.is_collision   = false;


	 /* 1) Apply Equation 1 to check position of COM */
	 if( dot( surface.normal, ( Pos_P - surface.centroid ) ) > 0.000000f )
	 {

          /* 3) Check distance between bounding sphere */
          if(  dot( surface.normal,( Pos_P - surface.centroid)) <= ParticleObject[P_type].radius )
          {
               /* 1. Check if a vertex is colliding */



        	    Result = Collision_Detection_Vertcies_Surface( ParticleObject[P_type].num_vertex,
        	    	                                      	   P_Vertex,
        	    		                                       surface_vertex, surface       );

        	    /* No moment for face face contact */
        	    if( Result.collision_class==Face_Face )
        	    {
        	    	Result.contact_point = Pos_P;
        	    }

        	    return Result;

           } /* End Bound Collision */

       }/* End which side we on  */

   return Result;
}
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*       (5) Checks if there is contact between 2 sets of edges
 *                       Normal points away from A                           */
/*---------------------------------------------------------------------------*/
__device__ CollData Collision_Detection_Polyhedra_Polyhedra_EdgesN( Edge   *EdgesA,   Edge *EdgesB,
		                                                            uint   NumEdgesA, uint  NumEdgesB,
                                                                    float3 RelPos,    float3 PosA,
                                                                    float3 PosB,      float3 RVel      )
{

	CollData Result;

	float  dA,dB;

	float3 pointA;
	float3 pointB;


	float min_edgedis = 0.50f;//ParticleObject[0].radius*0.10f;

	int num_edgeColl=0;

	float3 selected_PA;
	float3 selected_PB;

	Result.is_collision  = false;

	/* Step one find intersecting Edges */

	bool found_vaild_Edges = false;

	Result.contact_point = make_float3(0.0f);

	float3 first_normal = make_float3(0.0f);
	float avg_dis = 0.0f;
	bool is_many_normal =false;

	  /* Loop over all edges A */
	  for ( int i=0; i<NumEdgesA; i++ )
	  {
		     float DEAD = dot(EdgesA[i].dir,EdgesA[i].dir);

		     /* Loop over all edges B */
		     for ( int j=0; j < NumEdgesB; j++ )
		     {
		    	/* Normal Plane to Both Edges */
		        float3 EPlaneDir = cross( EdgesA[i].dir/length(EdgesA[i].dir), EdgesB[j].dir/length(EdgesB[j].dir) );
		        float mag    = dot(EPlaneDir,EPlaneDir);

		        /* 1. If not Parallel Check */
		        if( mag > 0.00000f )
		        {
		        	/* Choose Normal Such that it points in Direction of RelPos */
		        	EPlaneDir/=sqrtf(mag);
				    EPlaneDir = 1.0f*(dot (RelPos, EPlaneDir)/fabs( dot(RelPos, EPlaneDir) ) )*EPlaneDir;
				    EPlaneDir/=length(EPlaneDir);

                    float PA = dot( EPlaneDir,( PosA - (EdgesA[i].Point + EdgesB[j].Point)*0.50f) );
				    float PB = dot( EPlaneDir,( PosB - (EdgesA[i].Point + EdgesB[j].Point)*0.50f) );



				    /* 2. Apply heuristic that COM points must be on either side and the normal is sane */
				    if(  ((PA/fabs(PA))*(PB/fabs(PB)) < 0.0) &&  dot(EPlaneDir,RelPos)>0.70f )
				    {

				      //printf("passed plane herustics %d %d \n",i,j);

			          float3 DV      = EdgesA[i].Point - EdgesB[j].Point;
			          float  DEBD    = dot(EdgesB[j].dir,EdgesB[j].dir);
			          float  DEABD   = dot(EdgesA[i].dir,EdgesB[j].dir);
			          float  DEADAB  = dot(EdgesA[i].dir,DV);
			          float  DEBDAB  = dot(EdgesB[j].dir,DV);
			          float  detJinv = 1.00000f/(DEABD*DEABD - DEAD*DEBD);

		        	  dA = detJinv*(DEBD*DEADAB - DEABD*DEBDAB);

		              /* 3. Check for valid point on target Edge */
		              if( dA >= -0.0100f &&  dA <= 1.0100f )
		              {

		                /* 3.1 Get point on approaching object */
		                dB = detJinv*(DEABD*DEADAB - DEAD*DEBDAB);

		    	      /* Check both points are valid */
			          if( dB >= -0.0100f && dB <= 1.0100f )
			          {
					    /* 3. Compute the intersection points on the respective Edges */
					    pointA = EdgesA[i].Point + EdgesA[i].dir*dA;
					    pointB = EdgesB[j].Point + EdgesB[j].dir*dB;



					    /* Overlap distance between the edges */
					    float dis = fabs(dot( EPlaneDir, ( pointB - pointA ) ));

					    /* 4. Check that pen dis is sane  */
					    if( dis > 0.00000f )//&& (fabs(dis) < ParticleObject[0].radius*0.010f ))
					    {
						      /* Checks if B is inside A */
						      float3 B_PB = PosB - pointB;
						      float3 B_PA = PosB - pointA;

						      float mag_BPB = length(B_PB);
						      float mag_BPA = length(B_PA);

							  /* Checks if B is inside A */
							  float3 A_PB = PosA - pointA;
							  float3 A_PA = PosA - pointB;

							  float mag_APB = length(A_PB);
							  float mag_APA = length(A_PA);

						      /* 5. check there is actual penetration that points are close to each other */
					          if( (mag_BPB > mag_BPA) &&  dot( B_PB/mag_BPB, B_PA/mag_BPA )  > 0.99f && ( ( dot( A_PB/mag_APB, A_PA/mag_APA ) ) > 0.99f) )

					          {
					        	 // printf("Epair %d EL %d EP %d  angle %f :",num_edgeColl,i,j, dot(EPlaneDir,RelPos));
                                 // PrintVectorD(EPlaneDir);

					        	 /* Some logic to choose normal */

					    	     if ( num_edgeColl==0 )
					    	     {
					    		   first_normal = EPlaneDir;
					    	     }
					    	     else
					    	     {
					    		   if (length(first_normal-EPlaneDir)>0.10f)
					    		   {
					    			is_many_normal = true;
					    		   }
					    	     }


						         Result.contact_point += pointA;

						         avg_dis += dis;
						         num_edgeColl++;

						         found_vaild_Edges    = true;

						         if( fabs(dis) < min_edgedis )
						         {
						           min_edgedis          = fabs(dis);
						           Result.Normal        = EPlaneDir;

						           selected_PA          = pointA;
						           selected_PB          = pointB;
						         }
					        }
					      }

					     } /* End dis herustic */

			        }/* Valid points*/

		         } /* End dA valid */

		    } /* End Loop over Edges B */
		  }


	     } /* End Search */


	  if(found_vaild_Edges)
	  {

		//  printf("\n");

		Result.is_collision = true;

		//printf(" Num Edge %d diff normal %d \n",num_edgeColl,is_many_normal);

	    /* If the intersection point is inside the particle we have contact */
	    if(  !is_many_normal )
	    {
		    Result.Normal *= 1.0f/length(Result.Normal);

		  /* Just a single normal check that  */

	    	  /* Ray Points toward B in Direction of the contact vertex so Pos_A is the start */
                float3 Ray = selected_PA - PosA;
				         Ray/=length(Ray);

				  float dir = dot( Ray, Result.Normal);

				   /* check that there is a common component at-least */
			 if ( fabs(dir) > 0 )
			 {

			   Result.Normal = 1.0f*(dot (RelPos, Result.Normal)/fabs( dot(RelPos, Result.Normal) ) )*Result.Normal;
			   Result.Normal /=length(Result.Normal);

			   Result.dis    =  min(min_edgedis,fabs(dot( Result.Normal, ( selected_PB - selected_PA ) )));
			   Result.contact_point/=(float)num_edgeColl;
			  // printf("Single edge\n");

		     }
		  else/* If we choose the wrong normal close to an edge */
		  {
             //printf("using Rel pos Calc Normal: %f  %f  %f \n",Result.Normal.x,Result.Normal.y,Result.Normal.z);

			 //printf("Single wrong edge\n");
			 Result.Normal        = 1.0f*RelPos/length(RelPos);
		     Result.contact_point = PosA;
		     Result.dis           = fabs(dot( Result.Normal, ( selected_PB - selected_PA ) )) + 0.000010f;
		  }

    	  Result.collision_class = Edge_Edge;

     	  return Result;
	   }
	   else /* Use a default */
	   {
//		   if (ID==0)
//		   {
//		     printf("Multiple edge\n");
//		   }
			   Result.Normal = 1.0f*(dot (RelPos, Result.Normal)/fabs( dot(RelPos, Result.Normal) ) )*Result.Normal;
			   Result.Normal /=length(Result.Normal);
		     Result.contact_point /=(float)num_edgeColl;
		     Result.dis           = min(min_edgedis,ParticleObject[0].radius*0.010f);
	    	 Result.collision_class = Edge_Edge;
	     	 return Result;
	   }


	 }
	 else/* Must be a SP */
	 {
	      Result.collision_class = SP_PLANE;
	      Result.is_collision    = false;
	 }


	 /* If the intersection point is not inside the particle then no contact*/
 	 Result.collision_class = SP_PLANE;
 	 Result.is_collision    = false;
	 return Result;
}
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   (6) Checks if a vertex of A is in a face of B      */
/*---------------------------------------------------------------------------*/
__device__ CollData Collision_Detection_Polyhedara_Polyhedra_VertexFaces( float3  Pos_A,
		                                                             uint    num_vertex_A,
		                                                             float3 *VertexList_A,
		                                                             float3  Pos_B,
		                                                             uint    num_faces_B,
		                                                             Plane  *FaceList_B ,
		                                                             uint    num_vertex_B,
		                                                             float3 *VertexList_B,
	                                                                 POBJ    PObjectB,

	                                                                 int      P_typeB,
	                                                                 int indexA,
	                                                                 int indexB,
	                                                                 int *contact_face=0 )
{
     CollData Result; /* Function returned */

/*---------------------------------------------------------------------------*/
/*                           TYPE 1: VERTEX FACE                             */
/*---------------------------------------------------------------------------*/

     int   num_pen_vertex = 0;

     int         valid_vertex    [8];
     PointInPoly selected_vertex [8];
     PointInPoly temp;

	 Result.is_collision = false;
     num_pen_vertex=0;

     float max_pen = 0.050*ParticleObject[P_typeB].radius;

	 /* Loop over Vertexes of Particle A and Faces particle B */
	 for ( uint i=0; i<num_vertex_A; i++ )
	 {

	  		temp = Collision_Detection_Vertex_InPoly_3Faces( Pos_A,
                                                     VertexList_A[i],
	  				                                 num_faces_B,
	  				                                 FaceList_B,
	  				                                 Pos_B  );

	  		 if( temp.is_point_inside )
	  	     {

	  			 selected_vertex [num_pen_vertex] = temp;
	  	    	 valid_vertex    [num_pen_vertex] = i;
	  	    	 num_pen_vertex++;
	  	     }
	  	  } /* End of checking all planes of A */


	  	  if( num_pen_vertex > 0 )
	  	  {
	  		Result.is_collision = true;

		  	/* Single vertex face contact */
	  		if( num_pen_vertex==1 )
	  		{
	  		    Result.collision_class = Vertex_Face;

	           /* By default we use the face with the smallest distance */
		        Result.dis           = selected_vertex[0].distance[0];
		    	Result.Normal        = FaceList_B[selected_vertex[0].Close_Face_num[0]].normal;
		    	Result.contact_point = VertexList_A[valid_vertex[0]];


		    	/* only check if there is another face that is valid */
		    	if( selected_vertex[0].Close_Face_num[1]==-1 )
		    	{
		    	  if(SimParms.use_hist)
		    	  {
		    	   contact_face[0] = selected_vertex[0].Close_Face_num[0];
		    	  }
		    	  return Result;
		    	}
		    	else
		    	{
		    	  /* Ray Points toward B in Direction of the contact vertex so Pos_A is the start */
                  float3 Ray = VertexList_A[valid_vertex[0]] - Pos_A;
				       Ray/=length(Ray);

				  float dir = dot( Ray, FaceList_B[selected_vertex[0].Close_Face_num[0]].normal );

				   /* check that there is a common component at-least */
				   if ( fabs(dir) > 0 )
				   {
				    /* compute where the ray intersects the plane of the face  */
				    float d = dot(  (FaceList_B[selected_vertex[0].Close_Face_num[0]].centroid - Pos_A) ,
						    FaceList_B[selected_vertex[0].Close_Face_num[0]].normal)/dir;


				     if( is_PointOnSurface( PObjectB.face[selected_vertex[0].Close_Face_num[0]].num_vertex,
									  VertexList_B,
									  PObjectB.face[selected_vertex[0].Close_Face_num[0]].vertex_Order,
									  Pos_A + d*Ray,PObjectB.face[selected_vertex[0].Close_Face_num[0]].area,0.10f ))
				     {
				       if(SimParms.use_hist)
				       {
				    	 contact_face[0] = selected_vertex[0].Close_Face_num[0];
				       }
					   return Result;
				     }
				    }

		    	}



		    	/* check if 2nd face is valid */
		    	if( selected_vertex[0].Close_Face_num[1] > -1  )
		    	{
		    	  /* Ray Points toward B in Direction of the contact vertex so Pos_A is the start */
                  float3 Ray = VertexList_A[valid_vertex[0]] - Pos_A;
				       Ray/=length(Ray);

				  float dir = dot( Ray, FaceList_B[selected_vertex[0].Close_Face_num[1]].normal );

				   /* check that there is a common component at-least */
				   if ( fabs(dir) > 0 )
				   {
				    /* compute where the ray intersects the plane of the face  */
				    float d = dot(  (FaceList_B[selected_vertex[0].Close_Face_num[1]].centroid - Pos_A) ,
						    FaceList_B[selected_vertex[0].Close_Face_num[1]].normal)/dir;


				     if( is_PointOnSurface( PObjectB.face[selected_vertex[0].Close_Face_num[1]].num_vertex,
									  VertexList_B,
									  PObjectB.face[selected_vertex[0].Close_Face_num[1]].vertex_Order,
									  Pos_A + d*Ray,PObjectB.face[selected_vertex[0].Close_Face_num[1]].area,0.10f ))
				     {
				       Result.Normal = FaceList_B[selected_vertex[0].Close_Face_num[1]].normal;
				       Result.dis           = selected_vertex[0].distance[1];

				    	  if(SimParms.use_hist)
				    	  {
				       contact_face[0] = selected_vertex[0].Close_Face_num[1];
				    	  }
					   return Result;
				     }
				   }
		    	}


		    	/* check if 3rd face is valid */
		    	if( selected_vertex[0].Close_Face_num[2]>-1 )
		    	{
		    	  /* Ray Points toward B in Direction of the contact vertex so Pos_A is the start */
                  float3 Ray = VertexList_A[valid_vertex[0]] - Pos_A;
				       Ray/=length(Ray);

				  float dir = dot( Ray, FaceList_B[selected_vertex[0].Close_Face_num[2]].normal );

				   /* check that there is a common component at-least */
				   if ( fabs(dir) > 0 )
				   {
				    /* compute where the ray intersects the plane of the face  */
				    float d = dot(  (FaceList_B[selected_vertex[0].Close_Face_num[2]].centroid - Pos_A) ,
						    FaceList_B[selected_vertex[0].Close_Face_num[2]].normal)/dir;


				     if( is_PointOnSurface( PObjectB.face[selected_vertex[0].Close_Face_num[2]].num_vertex,
									  VertexList_B,
									  PObjectB.face[selected_vertex[0].Close_Face_num[2]].vertex_Order,
									  Pos_A + d*Ray,PObjectB.face[selected_vertex[0].Close_Face_num[2]].area,0.10f ) )
				     {
					       Result.Normal = FaceList_B[selected_vertex[0].Close_Face_num[2]].normal;
					       Result.dis    = selected_vertex[0].distance[2];

					    	  if(SimParms.use_hist)
					    	  {
					       contact_face[0] = selected_vertex[0].Close_Face_num[2];
					    	  }
					   return Result;
				     }
				   }
		    	}

				 //printf(" Warning VF Logic Error \n");
				 return Result;


	  		}


	  		/* EDGE FACE */
	  		if (num_pen_vertex==2)
	  		{
	  		  Result.is_collision    = true;
	  		  Result.collision_class = Edge_Face;

	  		  Result.contact_point   = (VertexList_A[valid_vertex[0]]
	  		                           + VertexList_A[valid_vertex[1]])*0.50;

	  		  Result.dis = fabs(dot(FaceList_B[selected_vertex[0].Close_Face_num[0]].normal,
	                     (Result.contact_point-FaceList_B[selected_vertex[0].Close_Face_num[0]].centroid)));

	  		  Result.Normal = FaceList_B[selected_vertex[0].Close_Face_num[0]].normal;



	    	/* only check if there is another face that is valid */
	    	if( selected_vertex[0].Close_Face_num[1]==-1 && selected_vertex[1].Close_Face_num[1]==-1 )
	    	{
		      if(SimParms.use_hist)
		      {
	    	    contact_face[0] = selected_vertex[0].Close_Face_num[0];
		      }
	    	  return Result;
	    	}
	    	else
	    	{
	    	  /* Ray Points toward B in Direction of the contact vertex so Pos_A is the start */
              float3 Ray = Result.contact_point - Pos_A;
			       Ray/=length(Ray);

			  float dir = dot( Ray, FaceList_B[selected_vertex[0].Close_Face_num[0]].normal );

			   /* check that there is a common component at-least */
			   if ( fabs(dir) > 0 )
			   {
			    /* compute where the ray intersects the plane of the face  */
			    float d = dot(  (FaceList_B[selected_vertex[0].Close_Face_num[0]].centroid - Pos_A) ,
					    FaceList_B[selected_vertex[0].Close_Face_num[0]].normal)/dir;


			     if( is_PointOnSurface( PObjectB.face[selected_vertex[0].Close_Face_num[0]].num_vertex,
								  VertexList_B,
								  PObjectB.face[selected_vertex[0].Close_Face_num[0]].vertex_Order,
								  Pos_A + d*Ray,PObjectB.face[selected_vertex[0].Close_Face_num[0]].area,0.10f ) )
			     {
			    	 /* first face is correct */
				        Result.dis           =  fabs(dot(FaceList_B[selected_vertex[0].Close_Face_num[0]].normal,
				        		                     (Result.contact_point-FaceList_B[selected_vertex[0].Close_Face_num[0]].centroid)));
				    	Result.Normal        = FaceList_B[selected_vertex[0].Close_Face_num[0]].normal;
				   return Result;
			     }

			    }

	    	}


	    	/* check if 2nd face is valid */
	    	if( selected_vertex[0].Close_Face_num[1]>-1  )
	    	{
		    	  /* Ray Points toward B in Direction of the contact vertex so Pos_A is the start */
	              float3 Ray = Result.contact_point - Pos_A;
				       Ray/=length(Ray);

				  float dir = dot( Ray, FaceList_B[selected_vertex[0].Close_Face_num[1]].normal );

				   /* check that there is a common component at-least */
				   if ( fabs(dir) > 0 )
				   {
				    /* compute where the ray intersects the plane of the face  */
				    float d = dot(  (FaceList_B[selected_vertex[0].Close_Face_num[1]].centroid - Pos_A) ,
						    FaceList_B[selected_vertex[0].Close_Face_num[1]].normal)/dir;


				     if( is_PointOnSurface( PObjectB.face[selected_vertex[0].Close_Face_num[1]].num_vertex,
									  VertexList_B,
									  PObjectB.face[selected_vertex[0].Close_Face_num[1]].vertex_Order,
									  Pos_A + d*Ray,PObjectB.face[selected_vertex[0].Close_Face_num[1]].area,0.10f) )
				     {
				    	 /* first face is correct */
					        Result.dis           =  fabs(dot(FaceList_B[selected_vertex[0].Close_Face_num[1]].normal,
					        		                     (Result.contact_point-FaceList_B[selected_vertex[0].Close_Face_num[1]].centroid)));
					    	Result.Normal        = FaceList_B[selected_vertex[0].Close_Face_num[1]].normal;

					    	  if(SimParms.use_hist)
					    	  {
					    	contact_face[0] = selected_vertex[0].Close_Face_num[1];
					    	  }
					   return Result;
				     }

				    }
	    	}


	    	/* check if 3rd face is valid */
	    	if(selected_vertex[0].Close_Face_num[2]>-1  )
	    	{
		    	  /* Ray Points toward B in Direction of the contact vertex so Pos_A is the start */
	              float3 Ray = Result.contact_point - Pos_A;
				       Ray/=length(Ray);

				  float dir = dot( Ray, FaceList_B[selected_vertex[0].Close_Face_num[2]].normal );

				   /* check that there is a common component at-least */
				   if ( fabs(dir) > 0 )
				   {
				    /* compute where the ray intersects the plane of the face  */
				    float d = dot(  (FaceList_B[selected_vertex[0].Close_Face_num[2]].centroid - Pos_A) ,
						    FaceList_B[selected_vertex[0].Close_Face_num[2]].normal)/dir;


				     if( is_PointOnSurface( PObjectB.face[selected_vertex[0].Close_Face_num[2]].num_vertex,
									  VertexList_B,
									  PObjectB.face[selected_vertex[0].Close_Face_num[2]].vertex_Order,
									  Pos_A + d*Ray,PObjectB.face[selected_vertex[0].Close_Face_num[2]].area,0.10f) )
				     {
				    	 /* first face is correct */
					        Result.dis           =  fabs(dot(FaceList_B[selected_vertex[0].Close_Face_num[2]].normal,
					        		                     (Result.contact_point-FaceList_B[selected_vertex[0].Close_Face_num[2]].centroid)));
					    	Result.Normal        = FaceList_B[selected_vertex[0].Close_Face_num[2]].normal;

					    	  if(SimParms.use_hist)
					    	  {
					    	contact_face[0] = selected_vertex[0].Close_Face_num[2];
					    	  }
					    	return Result;
				     }


				    }
	    	}



	  		   //printf(" Warning Edge-Face logic error \n");
	  		   return Result;

	  		}/* End Edge Face */


	  		/* FACE FACE: */

	  		if (num_pen_vertex>2)
	  		{
		  		Result.is_collision = true;

		  		Result.collision_class = Face_Face;
                Result.contact_point = make_float3(0.0f);

                /* Sum all penetrated vertex */
    		  	for (int i=0; i<num_pen_vertex; i++)
		  		{
		  		   Result.contact_point += VertexList_A[valid_vertex[i]];
		  		}
		  		Result.contact_point/=(float)num_pen_vertex;

				/* Ray Points toward B in Direction of the contact vertex so Pos_A is the start */
				float3 Ray = Result.contact_point - Pos_A;
					   Ray/=length(Ray);

				  float dir = dot( Ray, FaceList_B[selected_vertex[0].Close_Face_num[0]].normal );

				   /* check that there is a common component at-least */
				   if ( fabs(dir) > 0 )
				   {
					/* compute where the ray intersects the plane of the face  */
					float d = dot(  (FaceList_B[selected_vertex[0].Close_Face_num[0]].centroid - Pos_A) ,
							FaceList_B[selected_vertex[0].Close_Face_num[0]].normal)/dir;


					 if( is_PointOnSurface( PObjectB.face[selected_vertex[0].Close_Face_num[0]].num_vertex,
									  VertexList_B,
									  PObjectB.face[selected_vertex[0].Close_Face_num[0]].vertex_Order,
									  Pos_A + d*Ray,PObjectB.face[selected_vertex[0].Close_Face_num[0]].area,0.10f ) )
					 {
						 /* first face is correct */
							Result.dis           =  fabs(dot(FaceList_B[selected_vertex[0].Close_Face_num[0]].normal,
														 (Result.contact_point-FaceList_B[selected_vertex[0].Close_Face_num[0]].centroid)));
							Result.Normal        = FaceList_B[selected_vertex[0].Close_Face_num[0]].normal;

					    	  if(SimParms.use_hist)
					    	  {
							contact_face[0] = selected_vertex[0].Close_Face_num[0];
					    	  }
					        return Result;
					 }

					}


				/* check if 2nd face is valid */
				if(selected_vertex[0].Close_Face_num[1]>-1  )
				{
					  /* Ray Points toward B in Direction of the contact vertex so Pos_A is the start */
					  float3 Ray = Result.contact_point - Pos_A;
						   Ray/=length(Ray);

					  float dir = dot( Ray, FaceList_B[selected_vertex[0].Close_Face_num[1]].normal );

					   /* check that there is a common component at-least */
					   if ( fabs(dir) > 0 )
					   {
						/* compute where the ray intersects the plane of the face  */
						float d = dot(  (FaceList_B[selected_vertex[0].Close_Face_num[1]].centroid - Pos_A) ,
								FaceList_B[selected_vertex[0].Close_Face_num[1]].normal)/dir;


						 if( is_PointOnSurface( PObjectB.face[selected_vertex[0].Close_Face_num[1]].num_vertex,
										  VertexList_B,
										  PObjectB.face[selected_vertex[0].Close_Face_num[1]].vertex_Order,
										  Pos_A + d*Ray,PObjectB.face[selected_vertex[0].Close_Face_num[1]].area,0.010f))
						 {
							 /* second face is correct */
								Result.dis           =  fabs(dot(FaceList_B[selected_vertex[0].Close_Face_num[1]].normal,
															 (Result.contact_point-FaceList_B[selected_vertex[0].Close_Face_num[1]].centroid)));
								Result.Normal        = FaceList_B[selected_vertex[0].Close_Face_num[1]].normal;
						  		if (Result.dis>max_pen)
						  		{
						  			//printf("ERROR Second face to large \n");
						  		}
						    	  if(SimParms.use_hist)
						    	  {
						  		contact_face[0] = selected_vertex[0].Close_Face_num[1];
						    	  }
						        return Result;
						 }

						}
				}


				/* check if 3rd face is valid */
				if(selected_vertex[0].Close_Face_num[2]>-1  )
				{
					  /* Ray Points toward B in Direction of the contact vertex so Pos_A is the start */
					  float3 Ray = Result.contact_point - Pos_A;
						   Ray/=length(Ray);

					  float dir = dot( Ray, FaceList_B[selected_vertex[0].Close_Face_num[2]].normal );

					   /* check that there is a common component at-least */
					   if ( fabs(dir) > 0 )
					   {
						/* compute where the ray intersects the plane of the face  */
						float d = dot(  (FaceList_B[selected_vertex[0].Close_Face_num[2]].centroid - Pos_A) ,
								FaceList_B[selected_vertex[0].Close_Face_num[2]].normal)/dir;


						 if( is_PointOnSurface( PObjectB.face[selected_vertex[0].Close_Face_num[2]].num_vertex,
										  VertexList_B,
										  PObjectB.face[selected_vertex[0].Close_Face_num[2]].vertex_Order,
										  Pos_A + d*Ray,PObjectB.face[selected_vertex[0].Close_Face_num[2]].area,0.10f))
						 {
							 /* first face is correct */
								Result.dis           =  fabs(dot(FaceList_B[selected_vertex[0].Close_Face_num[2]].normal,
															 (Result.contact_point-FaceList_B[selected_vertex[0].Close_Face_num[2]].centroid)));
								Result.Normal        = FaceList_B[selected_vertex[0].Close_Face_num[2]].normal;

						    	if(SimParms.use_hist)
						    	{
								  contact_face[0] = selected_vertex[0].Close_Face_num[2];
						    	}

						  		if (Result.dis>max_pen)
						  		{
						  			//printf("ERROR Third face to large \n");
						  		}

						   return Result;
						 }


						}
				}

				Result.Normal = (Pos_A-Pos_B)/length((Pos_A-Pos_B));
			    //printf(" Warning PFace Logic Error %f \n",Result.dis);
		  		return Result;

	  		}/* End Face Face */


	  	  }


 return Result;

}
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*   (6) Checks if a vertex of A is in a face of B      */
/*---------------------------------------------------------------------------*/
__device__ CollData Collision_Detection_Polyhedara_VolumeObject_VertexFace( float3  Pos_P,
		                                                               uint    num_vertex_P,
		                                                               float3 *VertexList_P,

		                                                               uint    num_faces_DObject,
		                                                               Surface *FaceList_VObject ,
		                                                               uint    num_vertex_DObject,
		                                                               float3 *VertexList_DObject )
{
     CollData Result; /* Function returned */

/*---------------------------------------------------------------------------*/
/*                           TYPE 1: VERTEX FACE                             */
/*---------------------------------------------------------------------------*/

     int   num_pen_vertex = 0;

     int         valid_vertex    [8];
     PointInPoly selected_vertcies_penetrated_faces [8];
     PointInPoly temp;

	 Result.is_collision    = false;
	 Result.collision_class = UNDEF;

     num_pen_vertex=0;

     //float max_pen;

	 /* Loop over Vertexes of Particle and Faces of Volume Object */
	 for ( uint i=0; i<num_vertex_P; i++ )
	 {

	  		temp = Collision_Detection_Vertex_InVObject_3Faces( Pos_P,
                                                                VertexList_P[i],
	  				                                            num_faces_DObject,
	  				                                            FaceList_VObject   );

	  		 if( temp.is_point_inside )
	  	     {

	  			 selected_vertcies_penetrated_faces [num_pen_vertex] = temp;
	  	    	 valid_vertex    [num_pen_vertex] = i;
	  	    	 num_pen_vertex++;
	  	     }

	  } /* End of checking all Vertex of Particle */


	 /* We now have a list of vertices which have penetrated  */
	  	  if( num_pen_vertex > 0 )
	  	  {

	  		Result.is_collision = true;

		  	  /* Single vertex face contact */
	  		if( num_pen_vertex==1 )
	  		{

	  			//printf (" Single Vertex \n");

	  		    Result.collision_class = Vertex_Face;

	           /* By default we use the face with the smallest distance */
		        Result.dis           = selected_vertcies_penetrated_faces[0].distance[0];
		    	Result.Normal        = FaceList_VObject[ selected_vertcies_penetrated_faces[0].Close_Face_num[0] ].normal;

		    	Result.contact_point = VertexList_P[ valid_vertex[0] ];


		    	/* only check if there is another face that is valid */
		    	if( selected_vertcies_penetrated_faces[0].Close_Face_num[1] == -1 )
		    	{
		    		//printf (" Only one face valid  \n");
		    	  return Result;
		    	}
		    	else /* We are close to an edge so check which face is the most correct */
		    	{
		    	  /* Ray Points toward B in Direction of the contact vertex so Pos_A is the start */
                  float3 Ray = VertexList_P[valid_vertex[0]] - Pos_P;
				         Ray/=length(Ray);

				  float dir = dot( Ray, FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[0]].normal );

				   /* check that there is a common component at-least */
				   if ( fabs(dir) > 0 )
				   {
				    /* compute where the ray intersects the plane of the face  */
				    float d = dot(  (FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[0]].centroid - Pos_P) ,
				    		FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[0]].normal)/dir;


				     if( is_PointOnSurface( FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[0]].num_vertex,
									  VertexList_DObject,
									  FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[0]].vertex_Order,
									  Pos_P + d*Ray,FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[0]].area,0.1f ))
				     {
					   return Result;
				     }

				    }

		    	}
                /* Only if the ray between the COMs doesnt not intersect the closest face then we will look at the second face */


		    	/* check if 2nd face is valid */
		    	if(selected_vertcies_penetrated_faces[0].Close_Face_num[1]>-1  )
		    	{
		    	  /* Ray Points toward B in Direction of the contact vertex so Pos_A is the start */
                  float3 Ray = VertexList_P[valid_vertex[0]] - Pos_P;
				       Ray/=length(Ray);

				  float dir = dot( Ray, FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[1]].normal );

				   /* check that there is a common component at-least */
				   if ( fabs(dir) > 0 )
				   {
				    /* compute where the ray intersects the plane of the face  */
				    float d = dot(  (FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[1]].centroid - Pos_P) ,
				    		FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[1]].normal)/dir;


				     if( is_PointOnSurface( FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[1]].num_vertex,
									  VertexList_DObject,
									  FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[1]].vertex_Order,
									  Pos_P + d*Ray,FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[1]].area,0.10f ) )
				     {
				       Result.Normal = FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[1]].normal;
				       Result.dis    = selected_vertcies_penetrated_faces[0].distance[1];
					   return Result;
				     }
				   }
		    	}


		    	/* check if 3rd face is valid */
		    	if( selected_vertcies_penetrated_faces[0].Close_Face_num[2]>-1 )
		    	{
		    	  /* Ray Points toward B in Direction of the contact vertex so Pos_A is the start */
                  float3 Ray = VertexList_P[valid_vertex[0]] - Pos_P;
				       Ray/=length(Ray);

				  float dir = dot( Ray, FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[2]].normal );

				   /* check that there is a common component at-least */
				   if ( fabs(dir) > 0 )
				   {
				    /* compute where the ray intersects the plane of the face  */
				    float d = dot(  (FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[2]].centroid - Pos_P) ,
				    		FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[2]].normal)/dir;


				     if( is_PointOnSurface( FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[2]].num_vertex,
									  VertexList_DObject,
									  FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[2]].vertex_Order,
									  Pos_P + d*Ray,FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[2]].area,0.10f ) )
				     {
					       Result.Normal = FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[2]].normal;
					       Result.dis    = selected_vertcies_penetrated_faces[0].distance[2];
					   return Result;
				     }
				   }
		    	}

				 //printf(" VF Default \n");
				 return Result;


	  		}


	  		/* EDGE FACE we use the midpoint as the contact point */
	  		if ( num_pen_vertex==2 )
	  		{
	  		  Result.is_collision    = true;
	  		  Result.collision_class = Edge_Face;

	  		  Result.contact_point   = (VertexList_P[valid_vertex[0]]
	  		                           + VertexList_P[valid_vertex[1]])*0.50;

	  		  Result.Normal        = FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[0]].normal;

	  		  Result.dis = fabs( dot( FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[0]].normal,
	                            (Result.contact_point-FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[0]].centroid)) );



	    	/* only check if there is another face that is valid: we have */
	    	if( selected_vertcies_penetrated_faces[0].Close_Face_num[1]==-1 && selected_vertcies_penetrated_faces[1].Close_Face_num[1]==-1 )
	    	{
	    	  return Result;
	    	}
	    	else
	    	{
	    	  /* Ray Points toward B in Direction of the contact vertex so Pos_A is the start */
              float3 Ray = Result.contact_point - Pos_P;
			       Ray/=length(Ray);

			  float dir = dot( Ray, FaceList_VObject [selected_vertcies_penetrated_faces[0].Close_Face_num[0]].normal );

			   /* check that there is a common component at-least */
			   if ( fabs(dir) > 0 )
			   {
			    /* compute where the ray intersects the plane of the face  */
			    float d = dot(  (FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[0]].centroid - Pos_P) ,
			    		FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[0]].normal)/dir;


			     if( is_PointOnSurface( FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[0]].num_vertex,
								  VertexList_DObject,
								  FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[0]].vertex_Order,
								  Pos_P + d*Ray,FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[0]].area, 0.10f ) )
			     {
			    	 /* first face is correct */
				        Result.dis           =  fabs(dot(FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[0]].normal,
				        		                     (Result.contact_point-FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[0]].centroid)));
				    	Result.Normal        = FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[0]].normal;
				   return Result;
			     }

			    }

	    	}


	    	/* check if 2nd face is valid */
	    	if(selected_vertcies_penetrated_faces[0].Close_Face_num[1]>-1  )
	    	{
		    	  /* Ray Points toward B in Direction of the contact vertex so Pos_A is the start */
	              float3 Ray = Result.contact_point - Pos_P;
				       Ray/=length(Ray);

				  float dir = dot( Ray, FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[1]].normal );

				   /* check that there is a common component at-least */
				   if ( fabs(dir) > 0 )
				   {
				    /* compute where the ray intersects the plane of the face  */
				    float d = dot(  (FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[1]].centroid - Pos_P) ,
				    		FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[1]].normal)/dir;


				     if( is_PointOnSurface( FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[1]].num_vertex,
									  VertexList_DObject,
									  FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[1]].vertex_Order,
									  Pos_P + d*Ray,FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[1]].area,0.10f ) )
				     {
				    	 /* first face is correct */
					        Result.dis           =  fabs(dot(FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[1]].normal,
					        		                     (Result.contact_point-FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[1]].centroid)));
					    	Result.Normal        = FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[1]].normal;
					   return Result;
				     }

				    }
	    	}


	    	/* check if 3rd face is valid */
	    	if(selected_vertcies_penetrated_faces[0].Close_Face_num[2]>-1  )
	    	{
		    	  /* Ray Points toward B in Direction of the contact vertex so Pos_A is the start */
	              float3 Ray = Result.contact_point - Pos_P;
				       Ray/=length(Ray);

				  float dir = dot( Ray, FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[2]].normal );

				   /* check that there is a common component at-least */
				   if ( fabs(dir) > 0 )
				   {
				    /* compute where the ray intersects the plane of the face  */
				    float d = dot(  (FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[2]].centroid - Pos_P) ,
				    		FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[2]].normal)/dir;


				     if( is_PointOnSurface( FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[2]].num_vertex,
									  VertexList_DObject,
									  FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[2]].vertex_Order,
									  Pos_P + d*Ray, FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[2]].area, 0.10f ) )
				     {
				    	 /* first face is correct */
					        Result.dis           =  fabs(dot(FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[2]].normal,
					        		                     (Result.contact_point-FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[2]].centroid)));
					    	Result.Normal        = FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[2]].normal;
					   return Result;
				     }


				    }
	    	}
	  		   return Result;

	  		}/* End Edge Face */

	  	     float max_pen = 0.050f;
	  		/* FACE FACE: */

	  		  		if (num_pen_vertex>2)
	  		  		{
	  			  		Result.is_collision = true;
	  			  		Result.collision_class = Face_Face;
	  	                Result.contact_point = make_float3(0.0f);

	  	                /* Sum all penetrated vertex */
	  	    		  	for (int i=0; i<num_pen_vertex; i++)
	  			  		{
	  			  		   Result.contact_point += VertexList_P[valid_vertex[i]];
	  			  		}
	  			  		Result.contact_point/=(float)num_pen_vertex;

	  					/* Ray Points toward B in Direction of the contact vertex so Pos_P is the start */
	  					float3 Ray = Result.contact_point - Pos_P;
	  						   Ray/=length(Ray);

	  					  float dir = dot( Ray, FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[0]].normal );

	  					   /* check that there is a common component at-least */
	  					   if ( fabs(dir) > 0 )
	  					   {
	  						/* compute where the ray intersects the plane of the face  */
	  						float d = dot(  (FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[0]].centroid - Pos_P) ,
	  								FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[0]].normal)/dir;


	  						 if( is_PointOnSurface( FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[0]].num_vertex,
	  										  VertexList_DObject,
	  										FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[0]].vertex_Order,
	  										  Pos_P + d*Ray,FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[0]].area,0.10f) )
	  						 {
	  							 /* first face is correct */
	  								Result.dis           =  fabs(dot(FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[0]].normal,
	  															 (Result.contact_point-FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[0]].centroid)));
	  								Result.Normal        = FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[0]].normal;

//	  						    	  if(SimParms.use_hist)
//	  						    	  {
//	  								contact_face[0] = selected_vertcies_penetrated_faces[0].Close_Face_num[0];
//	  						    	  }
	  						        return Result;
	  						 }

	  						}


	  					/* check if 2nd face is valid */
	  					if(selected_vertcies_penetrated_faces[0].Close_Face_num[1]>-1  )
	  					{
	  						  /* Ray Points toward B in Direction of the contact vertex so Pos_P is the start */
	  						  float3 Ray = Result.contact_point - Pos_P;
	  							   Ray/=length(Ray);

	  						  float dir = dot( Ray, FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[1]].normal );

	  						   /* check that there is a common component at-least */
	  						   if ( fabs(dir) > 0 )
	  						   {
	  							/* compute where the ray intersects the plane of the face  */
	  							float d = dot(  (FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[1]].centroid - Pos_P) ,
	  									FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[1]].normal)/dir;


	  							 if( is_PointOnSurface( FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[1]].num_vertex,
	  											  VertexList_DObject,
	  											FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[1]].vertex_Order,
	  											  Pos_P + d*Ray,FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[1]].area,0.010f))
	  							 {
	  								 /* second face is correct */
	  									Result.dis           =  fabs(dot(FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[1]].normal,
	  																 (Result.contact_point-FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[1]].centroid)));
	  									Result.Normal        = FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[1]].normal;
	  							  		if (Result.dis>max_pen)
	  							  		{
	  							  			printf("ERROR Second face to large \n");
	  							  		}
//	  							    	  if(SimParms.use_hist)
//	  							    	  {
//	  							  		contact_face[0] = selected_vertcies_penetrated_faces[0].Close_Face_num[1];
//	  							    	  }
	  							        return Result;
	  							 }

	  							}
	  					}


	  					/* check if 3rd face is valid */
	  					if(selected_vertcies_penetrated_faces[0].Close_Face_num[2]>-1  )
	  					{
	  						  /* Ray Points toward B in Direction of the contact vertex so Pos_P is the start */
	  						  float3 Ray = Result.contact_point - Pos_P;
	  							   Ray/=length(Ray);

	  						  float dir = dot( Ray, FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[2]].normal );

	  						   /* check that there is a common component at-least */
	  						   if ( fabs(dir) > 0 )
	  						   {
	  							/* compute where the ray intersects the plane of the face  */
	  							float d = dot(  (FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[2]].centroid - Pos_P) ,
	  									FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[2]].normal)/dir;


	  							 if( is_PointOnSurface( FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[2]].num_vertex,
	  											  VertexList_DObject,
	  											FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[2]].vertex_Order,
	  											  Pos_P + d*Ray,FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[2]].area,0.10f))
	  							 {
	  								 /* first face is correct */
	  									Result.dis           =  fabs(dot(FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[2]].normal,
	  																 (Result.contact_point-FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[2]].centroid)));
	  									Result.Normal        = FaceList_VObject[selected_vertcies_penetrated_faces[0].Close_Face_num[2]].normal;

//	  							    	if(SimParms.use_hist)
//	  							    	{
//	  									  contact_face[0] = selected_vertcies_penetrated_faces[0].Close_Face_num[2];
//	  							    	}

	  							  		if (Result.dis>max_pen)
	  							  		{
	  							  			printf("ERROR Third face to large \n");
	  							  		}

	  							   return Result;
	  							 }


	  							}
	  					}

	  				    printf(" Warning Face Logic Error \n");
	  			  		return Result;

	  		  		}/* End Edge Face */


	  	  }


 return Result;

}
/*------------------------------------*/



/*___________________________________________________________________________*/
/*    (1)        Checks if there is Vertex Surface Contact                   */
/*___________________________________________________________________________*/
__device__ CollData Contact_Vertcies_VOsurface( uint    num_vertex,
		                                      float3   *P_VertexList,
		                                      float3   *surface_vertex,
		                                      Surface   surface         )
{
	CollData Result;

	float  dis     = 1.0; /* Max variable CM */
	int    num_vertexCount = 0;
	Result.collision_class = UNDEF; /* Hold for Valid Face */


	Result.contact_point = make_float3(0.0f); /* point on particle surface*/
	Result.is_collision  = false;





     /* Need to ensure that we dont select the wrong face */
     for( uint i=0; i<num_vertex; i++ )
     {

	    /* 5.1) Now get the _|_ intersection distance EQ 1*/

    	 Result.dis = dot( surface.normal, ( P_VertexList[i]  - surface.centroid ) );


        /* Check if EQ1 < 0 */
        if( Result.dis < 0.000000f )
        {

        	Result.is_collision = true;
            /* Ray Trace point to Surface Plane */

        	/* Check if the point is on the surface TOL= 1% */
        	if( is_PointOnSurface( surface.num_vertex, surface_vertex,
        	                               surface.vertex_Order,P_VertexList[i]  +fabs(Result.dis)*surface.normal,
        	                               surface.area,surface.Tol )       )

        	{

              /* Increment contact points */
                num_vertexCount++;

                /* Average contact points */
                Result.contact_point  += P_VertexList[i] ;

                /* Get the largest vertex penetration */
                if ( ( Result.dis < dis)  )
                {

  		           dis        = Result.dis;
  		           Result.collision_class = Vertex_Face;
  		        }

        	}

	    }/* End vertex plane penetration */

    }/* End loop over all vertex */

     /* If there is a SP then NO contact */
     if( !Result.is_collision )
     {
       Result.collision_class = SP_PLANE;
       return Result;
     }


    /* Return information for vertex-face contact */
    if( (num_vertexCount>0) && Result.collision_class==Vertex_Face )
    {
       Result.is_collision   = true;

       if(num_vertexCount==1)
       {
    	  Result.collision_class = Vertex_Face;
       }
       else if(num_vertexCount==2)
       {
    	   Result.collision_class = Edge_Face;
       }
       else
       {
         Result.is_collision   = true;
         Result.collision_class = Face_Face;
       }

       Result.Normal         = surface.normal;
       Result.dis            = fabs(dis);

       /* Contact point in world coordinates */
       Result.contact_point     /=((float)num_vertexCount);

       return Result;
    }

    /* Edge-Edge Candidate contact */
    Result.is_collision   = false;
    Result.collision_class = Edge_Edge;

    return Result;

}
/*___________________________________________________________________________*/


/*___________________________________________________________________________*/
/*   (5)  Checks for Polyhedra contacting with a Dynamic World Object        */
/*___________________________________________________________________________*/
__device__ CollData Collision_Detection_Polyhedra_DObject( uint     index,
														   uint     index_d,
														   uint     P_type,
														   float3   Pos_P,
														   POBJ_ROT PObjectRot,
														   float3   Vel_P )
{

     CollData Result;
     Result.is_collision = false;
     Result.collision_class = UNDEF;


     /* Check for possible collisions 1 surface at a time */
     for ( uint i=0; i < VolumeObject[index_d].num_faces; i++ )
     {

        Surface surface = 	VolumeObject[index_d].faces[i];

	   /* 1) Check COM pos is valid */
       if( dot( surface.normal, ( Pos_P - surface.centroid ) ) > 0.000000f )
	   {

		      /* 3) Check distance between bounding sphere and Plane */
			  if(  dot( surface.normal,( Pos_P - surface.centroid) )
			          <= ParticleObject[P_type].radius )
			  {
                 printf("checking if on surface \n");
			     /* 4) Check if a vertex is colliding with surface */
				 Result = Collision_Detection_Vertcies_Surface
						   ( ParticleObject[P_type].num_vertex,
								   PObjectRot.Vertex,
								   VolumeObject[index_d].vertex, surface );

			   } /* End Bound Collision */

			  if(Result.is_collision)
			  {
			    return Result;
			  }


		}/* End which side we on  */

     }/* End Loop over surfaces */

     Result.is_collision =false;
     return Result;

}
/*___________________________________________________________________________*/





/*___________________________________________________________________________*/
/*   (4.1) Determines the actual vertex in contact with a cylinder surface     */
/*___________________________________________________________________________*/
__device__ CollData Collision_Detection_Cylinder_PolyF( uint index,
		                                                     uint P_type,
		                                                     float3 Pos_P,
		                                                     float3 Vel_P,
		                                                     Macro_Cylinder Cyl,
		                                                     POBJ_ROT       POBJ_RT,
		                                                     int WobjNum          )
{

	CollData Result;
	Result.is_collision=false;

    float3 AXIS = Cyl.AXIS;
	float3 CapB = Cyl.center_bot_cap;

	float3 length_Axis = make_float3(1.0,1.0,1.0)-AXIS;

	if( (dot(Pos_P,length_Axis) > dot(Cyl.center_bot_cap,length_Axis) && (dot(Pos_P,length_Axis)) < dot(Cyl.center_top_cap,length_Axis)))
	{
	/* The first Cylinder dictates where the origin */
	 if (WobjNum>0 && SimParms.Simulation_Type==silo )
	 {
	   CapB = WorldObject[0].cyl_geo.center_bot_cap;
	   /*Only want to change the bottom cap along the length */
	   if(length_Axis.x==1)
	   {
	    CapB.x = Cyl.center_bot_cap.x;
	   }

	   if(length_Axis.y==1)
	   {
	    CapB.y = Cyl.center_bot_cap.y;
	   }

	   if(length_Axis.z==1)
	   {
	    CapB.z = Cyl.center_bot_cap.z;
	   }
	 }

	float3 Pvec = Pos_P*AXIS - CapB*AXIS;/* Relative to cent*/


	float Pvec_dis  = length(Pvec) + ParticleObject[P_type].radius;

    /* Check if we inside the cylinder */
    if (  Pvec_dis > Cyl.radius   )
    {

        /* Find colliding Vertex */
        for( uint i=0; i<ParticleObject[P_type].num_vertex; i++ )
        {

    	    /* Vertex position relative to cent of cylinder */
        	Pvec = ( POBJ_RT.Vertex[i] )*AXIS - CapB*AXIS  ;
            if (  length(Pvec) > Cyl.radius   )
            {

              /* get the normal to the surface */
    	      Result.Normal = -Pvec;
    	      Result.Normal *= (1.0f/length(Result.Normal));
              Result.contact_point = POBJ_RT.Vertex[i];
    	      Result.dis = fabs(length(Pvec) - Cyl.radius);
              Result.is_collision = true;
           }
        }
    }

	}/* Board pahse check */


   return Result;
}
/*___________________________________________________________________________*/



/*___________________________________________________________________________*/
/*   (4) Determines the actual vertex in contact with a cylinder surface     */
/*___________________________________________________________________________*/
__device__ CollData Collision_Detection_Polyhedra_Cylinder( uint           index,
		                                                    uint           P_type,
		                                                    float3         Pos_P,
		                                                    float3         Vel_P,
		                                                    Macro_Cylinder Cyl,
		                                                    POBJ_ROT       POBJ_RT    )
{

	CollData Result;
	Result.is_collision = false;

    float3 AXIS = Cyl.AXIS;
	float3 Pvec = Pos_P*AXIS - Cyl.center_bot_cap;/* Relative to cent*/


	float Pvec_dis  = length(Pvec) + ParticleObject[P_type].radius;

    /* Check if bounding sphere intersects cylinder */
    if (  Pvec_dis > Cyl.radius   )
    {

    	float3 v_point;       /* vertex of particle */
    	float  vertex_dis;    /* distance to coll   */

    	float  dis = 0.0;


    	int num_vertexCount = 0;
    	float3 s_point      = make_float3(0.0f);

        /* Find colliding Vertex */
        for( uint i=0; i<ParticleObject[P_type].num_vertex; i++ )
        {

    	    /* Vertex position relative to cent of cylinder */
    	    v_point = ( POBJ_RT.Vertex[i] )*AXIS - Cyl.center_bot_cap  ;

    	    float rc = length(v_point);

            /* Check if vertex point is intersecting */
           if( rc > Cyl.radius )
    	   {

         	  vertex_dis =  rc - Cyl.radius ;

        	     /* Increment contact points */
        	       num_vertexCount++;

        	       /* Average contact points */
        	       s_point  += POBJ_RT.Vertex[i];

        	            dis  += vertex_dis;

        	       /* Exit Face Contact */
        	       if(num_vertexCount>2)
        	       {
        	    	  dis*=1.0/((float)num_vertexCount);
        	         Result.is_collision   = true;
        	         Result.collision_class = Face_Face;

        	         Result.Normal  = -Pvec;//1.0f*( Cyl.normal_bot_cap - Result.body_point*Cyl.AXIS );
        	         Result.Normal *= (1.0f/length(Result.Normal));
        	         Result.dis            = min(0.00100f,fabs(dis));
                     /* no moment for face contact */
        	         Result.contact_point = Pos_P;

        	         return Result;

        	       }


    	   }

        }

        /* No Collision particle can be moved */
        if( num_vertexCount==0 )
        {
          Result.is_collision = false;
          return Result;
        }

	       if(num_vertexCount==2)
	       {

	    	   Result.is_collision   = true;
	    	   dis*=1.0/((float)num_vertexCount);
	         Result.is_collision   = true;
	         Result.collision_class = Edge_Face;
	         Result.contact_point     = s_point/((float)num_vertexCount);

	         Result.Normal  = -Pvec;//-1.0f*( Cyl.normal_bot_cap - Result.body_point*Cyl.AXIS );
	         Result.Normal *= (1.0f/length(Result.Normal));
	         Result.dis            = fabs(dis);
	         /* Contact point world coordinates */



	         return Result;

	       }

	       dis*=1.0/((float)num_vertexCount);

        /* get contact point */
        Result.contact_point     = s_point;
        /* get the normal to the surface */
        Result.Normal  = -Pvec;//-1.0f*( Cyl.normal_bot_cap - Result.body_point*Cyl.AXIS );
        Result.Normal *= (1.0f/length(Result.Normal));
        Result.dis     = dis;
        Result.is_collision   = true;
        Result.collision_class = Vertex_Face;


        return Result;


    }



   return Result;
}
/*___________________________________________________________________________*/



__device__ CollData Collision_Detection_Polyhedra_Polyhedra_EdgesV( Edge   *EdgesA,   Edge *EdgesB,
		                                                           uint   NumEdgesA, uint  NumEdgesB,
		                                                           float3 RelPos, float3 PosB   )
{

	CollData Result;

	float  dA,dB;
	float3 pointA;
	float3 pointB;


	float min_edgedis = 5.0f;

	int num_edgeColl=0;

	float3 selected_PA;
	float3 selected_PB;

	Result.is_collision = false;

	/* Step one find intersecting Edges */

	bool found_vaild_Edges=false;

	Result.contact_point=make_float3(0.0f);

	  /* Now loop over all edges A */
	  for ( int i=0;i<NumEdgesA;i++ )
	  {
		     float DEAD = dot(EdgesA[i].dir,EdgesA[i].dir);

		     /* Now loop over all edges B */
		     for ( int j=0;j<NumEdgesB;j++ )
		     {
		        float3 dir = cross(EdgesA[i].dir,EdgesB[j].dir);




		        float mag    = dot(dir,dir);


		        /* Parallel Check */
		        if( mag > 0.00000f )
		        {
			        float3 DV    = EdgesA[i].Point - EdgesB[j].Point;
			        float DEBD   = dot(EdgesB[j].dir,EdgesB[j].dir);

			        float DEABD  = dot(EdgesA[i].dir,EdgesB[j].dir);

			        float DEADAB = dot(EdgesA[i].dir,DV);
			        float DEBDAB = dot(EdgesB[j].dir,DV);

			        float detJinv = 1.00000f/(DEABD*DEABD - DEAD*DEBD);


		        	dA = detJinv*(DEBD*DEADAB - DEABD*DEBDAB);


		        /* Valid point on target Edge */
		        if( dA >= -0.010f &&  dA <= 1.0100f )
		        {
		           /* Get point on approaching object */
		          dB = detJinv*(DEABD*DEADAB - DEAD*DEBDAB);

		    	   /* Now both points are valid */
			       if( dB >= -0.0100f && dB <= 1.010f )
			       {

			    	 /* Get the points on the respective Edges */
			    	 pointA = EdgesA[i].Point + EdgesA[i].dir*dA;
			    	 pointB = EdgesB[j].Point + EdgesB[j].dir*dB;

			    	  /* Compute distance between the points */
			    	  DV = ( pointA - pointB );

			    	  float dis = dot(DV,DV);

			    	  /* If we have a non-zero distance */
			    	  if(dis > 0.00000f)
			    	  {

			    		dis  = sqrt(dis);
			    	    DV*=(1.0f/dis);




                           /* Find smallest edge */
			    	       if( dis < min_edgedis )
			    	       {

			    	    	   num_edgeColl++;
			    	    	   found_vaild_Edges    = true;
		    	               min_edgedis          = dis;

		    				   Result.contact_point = pointA;
		    				   Result.Normal = cross(EdgesA[i].dir, EdgesB[j].dir);
		    		           selected_PA =pointA;
		    		           selected_PB =pointB;
			    	       }

			    	  }

			        }/* Valid points*/

		         } /* End dA valid */

		    } /* End Loop over Edges B */
		  }


	     } /* End Search */


	  if(found_vaild_Edges)
	  {
		bool valid_Point = true;
	    Result.Normal *= 1.0f/length(Result.Normal);



	  if ( (length(PosB-selected_PB) < length(PosB-selected_PA))  )
	  {
		  valid_Point = false;
	  }

	  if(valid_Point)
	  {


		  Result.is_collision = true;

          /* TODO: Selecting Normal is not robust use RelPos as a FIX */
		  if( fabs( dot( RelPos, Result.Normal ) ) > 0.0f)
		  {
			   Result.Normal = (dot (RelPos, Result.Normal)/fabs( dot(RelPos, Result.Normal) ) )*Result.Normal;
			   Result.dis    = fabs(dot( Result.Normal, ( selected_PB - selected_PA ) ))+0.000010f;
		  }
		  else
		  {
			 Result.Normal        = RelPos/length(RelPos);
		     Result.contact_point = PosB;
		     Result.dis           = fabs(dot( Result.Normal, ( selected_PB - selected_PA ) ))+0.000010f;
		  }

    	  Result.collision_class = Edge_Edge;

     	  return Result;
	   }
	 }
	 else/* Must be a SP */
	 {
	      Result.collision_class = SP_PLANE;
	      Result.is_collision    = false;
	 }


 	 Result.collision_class = SP_PLANE;
 	 Result.is_collision    = false;
	 return Result;
}
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/


