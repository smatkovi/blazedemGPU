

/*---------------------------------------------------------------------------*/
/*     (1) Checks if a point is on the surface define by vertcies            */
/*---------------------------------------------------------------------------*/
__device__ bool is_PointOnSurface( int num_vertex, float3 *surface_vertex,
		                           int *v_order,float3 point, float SArea,
		                           float TOL  )
{
	float area = 0.000000f;

	/* Loop over all vertex
	 * Select pairs and compute the vector from the test point*/
	for( int i=0; i< (num_vertex-1); i++)
	{
	  /* Using the cross product to get the are between the two vector*/
	  area = area + 0.500000f*length(cross( (point - surface_vertex[ v_order[i]   ]), (point - surface_vertex[ v_order[i+1] ]) ));

	}
    /* Loop misses the area between first and last vertex*/
	area = area + 0.500000f*length( cross( (point - surface_vertex[ v_order[ num_vertex-1 ] ]), (point - surface_vertex[ v_order[ 0 ] ]) ) );

	if ( (fabs(area-SArea)/SArea)*100.0000f <= TOL)
	{
		return true;
	}

	return false;

}
/*---------------------------------------------------------------------------*/




/*---------------------------------------------------------------------------*/
         /*   (2) Checks if a point is bound by all planes        */
/*---------------------------------------------------------------------------*/
__device__ bool isPointIn_ConvexPoly( int     num_faces_A,
                                      float3  Pos_A,
                                      Plane   *faceListA,
                                      float3  point        )
{
   bool   is_collision_vp;
   float  vertex_dis;

/*---------------------------------------------------------------------------*/
/*                     TYPE 1: PLANE SEARCH FACE A                           */
/*---------------------------------------------------------------------------*/

      /* Loop over Planes of Particle A */
	  for ( uint i=0; i<num_faces_A; i++)
	  {

           is_collision_vp = false; /* check for intersection */

           /* _|_ distance between vertex and plane */
           vertex_dis = dot( faceListA[i].normal,( point - faceListA[i].centroid ) );

           /* if we past the plane report collision */
	       if( vertex_dis <= 0.000100f )
	       {
	          is_collision_vp = true; /* point inside*/

	       }

           if( !is_collision_vp )
           {
              return false;
           }
	   } /* End of checking all planes of A */

	  return true;
}





/*---------------------------------------------------------------------------*/
     /* (3) Quick Check for Vertex-Face SP O( F1V2 + F2V1 )= O(2FV) 6Bit Reg */
/*---------------------------------------------------------------------------*/
__device__ bool is_SP_FacesA_VertexB( uint num_faces_A,
		                              Plane *Faces_A,
		                              float3 Pos_B,
		                              uint num_vertex_B,
		                              float3 *Vertex_B)
{
   bool   is_contact;
   float  vertex_dis ;

/*---------------------------------------------------------------------------*/
/*                     TYPE 1: PLANE SEARCH FACE A                           */
/*---------------------------------------------------------------------------*/

      /* Loop over Planes of Particle A */
	  for ( uint i=0; i<num_faces_A; i++)
	  {
         is_contact = false; /* check for intersection */

         if( dot( Faces_A[i].normal, ( Pos_B - Faces_A[i].centroid ) )>0.0 )
         {
           /* Loop over Vertexes of Particle B */
           for ( uint j=0; j<num_vertex_B; j++ )
	       {

               /* _|_ distance between vertex and plane */
               vertex_dis = dot( Faces_A[i].normal,( Vertex_B[j] - Faces_A[i].centroid) );

               /* if we past the plane report collision */
	           if( vertex_dis <= 0.00000f )
	           {
	        	   is_contact = true;
	        	   break;
	           }

	        }/* End of Vertexes of Particle B */


             /* Found SP so exit */
             if( !is_contact )
             {
	           return true;
             }

          }/* end of COM check */

	  } /* End of checking all planes of A */


     return false;
}
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
  /* (4) Checks if single Vertex-Face SP O( F1V2 + F2V1 )= O(2FV) 6Bit Reg*/
/*---------------------------------------------------------------------------*/
__device__ bool is_SP_Face( Plane Faces_A,float3 Pos_B, uint num_vertex_B, float3 *Vertex_B)
{
   bool   is_contact;
   float  vertex_dis ;

/*---------------------------------------------------------------------------*/
/*                     TYPE 1: PLANE SEARCH FACE A                           */
/*---------------------------------------------------------------------------*/

         is_contact = false; /* check for intersection */

         if( dot( Faces_A.normal, ( Pos_B - Faces_A.centroid ) )>0.0 )
         {
           /* Loop over Vertexes of Particle B */
           for ( uint j=0; j<num_vertex_B; j++ )
	       {

               /* _|_ distance between vertex and plane */
               vertex_dis = dot( Faces_A.normal,( Vertex_B[j] - Faces_A.centroid) );

               /* if we past the plane report collision */
	           if( vertex_dis <= 0.00000f )
	           {
	        	   is_contact = true;
	        	   break;
	           }

	        }/* End of Vertexes of Particle B */


             /* Found SP so exit */
             if( !is_contact )
             {
	           return true;
             }

          }/* end of COM check */



     return false;
}
/*---------------------------------------------------------------------------*/




/*---------------------------------------------------------------------------*/
     /* (6) Quick Check for Vertex-Face SP O( F1V2 + F2V1 )= O(2FV)*/
/*---------------------------------------------------------------------------*/
__device__ bool is_SP_FacesVO_A_VertexB( uint num_faces_DObject,
		                                 Surface *Faces_DObject,
		                                 float3 Pos_B,
		                                 uint num_vertex_B,
		                                 float3 *Vertex_B        )
{
   bool is_contact;
   float  vertex_dis ;

/*---------------------------------------------------------------------------*/
/*                     TYPE 1: PLANE SEARCH FACE A                           */
/*---------------------------------------------------------------------------*/

      /* Loop over Planes of Particle A */
	  for ( uint i=0; i<num_faces_DObject; i++)
	  {
         is_contact = false; /* check for intersection */

         if( dot( Faces_DObject[i].normal, ( Pos_B - Faces_DObject[i].centroid ) )>0.0 )
         {
           /* Loop over Vertexes of Particle B */
           for ( uint j=0; j<num_vertex_B; j++ )
	       {

               /* _|_ distance between vertex and plane */
               vertex_dis = dot( Faces_DObject[i].normal,( Vertex_B[j] - Faces_DObject[i].centroid) );

               /* if we past the plane report collision */
	           if( vertex_dis <= 0.00000f )
	           {
	        	   is_contact =true;
	        	   break;
	           }

	        }/* End of Vertexes of Particle B */


             /* Found SP so exit */
             if( !is_contact )
             {
	           return true;
             }

          }/* end of COM check */

	  } /* End of checking all planes of A */


     return false;
}
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
/*            (7) Checks if there is Vertex Plane contact                    */
/*---------------------------------------------------------------------------*/
__device__ CollData Collision_Detection_Vertcies_Plane( uint     num_vertex,
		                                                float3  *Vertex_List,
		                                                float3   plane_normal,
		                                                float3   plane_centroid )
{
	CollData Result;

	Result.contact_point = make_float3(0.0f); /* Point of contact */

	float  vertex_dis;    /* _|_ distance to the surface*/

    int    num_vertexCount  = 0;
	float  dis_max          = 0.0; /* Max variable */


     /*  Loop  over all vertex of the particle */
     for( uint i=0; i< num_vertex; i++ )
     {
	    /* 5.1) Now get the _|_ intersection distance EQ 1*/
        vertex_dis = dot( plane_normal, ( Vertex_List[i]  - plane_centroid ) );

        /* Check if EQ1 < 0 */
        if( vertex_dis < 0.000000f )
        {
           num_vertexCount++;

           Result.contact_point  += Vertex_List[i] ; /* Average surface point */

           /* get largest  vertex penetration */
           if ( ( fabs(vertex_dis) > dis_max) )
           {
  		     dis_max     = fabs(vertex_dis);
  		   }

	    }

    }/* End all vertex */


    /* No Collision Exit */
    if( num_vertexCount==0 )
    {
      Result.is_collision    = false;
      Result.collision_class = SP_PLANE;
      return Result;
    }

    /* We have contact */
    Result.is_collision    = true;


    if( num_vertexCount==1 )
    {
        Result.collision_class = Vertex_Face;
    }
    if (num_vertexCount==2 )
    {
    	Result.collision_class = Edge_Face;
    }
    else
    {
    	Result.collision_class = Face_Face;
    }

    Result.contact_point /= (float)num_vertexCount;
    Result.Normal         =  plane_normal;
	Result.dis            =  dis_max;

    return Result;
}
/*---------------------------------------------------------------------------*/







/*---------------------------------------------------------------------------*/
/*    (8)        Checks if there is Vertex Surface Contact
 *               returns average contact point                  */
/*---------------------------------------------------------------------------*/
__device__ CollData Collision_Detection_Vertcies_Surface( uint     num_vertex,
		                                                  float3  *Vertex_List,
		                                                  float3  *surface_vertex,
		                                                  Surface  surface         )
{
	CollData Result;
	Result.is_collision  = false;
	Result.contact_point = make_float3(0.0f);

	float3  hit_point;   /* plane intersection point */
	float   vertex_dis;  /* _|_ distance to the surface*/

	float dis_max        = 0.0; /* Max variable CM */
    uint  num_vertexCount = 0;

    for( uint i=0; i<num_vertex; i++ ) /* Loop over Particle vertex */
    {

	    /* distance between vertex and surface */
        vertex_dis = dot( surface.normal, ( Vertex_List[i] - surface.centroid ) );

        /* if this is negative then possible contact */
        if( vertex_dis < 0.000000f )
        {

            /* Ray Trace point to Surface Plane normal points out of surface */
        	hit_point  = Vertex_List[i] + fabs(vertex_dis)*surface.normal;



        	/* Check if the point is on the surface TOL= 1% */
        	if( is_PointOnSurface( surface.num_vertex, surface_vertex,
        	                       surface.vertex_Order, hit_point,
        	                       surface.area,surface.Tol)     )
        	{


                /* Increment contact points */
                num_vertexCount++;

                /* Average contact points */
                Result.contact_point  += hit_point;

                /* Get the largest  vertex penetration */
                if ( ( fabs(vertex_dis) > dis_max ) )
                {
                   Result.is_collision   = true;
  		           dis_max               = fabs(vertex_dis);
  		        }

                /* face face contact so no moment */
                if( num_vertexCount > 2 )
                {

                  Result.collision_class = Face_Face;
                  Result.Normal          = surface.normal;
                  Result.dis             = dis_max;
                  /* Contact point world coordinates */
                  Result.contact_point   = make_float3(0.0f);

                  return Result;
                }
        	}

	    }/* End vertex plane penetration */

    }/* End loop over all vertex */

     /* If there is a SP then NO contact */
     if( Result.is_collision == false )
     {
       Result.collision_class = SP_PLANE;
       return Result;
     }

       Result.is_collision   = true;

       if( num_vertexCount==1 )
       {
           Result.collision_class = Vertex_Face;
       }
       else
       {
       	Result.collision_class = Edge_Face;
       }

       Result.Normal         = surface.normal;
       Result.dis            = dis_max;


       /* Contact point in world coordinates */
       Result.contact_point /= ((float)num_vertexCount);


       return Result;

}
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*    (3)        Checks if there is Vertex Surface Contact
 *               *NB: check tol for highly faceted objects                   */
/*---------------------------------------------------------------------------*/
__device__ CollData Collision_Detection_Vertcies_Face( uint     num_vertex,
		                                               float3  *Vertex_List,
		                                               float3  *face_vertex,
		                                               Plane    Face,
		                                               Surface  Face_Data)
{
	CollData Result;
	Result.is_collision  = false;
	Result.contact_point = make_float3(0.0f);
	Result.collision_class = UNDEF;


	float3  hit_point;   /* plane intersection point */
	float   vertex_dis;  /* _|_ distance to the Face*/

	float dis_max        = 0.0; /* Max variable CM */
    uint  num_vertexCount = 0;


    bool neg=false;
    for( uint i=0; i<num_vertex; i++ ) /* Loop over Particle vertex */
    {

	    /* distance between vertex and Face */
        vertex_dis = dot( Face.normal, ( Vertex_List[i] - Face.centroid ) );

        /* if this is negative then possible contact */
        if( vertex_dis < 0.000000f && vertex_dis > -0.000100f )
        {
           neg= true;


            /* Ray Trace point to Surface Plane normal points out of Face */
        	hit_point  = Vertex_List[i] + fabs(vertex_dis)*Face.normal;

        	/* Check if the point is on the surface TOL= 1% */
        	if( is_PointOnSurface( Face_Data.num_vertex, face_vertex,
        			               Face_Data.vertex_Order, hit_point,
        			               Face_Data.area,0.10f)     )
        	{

        		Result.is_collision   = true;

                /* Increment contact points */
                num_vertexCount++;

                /* Average contact points */
                Result.contact_point  += hit_point;

                /* Get the largest  vertex penetration */
                if ( ( fabs(vertex_dis) > dis_max ) )
                {
  		           dis_max               = fabs(vertex_dis);
  		        }

                /* face face contact so no moment */
                if( num_vertexCount > 2 )
                {

                  Result.collision_class = Face_Face;
                  Result.Normal          = Face.normal;
                  Result.dis             = dis_max;
                  /* Contact point world coordinates */
                  Result.contact_point   = make_float3(0.0f);

                  return Result;
                }
        	}


	    }/* End vertex plane penetration */

    }/* End loop over all vertex */

     /* If there is a SP then NO contact */
     if( !neg )
     {
       Result.collision_class = SP_PLANE;
       return Result;
     }

       if (Result.is_collision)
       {


       if( num_vertexCount==1 )
       {
           Result.collision_class = Vertex_Face;
       }
       else
       {
       	Result.collision_class = Edge_Face;
       }

       Result.Normal         = Face.normal;
       Result.dis            = dis_max;

       /* Contact point in world coordinates */
       Result.contact_point /= ((float)num_vertexCount);
       }

       return Result;

}
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/*    (3)        Checks if there is Vertex Surface Contact
 *               *NB: check tol for highly faceted objects                   */
/*---------------------------------------------------------------------------*/
__device__ CollData Collision_Detection_Vertcies_VFace( uint     num_vertex,
		                                               float3  *Vertex_List,
		                                               float3  *face_vertex,
		                                               Surface  Face_Data,
		                                               float3    Pos_P)
{
	CollData Result;
	Result.is_collision  = false;
	Result.contact_point = make_float3(0.0f);
	Result.collision_class = UNDEF;


	float3  hit_point;   /* plane intersection point */
	float   vertex_dis;  /* _|_ distance to the Face*/

	float dis_max        = 0.0; /* Max variable CM */
    uint  num_vertexCount = 0;

    bool neg=false;

    if( dot( Face_Data.normal,( Pos_P - Face_Data.centroid ) ) > 0.0 )
    {

    for( uint i=0; i<num_vertex; i++ ) /* Loop over Particle vertex */
    {

	    /* distance between vertex and Face */
        vertex_dis = dot( Face_Data.normal, ( Vertex_List[i] - Face_Data.centroid ) );

        /* if this is negative then possible contact */
        if( vertex_dis < 0.000000f && vertex_dis > -0.000100f)
        {
           neg= true;


            /* Ray Trace point to Surface Plane normal points out of Face */
        	hit_point  = Vertex_List[i] + fabs(vertex_dis)*Face_Data.normal;

        	/* Check if the point is on the surface TOL= 1% */
        	if( is_PointOnSurface( Face_Data.num_vertex, face_vertex,
        			               Face_Data.vertex_Order, hit_point,
        			               Face_Data.area,Face_Data.Tol)     )
        	{

        		Result.is_collision   = true;

                /* Increment contact points */
                num_vertexCount++;

                /* Average contact points */
                Result.contact_point  += hit_point;

                /* Get the largest  vertex penetration */
                if ( ( fabs(vertex_dis) > dis_max ) )
                {
  		           dis_max               = fabs(vertex_dis);
  		        }

                /* face face contact so no moment */
                if( num_vertexCount > 2 )
                {

                  Result.collision_class = Face_Face;
                  Result.Normal          = Face_Data.normal;
                  Result.dis             = dis_max;
                  /* Contact point world coordinates */
                  Result.contact_point   = make_float3(0.0f);

                  return Result;
                }
        	}


	    }/* End vertex plane penetration */

    }/* End loop over all vertex */

     /* If there is a SP then NO contact */
     if( !neg )
     {
       Result.collision_class = SP_PLANE;
       return Result;
     }

       if (Result.is_collision)
       {


       if( num_vertexCount==1 )
       {
           Result.collision_class = Vertex_Face;
       }
       else
       {
       	Result.collision_class = Edge_Face;
       }

       Result.Normal         = Face_Data.normal;
       Result.dis            = dis_max;

       /* Contact point in world coordinates */
       Result.contact_point /= ((float)num_vertexCount);
       }
    }

       return Result;

}
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
/*   (4) Returns the Closest 3 faces to a vertex inside a polyhedra
 *               bool is_point inside                                        */
/*---------------------------------------------------------------------------*/
__device__ PointInPoly Collision_Detection_Vertex_InPoly_3Faces( float3  Pos_P,
		                                                         float3  P_vertex,
                                                                 uint    num_faces,
                                                                 Plane  *Faces,
                                                                 float3  Pos_Faces )
{

      bool   is_inside_plane;
	  float  vertex_dis ;

	  PointInPoly res;

      res.is_point_inside = false;

	  res.Close_Face_num[0] = -1;
	  res.Close_Face_num[1] = -1;
	  res.Close_Face_num[2] = -1;

	  res.distance[0] = 10.0f;/* should use pdata.bradius*2 */
	  res.distance[1] = 10.0f;
	  res.distance[2] = 10.0f;

      /* Loop over Planes of Particle A */
	  for ( uint i=0; i< num_faces; i++)
	  {
         is_inside_plane = false;

         /* _|_ distance between vertex and face */
         vertex_dis = dot(Faces[i].normal,( P_vertex - Faces[i].centroid ));

         /* if we past the plane report collision */
	     if( vertex_dis < -0.0000100f )
	     {
	        is_inside_plane     = true; /* point inside*/

	        /* We only want faces that the COM of the particle is on the correct side of */
	        if( dot( Faces[i].normal,( Pos_P - Faces[i].centroid ) ) > 0.0 )
	        {
	            if(fabs( vertex_dis) < res.distance[0] )
	            {
	               /* update the other 2 distances */
	               if( res.Close_Face_num[0] > -1 )
	               {
		              if( res.Close_Face_num[1] >-1 )
	                  {
		            	 res.distance  [2] = res.distance  [1];
		            	 res.Close_Face_num[2] = res.Close_Face_num[1];
	                  }
	            	  res.distance  [1] = res.distance  [0];
	            	  res.Close_Face_num[1] = res.Close_Face_num[0];
	               }

	                res.distance  [0]   = fabs(vertex_dis);
	            	res.Close_Face_num[0] = i;

	              }
	              else if( fabs(vertex_dis) < res.distance[1] )
	              {
	            	  /* update 3rd close face */
	            	  if( res.Close_Face_num[1] > -1 )
                      {
	            	     res.distance  [2] = res.distance  [1];
	            	 	 res.Close_Face_num[2] = res.Close_Face_num[1];
                      }

	            	  res.distance  [1] = fabs(vertex_dis);
	            	  res.Close_Face_num[1] = i;
	              }
	              else if( fabs(vertex_dis) < res.distance[2] )
	              {
	            	  res.distance  [2] = fabs(vertex_dis);
	            	  res.Close_Face_num[2] = i;
	              }

	           }/* end of checking COM_A */

	         }/* End if we past a face */

             /* early exit We must be past all faces for there to be contact */
             if( !is_inside_plane )
             {
            	return res;
             }

	     } /* End of checking all faces */

	  res.is_point_inside = true;
	  return res;
}



/*---------------------------------------------------------------------------*/
/*   (4.1) Returns the Closest 3 faces to a vertex inside a Dpolyhedra
 *               bool is_point inside                                        */
/*---------------------------------------------------------------------------*/
__device__ PointInPoly Collision_Detection_Vertex_InVObject_3Faces(
		                                                 float3  Pos_P,
		                                                 float3  P_vertex,
                                                         uint    num_faces,
                                                         Surface  *Faces)
{

      bool   is_inside_plane;
	  float  vertex_dis ;

	  PointInPoly res;

      res.is_point_inside = false;

	  res.Close_Face_num[0] = -1;
	  res.Close_Face_num[1] = -1;
	  res.Close_Face_num[2] = -1;

	  res.distance[0] = 10.0f;
	  res.distance[1] = 10.0f;
	  res.distance[2] = 10.0f;


      /* Loop over Planes of Dynamic Objects */
	  for ( uint i=0; i< num_faces; i++ )
	  {
         is_inside_plane = false;

         /* _|_ distance between vertex and face */
         vertex_dis = dot( Faces[i].normal,( P_vertex - Faces[i].centroid ));



         /* if we past the plane report collision */
	     if( vertex_dis <= 0.000010f )
	     {
	        is_inside_plane     = true; /* point inside*/


	        /* We only want faces that the COM of the particle is on the correct side of */
	        if( dot( Faces[i].normal,( Pos_P - Faces[i].centroid ) ) > 0.0f )
	        {

	            if( fabs( vertex_dis) < res.distance[0] )
	            {
	               /* update the other 2 distances */
	               if( res.Close_Face_num[0] > -1 )
	               {
		              if( res.Close_Face_num[1] >-1 )
	                  {
		            	 res.distance  [2] = res.distance  [1];
		            	 res.Close_Face_num[2] = res.Close_Face_num[1];
	                  }
	            	  res.distance  [1] = res.distance  [0];
	            	  res.Close_Face_num[1] = res.Close_Face_num[0];
	               }

	                res.distance  [0]   = fabs(vertex_dis);
	            	res.Close_Face_num[0] = i;

	              }
	              else if( fabs(vertex_dis) < res.distance[1] )
	              {
	            	  /* update 3rd close face */
	            	  if( res.Close_Face_num[1] > -1 )
                      {
	            	     res.distance  [2] = res.distance  [1];
	            	 	 res.Close_Face_num[2] = res.Close_Face_num[1];
                      }

	            	  res.distance  [1] = fabs(vertex_dis);
	            	  res.Close_Face_num[1] = i;
	              }
	              else if( fabs(vertex_dis) < res.distance[2] )
	              {
	            	  res.distance  [2] = fabs(vertex_dis);
	            	  res.Close_Face_num[2] = i;
	              }

	           }/* end of checking COM_A */

	         }/* End if we past a face */

             /* early exit We must be past all faces for there to be contact */
             if( !is_inside_plane )
             {
            	return res;

             }

	     } /* End of checking all faces */

	  res.is_point_inside = true;
	  return res;
}



/*---------------------------------------------------------------------------*/
   /*   (3) Returns all the average of vertcies that are inside the convex poly
          * and returns the average contact point and number of vertcies
          *    26 Bit registers */
/*---------------------------------------------------------------------------*/
__device__ float4  VertciesIn_ConvexPoly( int     num_faces_A,
                                          Plane   *faceListA,
                                          int     num_vertex,
                                          float3  *vertexList   )
{

  bool is_collision_vp;
  float  vertex_dis ;

/*---------------------------------------------------------------------------*/
/*                     TYPE 1: PLANE SEARCH FACE A                           */
/*---------------------------------------------------------------------------*/

   int    num_faces_in;/* Counter number of passed faces */
   float3 c_point=make_float3(0.0f);
   int    num_pen=0;

   for ( uint k=0; k<num_vertex; k++)
   {
      /* Loop over Planes of Particle A */
	  for ( uint i=0; i<num_faces_A; i++)
	  {

           is_collision_vp = false; /* check for intersection */

           /* _|_ distance between vertex and plane */
           vertex_dis = dot( faceListA[i].normal,( vertexList[k] - faceListA[i].centroid ) );

           /* if we past the plane report collision */
	       if( vertex_dis < 0.00000f )
	       {
	          is_collision_vp = true; /* point inside*/
	          num_faces_in++;
	       }

           if( !is_collision_vp )
           {
              break;
           }
	   } /* End of checking all faces of A */

	  /* Only if all faces passed then store as a contact point */
	   if (num_faces_in==num_faces_A)
	   {
		   c_point += vertexList[k];
		   num_pen++;
	   }

   	}/* End Checking All vertcies */

   if(num_pen>1)
   {
	   c_point/=(float)num_pen;
   }

	 return make_float4(c_point,(float)num_pen);
}
/*---------------------------------------------------------------------------*/

