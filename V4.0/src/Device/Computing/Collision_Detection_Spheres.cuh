/* govender.nicolin@gmail.com */

/*___________________________________________________________________________*/
/*     (1)    Checks if a sphere is in contact with a surface                */
/*___________________________________________________________________________*/
__device__ CollData_Sphere Collision_Detection_Sphere_Surface(
		                                                uint P_type,
												        float3  Pos_P,
												        Surface surface,
												        float3  *vertex    )
{

    CollData_Sphere Result;
    Result.is_collision = false;


      /* Check distance between sphere Pt and surface  */
      Result.dis = dot( surface.normal,( (Pos_P - ParticleObject[P_type].radius*surface.normal) - surface.centroid));

      /* Check if penetration is valid */
      if(  Result.dis < 0.000000f && Result.dis > -ParticleObject[P_type].radius )
      {

  	    Result.dis = fabs(Result.dis);

        /* Now check if the point is on the surface*/


    /* 1. Sphere Surface contact */
        if( is_PointOnSurface( surface.num_vertex,
                vertex,
               surface.vertex_Order, (Pos_P - ParticleObject[P_type].radius*surface.normal) + surface.normal*Result.dis,
               surface.area,surface.Tol                 ) )
        {
      	   Result.is_collision   = true;
      	   Result.contact_point  = (Pos_P - ParticleObject[P_type].radius*surface.normal);
      	   Result.Normal         = surface.normal;
           return Result;
         }

      }

   return Result;
}
/*___________________________________________________________________________*/




/*___________________________________________________________________________*/
/*       (2) Checks if a sphere is in contact with a infinite Plane          */
/*___________________________________________________________________________*/
__device__ CollData_Sphere Collision_Detection_Sphere_Plane (
		                                     uint P_type,
		                                     float3 Pos_P,
		                                     float3 surface_normal,
		                                     float3 surface_centroid )
{

	 CollData_Sphere Result;
	    Result.is_collision = false;

	      /* Get the point on the surface of sphere */
	      float3  Pt =  Pos_P - ParticleObject[P_type].radius*surface_normal ;

	      /* Check distance between sphere and surface  */
	      Result.dis = dot( surface_normal,( Pt - surface_centroid));

	      /* Check if penetration is valid */
	      if(  Result.dis < 0.000000f && Result.dis > -ParticleObject[P_type].radius )
	      {
	  	    Result.dis = fabs(Result.dis);

	        /* Get the point on the plane */
	        //float3 point  = Pt + surface_normal*Result.dis;


	      	  Result.is_collision   = true;
	      	  Result.contact_point     = Pt;
	      	  Result.Normal         = surface_normal;
	          return Result;

	      }


	   return Result;
}






/*___________________________________________________________________________*/
__device__ CollData Collision_Detection_Sphere_DObject ( uint   index,
		                                                 uint   index_d,
		                                                 float3 Pos_P,
                        		                         float3 Vel_P,
                        		                         Edge   *EdgesD,
                        		                         uint    num_edgeD,
                        		                         float radius,
                        		                         Contact_Info Lifter_Contact_Hist)
{

    CollData Result;
    Result.is_collision = false;

    if( Lifter_Contact_Hist.cont_type!=2 )
    {
     /* Check for possible collisions 1 surface at a time */
     for ( int j=0; j < VolumeObject[index_d].num_faces; j++ )
     {
        Surface surface = VolumeObject[index_d].faces[j];

        /* Get the point on the surface of sphere */
        float3  Pt =  Pos_P - radius*surface.normal ;

        /* Check distance between sphere and surface  */
        Result.dis = dot( surface.normal,( Pt - surface.centroid));

        /* Check if penetration is valid */
        if(  Result.dis < 0.000000f && Result.dis > -radius*0.50f )
        {
    	  Result.dis = fabs(Result.dis);

          /* Get the point on the plane */
          float3 point  = Pt + surface.normal*Result.dis;

          /* Now check if the point is on the surface*/
          if ( is_PointOnSurface( surface.num_vertex,
                  		                        VolumeObject[index_d].vertex,
                                                  surface.vertex_Order, point,
                                                  surface.area, surface.Tol)       )
          {
        	  Result.is_collision   = true;
        	  Result.contact_point  = Pt;
        	  Result.Normal         = surface.normal;
        	  Lifter_Contact_Hist.snum = j;
        	  Lifter_Contact_Hist.cont_type = 1;
              return Result;
           }


       } /* End Potential Collision */

     }/* End check over all surfaces */

    }
    /* End surface check */


     /* 2. Check Edges */
	   for ( uint i=0; i<num_edgeD; i++ )
	   {

		 float3 Edir = EdgesD[i].dir/length(EdgesD[i].dir);


	     /* Vector between Sphere and Edge Start  */
	     float3 vpc = Pos_P - EdgesD[i].Point ;
         float dVec = dot(Edir,vpc);


	     /* Check if sphere projects onto the Edge */
	     if(  dVec > 0 )
	     {

	      /* Project the Sphere onto the Edge */
	      float3 P2 = EdgesD[i].Point + dVec*Edir;

	      float3 d2 = Pos_P - P2;
	      float  d  = length(d2) ;

	      /* We have a collision */
	     if( d < radius )
	     {

           /* Penetration distance */
	    	Result.dis = fabs(d -radius);

	    	/* Limit of Half penetration */
	        if( Result.dis < radius )
	        {

	           float3 COM_Poly = VolumeObject[index_d].COM;
	           COM_Poly.z = Pos_P.z;

               float3  NPoly   = P2 - COM_Poly;
               float3  NSphere = Pos_P -P2 ;

               Result.Normal = (NSphere);
               Result.Normal*=(1.0f/length(Result.Normal));

               Result.contact_point =  make_float3(0.0f);

               /* No spin for edge contact */
               Lifter_Contact_Hist.cont_type = 2;
               Result.collision_class = Edge_Face;
	    	   Result.is_collision  = true;

	    	   return Result;
	         }

	       }

	     }

	   }


   return Result;
}
/*___________________________________________________________________________*/




/*___________________________________________________________________________*/
/*   (4.1) Determines the actual vertex in contact with a cylinder surface     */
/*___________________________________________________________________________*/
__device__ CollData_Sphere Collision_Detection_Sphere_Cylinder( uint index, uint P_type,
		                                       float3 Pos_P, float3 Vel_P,
		                                       Macro_Cylinder Cyl, int WobjNum          )
{

	CollData_Sphere Result;
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
        /* get the normal to the surface */
    	Result.Normal = -Pvec;
    	Result.Normal *= (1.0f/length(Result.Normal));
        Result.contact_point = Pos_P -Result.Normal*ParticleObject[P_type].radius;
    	Result.dis = fabs(Pvec_dis -Cyl.radius);
        Result.is_collision = true;

		//printf("contact %f cylrad %f \n",Pvec_dis, Cyl.radius);
		//PrintVectorD(AXIS);
		//PrintVectorD(Cyl.center_bot_cap);
        /* Ray trace so particle is on the surface */
        //position_com[index]+=Result.dis*Result.Normal;
    }

	}/* Board pahse check */


   return Result;
}
/*___________________________________________________________________________*/




__device__ CollData_Sphere Collision_Detection_Sphere_ConeLL( uint index, uint P_type,
		                                       float3 Pos_P, float3 Vel_P,
		                                       Macro_Cylinder Cone          )
{
	CollData_Sphere Result;
	Result.is_collision=false;


	if ( Pos_P.y >= Cone.center_bot_cap.y && Pos_P.y < Cone.center_top_cap.y)
	{
	float alpha = atan( ((Cone.radius-Cone.radius_bot)/Cone.height));

	Cone.AXIS = Cone.center_top_cap - Cone.center_bot_cap;
	Cone.AXIS/=length(Cone.AXIS);
	
	/* Ray  cone  bot to particle pos */
	float3 BP = Pos_P - Cone.center_bot_cap; 

	/* 1) Project Ray onto the cone axis */
	float3 BP_Proj = dot(BP,Cone.AXIS)*Cone.AXIS;
	float h_BP_Proj = length(BP_Proj);

	/* 2) Shortest distance between pos and  cone axis */
	float3 P_Pp = cross(BP,Cone.AXIS);
	float  dist_P_Pp = length(P_Pp);

	/* 3) Project Pos onto axis but in direction of normal from cone surface */
	float h_BPp_Proj = ((Cone.radius-Cone.radius_bot)/Cone.height)*dist_P_Pp + h_BP_Proj;

	/* 4) Shortest distance between Pos and Cone surface _|_ */
	float dist_P_cone = cos(alpha) * ( (Cone.radius - Cone.radius_bot) / Cone.height*h_BP_Proj + Cone.radius_bot - dist_P_Pp );

	if( dist_P_cone < ParticleObject[P_type].radius )
	{
		Result.is_collision = true;

	    Result.dis = fabs(dist_P_cone -ParticleObject[P_type].radius);
         
         // Compute Normal at the contact
     
		float3 PostPp = (Cone.center_bot_cap + h_BPp_Proj*Cone.AXIS);
		Result.Normal = Pos_P - PostPp;         
        Result.Normal*=(-1.0f/length(Result.Normal));
         // And contact point
		Result.contact_point =  Pos_P -1.0f*dist_P_cone*Result.Normal; 
		
		//PrintVectorND(Pos_P);
		//PrintVectorND(Result.contact_point);
		//PrintVectorND(Result.Normal);
		//printf(" Pen %f \n",Result.dis);
	}
	}
	return Result;

}

/*___________________________________________________________________________*/
/*   (4.2) Determines the actual vertex in contact with a cylinder surface   */
/*___________________________________________________________________________*/
__device__ CollData isContact_Sphere_CylinderR( uint index, uint P_type,
		                                       float3 Pos_P, float3 Vel_P,
		                                       Macro_Cylinder Cyl          )
{

	CollData Result;
	Result.is_collision=false;

   float3 pos= Pos_P - make_float3(Cyl.radius,Cyl.radius,Pos_P.z);

	float distToCenter = length(pos) ;/* Relative to cent*/


	float dis  = fabs(distToCenter - (Cyl.radius -ParticleObject[P_type].radius));

    /* Check if we inside the cylinder */
    if (  dis<0.05*ParticleObject[P_type].radius   )
    {
        /* get the normal to the surface */
    	Result.Normal = -pos/distToCenter;

    	Result.contact_point = Pos_P -Result.Normal*ParticleObject[P_type].radius;
    	Result.dis = dis;
        Result.is_collision = true;
    }



   return Result;
}
/*___________________________________________________________________________*/
