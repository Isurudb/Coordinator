#pragma once
#include "coordinator/primary_nodelet.h"
//#include "std_msgs/String.h"



/************************************************************************/
void PrimaryNodelet::RunTest0(ros::NodeHandle *nh){
    int system_ret;
    std::string undock_command;
     undock_command = "rosrun dock dock_tool -undock";
    NODELET_INFO_STREAM("[PRIMARY_COORD]: Congratulations, you have passed quick checkout. " 
    "May your days be blessed with only warnings and no errors.");
    

    
    ros::Duration(5.0).sleep();
    /* ROS_INFO("Undocking the Astrobee ");
    NODELET_INFO_STREAM("Calling " << undock_command);
    system_ret = system(undock_command.c_str()); */

    if(system_ret != 0){
        NODELET_ERROR_STREAM("[PRIMARY/DMPC] Failed to Launch DMPC nodes.");
    }
    ROS_INFO("Initializing the position data....");
    position_ref.x = 11.3;
    position_ref.y = -5.4;
    position_ref.z =  4.4;

    //debug quaternion ambiguity
    /*q 0_x = attitude.x;
    q0_y = attitude.y;
    q0_z = attitude.z; */

   /*  z_nominal[0]=x0[0];
      z_nominal[1]=x0[1];
      z_nominal[2]=x0[2];
      z_nominal[3]=x0[3];
      z_nominal[4]=x0[4];
      z_nominal[5]=x0[5]; */
    initialzation=true;

    ROS_INFO("Position data successfully initialized!");

    ROS_INFO("End of initialization............. <<< Test 0 >>> ..............");

   // disable_default_ctl();
    //check_regulate();  // check regulation until satisfied
    //ROS_INFO("Setting up the publisher ");

    // pub_ctl_=nh->advertise<ff_msgs::FamCommand>(TOPIC_GNC_CTL_CMD,1);
    

   // RunTest2(nh);

    NODELET_DEBUG_STREAM("[PRIMARY COORD]: ...test complete!");
    base_status_.test_finished = false;
};


/************************************************************************/
void PrimaryNodelet::RunTest5(ros::NodeHandle *nh){
    /* RATTLE test: hand off control to RATTLE coordinator
    */
ROS_INFO("Runnig Test 1 -- Worst Estimate -- TRMPC");
Estimate_status="Worst";
RunTest0(nh);
    ROS_INFO("Runnig Test 1 now ");
primary_status_.control_mode = "regulate";
    ros::Duration(0.4).sleep(); // make sure controller gets the regulate settings before disabling default controller.
    // geometry_msgs::Vector3 torque;
    // double r=0, p=0, y=3.14159;  // Rotate the previous pose by 180* about Z

    //     q_ref.setRPY(r, p, y);
    //     tf2::convert(attitude,attitude_);
    //     q_ref_inv=q_ref;//.inverse();
    NODELET_DEBUG_STREAM("[PRIMARY COORD]: Disabling default controller...");
    disable_default_ctl();
 ROS_INFO("Initiating the Quaternion Feedback Controller");
    ros::Rate loop_rate(62.5);
 ROS_INFO("Setting up the publisher ");
 int t=0;
    while(ros::ok()){
        
        // Rotation matrix
        // q_e= q_ref_inv*attitude_;  // Calculate the new orientation
        // q_e.normalize();
        float R_11 = 2*(attitude.x*attitude.x + attitude.w*attitude.w)-1;
        float R_12 = 2*(attitude.x*attitude.y - attitude.w*attitude.z);
        float R_13 = 2*(attitude.x*attitude.z + attitude.w*attitude.y); 
        float R_21 = 2*(attitude.x*attitude.y + attitude.w*attitude.z);
        float R_22 = 2*(attitude.y*attitude.y + attitude.w*attitude.w)-1;
        float R_23 = 2*(attitude.y*attitude.z - attitude.w*attitude.x);
        float R_31 = 2*(attitude.x*attitude.z - attitude.w*attitude.y); 
        float R_32 = 2*(attitude.y*attitude.z + attitude.w*attitude.x);
        float R_33 = 2*(attitude.z*attitude.z + attitude.w*attitude.w)-1;

        float u_x = kN[0];//Fx;//arg_fx;//-13.5*velocity_.x -0.85*position_error.x;
        float u_y = kN[1];//]Fy;//arg_fy;//-13.5*velocity_.y -0.85*position_error.y;
        float u_z = kN[2];//Fz;//arg_fz;//-1.0*velocity_.z -0.1*position_error.z;

        

        

        float f_max=0.6;
        float f_min=-0.6;

        if (u_x>f_max){
            u_x=f_max;
        }

         if (u_y>f_max){
            u_y=f_max;
        }

         if (u_z>f_max){
            u_z=f_max;
        }

        if (u_x<f_min){
            u_x=f_min;
        }

        if (u_y<f_min){
            u_y=f_min;
        }
        if (u_z<f_min){
            u_z=f_min;
        }


            
            ctl_input.force.x = u_x*R_11 + u_y*R_21 + u_z*R_31;//-0.05*velocity_.x +0.005*position_error.x;
            ctl_input.force.y = u_x*R_12 + u_y*R_22 + u_z*R_32;//-0.05*velocity_.y -0.005*position_error.y;
            ctl_input.force.z = u_x*R_13 + u_y*R_23 + u_z*R_33;//-0.05*velocity_.z +0.005*position_error.z;
            
        
        float ex =position_error.x;
        float ey =position_error.y;
        float ez =position_error.z;
       
         
         if(t==60){ 
         if(sqrt(ex*ex+ey*ey+ez*ez)<0.1)
            {
            ROS_INFO(" ---------------------------------------Goal Position arrived--------------------------------");
        }
        else{
        ROS_INFO(" Deploying TRMPC for transverse motion  ex: [%f]  ey: [%f] ez: [%f]\n Fx: [%f] Fy: [%f] Fz: [%f] ",
            position_error.x, position_error.y, position_error.z,ctl_input.force.x,ctl_input.force.y,ctl_input.force.z);
           
        ROS_INFO("qx: [%f]  qy: [%f] qz: [%f] qw: [%f]", q_e.getX()*q_e.getX(),q_e.getY()*q_e.getY(),q_e.getZ()*q_e.getZ(),q_e.getW());
        }
         t=0;
        
         }

        gnc_setpoint.header.frame_id="body";
        gnc_setpoint.header.stamp=ros::Time::now();
        gnc_setpoint.wrench=ctl_input;
        gnc_setpoint.status=3;
        gnc_setpoint.control_mode=2;

        
        ctl_input.torque.x=arg_tau_x;//-0.02*q_e.getX()-0.2*omega.x;
        ctl_input.torque.y=arg_tau_y;//-0.02*q_e.getY()-0.2*omega.y;
        ctl_input.torque.z=arg_tau_z;//-0.02*q_e.getZ()-0.2*omega.z;
  
        

        pub_ctl_.publish(gnc_setpoint);
        t+=1;
 
        loop_rate.sleep();

        ros::spinOnce();



    };
    //****************************************************************************************************
    NODELET_DEBUG_STREAM("[PRIMARY COORD]: ...test complete!");
    base_status_.test_finished = true;
}

    
// ***************************************************************************************************
    // Additional test commands go here
    // Test commands can be anything you want! Talk to as many custom nodes as desired.
//******************************************************************************************************
void PrimaryNodelet::RunTest6(ros::NodeHandle *nh){
    /* RATTLE test: hand off control to RATTLE coordinator
    */
ROS_INFO("Test 2 -- Worst Estimate -- MPC");
Estimate_status="Worst";
RunTest0(nh);
primary_status_.control_mode = "regulate";
    ros::Duration(0.4).sleep(); // make sure controller gets the regulate settings before disabling default controller.
    NODELET_DEBUG_STREAM("[PRIMARY COORD]: Disabling default controller...");
    disable_default_ctl();
 ROS_INFO("Initiating the Quaternion Feedback Controller");
    ros::Rate loop_rate(62.5);
 ROS_INFO("Setting up the publisher ");
 int t=0;
    while(ros::ok()){
        // Rotation matrix
        // q_e= q_ref_inv*attitude_;  // Calculate the new orientation
        // q_e.normalize();
        float R_11 = 2*(attitude.x*attitude.x + attitude.w*attitude.w)-1;
        float R_12 = 2*(attitude.x*attitude.y - attitude.w*attitude.z);
        float R_13 = 2*(attitude.x*attitude.z + attitude.w*attitude.y); 
        float R_21 = 2*(attitude.x*attitude.y + attitude.w*attitude.z);
        float R_22 = 2*(attitude.y*attitude.y + attitude.w*attitude.w)-1;
        float R_23 = 2*(attitude.y*attitude.z - attitude.w*attitude.x);
        float R_31 = 2*(attitude.x*attitude.z - attitude.w*attitude.y); 
        float R_32 = 2*(attitude.y*attitude.z + attitude.w*attitude.x);
        float R_33 = 2*(attitude.z*attitude.z + attitude.w*attitude.w)-1;

        float u_x = X_QP[0];//Fx;//arg_fx;//-13.5*velocity_.x -0.85*position_error.x;
        float u_y = X_QP[1];//]Fy;//arg_fy;//-13.5*velocity_.y -0.85*position_error.y;
        float u_z = X_QP[2];//Fz;//arg_fz;//-1.0*velocity_.z -0.1*position_error.z;

               
        ctl_input.force.x = u_x*R_11 + u_y*R_21 + u_z*R_31;//-0.05*velocity_.x +0.005*position_error.x;
        ctl_input.force.y = u_x*R_12 + u_y*R_22 + u_z*R_32;//-0.05*velocity_.y -0.005*position_error.y;
        ctl_input.force.z = u_x*R_13 + u_y*R_23 + u_z*R_33;//-0.05*velocity_.z +0.005*position_error.z;
            
        float ex =position_error.x;
        float ey =position_error.y;
        float ez =position_error.z;
        
       
         
         if(t==60){ 
         if(sqrt(ex*ex+ey*ey+ez*ez)<0.1)
            {
            ROS_INFO(" ---------------------------------------Goal Position arrived--------------------------------");
        }
         else{  
        ROS_INFO(" Deploying MPC for transverse motion  ex: [%f]  ey: [%f] ez: [%f]\n Fx: [%f] Fy: [%f] Fz: [%f] ",
            position_error.x, position_error.y, position_error.z,ctl_input.force.x,ctl_input.force.y,ctl_input.force.z);
           
        ROS_INFO("qx: [%f]  qy: [%f] qz: [%f] qw: [%f]", q_e.getX()*q_e.getX(),q_e.getY()*q_e.getY(),q_e.getZ()*q_e.getZ(),q_e.getW());
         }
         t=0;
         }

        gnc_setpoint.header.frame_id="body";
        gnc_setpoint.header.stamp=ros::Time::now();
        gnc_setpoint.wrench=ctl_input;
        gnc_setpoint.status=3;
        gnc_setpoint.control_mode=2;

        
        ctl_input.torque.x=arg_tau_x;//-0.02*q_e.getX()-0.2*omega.x;
        ctl_input.torque.y=arg_tau_y;//-0.02*q_e.getY()-0.2*omega.y;
        ctl_input.torque.z=arg_tau_z;//-0.02*q_e.getZ()-0.2*omega.z;
  
        

        pub_ctl_.publish(gnc_setpoint);
        

        t+=1;

        
        loop_rate.sleep();

        ros::spinOnce();



    };
    //****************************************************************************************************
    NODELET_DEBUG_STREAM("[PRIMARY COORD]: ...test complete!");
    base_status_.test_finished = true;
}
/**************************************************************************************************/
//Test 3
/**************************************************************************************************/
void PrimaryNodelet::RunTest3(ros::NodeHandle *nh){
    /*  Good estimate TRMPC
    */
ROS_INFO("Runnig test 3 ---- TRMPC -- Good Estimate ");
RunTest0(nh);
    ROS_INFO("Runnig Test 1 now ");
primary_status_.control_mode = "regulate";
    ros::Duration(0.4).sleep(); // make sure controller gets the regulate settings before disabling default controller.
   
    NODELET_DEBUG_STREAM("[PRIMARY COORD]: Disabling default controller...");
    disable_default_ctl();
 ROS_INFO("Initiating the Quaternion Feedback Controller");
    ros::Rate loop_rate(62.5);
 ROS_INFO("Setting up the publisher ");
 int t=0;
    while(ros::ok()){
        
         
 // Rotation matrix
        // q_e= q_ref_inv*attitude_;  // Calculate the new orientation
        // q_e.normalize();
        float R_11 = 2*(attitude.x*attitude.x + attitude.w*attitude.w)-1;
        float R_12 = 2*(attitude.x*attitude.y - attitude.w*attitude.z);
        float R_13 = 2*(attitude.x*attitude.z + attitude.w*attitude.y); 
        float R_21 = 2*(attitude.x*attitude.y + attitude.w*attitude.z);
        float R_22 = 2*(attitude.y*attitude.y + attitude.w*attitude.w)-1;
        float R_23 = 2*(attitude.y*attitude.z - attitude.w*attitude.x);
        float R_31 = 2*(attitude.x*attitude.z - attitude.w*attitude.y); 
        float R_32 = 2*(attitude.y*attitude.z + attitude.w*attitude.x);
        float R_33 = 2*(attitude.z*attitude.z + attitude.w*attitude.w)-1;

        float u_x = kN[0];//Fx;//arg_fx;//-13.5*velocity_.x -0.85*position_error.x;
        float u_y = kN[1];//]Fy;//arg_fy;//-13.5*velocity_.y -0.85*position_error.y;
        float u_z = kN[2];//Fz;//arg_fz;//-1.0*velocity_.z -0.1*position_error.z;

        float f_max=0.6;
        float f_min=-0.6;

        if (u_x>f_max){
            u_x=f_max;
        }

         if (u_y>f_max){
            u_y=f_max;
        }

         if (u_z>f_max){
            u_z=f_max;
        }

        if (u_x<f_min){
            u_x=f_min;
        }

        if (u_y<f_min){
            u_y=f_min;
        }
        if (u_z<f_min){
            u_z=f_min;
        }


            
            ctl_input.force.x = u_x*R_11 + u_y*R_21 + u_z*R_31;//-0.05*velocity_.x +0.005*position_error.x;
            ctl_input.force.y = u_x*R_12 + u_y*R_22 + u_z*R_32;//-0.05*velocity_.y -0.005*position_error.y;
            ctl_input.force.z = u_x*R_13 + u_y*R_23 + u_z*R_33;//-0.05*velocity_.z +0.005*position_error.z;
            
        
        float ex =position_error.x;
        float ey =position_error.y;
        float ez =position_error.z;
       
         
         if(t==60){ 
        if(sqrt(ex*ex+ey*ey+ez*ez)<0.1)
            {
            ROS_INFO(" ---------------------------------------Goal Position arrived--------------------------------");
        }
        else{    
        ROS_INFO(" Deploying TRMPC for transverse motion  ex: [%f]  ey: [%f] ez: [%f]\n Fx: [%f] Fy: [%f] Fz: [%f] ",
            position_error.x, position_error.y, position_error.z,ctl_input.force.x,ctl_input.force.y,ctl_input.force.z);
           
        ROS_INFO("qx: [%f]  qy: [%f] qz: [%f] qw: [%f]", q_e.getX()*q_e.getX(),q_e.getY()*q_e.getY(),q_e.getZ()*q_e.getZ(),q_e.getW());
        }
         t=0;
         }

        gnc_setpoint.header.frame_id="body";
        gnc_setpoint.header.stamp=ros::Time::now();
        gnc_setpoint.wrench=ctl_input;
        gnc_setpoint.status=3;
        gnc_setpoint.control_mode=2;

        
        ctl_input.torque.x=arg_tau_x;//-0.02*q_e.getX()-0.2*omega.x;
        ctl_input.torque.y=arg_tau_y;//-0.02*q_e.getY()-0.2*omega.y;
        ctl_input.torque.z=arg_tau_z;//-0.02*q_e.getZ()-0.2*omega.z;
  
        

        pub_ctl_.publish(gnc_setpoint);
        t+=1;
   

        loop_rate.sleep();

        ros::spinOnce();

    };
    //****************************************************************************************************
    NODELET_DEBUG_STREAM("[PRIMARY COORD]: ...test complete!");
    base_status_.test_finished = true;
}


void PrimaryNodelet::RunTest4(ros::NodeHandle *nh){
    /*  Good estimate MPC
    */
 ROS_INFO("Runnig test 4 ---- MPC -- Good Estimate ");
 Estimate_status=="Good";
RunTest0(nh);
    ROS_INFO("Runnig Test 1 now ");
primary_status_.control_mode = "regulate";
    ros::Duration(0.4).sleep(); // make sure controller gets the regulate settings before disabling default controller.
   
    NODELET_DEBUG_STREAM("[PRIMARY COORD]: Disabling default controller...");
    disable_default_ctl();
 ROS_INFO("Initiating the Quaternion Feedback Controller");
    ros::Rate loop_rate(62.5);
 ROS_INFO("Setting up the publisher ");
 int t=0;
    while(ros::ok()){
        
        
        // Rotation matrix
        // q_e= q_ref_inv*attitude_;  // Calculate the new orientation
        // q_e.normalize();
        float R_11 = 2*(attitude.x*attitude.x + attitude.w*attitude.w)-1;
        float R_12 = 2*(attitude.x*attitude.y - attitude.w*attitude.z);
        float R_13 = 2*(attitude.x*attitude.z + attitude.w*attitude.y); 
        float R_21 = 2*(attitude.x*attitude.y + attitude.w*attitude.z);
        float R_22 = 2*(attitude.y*attitude.y + attitude.w*attitude.w)-1;
        float R_23 = 2*(attitude.y*attitude.z - attitude.w*attitude.x);
        float R_31 = 2*(attitude.x*attitude.z - attitude.w*attitude.y); 
        float R_32 = 2*(attitude.y*attitude.z + attitude.w*attitude.x);
        float R_33 = 2*(attitude.z*attitude.z + attitude.w*attitude.w)-1;

        float u_x = X_QP[0];//Fx;//arg_fx;//-13.5*velocity_.x -0.85*position_error.x;
        float u_y = X_QP[1];//]Fy;//arg_fy;//-13.5*velocity_.y -0.85*position_error.y;
        float u_z = X_QP[2];//Fz;//arg_fz;//-1.0*velocity_.z -0.1*position_error.z;

               
        ctl_input.force.x = u_x*R_11 + u_y*R_21 + u_z*R_31;//-0.05*velocity_.x +0.005*position_error.x;
        ctl_input.force.y = u_x*R_12 + u_y*R_22 + u_z*R_32;//-0.05*velocity_.y -0.005*position_error.y;
        ctl_input.force.z = u_x*R_13 + u_y*R_23 + u_z*R_33;//-0.05*velocity_.z +0.005*position_error.z;
            
        float ex =position_error.x;
        float ey =position_error.y;
        float ez =position_error.z;
        
       
         
         if(t==60){ 
             if(sqrt(ex*ex+ey*ey+ez*ez)<0.1)
            {
            ROS_INFO(" ---------------------------------------Goal Position arrived--------------------------------");
        }
        else{
        ROS_INFO(" Deploying MPC for transverse motion  ex: [%f]  ey: [%f] ez: [%f]\n Fx: [%f] Fy: [%f] Fz: [%f] ",
            position_error.x, position_error.y, position_error.z,ctl_input.force.x,ctl_input.force.y,ctl_input.force.z);
           
        ROS_INFO("qx: [%f]  qy: [%f] qz: [%f] qw: [%f]", q_e.getX()*q_e.getX(),q_e.getY()*q_e.getY(),q_e.getZ()*q_e.getZ(),q_e.getW());
        }
         t=0;
         }

        gnc_setpoint.header.frame_id="body";
        gnc_setpoint.header.stamp=ros::Time::now();
        gnc_setpoint.wrench=ctl_input;
        gnc_setpoint.status=3;
        gnc_setpoint.control_mode=2;

        
        ctl_input.torque.x=arg_tau_x;//-0.02*q_e.getX()-0.2*omega.x;
        ctl_input.torque.y=arg_tau_y;//-0.02*q_e.getY()-0.2*omega.y;
        ctl_input.torque.z=arg_tau_z;//-0.02*q_e.getZ()-0.2*omega.z;
  
        

        pub_ctl_.publish(gnc_setpoint);
        

        t+=1;

        loop_rate.sleep();

        ros::spinOnce();

    };
    //****************************************************************************************************
    NODELET_DEBUG_STREAM("[PRIMARY COORD]: ...test complete!");
    base_status_.test_finished = true;
}

void PrimaryNodelet::RunTest1(ros::NodeHandle *nh){
    /*  Best estimate TRMPC
    */
 ROS_INFO("Test 5  TRMPC -- Best Estimate -- status -- <<< Initiated >>> ");
RunTest0(nh);
    ROS_INFO("Runnig Test 1 now ");
primary_status_.control_mode = "regulate";
    ros::Duration(0.4).sleep(); // make sure controller gets the regulate settings before disabling default controller.
   
    NODELET_DEBUG_STREAM("[PRIMARY COORD]: Disabling default controller...");
    disable_default_ctl();
    Estimate_status=="Best";
 ROS_INFO("Initiating the Quaternion Feedback Controller");
    ros::Rate loop_rate(62.5);
 ROS_INFO("Setting up the publisher ");
 int t=0;
    while(ros::ok()){
        
        
    // Rotation matrix
        // q_e= q_ref_inv*attitude_;  // Calculate the new orientation
        // q_e.normalize();
        float R_11 = 2*(attitude.x*attitude.x + attitude.w*attitude.w)-1;
        float R_12 = 2*(attitude.x*attitude.y - attitude.w*attitude.z);
        float R_13 = 2*(attitude.x*attitude.z + attitude.w*attitude.y); 
        float R_21 = 2*(attitude.x*attitude.y + attitude.w*attitude.z);
        float R_22 = 2*(attitude.y*attitude.y + attitude.w*attitude.w)-1;
        float R_23 = 2*(attitude.y*attitude.z - attitude.w*attitude.x);
        float R_31 = 2*(attitude.x*attitude.z - attitude.w*attitude.y); 
        float R_32 = 2*(attitude.y*attitude.z + attitude.w*attitude.x);
        float R_33 = 2*(attitude.z*attitude.z + attitude.w*attitude.w)-1;

        float u_x = kN[0];//Fx;//arg_fx;//-13.5*velocity_.x -0.85*position_error.x;
        float u_y = kN[1];//]Fy;//arg_fy;//-13.5*velocity_.y -0.85*position_error.y;
        float u_z = kN[2];//Fz;//arg_fz;//-1.0*velocity_.z -0.1*position_error.z;

        float f_max=0.6;
        float f_min=-0.6;

        if (u_x>f_max){
            u_x=f_max;
        }

         if (u_y>f_max){
            u_y=f_max;
        }

         if (u_z>f_max){
            u_z=f_max;
        }

        if (u_x<f_min){
            u_x=f_min;
        }

        if (u_y<f_min){
            u_y=f_min;
        }
        if (u_z<f_min){
            u_z=f_min;
        }


            
            ctl_input.force.x = u_x*R_11 + u_y*R_21 + u_z*R_31;//-0.05*velocity_.x +0.005*position_error.x;
            ctl_input.force.y = u_x*R_12 + u_y*R_22 + u_z*R_32;//-0.05*velocity_.y -0.005*position_error.y;
            ctl_input.force.z = u_x*R_13 + u_y*R_23 + u_z*R_33;//-0.05*velocity_.z +0.005*position_error.z;
            
        
        
       
         
        float ex =position_error.x;
        float ey =position_error.y;
        float ez =position_error.z;
       
         
         if(t==60){ 
            if(sqrt(ex*ex+ey*ey+ez*ez)<0.1)
            {
            ROS_INFO(" ---------------------------------------Goal Position arrived--------------------------------");
        }
        else{
        ROS_INFO(" Deploying TRMPC for transverse motion  ex: [%f]  ey: [%f] ez: [%f]\n Fx: [%f] Fy: [%f] Fz: [%f] ",
            position_error.x, position_error.y, position_error.z,ctl_input.force.x,ctl_input.force.y,ctl_input.force.z);
           
        ROS_INFO("qx: [%f]  qy: [%f] qz: [%f] qw: [%f]", q_e.getX()*q_e.getX(),q_e.getY()*q_e.getY(),q_e.getZ()*q_e.getZ(),q_e.getW());
        }
         t=0;
        
         }

        gnc_setpoint.header.frame_id="body";
        gnc_setpoint.header.stamp=ros::Time::now();
        gnc_setpoint.wrench=ctl_input;
        gnc_setpoint.status=3;
        gnc_setpoint.control_mode=2;

        
        ctl_input.torque.x=arg_tau_x;//-0.02*q_e.getX()-0.2*omega.x;
        ctl_input.torque.y=arg_tau_y;//-0.02*q_e.getY()-0.2*omega.y;
        ctl_input.torque.z=arg_tau_z;//-0.02*q_e.getZ()-0.2*omega.z;
  
        

        pub_ctl_.publish(gnc_setpoint);
        loop_rate.sleep();

        ros::spinOnce();

        t+=1;

    };
    //****************************************************************************************************
    NODELET_DEBUG_STREAM("[PRIMARY COORD]: ...test complete!");
    base_status_.test_finished = true;
}


void PrimaryNodelet::RunTest2(ros::NodeHandle *nh){
    /*  Best estimate MPC
    */
ROS_INFO("Test 6  MPC -- Best Estimate -- status -- <<< Initiated >>>");
RunTest0(nh);
    ROS_INFO("Runnig Test 1 now ");
primary_status_.control_mode = "regulate";
    ros::Duration(0.4).sleep(); // make sure controller gets the regulate settings before disabling default controller.
   
    NODELET_DEBUG_STREAM("[PRIMARY COORD]: Disabling default controller...");
    disable_default_ctl();
 Estimate_status=="Best";
 ROS_INFO("Initiating the Quaternion Feedback Controller");
    ros::Rate loop_rate(62.5);
 ROS_INFO("Setting up the publisher ");
 int t=0;
    while(ros::ok()){
        // Rotation matrix
        // q_e= q_ref_inv*attitude_;  // Calculate the new orientation
        // q_e.normalize();
        float R_11 = 2*(attitude.x*attitude.x + attitude.w*attitude.w)-1;
        float R_12 = 2*(attitude.x*attitude.y - attitude.w*attitude.z);
        float R_13 = 2*(attitude.x*attitude.z + attitude.w*attitude.y); 
        float R_21 = 2*(attitude.x*attitude.y + attitude.w*attitude.z);
        float R_22 = 2*(attitude.y*attitude.y + attitude.w*attitude.w)-1;
        float R_23 = 2*(attitude.y*attitude.z - attitude.w*attitude.x);
        float R_31 = 2*(attitude.x*attitude.z - attitude.w*attitude.y); 
        float R_32 = 2*(attitude.y*attitude.z + attitude.w*attitude.x);
        float R_33 = 2*(attitude.z*attitude.z + attitude.w*attitude.w)-1;

        float u_x = X_QP[0];//Fx;//arg_fx;//-13.5*velocity_.x -0.85*position_error.x;
        float u_y = X_QP[1];//]Fy;//arg_fy;//-13.5*velocity_.y -0.85*position_error.y;
        float u_z = X_QP[2];//Fz;//arg_fz;//-1.0*velocity_.z -0.1*position_error.z;

               
        ctl_input.force.x = u_x*R_11 + u_y*R_21 + u_z*R_31;//-0.05*velocity_.x +0.005*position_error.x;
        ctl_input.force.y = u_x*R_12 + u_y*R_22 + u_z*R_32;//-0.05*velocity_.y -0.005*position_error.y;
        ctl_input.force.z = u_x*R_13 + u_y*R_23 + u_z*R_33;//-0.05*velocity_.z +0.005*position_error.z;
            
        
        float ex =position_error.x;
        float ey =position_error.y;
        float ez =position_error.z;
       
        
       
         
         if(t==60){ 
             if(sqrt(ex*ex+ey*ey+ez*ez)<0.1)
            {
            ROS_INFO(" ---------------------------------------Goal Position arrived--------------------------------");
        }
        else{
        ROS_INFO(" Deploying MPC for transverse motion  ex: [%f]  ey: [%f] ez: [%f]\n Fx: [%f] Fy: [%f] Fz: [%f] ",
            position_error.x, position_error.y, position_error.z,ctl_input.force.x,ctl_input.force.y,ctl_input.force.z);
           
        ROS_INFO("qx: [%f]  qy: [%f] qz: [%f] qw: [%f]", q_e.getX()*q_e.getX(),q_e.getY()*q_e.getY(),q_e.getZ()*q_e.getZ(),q_e.getW());
        }
         t=0;
         }

        gnc_setpoint.header.frame_id="body";
        gnc_setpoint.header.stamp=ros::Time::now();
        gnc_setpoint.wrench=ctl_input;
        gnc_setpoint.status=3;
        gnc_setpoint.control_mode=2;

        
        ctl_input.torque.x=arg_tau_x;//-0.02*q_e.getX()-0.2*omega.x;
        ctl_input.torque.y=arg_tau_y;//-0.02*q_e.getY()-0.2*omega.y;
        ctl_input.torque.z=arg_tau_z;//-0.02*q_e.getZ()-0.2*omega.z;
  
        

        pub_ctl_.publish(gnc_setpoint);
        

        t+=1;

        loop_rate.sleep();

        ros::spinOnce();

    };
    //****************************************************************************************************
    NODELET_DEBUG_STREAM("[PRIMARY COORD]: ...test complete!");
    base_status_.test_finished = true;
}



/* ************************************************************************** */
void PrimaryNodelet::control_mode_callback(const std_msgs::String::ConstPtr msg) {
    /* Update control_mode form an external node.
    */
    primary_status_.control_mode = msg->data;
}
