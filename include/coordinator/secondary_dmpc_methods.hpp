#pragma once
#include "coordinator/secondary_nodelet.h"
//#include "std_msgs/String.h"



/************************************************************************/
void SecondaryNodelet::RunTest0(ros::NodeHandle *nh){
    int system_ret;
    std::string undock_command;
     undock_command = "rosrun executive teleop_tool -move -att '1.5 0 0 1' -ns 'bumble'";//"rosrun dock dock_tool -undock -ns 'queen'";
    NODELET_INFO_STREAM("[SECONDARY_COORD]: Congratulations, you have passed quick checkout. " 
    "May your days be blessed with only warnings and no errors.");
    

    
    ros::Duration(5.0).sleep();
    ROS_INFO("Undocking the Astrobee ");
    NODELET_INFO_STREAM("Calling " << undock_command);
    system_ret = system(undock_command.c_str());
 
    if(system_ret != 0){
        NODELET_ERROR_STREAM("[SECONDARY/DMPC] Failed to Launch DMPC nodes.");
    }
    ROS_INFO("Rotate the previous pose by 180* about Z ....");

    //disable_default_ctl();
    //check_regulate();  // check regulation until satisfied
    //ROS_INFO("Setting up the publisher ");

    // pub_ctl_=nh->advertise<ff_msgs::FamCommand>(TOPIC_GNC_CTL_CMD,1);
    
     robot = "Secondary";
    //RunTest1(nh);
    position_ref.x = position_.x + x0_(0);
    position_ref.y = position_.y + x0_(1);
    position_ref.z = position_.z; + x0_(2);

     //double L0 L;
    L0= sqrt( (pos_ref2.x - position_.x)*(pos_ref2.x - position_.x) + (pos_ref2.y - position_.y)*(pos_ref2.y - position_.y) +  (pos_ref2.z - position_.z)*(pos_ref2.z - position_.z) );
    L=L0;
    for (int i = 0; i < 50; i++) 
    {
        L0 = sqrt( (pos_ref2.x - position_.x)*(pos_ref2.x - position_.x) + (pos_ref2.y - position_.y)*(pos_ref2.y - position_.y) +  (pos_ref2.z - position_.z)*(pos_ref2.z - position_.z) );
        L=L+0.01*(L-L0);
        //ROS_INFO("Esitmated L is L0: %f  L: %f",L0,L); 
    
    }



     //run_test_0=true;
    NODELET_DEBUG_STREAM("[PRIMARY COORD]: ...test complete!");
    //ROS_INFO("New Goal positions are x: %f y: %f z: %f",position_ref.x,position_ref.y,position_ref.z); 
      ROS_INFO("Esitmated L is : %f ",L); 
    base_status_.test_finished = false;
};


/************************************************************************/
void SecondaryNodelet::RunTest1(ros::NodeHandle *nh){
    /* RATTLE test: hand off control to RATTLE coordinator
    */
    ROS_INFO("Test 2 -- Worst Estimate -- MPC");
Estimate_status="Worst";
RunTest0(nh);
secondary_status_.control_mode = "regulate";
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
            
        float ex =position_error_2.x;
        float ey =position_error_2.y;
        float ez =position_error_2.z;
        
       
         
         if(t==60){ 
         if(sqrt(ex*ex+ey*ey+ez*ez)<0.01)
            {
            ROS_INFO(" ---------------------------------------Goal Position arrived--------------------------------");
            ROS_INFO(" Deploying MPC for transverse motion  ref_x: [%f]  ref_y: [%f] ref_z: [%f]\n pose_x: [%f] pose_y: [%f] Pose_z: [%f] \n ",
            pos_ref2.x, pos_ref2.y, pos_ref2.z,position_.x,position_.y,position_.z);

            ROS_INFO(" Deploying MPC for transverse motion  ex: [%f]  ey: [%f] ez: [%f]\n Fx: [%f] Fy: [%f] Fz: [%f] \n",
            position_error_2.x, position_error_2.y, position_error_2.z,x0[0],x0[1],x0[2]);
          
            // ROS_INFO(" Deploying MPC for transverse motion  ex: [%f]  ey: [%f] ez: [%f]\n Fx: [%f] Fy: [%f] Fz: [%f] \n",
            // position_error_2.x, position_error_2.y, position_error_2.z,ctl_input.force.x,ctl_input.force.y,ctl_input.force.z);
          
          
        }
         else{  
        ROS_INFO(" Deploying MPC for transverse motion  ex: [%f]  ey: [%f] ez: [%f]\n Fx: [%f] Fy: [%f] Fz: [%f] ",
            pos_ref2.x, pos_ref2.y, pos_ref2.z,ctl_input.force.x,ctl_input.force.y,ctl_input.force.z);
           
        ROS_INFO("qx: [%f]  qy: [%f] qz: [%f] qw: [%f]", q_e.getX()*q_e.getX(),q_e.getY()*q_e.getY(),q_e.getZ()*q_e.getZ(),q_e.getW());
         }
         t=0;
         }
            
           
        if ( (arg_tau_x>0.01))
        {
            ctl_input.torque.x=0.01;
        }
        else if (arg_tau_x<-0.01)
        {
            ctl_input.torque.x=-0.01;

        }

        else
        {

            ctl_input.torque.x=arg_tau_x;
        }
       // -------------------------------------
        
        if ( (arg_tau_y>0.01))
        {
            ctl_input.torque.y=0.01;
        }
        else if (arg_tau_y<-0.01)
        {
            ctl_input.torque.y=-0.01;

        }

        else
        {

            ctl_input.torque.y=arg_tau_y;
        }

        //-----------------------------

        if ( (arg_tau_z>0.01))
        {
            ctl_input.torque.z=0.01;
        }
        else if (arg_tau_z<-0.01)
        {
            ctl_input.torque.z=-0.01;

        }

        else
        {

            ctl_input.torque.z=arg_tau_z;
        }
        
        
        gnc_setpoint.header.frame_id="body";
        gnc_setpoint.header.stamp=ros::Time::now();
        gnc_setpoint.wrench=ctl_input;
        gnc_setpoint.status=3;
        gnc_setpoint.control_mode=2;
        
        
        

        
        //ctl_input.torque.x=arg_tau_x;//-0.02*q_e.getX()-0.2*omega.x;
        //ctl_input.torque.y=arg_tau_y;//-0.02*q_e.getY()-0.2*omega.y;
        //ctl_input.torque.z=arg_tau_z;//-0.02*q_e.getZ()-0.2*omega.z;
  
        

        pub_ctl_.publish(gnc_setpoint);
        //VL_status.publish(mpc_pred);
        

        t+=1;

        
        loop_rate.sleep();

        ros::spinOnce();



    };
    //****************************************************************************************************
    NODELET_DEBUG_STREAM("[PRIMARY COORD]: ...test complete!");
    base_status_.test_finished = true;
}


/* ************************************************************************** */
void SecondaryNodelet::control_mode_callback(const std_msgs::String::ConstPtr msg) {
    /* Update control_mode form an external node.
    */
    secondary_status_.control_mode = msg->data;
}