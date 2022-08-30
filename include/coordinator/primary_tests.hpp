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
    ROS_INFO("Undocking the Astrobee ");
    NODELET_INFO_STREAM("Calling " << undock_command);
    system_ret = system(undock_command.c_str());

    if(system_ret != 0){
        NODELET_ERROR_STREAM("[PRIMARY/DMPC] Failed to Launch DMPC nodes.");
    }
    ROS_INFO("Spinnig the Astrobee ....");

    //disable_default_ctl();
    //check_regulate();  // check regulation until satisfied
    //ROS_INFO("Setting up the publisher ");

   // pub_ctl_=nh->advertise<ff_msgs::FamCommand>(TOPIC_GNC_CTL_CMD,1);

    RunTest1(nh);

    NODELET_DEBUG_STREAM("[PRIMARY COORD]: ...test complete!");
    base_status_.test_finished = true;
};


/************************************************************************/
void PrimaryNodelet::RunTest1(ros::NodeHandle *nh){
    /* RATTLE test: hand off control to RATTLE coordinator
    */

    ROS_INFO("Runnig Test 1 now ");
primary_status_.control_mode = "regulate";
    ros::Duration(0.4).sleep(); // make sure controller gets the regulate settings before disabling default controller.

    NODELET_DEBUG_STREAM("[PRIMARY COORD]: Disabling default controller...");
    disable_default_ctl();

    ros::Rate loop_rate(62.5);
 ROS_INFO("Setting up the publisher ");
    while(ros::ok()){

        gnc_setpoint.header.frame_id="body";
        gnc_setpoint.header.stamp=ros::Time::now();
        gnc_setpoint.wrench=ctl_input;
        gnc_setpoint.status=3;
        gnc_setpoint.control_mode=2;

        ctl_input.force.x=0.0;
        ctl_input.force.y=0.0;
        ctl_input.force.z=0.0;
        ctl_input.torque.x=0.0;
        ctl_input.torque.y=0.0;
        ctl_input.torque.z=0.1;

        pub_ctl_.publish(gnc_setpoint);
        loop_rate.sleep();

        ros::spinOnce();


    };
    

    // Additional test commands go here
    // Test commands can be anything you want! Talk to as many custom nodes as desired.

    NODELET_DEBUG_STREAM("[PRIMARY COORD]: ...test complete!");
    base_status_.test_finished = true;
}


/* ************************************************************************** */
void PrimaryNodelet::control_mode_callback(const std_msgs::String::ConstPtr msg) {
    /* Update control_mode form an external node.
    */
    primary_status_.control_mode = msg->data;
}
