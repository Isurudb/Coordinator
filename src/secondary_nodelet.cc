/* primary_nodelet.cc

The primary coordinator, which derives from CoorindatorBase and adds methods for desired tests.
*/

#include "coordinator/secondary_nodelet.h"
#include "coordinator/secondary_dmpc_methods.hpp"
#include "ros/ros.h"
//#include "std_msgs/Int32.h"
/* ************************************************************************* */
void SecondaryNodelet::Initialize(ros::NodeHandle* nh) {
  /**
  * @brief This is called when the nodelet is loaded into the nodelet manager
  * 
  */
  // Create Multi-threaded NH
  MTNH = getMTNodeHandle();

  // Load Params
  load_params();
  std::cout << "[SECONDARY_COORD] Parameters loaded sucessfully." << std::endl;
  // publishers
  std::cout << "[SECONDARY_COORD] Assigning publishers ..." << std::endl;
  pub_flight_mode_ = nh->advertise<ff_msgs::FlightMode>(TOPIC_MOBILITY_FLIGHT_MODE_, 1, true);
 // std::cout << "[SECONDARY_COORD]     1. TOPIC_MOBILITY_FLIGHT_MODE_ : assigned..." << std::endl;  // FlightMode
  pub_status_ = nh->advertise<coordinator::StatusSecondary>(TOPIC_ASAP_STATUS, 5, true);
 // std::cout << "[SECONDARY_COORD]     2. TOPIC_ASAP_STATUS : assigned..." << std::endl;
  pub_ctl_=nh->advertise<ff_msgs::FamCommand>(TOPIC_GNC_CTL_CMD,1);
 // std::cout << "[SECONDARY_COORD]     3. TOPIC_GNC_CTL_CMD : assigned..." << std::endl;
  
  // subscribers
  std::cout << "[SECONDARY_COORD] Assigning subscribers ..." << std::endl;
  sub_flight_mode_= nh->subscribe<ff_msgs::FlightMode>(TOPIC_MOBILITY_FLIGHT_MODE_, 5,
    boost::bind(&SecondaryNodelet::flight_mode_callback, this, _1));  // flight mode getter
 // std::cout << "[SECONDARY_COORD]     1. TOPIC_MOBILITY_FLIGHT_MODE_ : assigned..." << std::endl;
  // sub_ekf_ = nh->subscribe<ff_msgs::EkfState>("gnc/ekf", 5,
  //   boost::bind(&PrimaryNodelet::ekf_callback, this, _1));;
  sub_test_number_ = nh->subscribe<coordinator::TestNumber>(TOPIC_ASAP_TEST_NUMBER, 5,
    boost::bind(&SecondaryNodelet::test_num_callback, this, _1));
 // std::cout << "[SECONDARY_COORD]     2. TOPIC_ASAP_TEST_NUMBER : assigned..." << std::endl;
  sub_flight_mode_= nh->subscribe<ff_msgs::FlightMode>(TOPIC_MOBILITY_FLIGHT_MODE_, 5,
    boost::bind(&SecondaryNodelet::flight_mode_callback, this, _1));  // flight mode setter
 // std::cout << "[SECONDARY_COORD]     3. TOPIC_MOBILITY_FLIGHT_MODE_ : assigned..." << std::endl;

   sub_ekf_ = nh->subscribe<ff_msgs::EkfState>(TOPIC_GNC_EKF_ , 3,
    boost::bind(&SecondaryNodelet::ekf_callback, this, _1)); //TOPIC_GNC_EKF_ TOPIC_GNC_EKF_ "/bumble/gnc/ekf"
 // std::cout << "[SECONDARY_COORD]     4. TOPIC_GNC_EKF_ : assigned..." << std::endl;

  sub_ekf_VL = nh->subscribe<ff_msgs::EkfState>(VIRTUAL_LEADER_TOPIC, 3,
    boost::bind(&SecondaryNodelet::VL_callback, this, _1));
  
 /*  sub_ekf_VL = nh->subscribe<ff_msgs::EkfState>(VIRTUAL_LEADER_TOPIC, 3,
    boost::bind(&SecondaryNodelet::VL_callback, this, _1));  */
 // std::cout << "[SECONDARY_COORD]     5. VIRTUAL_LEADER_TOPIC : assigned..." << std::endl;

  // services SERVICE_GNC_CTL_ENABLE
  serv_ctl_enable_ = nh->serviceClient<std_srvs::SetBool>(SERVICE_GNC_CTL_ENABLE_);
  std::cout << "[SECONDARY_COORD] SERVICE_GNC_CTL_ENABLE_ service assigned." << std::endl;
  // tracking points
  try{
    std::vector<double> param_data;
    nh->getParam("/asap/secondary/point_a_granite", param_data);
    POINT_A_GRANITE = Eigen::Matrix<double, 7, 1>(param_data.data());  // init from std::vector param
    nh->getParam("/asap/secondary/point_a_iss", param_data);
    POINT_A_ISS = Eigen::Matrix<double, 7, 1>(param_data.data());  
    nh->getParam("/asap/secondary/point_b_granite", param_data);
    POINT_B_GRANITE = Eigen::Matrix<double, 7, 1>(param_data.data());  
    nh->getParam("/asap/secondary/point_b_iss", param_data);
    POINT_B_ISS= Eigen::Matrix<double, 7, 1>(param_data.data());  
    nh->getParam("/asap/secondary/point_c_granite", param_data);
    POINT_C_GRANITE = Eigen::Matrix<double, 7, 1>(param_data.data());  
    nh->getParam("/asap/secondary/point_c_iss", param_data);
    POINT_C_ISS = Eigen::Matrix<double, 7, 1>(param_data.data()); 
    std::cout << "[SECONDARY_COORD] Matrix parameters loaded" << std::endl;
  }
  catch (const std::exception &exc) {
    std::cerr << exc.what() << std::endl;
    std::cout << "[SECONDARY_COORD]: Input parameters are invalid!" << std::endl;
  }

  // Pass control to Run method (activate timers and spin)
  std::cout << "[SECONDARY_COORD] Initialized." << std::endl;
  thread_.reset(new std::thread(&CoordinatorBase::Run, this, nh));
}



/* ************************************************************************** */
void SecondaryNodelet::get_status_msg(coordinator::StatusSecondary& msg){
  /**
   * @brief Fills the StatusPrimary message with the state of the private
   * variables. Published by a coordinator Timer.
   * Inputs: base_status_ and primary_status_
   * 
   */
  msg.stamp = ros::Time::now();
  msg.test_number = base_status_.test_number;
  msg.default_control = base_status_.default_control;
  msg.flight_mode = base_status_.flight_mode;
  msg.test_finished = base_status_.test_finished;
  msg.coord_ok = base_status_.coord_ok;
  msg.regulate_finished = base_status_.regulate_finished;

  msg.control_mode = secondary_status_.control_mode;
  msg.description = secondary_status_.description;
}

/* ************************************************************************** */
void SecondaryNodelet::load_params(){
  // Get sim and ground flags
  std::string sim_str, ground_str,coupled_str;
  ros::param::get("/asap/sim", sim_str);
  sim_ = !std::strcmp(sim_str.c_str(), "true"); // convert to bool
  ros::param::get("/asap/ground", ground_str);
  ground_ = !std::strcmp(ground_str.c_str(), "true");  // convert to bool, 1 if it's true

  ros::param::get("/asap/coupled", coupled_str);
  coupled = !std::strcmp(ground_str.c_str(), "true");  // convert to bool, 1 if it's true
  if(coupled)
  {
    std::cout << "[SECONDARY_COORD] Coupled mode is activated ................" << std::endl;
  }
  else
  {
    std::cout << "[SECONDARY_COORD] Coupled mode is deactivated ................" << std::endl;
  }
  


    // get the robot name

  std::string  secondary_robot_ns;
  ros::param::get("/asap/secondary_robot_name",secondary_robot_ns); //secondary_robot_name  

  // get the follower robot name

  std::string  robot_ns;
  ros::param::get("/asap/primary_robot_name",robot_ns); //primary_robot_name ie relatively secondary


  std::cout << "[SECONDARY_COORD] name spaced topics ................" << std::endl;
  
  TOPIC_ASAP_STATUS = robot_ns + TOPIC_ASAP_STATUS ;
  std::cout << TOPIC_ASAP_STATUS<< std::endl;

  TOPIC_ASAP_TEST_NUMBER = robot_ns + TOPIC_ASAP_TEST_NUMBER;
  std::cout << TOPIC_ASAP_TEST_NUMBER<< std::endl;

  TOPIC_GNC_CTL_CMD = robot_ns + TOPIC_GNC_CTL_CMD;
  std::cout << TOPIC_GNC_CTL_CMD<< std::endl;

  SERVICE_GNC_CTL_ENABLE_= robot_ns + SERVICE_GNC_CTL_ENABLE_;
  std::cout << SERVICE_GNC_CTL_ENABLE_<< std::endl;

  TOPIC_GNC_EKF_ = robot_ns + TOPIC_GNC_EKF_;
  std::cout << TOPIC_GNC_EKF_<< std::endl;

  VIRTUAL_LEADER_TOPIC = secondary_robot_ns + VIRTUAL_LEADER_TOPIC;
  std::cout << VIRTUAL_LEADER_TOPIC<< std::endl;

  TOPIC_MOBILITY_FLIGHT_MODE_ = robot_ns + TOPIC_MOBILITY_FLIGHT_MODE_;
   std::cout << TOPIC_MOBILITY_FLIGHT_MODE_<< std::endl;

  std::cout << "[SECONDARY_COORD] ................" << std::endl;
  
  // regulation
  ros::param::getCached("/asap/secondary/reg_time", reg_time_);
  ros::param::getCached("/asap/secondary/x_start", x0_(0));
  ros::param::getCached("/asap/secondary/y_start", x0_(1));
  ros::param::getCached("/asap/secondary/z_start", x0_(2));
  ros::param::getCached("/asap/secondary/qx_start", a0_(0));
  ros::param::getCached("/asap/secondary/qy_start", a0_(1));
  ros::param::getCached("/asap/secondary/qz_start", a0_(2));
  ros::param::getCached("/asap/secondary/qw_start", a0_(3));
  ros::param::getCached("/asap/secondary/pos_reg_thresh", pos_reg_thresh_);
  ros::param::getCached("/asap/secondary/vel_reg_thresh", vel_reg_thresh_);
  ros::param::getCached("/asap/secondary/att_reg_thresh", att_reg_thresh_);
  ros::param::getCached("/asap/secondary/omega_reg_thresh", omega_reg_thresh_);
  ros::param::getCached("/asap/secondary/roll", roll);
  ros::param::getCached("/asap/secondary/pitch", pitch);
  ros::param::getCached("/asap/secondary/yaw", yaw);
  

  ROS_INFO("[SECONDARY_COORD]....Goal position: X: %f Y: %f Z: %f ",x0_(0),x0_(1),x0_(2));
}
