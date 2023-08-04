/* primary_nodelet.cc

The primary coordinator, which derives from CoorindatorBase and adds methods for desired tests.
*/

#include "coordinator/primary_nodelet.h"
#include "coordinator/primary_tests.hpp"
#include "ros/ros.h"
//#include "std_msgs/Int32.h"
/* ************************************************************************* */
void PrimaryNodelet::Initialize(ros::NodeHandle* nh) {
  /**
  * @brief This is called when the nodelet is loaded into the nodelet manager
  * 
  */
  // Create Multi-threaded NH
  MTNH = getMTNodeHandle();

  // Load Params
  load_params();

  // publishers
  pub_flight_mode_ = nh->advertise<ff_msgs::FlightMode>(TOPIC_MOBILITY_FLIGHT_MODE, 1, true);  // FlightMode
  pub_status_ = nh->advertise<coordinator::StatusPrimary>(TOPIC_ASAP_STATUS, 5, true);
  pub_ctl_=nh->advertise<ff_msgs::FamCommand>(TOPIC_GNC_CTL_CMD,1);
  VL_status=nh->advertise<coordinator::Prediction>(VIRTUAL_LEADER_TOPIC,1);
  
  // subscribers
  sub_flight_mode_= nh->subscribe<ff_msgs::FlightMode>(TOPIC_MOBILITY_FLIGHT_MODE, 5,
    boost::bind(&PrimaryNodelet::flight_mode_callback, this, _1));  // flight mode getter
  // sub_ekf_ = nh->subscribe<ff_msgs::EkfState>("gnc/ekf", 5,
  //   boost::bind(&PrimaryNodelet::ekf_callback, this, _1));;
  sub_test_number_ = nh->subscribe<coordinator::TestNumber>(TOPIC_ASAP_TEST_NUMBER, 5,
    boost::bind(&PrimaryNodelet::test_num_callback, this, _1));
  sub_flight_mode_= nh->subscribe<ff_msgs::FlightMode>(TOPIC_MOBILITY_FLIGHT_MODE, 5,
    boost::bind(&PrimaryNodelet::flight_mode_callback, this, _1));  // flight mode setter
  sub_ekf_ = nh->subscribe<ff_msgs::EkfState>(TOPIC_GNC_EKF_,3,
    boost::bind(&PrimaryNodelet::ekf_callback, this, _1));
  sub_ekf_VF = nh->subscribe<ff_msgs::EkfState>(VIRTUAL_FOLLOWER_TOPIC,3,
    boost::bind(&PrimaryNodelet::VF_callback, this, _1));
  // sub_VL_status= nh->subscribe<coordinator::Prediction>(VIRTUAL_LEADER_TOPIC, 5,
  //   boost::bind(&PrimaryNodelet::VL_callback, this, _1));


  // services
  serv_ctl_enable_ = nh->serviceClient<std_srvs::SetBool>(SERVICE_GNC_CTL_ENABLE_); //"/gnc/ctl/eanble" /queen/gnc/ctl/enable
  //serv_ctl_enable_ = nh->serviceClient<std_srvs::SetBool>("/queen/gnc/ctl/enable"); //


  // tracking points
  try{
    std::vector<double> param_data;
    nh->getParam("/asap/primary/point_a_granite", param_data);
    POINT_A_GRANITE = Eigen::Matrix<double, 7, 1>(param_data.data());  // init from std::vector param
    nh->getParam("/asap/primary/point_a_iss", param_data);
    POINT_A_ISS = Eigen::Matrix<double, 7, 1>(param_data.data());  
    nh->getParam("/asap/primary/point_b_granite", param_data);
    POINT_B_GRANITE = Eigen::Matrix<double, 7, 1>(param_data.data());  
    nh->getParam("/asap/primary/point_b_iss", param_data);
    POINT_B_ISS= Eigen::Matrix<double, 7, 1>(param_data.data());  
    nh->getParam("/asap/primary/point_c_granite", param_data);
    POINT_C_GRANITE = Eigen::Matrix<double, 7, 1>(param_data.data());  
    nh->getParam("/asap/primary/point_c_iss", param_data);
    POINT_C_ISS = Eigen::Matrix<double, 7, 1>(param_data.data()); 
  }
  catch (const std::exception &exc) {
    std::cerr << exc.what() << std::endl;
    std::cout << "[PRIMARY_COORD]: Input parameters are invalid!" << std::endl;
  }

  // Pass control to Run method (activate timers and spin)
  std::cout << "[PRIMARY_COORD] Initialized." << std::endl;
  thread_.reset(new std::thread(&CoordinatorBase::Run, this, nh));
}



/* ************************************************************************** */
void PrimaryNodelet::get_status_msg(coordinator::StatusPrimary& msg){
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

  msg.control_mode = primary_status_.control_mode;
  msg.description = primary_status_.description;
}

/* ************************************************************************** */
void PrimaryNodelet::load_params(){
  // Get sim and ground flags
  std::string sim_str, ground_str;
  ros::param::get("/asap/sim", sim_str);
  sim_ = !std::strcmp(sim_str.c_str(), "true"); // convert to bool
  ros::param::get("/asap/ground", ground_str);
  ground_ = !std::strcmp(ground_str.c_str(), "true");  // convert to bool, 1 if it's true
  
  // get the robot name

  std::string robot_ns ;
  ros::param::get("/asap/primary_robot_name", robot_ns);

  // get the follower robot name

  std::string secondary_robot_ns ;
  ros::param::get("/asap/secondary_robot_name", secondary_robot_ns);


  std::cout << "[PRIMARY_COORD] name spaced topics ................" << std::endl;
  
  TOPIC_ASAP_STATUS = robot_ns + TOPIC_ASAP_STATUS ;
  std::cout << TOPIC_ASAP_STATUS<< std::endl;

  TOPIC_ASAP_TEST_NUMBER = robot_ns + TOPIC_ASAP_TEST_NUMBER;
  std::cout << TOPIC_ASAP_TEST_NUMBER<< std::endl;

  TOPIC_GNC_CTL_CMD = robot_ns + TOPIC_GNC_CTL_CMD;
  std::cout << TOPIC_GNC_CTL_CMD<< std::endl;

  SERVICE_GNC_CTL_ENABLE_ = robot_ns + SERVICE_GNC_CTL_ENABLE_;
  std::cout << SERVICE_GNC_CTL_ENABLE_<< std::endl;

  TOPIC_GNC_EKF_ = robot_ns + TOPIC_GNC_EKF_;
  std::cout << TOPIC_GNC_EKF_<< std::endl;

  VIRTUAL_FOLLOWER_TOPIC = secondary_robot_ns + VIRTUAL_FOLLOWER_TOPIC;
  std::cout << VIRTUAL_FOLLOWER_TOPIC<< std::endl;

  std::cout << "[PRIMARY_COORD] ................" << std::endl;

/* static std::string TOPIC_ASAP_STATUS = "/queen/asap/status";
static std::string TOPIC_ASAP_TEST_NUMBER = "/queen/asap/test_number";
static std::string TOPIC_GNC_CTL_CMD = "/queen/gnc/ctl/command";
 */  
  // regulation
  ros::param::getCached("/asap/primary/reg_time", reg_time_);
  ros::param::getCached("/asap/primary/x_start", x0_(0));
  ros::param::getCached("/asap/primary/y_start", x0_(1));
  ros::param::getCached("/asap/primary/z_start", x0_(2));
  ros::param::getCached("/asap/primary/qx_start", a0_(0));
  ros::param::getCached("/asap/primary/qy_start", a0_(1));
  ros::param::getCached("/asap/primary/qz_start", a0_(2));
  ros::param::getCached("/asap/primary/qw_start", a0_(3));
  ros::param::getCached("/asap/primary/pos_reg_thresh", pos_reg_thresh_);
  ros::param::getCached("/asap/primary/vel_reg_thresh", vel_reg_thresh_);
  ros::param::getCached("/asap/primary/att_reg_thresh", att_reg_thresh_);
  ros::param::getCached("/asap/primary/omega_reg_thresh", omega_reg_thresh_);

  ROS_INFO("[PRIMARY_COORD]....Goal position: X: %f Y: %f Z: %f ",x0_(0),x0_(1),x0_(2));
}
