#pragma once

/*
# coordinator

High-level logic coordinating test operation.

Every test has a test#() function available in case it is needed by asap.py
*/

// ROS
#include <ros/ros.h>
#include <nodelet/nodelet.h>
#include <pluginlib/class_list_macros.h>
#include <ros/package.h>
#include "eigen_conversions/eigen_msg.h"
#include <tf2/LinearMath/Quaternion.h>
#include <tf2/LinearMath/Transform.h>
#include <tf2_geometry_msgs/tf2_geometry_msgs.h>



// FSW includes
#include <ff_util/ff_nodelet.h>
#include <ff_util/ff_names.h>
#include <ff_util/ff_flight.h>
#include <ff_util/ff_action.h>

// msg includes
#include <ff_msgs/ControlState.h>
#include <ff_msgs/FlightMode.h>
#include <std_msgs/Float64MultiArray.h>
#include <std_msgs/String.h>
#include <msg_conversions/msg_conversions.h>
#include <ff_msgs/FamCommand.h>
#include <geometry_msgs/PoseWithCovariance.h>
#include <geometry_msgs/TwistWithCovariance.h>
#include <geometry_msgs/Inertia.h>
#include <geometry_msgs/Wrench.h>
#include <coordinator/StatusPrimary.h>
#include <coordinator/StatusSecondary.h>
#include <coordinator/TestNumber.h>
//#include <std_msgs/Int32.h>
// Actions
#include <ff_msgs/ControlAction.h>

// Service message
#include <std_srvs/SetBool.h>

// C++
#include <Eigen/Dense>
#include <vector>
#include <fstream>
#include <chrono>
#include <string.h>
#include <sstream>
#include <math.h>

//mathlab code generation
// Include Files
//#include "MPC_Guidance_v3_sand.h"
#include "mldivide.h"
#include "rt_nonfinite.h"
#include <cmath>
#include <cstring>
#include "rtwtypes.h"
#include <cstddef>
#include <cstdlib>





static std::string TOPIC_ASAP_STATUS = "asap/status";
static std::string TOPIC_ASAP_TEST_NUMBER = "asap/test_number";
static std::string TOPIC_GNC_CTL_CMD = "gnc/ctl/command";




// base status struct (key information)
struct BaseStatus {
  int test_number = -2;
  std::string flight_mode = "nominal";
  bool test_finished = false;

  bool coord_ok = true;
  bool regulate_finished = false;

  bool default_control = true;  // {true, false}: allow the default controller to run?
};


template <typename T>  // for PrimaryStatus or SecondaryStatus 
class CoordinatorBase 
{
 public:
  CoordinatorBase() {}; // don't do anything ROS-related in the constructor!
  ~CoordinatorBase() {};
  
  // Base main functions
  void Run(ros::NodeHandle *nh);
 protected:
  BaseStatus base_status_;

  ros::NodeHandle MTNH;
  std::shared_ptr<std::thread> thread_;

  ros::Publisher pub_flight_mode_;
  ros::Publisher pub_status_;
  ros::Publisher pub_ctl_;

  ros::Subscriber sub_flight_mode_;
  ros::Subscriber sub_ekf_;
  ros::Subscriber sub_test_number_;

  ros::ServiceClient serv_ctl_enable_;

  ros::Timer status_timer_;
  ros::Timer ctl_disable_timer_;

  ff_msgs::FlightMode flight_mode_;
  ff_msgs::EkfState ekf_state_;
  ff_msgs::FamCommand gnc_setpoint;

  geometry_msgs::Wrench ctl_input;
  geometry_msgs::Quaternion attitude;
  geometry_msgs::Vector3 omega,velocity_, position_, position_error, position_ref;
  tf2::Quaternion attitude_,q_ref,q_e,q_ref_inv;

  // Parameters
  bool ground_ = false;  // whether or not this is a ground test
  bool sim_ = false;
  std::string flight_mode_name_;

  // Stored status parameters
  std::string stored_control_mode_ = "track";  // stored control_mode, set by parameter inputs

  // Ekf state
  Eigen::Matrix<double, 16, 1> x_real_complete_;

  // Test number processing
  void process_test_number();
  void get_flight_mode();
  
  void publish_status(const ros::TimerEvent&);  // templated for Primary or Secondary status
  virtual void get_status_msg(coordinator::StatusPrimary& msg) {};
  virtual void get_status_msg(coordinator::StatusSecondary& msg) {};

  void test_num_callback(const coordinator::TestNumber::ConstPtr msg);
  void flight_mode_callback(const ff_msgs::FlightMode::ConstPtr msg);
  void ekf_callback(const ff_msgs::EkfState::ConstPtr msg);

  void debug();

  // Astrobee GNC interface
  void disable_default_ctl_callback(const ros::TimerEvent&);
  void disable_default_ctl();
  void enable_default_ctl();
  

  // Virtual test list: to be replaced on each derived coordinator
  virtual void RunTest0(ros::NodeHandle *nh) {};
  virtual void RunTest1(ros::NodeHandle *nh) {};
  virtual void RunTest2(ros::NodeHandle *nh) {};


  //controller1ModelClass controller1_Obj;// Instance of model class

// '<Root>/x_e'
float arg_x_e = 0.0;

// '<Root>/y_e'
 float arg_y_e = 0.0;

// '<Root>/z_e'
 float arg_z_e = 0.0;

// '<Root>/vx'
 float arg_vx = 0.0;

// '<Root>/vy'
 float arg_vy = 0.0;

// '<Root>/vz'
 float arg_vz = 0.0;

// '<Root>/qx'
 float arg_qx = 0.0;

// '<Root>/qy'
 float arg_qy = 0.0;

// '<Root>/qz'
 float arg_qz = 0.0;

// '<Root>/qw'
 float arg_qw = 0.0;

// '<Root>/omegax'
 float arg_omegax = 0.0;

// '<Root>/omegay'
 float arg_omegay = 0.0;

// '<Root>/omegaz'
 float arg_omegaz = 0.0;

// '<Root>/fx'
 float arg_fx;

// '<Root>/fy'
 float arg_fy;

// '<Root>/fz'
 float arg_fz;

// '<Root>/tau_x'
 float arg_tau_x;

// '<Root>/tau_y'
 float arg_tau_y;

// '<Root>/tau_z'
 float arg_tau_z;

double x0[6];
double Fx;
double Fy;
double Fz;
double target_state[6];
double dock_flag;
double CollAvoid_flag;
double dock_complete;
double num_iter;
double X_QP[60]={0};
double pt_sel;
double dr;

double z_nominal[6];
double zp_nextNominal[6];
double v_mpc[3];
bool initial_run = true;
bool initialzation = false;

double kn_tilda[3];
double kN[3];
bool rotation_done = false;
void step_PID();

// Function Declarations
void main_MPC_Guidance_v3_sand();
void MPC_Guidance_v3_sand();
void mldivide(double A[3600], double B[60]);
bool rtIsNaN(double value);
void nominal_dynamics();
void tubing_mpc();
  // add test definitions as necessary here
  // you can add more tests as desired in primary.h and secondary.h
};


/* ************************************************************************** */
// Coordinator template implementation
/* ************************************************************************** */

/* ************************************************************************** */
template<typename T>
void CoordinatorBase<T>::Run(ros::NodeHandle *nh) {
  /**
   * @brief Interpret the `/test_number` topic and run tests. High-level test management
   * takes place here, and is delegated out to other nodes. 
   * Each test function is intended to be run just ONCE per test number received.
   * This is the place to add new test numbers!
   * 
   */
  x_real_complete_ << 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 1.0,  0.0, 0.0, 0.0,  0.0, 0.0, 0.0,  0.0, 0.0, 0.0;

  // Status publisher in separate thread
  status_timer_ = MTNH.createTimer(ros::Duration(0.2),
    boost::bind(&CoordinatorBase::publish_status, this, _1));  // send commands (5 Hz)

  // Controller disabler in separate thread
  // ctl_disable_timer_ = MTNH.createTimer(ros::Duration(5.0),
  //   boost::bind(&CoordinatorBase::disable_default_ctl_callback, this, _1));  // send disable commands (0.2 Hz)

  // Check updates at 10Hz
  ros::Rate sleep_rate(10.0);

  // (1) Start up and wait for test_number_
  while (ros::ok() && base_status_.test_number == -2) {  // startup test number
    ros::spinOnce();
    sleep_rate.sleep();
  }

  // (2) Publish default flight mode so FAM will actually perform actuation
  get_flight_mode();
  pub_flight_mode_.publish(flight_mode_);  // TODO: should this be a more aggressive flight mode?
  ros::Duration(2.0).sleep();  // Pause so flight mode actually gets registered

  // (3) Execute test_number_ logic
  while (ros::ok()) { 
    if (!base_status_.test_finished) {
      // Tests go below...
      if (base_status_.test_number == 0) {
        RunTest0(nh);
      }
      else if (base_status_.test_number  == 1) {
        RunTest1(nh);
      }
      else if(base_status_.test_number  == 2) {
        RunTest2(nh);
      }
      // add additional test checks here

      base_status_.test_finished = true;
    }
    ros::spinOnce();
    sleep_rate.sleep();
  }
}


/* ************************************************************************** */
template <typename T>
void CoordinatorBase<T>::get_flight_mode() {
  /* Get a nominal flight mode for desired test.
  */
  if (base_status_.test_number != -1 ) {  // NOT shutdown test number or checkout test {and base_status_.test_number != 0}
    // create a nominal FlightMode
    if (!ff_util::FlightUtil::GetFlightMode(flight_mode_, base_status_.flight_mode)) {
      return;
    } 
  }
  else { // -1 shutdown test number
    // Shutdown Astrobee (turn off impellers)
    if (!ff_util::FlightUtil::GetFlightMode(flight_mode_, "off")) {
      return;
    }  
  }
}


/* ************************************************************************** */
template <typename T>
void CoordinatorBase<T>::publish_status(const ros::TimerEvent&) {
  /**
   * @brief Main coordinator of Base logic. Relies on get_status_msg, defined in derived class.
   * Uses either Primary or Secondary status logic.
   * 
   */
  T msg;
  get_status_msg(msg);
  pub_status_.publish(msg);
}


/* ************************************************************************** */
template<typename T>
void CoordinatorBase<T>::process_test_number() {
  /**
   * @brief Process test number logic for parameters
   * An example is provided here for processing individual digits of a test greater than 100.
   * 
   */

  if (base_status_.test_number > 100) {
    std::string test_number_str = std::to_string(base_status_.test_number);

    // Parameter settings xx(xxxxxx)
    // controller
    if (test_number_str[2] == '2') {  // standard MPC
      stored_control_mode_ = "track";
    }
    else if (test_number_str[2] == '3') {  // tube MPC
      stored_control_mode_ = "track_tube";
    }
  }
}


/* ************************************************************************** */
template<typename T>
void CoordinatorBase<T>::test_num_callback(const coordinator::TestNumber::ConstPtr msg) {
  /**
   * @brief Updates test numbers received from exec_asap
   * 
   */
  base_status_.test_number = msg->test_number;
  if (base_status_.test_number == -1) {
    // Re-enable default controller
    enable_default_ctl();

    // Set flight mode to off
    base_status_.flight_mode = "off";
    if (!ff_util::FlightUtil::GetFlightMode(flight_mode_, base_status_.flight_mode)) {
        return;
    }
    pub_flight_mode_.publish(flight_mode_);
  }
  // ROS_INFO("It is working bro dont be sad ;) ");
}


/* ************************************************************************** */
template<typename T>
void CoordinatorBase<T>::flight_mode_callback(const ff_msgs::FlightMode::ConstPtr msg) {
  /**
   * @brief Callback for flight mode.
   * 
   */
  flight_mode_name_ = msg->name;

  // kill ctl if it tries to turn on
  if (base_status_.default_control == false){
    auto serv_start = std::chrono::high_resolution_clock::now();

    // Disable the default controller so custom controller can run
    std_srvs::SetBool srv;
    srv.request.data = false;
    serv_ctl_enable_.call(srv);

    auto serv_finish = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> serv_elapsed = serv_finish - serv_start;
    std::string response = srv.response.message;

    std::cout << "[COORDINATOR]: Controller disable service time: " << serv_elapsed.count() << " seconds."<< std::endl;
    std::cout << "[COORDINATOR]: Controller disable service result: " << response << std::endl;
  }
}


/* ************************************************************************** */
template<typename T>
void CoordinatorBase<T>::ekf_callback(const ff_msgs::EkfState::ConstPtr msg) {
  /**
   * @brief The `gnc/ekf` subscriber callback. Called at 62.5 Hz.
   * Used to check if regulation is finished.
   * 
   */
  float qx = msg->pose.orientation.x;
  float qy = msg->pose.orientation.y;
  float qz = msg->pose.orientation.z;
  float qw = msg->pose.orientation.w;

  float px = msg->pose.position.x;
  float py = msg->pose.position.y;
  float pz = msg->pose.position.z;

  float vx = msg->velocity.x;
  float vy = msg->velocity.y;
  float vz = msg->velocity.z;


  float wx = msg->omega.x;
  float wy = msg->omega.y;
  float wz = msg->omega.z;

  if (qx != 0 || qy != 0 || qz != 0 || qw != 0) {
    x_real_complete_(0) = msg->pose.position.x;
    x_real_complete_(1) = msg->pose.position.y;
    x_real_complete_(2) = msg->pose.position.z;
    x_real_complete_(3) = msg->pose.orientation.x;
    x_real_complete_(4) = msg->pose.orientation.y;
    x_real_complete_(5) = msg->pose.orientation.z;
    x_real_complete_(6) = msg->pose.orientation.w;
    x_real_complete_(7) = msg->velocity.x;
    x_real_complete_(8) = msg->velocity.y;
    x_real_complete_(9) = msg->velocity.z;
    x_real_complete_(10) = msg->omega.x;
    x_real_complete_(11) = msg->omega.y;
    x_real_complete_(12) = msg->omega.z;
    x_real_complete_(13) = 0.0;
    x_real_complete_(14) = 0.0;
    x_real_complete_(15) = 0.0;
    }
    attitude.x=qx;
    attitude.y=qy;
    attitude.z=qz;
    attitude.w=qw;
    omega.x=wx;
    omega.y=wy;
    omega.z=wz;
    geometry_msgs::Vector3 torque;
    double r=0, p=0, y=-3.1459/4;  // Rotate the previous pose by 45* about Z

        q_ref.setRPY(r, p, y);
        tf2::convert(attitude,attitude_);
        q_ref_inv=q_ref.inverse();//
  q_e= q_ref_inv*attitude_;  // Calculate the new orientation
  q_e.normalize();
        float R_11 = 2*(attitude.x*attitude.x + attitude.w*attitude.w)-1;
        float R_12 = 2*(attitude.x*attitude.y - attitude.w*attitude.z);
        float R_13 = 2*(attitude.x*attitude.z + attitude.w*attitude.y); 
        float R_21 = 2*(attitude.x*attitude.y + attitude.w*attitude.z);
        float R_22 = 2*(attitude.y*attitude.y + attitude.w*attitude.w)-1;
        float R_23 = 2*(attitude.y*attitude.z - attitude.w*attitude.x);
        float R_31 = 2*(attitude.x*attitude.z - attitude.w*attitude.y); 
        float R_32 = 2*(attitude.y*attitude.z + attitude.w*attitude.x);
        float R_33 = 2*(attitude.z*attitude.z + attitude.w*attitude.w)-1;
        
        float cm_x =-0.3;
        float cm_y =0.0;
        float cm_z =0.0;


    position_.x = px ;//+ (cm_x*R_11 + cm_y*R_21 + cm_z*R_31);
    position_.y = py ;//+ (cm_x*R_12 + cm_y*R_22 + cm_z*R_32);
    position_.z = pz ;//+ (cm_x*R_13 + cm_y*R_23 + cm_z*R_33);

    /* position_ref.x = 10.8333388725;
    position_ref.y = -9.41988714508+0.5;
    position_ref.z = 4.20110343832;  */

  if(initialzation){
    position_error.x = position_.x - position_ref.x;
    position_error.y = position_.y - position_ref.y;
    position_error.z = position_.z - position_ref.z;

    velocity_.x=vx;
    velocity_.y=vy;
    velocity_.z=vz;

  

     arg_x_e = position_error.x;
     arg_y_e = position_error.y;
     arg_z_e = position_error.z;
     arg_vx  = vx;
     arg_vy  = vy;
     arg_vz  = vz;
     arg_qx = q_e.getX();
     arg_qy = q_e.getY();
     arg_qz = q_e.getZ();
     arg_qw = q_e.getW();
     arg_omegax = wx;
     arg_omegay = wy;
     arg_omegaz = wz;
    /* double *fx;
    double *fy;
    double *fz;
    double *taux;
    double *tauy;
    double *tauz; */
    step_PID();
    x0[0]=position_error.x;
    x0[1]=position_error.y;
    x0[2]=position_error.z;
    x0[3]=vx;
    x0[4]=vy;
    x0[5]=vz;
    main_MPC_Guidance_v3_sand();
    if (sqrt(q_e.getX()*q_e.getX()+q_e.getY()*q_e.getY()+q_e.getZ()*q_e.getZ())>0.005){
      z_nominal[0]=x0[0];
      z_nominal[1]=x0[1];
      z_nominal[2]=x0[2];
      z_nominal[3]=x0[3];
      z_nominal[4]=x0[4];
      z_nominal[5]=x0[5];

      initial_run=false;


    }
    else{
      z_nominal[0]=zp_nextNominal[0];
      z_nominal[1]=zp_nextNominal[1];
      z_nominal[2]=zp_nextNominal[2];
      z_nominal[3]=zp_nextNominal[3];
      z_nominal[4]=zp_nextNominal[4];
      z_nominal[5]=zp_nextNominal[5];


    }

    v_mpc[0]=Fx;
    v_mpc[1]=Fy;
    v_mpc[2]=Fz;
    nominal_dynamics();
    kn_tilda[0]=Fx;
    kn_tilda[1]=Fy;
    kn_tilda[2]=Fz;

      double sx=x0[0]-zp_nextNominal[0];
      double sy=x0[1]-zp_nextNominal[1];
      double sz=x0[2]-zp_nextNominal[2];
      double svx=x0[3]-zp_nextNominal[3];
      double svy=x0[4]-zp_nextNominal[4];
      double svz=x0[5]-zp_nextNominal[5];

    tubing_mpc();
    //X_QP=X_Qp
   // rt_OneStep();
   //ROS_INFO("ex: [%f]  ey: [%f] ez: [%f] ev_x: [%f] ev_y: [%f] ev_z: [%f]", sx,sy,sz,svx,svy,svz);

   // ROS_INFO("fx: [%f]  fy: [%f] fz: [%f] tau_x: [%f] tau_y: [%f] tau_y: [%f]", Fx,Fy,Fz, arg_tau_x,arg_tau_y,arg_tau_z);
  }
}


/* ************************************************************************** */
template<typename T>
void CoordinatorBase<T>::debug(){
  /**
   * @brief debug function call to test compatibility with other Bases
   * 
   */

}
/* *************************************************************************** */
template<typename T>
void CoordinatorBase<T>::tubing_mpc(){
  double a[18] = { 2.0202, 0.0, 0.0, 0.0, 2.0202, 0.0, 0.0, 0.0,
    2.0202, 4.0895, 0.0, 0.0, 0.0, 4.0895, 0.0, 0.0, 0.0, 4.0895 };
  double b_x[6];
  double d;

  if (((kn_tilda[0] == 0.0) && (kn_tilda[1] == 0.0)) && (kn_tilda[2] == 0.0)) {
    kN[0] = 0.0;
    kN[1] = 0.0;
    kN[2] = 0.0;
  } else {
    int i;
    for (i = 0; i < 6; i++) {
      b_x[i] = x0[i] - z_nominal[i];
    }

    for (i = 0; i < 3; i++) {
      d = 0.0;
      for (int i1 = 0; i1 < 6; i1++) {
        d += a[i + (3 * i1)] * b_x[i1];
      }

      kN[i] = kn_tilda[i] - d;
    }
  }



}
/* **************************************************************************** */
template<typename T>
void CoordinatorBase<T>::nominal_dynamics(){
   double b_a[36] = { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 2.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 1.0 };

   double a[18] = { 0.086693, 0.0, 0.0, 0.086693, 0.0, 0.0, 0.0,
    0.086693, 0.0, 0.0, 0.086693, 0.0, 0.0, 0.0, 0.086693, 0.0, 0.0,
    0.086693 };

   double d;

  // MPCParams.A;
  // MPCParams.B;
  for (int i = 0; i < 6; i++) {
    d = 0.0;
    for (int i1 = 0; i1 < 6; i1++) {
      d += b_a[i + (6 * i1)] * z_nominal[i1];
    }

    zp_nextNominal[i] = d + (((a[i] * v_mpc[0]) + (a[i + 6] * v_mpc[1])) + (a[i + 12] * v_mpc[2]));
  }



}

template<typename T>
void CoordinatorBase<T>::step_PID(){
  /**
   * @brief debug function call to test compatibility with other Bases
   * 
   */

    float a[18] = { -0.3832000000000001, -0.0, -0.0, -0.0,
    -0.3832000000000001, -0.0, -0.0, -0.0, -0.3832000000000001,
    -3.8281680000000002, -0.0, -0.0, -0.0, -3.8281680000000002, -0.0, -0.0, -0.0,
    -3.8281680000000002 };
    
    /* { -4.0, -0.0, -0.0, -0.0, -4.0, -0.0, -0.0, -0.0,
    -4.0, -2.828, -0.0, -0.0, -0.0, -2.828, -0.0, -0.0, -0.0, -2.828 }; */

   float a_0[18] ={ -0.0095625, -0.0, -0.0, -0.0, -0.0089375, -0.0,
    -0.0, -0.0, -0.0101875, -0.0764235, -0.0, -0.0, -0.0, -0.071428499999999992,
    -0.0, -0.0, -0.0, -0.0814185 };
   
    /* { -0.612, -0.0, -0.0, -0.0, -0.572, -0.0, -0.0,
    -0.0, -0.652, -0.43268399999999996, -0.0, -0.0, -0.0, -0.40440399999999993,
    -0.0, -0.0, -0.0, -0.460964 }; */

  float tmp[6];
  float u[3];
  float tau[3];
  int i;
  int i_0;
  //UNUSED_PARAMETER(arg_qw);

  // MATLAB Function: '<Root>/Position Controller' incorporates:
  //   Inport: '<Root>/vx'
  //   Inport: '<Root>/vy'
  //   Inport: '<Root>/vz'
  //   Inport: '<Root>/x_e'
  //   Inport: '<Root>/y_e'
  //   Inport: '<Root>/z_e'

  tmp[0] = arg_x_e;
  tmp[1] = arg_y_e;
  tmp[2] = arg_z_e;
  tmp[3] = arg_vx;
  tmp[4] = arg_vy;
  tmp[5] = arg_vz;
  for (i = 0; i < 3; i++) {
    u[i] = 0.0;
    for (i_0 = 0; i_0 < 6; i_0++) {
      u[i] += a[3 * i_0 + i] * tmp[i_0];
    }
  }


 /*  tmp[0] = arg_qx;
  tmp[1] = arg_qy;
  tmp[2] = arg_qz;
  tmp[3] = arg_omegax;
  tmp[4] = arg_omegay;
  tmp[5] = arg_omegaz;
  for (i = 0; i < 3; i++) {
    tau[i] = 0.0;
    for (i_0 = 0; i_0 < 6; i_0++) {
      tau[i] += a_0[3 * i_0 + i] * tmp[i_0];
    }
  } */
 arg_fx=u[0];
 arg_fy=u[1];
 arg_fz=u[2];


 double qw = arg_qx;
 double qx = arg_qx;
 double qy = arg_qy;
 double qz = arg_qz;
 double fx = Fx; 
 double fy = Fy;
 double fz = Fz;
 double Omega1 = arg_omegax;
 double Omega2 = arg_omegay;
 double Omega3 = arg_omegaz;

float a1[18] = { -0.0053271706630077078, -0.0, -0.0, -0.0,
    -0.020391594065038533, -0.0, -0.0, -0.0, -0.016794590442030823,
    -0.26635853315038532, -0.0, -0.0, -0.0, -1.0195797032519265, -0.0, -0.0,
    -0.0, -0.83972952210154106 };

  static const real_T b_a[9] = { 0.0, -0.071058772687986174, 0.0,
    0.071058772687986174, 0.0, 0.14211754537597235, -0.0, -0.14211754537597235,
    0.0 };

  real_T dv[9];
  real_T b_qx[6];
  real_T U[3];
  real_T dv1[3];
  real_T d;
  real_T d1;
  real_T d2;
  real_T d3;
  real_T d4;
  //int32_T i;

  // best case
  // 5.58;%9.58;%=2*9.58;
  // 9.58;
  //  d=0.6; % Distance between two Astrobee com
  // d/5;
  //  M=2;%linear density of arm kg/m
  // J1=M*l^3/4
  // J2=M*l^3/4
  // 0.453105990000000;% =2*0.153427995+13/16 Ml^2 Ixx of compund object
  // 0.540000000000000;%3*M*l^2;
  // 10.236728099999999 ;%=190/16M*l^2+2*(J1+J2+0.14271405)+9*9.58*l^2 Iyy of compund object 
  // 9.960905517999999 ;%=162/16*M*l^2+2*(J1+J2+0.162302759)+9*9.58*l^2 Izz of compund object 
  //  Ix=0.306855990;% =2*0.153427995 Ixx of compund object
  //  Iy=7.183028100000000 ;%=2*0.14271405+m*d^2 Iyy of compund object
  //  Iz=7.222205518000000 ;%=2*0.162302759+m*d^2 Izz of compund object
  //  PD gains
  // 1
  b_qx[0] = qx;
  b_qx[1] = qy;
  b_qx[2] = qz;
  b_qx[3] = Omega1;
  b_qx[4] = Omega2;
  b_qx[5] = Omega3;
  d = qw * qw;
  dv[0] = (2.0 * (d + (qx * qx))) - 1.0;
  d1 = qx * qy;
  d2 = qw * qz;
  dv[3] = 2.0 * (d1 - d2);
  d3 = qx * qz;
  d4 = qw * qy;
  dv[6] = 2.0 * (d3 + d4);
  dv[1] = 2.0 * (d1 + d2);
  dv[4] = (2.0 * (d + (qy * qy))) - 1.0;
  d1 = qy * qz;
  d2 = qw * qx;
  dv[7] = 2.0 * (d1 - d2);
  dv[2] = 2.0 * (d3 - d4);
  dv[5] = 2.0 * (d1 + d2);
  dv[8] = (2.0 * (d + (qz * qz))) - 1.0;
  for (i = 0; i < 3; i++) {
    dv1[i] = ((dv[i] * Fx) + (dv[i + 3] * Fy)) + (dv[i + 6] * Fz);
    d = 0.0;
    for (int32_T i1 = 0; i1 < 6; i1++) {
      d += a[i + (3 * i1)] * b_qx[i1];
    }

    U[i] = d;
  }

  d = dv1[0];
  d1 = dv1[1];
  d2 = dv1[2];
  for (i = 0; i < 3; i++) {
    U[i] -= ((b_a[i] * d) + (b_a[i + 3] * d1)) + (b_a[i + 6] * d2);
  }

  arg_tau_x = U[0];
  arg_tau_y = U[1];
  arg_tau_z = U[2];


}


//void rt_OneStep(void);
//void rt_OneStep(void)
/* template<typename T>
void CoordinatorBase<T>::rt_OneStep()
{
  static boolean_T OverrunFlag = false;

  // Disable interrupts here

  // Check for overrun
  if (OverrunFlag) {
    rtmSetErrorStatus(controller1ModelClass::getRTM(), "Overrun");
    return;
  }

  OverrunFlag = true;

  // Save FPU context here (if necessary)
  // Re-enable timer or interrupt here
  // Set model inputs here

  // Step the model
  controller1ModelClass::step(arg_x_e, arg_y_e, arg_z_e, arg_vx, arg_vy, arg_vz, arg_qx,
                       arg_qy, arg_qz, arg_qw, arg_omegax, arg_omegay,
                       arg_omegaz, &arg_fx, &arg_fy, &arg_fz, &arg_tau_x,
                       &arg_tau_y, &arg_tau_z);

  // Get model outputs here

  // Indicate task complete
  OverrunFlag = false;

  // Disable interrupts here
  // Restore FPU context here (if necessary)
  // Enable interrupts here
} */


/* ************************************************************************** */
template<typename T>
void CoordinatorBase<T>::disable_default_ctl_callback(const ros::TimerEvent&) {
  /**
   * @brief Switch default controller off, repeatedly.
   * @param base_status.default_control is monitored for activation
   */
  if (base_status_.default_control == false){
    auto serv_start = std::chrono::high_resolution_clock::now();

    // Disable the default controller so custom controller can run
    std_srvs::SetBool srv;
    srv.request.data = false;
    serv_ctl_enable_.call(srv);

    auto serv_finish = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> serv_elapsed = serv_finish - serv_start;
    std::string response = srv.response.message;

    std::cout << "[COORDINATOR]: Controller disable service time: " << serv_elapsed.count() << " seconds."<< std::endl;
    std::cout << "[COORDINATOR]: Controller disable service result: " << response << std::endl;
  }
}


/* ************************************************************************** */
template<typename T>
void CoordinatorBase<T>::disable_default_ctl() {
  /**
   * @brief Switch default controller off.
   * @param base_status.default_control is monitored for activation
   */
  base_status_.default_control = false;

  // Disable the default controller so custom controller can run
  std_srvs::SetBool srv;
  srv.request.data = false;
  serv_ctl_enable_.call(srv);
}


/* ************************************************************************** */
template<typename T>
void CoordinatorBase<T>::enable_default_ctl() {
  /**
   * @brief Switch default controller on.
   * 
   */
  ROS_DEBUG_STREAM("[COORDINATOR]: Enabling default controller...");

  // Disable the default controller so tube-MPC can run
  base_status_.default_control = true;

  std_srvs::SetBool srv;
  srv.request.data = true;
  serv_ctl_enable_.call(srv);

  ROS_DEBUG_STREAM("[COORDINATOR]: Ctl enable service result: " << srv.response.message);
}


template<typename T>
void CoordinatorBase<T>::main_MPC_Guidance_v3_sand()
{

  // Initialize function 'MPC_Guidance_v3_sand' input arguments.
  // Initialize function input argument 'x0'.
  // Call the entry-point 'MPC_Guidance_v3_sand'.
  
  MPC_Guidance_v3_sand();
}


template<typename T>
bool  CoordinatorBase<T>::rtIsNaN(double value){
  return ((value!=value) ? true : false);
}

template<typename T>
void CoordinatorBase<T>::MPC_Guidance_v3_sand()
{
   static const real_T b[14400] = { 70.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 70.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 70.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100000.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100000.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100000.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    70.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 70.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 70.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 100000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 100000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 70.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 70.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 70.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100000.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    100000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 100000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 70.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 70.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 70.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100000.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100000.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 70.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 70.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    70.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 100000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 100000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 100000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 70.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 70.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 70.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100000.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100000.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    100000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 70.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 70.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 70.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 100000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100000.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 70.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 70.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 70.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    100000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 100000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 100000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 70.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 70.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 70.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100000.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100000.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100000.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    70.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 70.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 70.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 100000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 100000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 70.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 70.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 70.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100000.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    100000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 100000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 70.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 70.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 70.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100000.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100000.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 70.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 70.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    70.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 100000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 100000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 100000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 70.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 70.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 70.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100000.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100000.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    100000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 70.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 70.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 70.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 100000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100000.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 70.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 70.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 70.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    100000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 100000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 100000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 70.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 70.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 70.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100000.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100000.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100000.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    70.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 70.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 70.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 100000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 100000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 70.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 70.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 70.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100000.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    100000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 100000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1732.0163, 1.8E-12, 2.3679E-11, 16123.4203, 1.3401E-10,
    2.0746E-10, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.8E-12, 1732.0163,
    -2.7364E-12, 7.5606E-11, 16123.4203, -1.4753E-10, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 2.3679E-11, -2.7364E-12, 1732.0163, 2.2982E-10,
    8.8071E-11, 16123.4203, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    16123.4203, 7.5606E-11, 2.2982E-10, 831763.0736, 1.8698E-9, 2.4387E-9, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.3401E-10, 16123.4203, 8.8071E-11,
    1.8698E-9, 831763.0736, 5.2194E-10, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 2.0746E-10, -1.4753E-10, 16123.4203, 2.4387E-9, 5.2194E-10, 831763.0736
  };

  static const real_T b_b[7200] = { 0.052067, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.1562, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.26033, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.36447, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.4686, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.57274, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.67687, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.78101, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.88514, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.98927, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.0934, 0.0, 0.0, 0.052067, 0.0, 0.0,
    1.1975, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.3017, 0.0, 0.0, 0.052067, 0.0, 0.0,
    1.4058, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.5099, 0.0, 0.0, 0.052067, 0.0, 0.0,
    1.6141, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.7182, 0.0, 0.0, 0.052067, 0.0, 0.0,
    1.8223, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.9265, 0.0, 0.0, 0.052067, 0.0, 0.0,
    2.0306, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.052067, 0.0,
    0.0, 0.1562, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.26033, 0.0, 0.0, 0.052067, 0.0,
    0.0, 0.36447, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.4686, 0.0, 0.0, 0.052067, 0.0,
    0.0, 0.57274, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.67687, 0.0, 0.0, 0.052067, 0.0,
    0.0, 0.78101, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.88514, 0.0, 0.0, 0.052067, 0.0,
    0.0, 0.98927, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.0934, 0.0, 0.0, 0.052067, 0.0,
    0.0, 1.1975, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.3017, 0.0, 0.0, 0.052067, 0.0,
    0.0, 1.4058, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.5099, 0.0, 0.0, 0.052067, 0.0,
    0.0, 1.6141, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.7182, 0.0, 0.0, 0.052067, 0.0,
    0.0, 1.8223, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.9265, 0.0, 0.0, 0.052067, 0.0,
    0.0, 2.0306, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.052067,
    0.0, 0.0, 0.1562, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.26033, 0.0, 0.0, 0.052067,
    0.0, 0.0, 0.36447, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.4686, 0.0, 0.0, 0.052067,
    0.0, 0.0, 0.57274, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.67687, 0.0, 0.0, 0.052067,
    0.0, 0.0, 0.78101, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.88514, 0.0, 0.0, 0.052067,
    0.0, 0.0, 0.98927, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.0934, 0.0, 0.0, 0.052067,
    0.0, 0.0, 1.1975, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.3017, 0.0, 0.0, 0.052067,
    0.0, 0.0, 1.4058, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.5099, 0.0, 0.0, 0.052067,
    0.0, 0.0, 1.6141, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.7182, 0.0, 0.0, 0.052067,
    0.0, 0.0, 1.8223, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.9265, 0.0, 0.0, 0.052067,
    0.0, 0.0, 2.0306, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.052067,
    0.0, 0.0, 0.052067, 0.0, 0.0, 0.1562, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.26033,
    0.0, 0.0, 0.052067, 0.0, 0.0, 0.36447, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.4686,
    0.0, 0.0, 0.052067, 0.0, 0.0, 0.57274, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.67687,
    0.0, 0.0, 0.052067, 0.0, 0.0, 0.78101, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.88514,
    0.0, 0.0, 0.052067, 0.0, 0.0, 0.98927, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.0934,
    0.0, 0.0, 0.052067, 0.0, 0.0, 1.1975, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.3017,
    0.0, 0.0, 0.052067, 0.0, 0.0, 1.4058, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.5099,
    0.0, 0.0, 0.052067, 0.0, 0.0, 1.6141, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.7182,
    0.0, 0.0, 0.052067, 0.0, 0.0, 1.8223, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.9265,
    0.0, 0.0, 0.052067, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.052067,
    0.0, 0.0, 0.052067, 0.0, 0.0, 0.1562, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.26033,
    0.0, 0.0, 0.052067, 0.0, 0.0, 0.36447, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.4686,
    0.0, 0.0, 0.052067, 0.0, 0.0, 0.57274, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.67687,
    0.0, 0.0, 0.052067, 0.0, 0.0, 0.78101, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.88514,
    0.0, 0.0, 0.052067, 0.0, 0.0, 0.98927, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.0934,
    0.0, 0.0, 0.052067, 0.0, 0.0, 1.1975, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.3017,
    0.0, 0.0, 0.052067, 0.0, 0.0, 1.4058, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.5099,
    0.0, 0.0, 0.052067, 0.0, 0.0, 1.6141, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.7182,
    0.0, 0.0, 0.052067, 0.0, 0.0, 1.8223, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.9265,
    0.0, 0.0, 0.052067, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.052067,
    0.0, 0.0, 0.052067, 0.0, 0.0, 0.1562, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.26033,
    0.0, 0.0, 0.052067, 0.0, 0.0, 0.36447, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.4686,
    0.0, 0.0, 0.052067, 0.0, 0.0, 0.57274, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.67687,
    0.0, 0.0, 0.052067, 0.0, 0.0, 0.78101, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.88514,
    0.0, 0.0, 0.052067, 0.0, 0.0, 0.98927, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.0934,
    0.0, 0.0, 0.052067, 0.0, 0.0, 1.1975, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.3017,
    0.0, 0.0, 0.052067, 0.0, 0.0, 1.4058, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.5099,
    0.0, 0.0, 0.052067, 0.0, 0.0, 1.6141, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.7182,
    0.0, 0.0, 0.052067, 0.0, 0.0, 1.8223, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.9265,
    0.0, 0.0, 0.052067, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.052067, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.1562, 0.0, 0.0, 0.052067, 0.0,
    0.0, 0.26033, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.36447, 0.0, 0.0, 0.052067, 0.0,
    0.0, 0.4686, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.57274, 0.0, 0.0, 0.052067, 0.0,
    0.0, 0.67687, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.78101, 0.0, 0.0, 0.052067, 0.0,
    0.0, 0.88514, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.98927, 0.0, 0.0, 0.052067, 0.0,
    0.0, 1.0934, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.1975, 0.0, 0.0, 0.052067, 0.0,
    0.0, 1.3017, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.4058, 0.0, 0.0, 0.052067, 0.0,
    0.0, 1.5099, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.6141, 0.0, 0.0, 0.052067, 0.0,
    0.0, 1.7182, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.8223, 0.0, 0.0, 0.052067, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.1562, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.26033, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.36447, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.4686, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.57274, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.67687, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.78101, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.88514, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.98927, 0.0, 0.0, 0.052067, 0.0, 0.0,
    1.0934, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.1975, 0.0, 0.0, 0.052067, 0.0, 0.0,
    1.3017, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.4058, 0.0, 0.0, 0.052067, 0.0, 0.0,
    1.5099, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.6141, 0.0, 0.0, 0.052067, 0.0, 0.0,
    1.7182, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.8223, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.052067,
    0.0, 0.0, 0.052067, 0.0, 0.0, 0.1562, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.26033,
    0.0, 0.0, 0.052067, 0.0, 0.0, 0.36447, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.4686,
    0.0, 0.0, 0.052067, 0.0, 0.0, 0.57274, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.67687,
    0.0, 0.0, 0.052067, 0.0, 0.0, 0.78101, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.88514,
    0.0, 0.0, 0.052067, 0.0, 0.0, 0.98927, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.0934,
    0.0, 0.0, 0.052067, 0.0, 0.0, 1.1975, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.3017,
    0.0, 0.0, 0.052067, 0.0, 0.0, 1.4058, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.5099,
    0.0, 0.0, 0.052067, 0.0, 0.0, 1.6141, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.7182,
    0.0, 0.0, 0.052067, 0.0, 0.0, 1.8223, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.1562, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.26033, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.36447, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.4686, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.57274, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.67687, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.78101, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.88514, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.98927, 0.0, 0.0, 0.052067, 0.0, 0.0,
    1.0934, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.1975, 0.0, 0.0, 0.052067, 0.0, 0.0,
    1.3017, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.4058, 0.0, 0.0, 0.052067, 0.0, 0.0,
    1.5099, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.6141, 0.0, 0.0, 0.052067, 0.0, 0.0,
    1.7182, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.1562, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.26033, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.36447, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.4686, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.57274, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.67687, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.78101, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.88514, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.98927, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.0934, 0.0, 0.0,
    0.052067, 0.0, 0.0, 1.1975, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.3017, 0.0, 0.0,
    0.052067, 0.0, 0.0, 1.4058, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.5099, 0.0, 0.0,
    0.052067, 0.0, 0.0, 1.6141, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.7182, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.052067, 0.0,
    0.0, 0.1562, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.26033, 0.0, 0.0, 0.052067, 0.0,
    0.0, 0.36447, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.4686, 0.0, 0.0, 0.052067, 0.0,
    0.0, 0.57274, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.67687, 0.0, 0.0, 0.052067, 0.0,
    0.0, 0.78101, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.88514, 0.0, 0.0, 0.052067, 0.0,
    0.0, 0.98927, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.0934, 0.0, 0.0, 0.052067, 0.0,
    0.0, 1.1975, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.3017, 0.0, 0.0, 0.052067, 0.0,
    0.0, 1.4058, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.5099, 0.0, 0.0, 0.052067, 0.0,
    0.0, 1.6141, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.7182, 0.0, 0.0, 0.052067, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.052067, 0.0,
    0.0, 0.1562, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.26033, 0.0, 0.0, 0.052067, 0.0,
    0.0, 0.36447, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.4686, 0.0, 0.0, 0.052067, 0.0,
    0.0, 0.57274, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.67687, 0.0, 0.0, 0.052067, 0.0,
    0.0, 0.78101, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.88514, 0.0, 0.0, 0.052067, 0.0,
    0.0, 0.98927, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.0934, 0.0, 0.0, 0.052067, 0.0,
    0.0, 1.1975, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.3017, 0.0, 0.0, 0.052067, 0.0,
    0.0, 1.4058, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.5099, 0.0, 0.0, 0.052067, 0.0,
    0.0, 1.6141, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.1562, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.26033, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.36447, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.4686, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.57274, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.67687, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.78101, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.88514, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.98927, 0.0, 0.0,
    0.052067, 0.0, 0.0, 1.0934, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.1975, 0.0, 0.0,
    0.052067, 0.0, 0.0, 1.3017, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.4058, 0.0, 0.0,
    0.052067, 0.0, 0.0, 1.5099, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.6141, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.1562, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.26033, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.36447, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.4686, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.57274, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.67687, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.78101, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.88514, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.98927, 0.0, 0.0, 0.052067, 0.0, 0.0,
    1.0934, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.1975, 0.0, 0.0, 0.052067, 0.0, 0.0,
    1.3017, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.4058, 0.0, 0.0, 0.052067, 0.0, 0.0,
    1.5099, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.6141, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.052067,
    0.0, 0.0, 0.052067, 0.0, 0.0, 0.1562, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.26033,
    0.0, 0.0, 0.052067, 0.0, 0.0, 0.36447, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.4686,
    0.0, 0.0, 0.052067, 0.0, 0.0, 0.57274, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.67687,
    0.0, 0.0, 0.052067, 0.0, 0.0, 0.78101, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.88514,
    0.0, 0.0, 0.052067, 0.0, 0.0, 0.98927, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.0934,
    0.0, 0.0, 0.052067, 0.0, 0.0, 1.1975, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.3017,
    0.0, 0.0, 0.052067, 0.0, 0.0, 1.4058, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.5099,
    0.0, 0.0, 0.052067, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.1562, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.26033, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.36447, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.4686, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.57274, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.67687, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.78101, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.88514, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.98927, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.0934, 0.0, 0.0, 0.052067, 0.0, 0.0,
    1.1975, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.3017, 0.0, 0.0, 0.052067, 0.0, 0.0,
    1.4058, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.5099, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.052067, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.1562, 0.0, 0.0, 0.052067, 0.0,
    0.0, 0.26033, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.36447, 0.0, 0.0, 0.052067, 0.0,
    0.0, 0.4686, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.57274, 0.0, 0.0, 0.052067, 0.0,
    0.0, 0.67687, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.78101, 0.0, 0.0, 0.052067, 0.0,
    0.0, 0.88514, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.98927, 0.0, 0.0, 0.052067, 0.0,
    0.0, 1.0934, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.1975, 0.0, 0.0, 0.052067, 0.0,
    0.0, 1.3017, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.4058, 0.0, 0.0, 0.052067, 0.0,
    0.0, 1.5099, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.052067, 0.0,
    0.0, 0.052067, 0.0, 0.0, 0.1562, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.26033, 0.0,
    0.0, 0.052067, 0.0, 0.0, 0.36447, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.4686, 0.0,
    0.0, 0.052067, 0.0, 0.0, 0.57274, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.67687, 0.0,
    0.0, 0.052067, 0.0, 0.0, 0.78101, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.88514, 0.0,
    0.0, 0.052067, 0.0, 0.0, 0.98927, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.0934, 0.0,
    0.0, 0.052067, 0.0, 0.0, 1.1975, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.3017, 0.0,
    0.0, 0.052067, 0.0, 0.0, 1.4058, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.1562, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.26033, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.36447, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.4686, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.57274, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.67687, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.78101, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.88514, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.98927, 0.0, 0.0,
    0.052067, 0.0, 0.0, 1.0934, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.1975, 0.0, 0.0,
    0.052067, 0.0, 0.0, 1.3017, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.4058, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.1562, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.26033, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.36447, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.4686, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.57274, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.67687, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.78101, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.88514, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.98927, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.0934, 0.0, 0.0,
    0.052067, 0.0, 0.0, 1.1975, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.3017, 0.0, 0.0,
    0.052067, 0.0, 0.0, 1.4058, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.1562,
    0.0, 0.0, 0.052067, 0.0, 0.0, 0.26033, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.36447,
    0.0, 0.0, 0.052067, 0.0, 0.0, 0.4686, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.57274,
    0.0, 0.0, 0.052067, 0.0, 0.0, 0.67687, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.78101,
    0.0, 0.0, 0.052067, 0.0, 0.0, 0.88514, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.98927,
    0.0, 0.0, 0.052067, 0.0, 0.0, 1.0934, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.1975,
    0.0, 0.0, 0.052067, 0.0, 0.0, 1.3017, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.1562, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.26033, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.36447, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.4686, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.57274, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.67687, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.78101, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.88514, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.98927, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.0934, 0.0, 0.0,
    0.052067, 0.0, 0.0, 1.1975, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.3017, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.052067, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.1562, 0.0, 0.0, 0.052067,
    0.0, 0.0, 0.26033, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.36447, 0.0, 0.0, 0.052067,
    0.0, 0.0, 0.4686, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.57274, 0.0, 0.0, 0.052067,
    0.0, 0.0, 0.67687, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.78101, 0.0, 0.0, 0.052067,
    0.0, 0.0, 0.88514, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.98927, 0.0, 0.0, 0.052067,
    0.0, 0.0, 1.0934, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.1975, 0.0, 0.0, 0.052067,
    0.0, 0.0, 1.3017, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.1562, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.26033, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.36447, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.4686, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.57274, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.67687, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.78101, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.88514, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.98927, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.0934, 0.0, 0.0,
    0.052067, 0.0, 0.0, 1.1975, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.1562, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.26033, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.36447, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.4686, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.57274, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.67687, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.78101, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.88514, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.98927, 0.0, 0.0, 0.052067, 0.0, 0.0,
    1.0934, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.1975, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.1562, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.26033, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.36447, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.4686, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.57274, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.67687, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.78101, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.88514, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.98927, 0.0, 0.0,
    0.052067, 0.0, 0.0, 1.0934, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.1975, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.1562, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.26033, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.36447, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.4686, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.57274, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.67687, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.78101, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.88514, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.98927, 0.0, 0.0, 0.052067, 0.0, 0.0, 1.0934, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.1562, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.26033, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.36447, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.4686, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.57274, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.67687, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.78101, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.88514, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.98927, 0.0, 0.0, 0.052067, 0.0, 0.0,
    1.0934, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.052067, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.1562, 0.0, 0.0, 0.052067,
    0.0, 0.0, 0.26033, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.36447, 0.0, 0.0, 0.052067,
    0.0, 0.0, 0.4686, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.57274, 0.0, 0.0, 0.052067,
    0.0, 0.0, 0.67687, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.78101, 0.0, 0.0, 0.052067,
    0.0, 0.0, 0.88514, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.98927, 0.0, 0.0, 0.052067,
    0.0, 0.0, 1.0934, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.1562, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.26033, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.36447, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.4686, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.57274, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.67687, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.78101, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.88514, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.98927, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.052067,
    0.0, 0.0, 0.1562, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.26033, 0.0, 0.0, 0.052067,
    0.0, 0.0, 0.36447, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.4686, 0.0, 0.0, 0.052067,
    0.0, 0.0, 0.57274, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.67687, 0.0, 0.0, 0.052067,
    0.0, 0.0, 0.78101, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.88514, 0.0, 0.0, 0.052067,
    0.0, 0.0, 0.98927, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.1562, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.26033, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.36447, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.4686, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.57274, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.67687, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.78101, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.88514, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.98927, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.052067, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.1562, 0.0, 0.0, 0.052067, 0.0,
    0.0, 0.26033, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.36447, 0.0, 0.0, 0.052067, 0.0,
    0.0, 0.4686, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.57274, 0.0, 0.0, 0.052067, 0.0,
    0.0, 0.67687, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.78101, 0.0, 0.0, 0.052067, 0.0,
    0.0, 0.88514, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.052067, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.1562, 0.0, 0.0, 0.052067, 0.0,
    0.0, 0.26033, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.36447, 0.0, 0.0, 0.052067, 0.0,
    0.0, 0.4686, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.57274, 0.0, 0.0, 0.052067, 0.0,
    0.0, 0.67687, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.78101, 0.0, 0.0, 0.052067, 0.0,
    0.0, 0.88514, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.052067, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.1562, 0.0, 0.0, 0.052067, 0.0,
    0.0, 0.26033, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.36447, 0.0, 0.0, 0.052067, 0.0,
    0.0, 0.4686, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.57274, 0.0, 0.0, 0.052067, 0.0,
    0.0, 0.67687, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.78101, 0.0, 0.0, 0.052067, 0.0,
    0.0, 0.88514, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.1562, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.26033, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.36447, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.4686, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.57274, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.67687, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.78101, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.052067, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.1562, 0.0, 0.0, 0.052067,
    0.0, 0.0, 0.26033, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.36447, 0.0, 0.0, 0.052067,
    0.0, 0.0, 0.4686, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.57274, 0.0, 0.0, 0.052067,
    0.0, 0.0, 0.67687, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.78101, 0.0, 0.0, 0.052067,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.1562, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.26033, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.36447, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.4686, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.57274, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.67687, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.78101, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.052067, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.1562, 0.0, 0.0, 0.052067, 0.0,
    0.0, 0.26033, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.36447, 0.0, 0.0, 0.052067, 0.0,
    0.0, 0.4686, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.57274, 0.0, 0.0, 0.052067, 0.0,
    0.0, 0.67687, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.052067,
    0.0, 0.0, 0.052067, 0.0, 0.0, 0.1562, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.26033,
    0.0, 0.0, 0.052067, 0.0, 0.0, 0.36447, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.4686,
    0.0, 0.0, 0.052067, 0.0, 0.0, 0.57274, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.67687,
    0.0, 0.0, 0.052067, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.1562, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.26033, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.36447, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.4686, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.57274, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.67687, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.1562, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.26033, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.36447, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.4686, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.57274, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.1562, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.26033, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.36447, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.4686, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.57274, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.1562, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.26033, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.36447, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.4686, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.57274, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.052067, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.1562, 0.0, 0.0, 0.052067,
    0.0, 0.0, 0.26033, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.36447, 0.0, 0.0, 0.052067,
    0.0, 0.0, 0.4686, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.1562, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.26033, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.36447, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.4686, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.1562, 0.0,
    0.0, 0.052067, 0.0, 0.0, 0.26033, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.36447, 0.0,
    0.0, 0.052067, 0.0, 0.0, 0.4686, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.052067, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.1562, 0.0, 0.0, 0.052067, 0.0,
    0.0, 0.26033, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.36447, 0.0, 0.0, 0.052067, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.052067, 0.0,
    0.0, 0.1562, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.26033, 0.0, 0.0, 0.052067, 0.0,
    0.0, 0.36447, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.052067, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.1562, 0.0, 0.0, 0.052067, 0.0,
    0.0, 0.26033, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.36447, 0.0, 0.0, 0.052067, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.1562, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.26033, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.052067, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.1562, 0.0, 0.0, 0.052067,
    0.0, 0.0, 0.26033, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.052067, 0.0,
    0.0, 0.1562, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.26033, 0.0, 0.0, 0.052067, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.052067, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.1562, 0.0, 0.0, 0.052067,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.1562,
    0.0, 0.0, 0.052067, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.052067, 0.0, 0.0,
    0.052067, 0.0, 0.0, 0.1562, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.052067, 0.0, 0.0, 0.052067 };

  static const real_T b_H[3600] = { 59382.3876, 6.9307E-11, 3.0096E-10,
    37659.217, 6.7726E-11, 2.8846E-10, 35937.5645, 6.6144E-11, 2.7595E-10,
    34218.9483, 6.4563E-11, 2.6344E-10, 32504.8866, 6.2982E-11, 2.5094E-10,
    30796.8974, 6.1401E-11, 2.3843E-10, 29096.4989, 5.982E-11, 2.2593E-10,
    27405.2093, 5.8239E-11, 2.1342E-10, 25724.5467, 5.6658E-11, 2.0091E-10,
    24056.0293, 5.5077E-11, 1.8841E-10, 22401.1752, 5.3496E-11, 1.759E-10,
    20761.5025, 5.1915E-11, 1.6339E-10, 19138.5294, 5.0333E-11, 1.5089E-10,
    17533.7741, 4.8752E-11, 1.3838E-10, 15948.7546, 4.7171E-11, 1.2588E-10,
    14384.9892, 4.559E-11, 1.1337E-10, 12843.9959, 4.4009E-11, 1.0086E-10,
    11327.2929, 4.2428E-11, 8.8357E-11, 9836.3984, 4.0847E-11, 7.5851E-11,
    8372.8305, 3.9266E-11, 6.3344E-11, 6.9307E-11, 59382.3876, -3.231E-11,
    6.7092E-11, 37659.217, -3.2107E-11, 6.4878E-11, 35937.5645, -3.1905E-11,
    6.2663E-11, 34218.9483, -3.1703E-11, 6.0449E-11, 32504.8866, -3.1501E-11,
    5.8235E-11, 30796.8974, -3.1298E-11, 5.602E-11, 29096.4989, -3.1096E-11,
    5.3806E-11, 27405.2093, -3.0894E-11, 5.1591E-11, 25724.5467, -3.0692E-11,
    4.9377E-11, 24056.0293, -3.0489E-11, 4.7162E-11, 22401.1752, -3.0287E-11,
    4.4948E-11, 20761.5025, -3.0085E-11, 4.2734E-11, 19138.5294, -2.9883E-11,
    4.0519E-11, 17533.7741, -2.9681E-11, 3.8305E-11, 15948.7546, -2.9478E-11,
    3.609E-11, 14384.9892, -2.9276E-11, 3.3876E-11, 12843.9959, -2.9074E-11,
    3.1661E-11, 11327.2929, -2.8872E-11, 2.9447E-11, 9836.3984, -2.8669E-11,
    2.7232E-11, 8372.8305, -2.8467E-11, 3.0096E-10, -3.231E-11, 59382.3876,
    2.887E-10, -2.9552E-11, 37659.217, 2.7644E-10, -2.6795E-11, 35937.5645,
    2.6417E-10, -2.4038E-11, 34218.9483, 2.5191E-10, -2.1281E-11, 32504.8866,
    2.3964E-10, -1.8524E-11, 30796.8974, 2.2738E-10, -1.5767E-11, 29096.4989,
    2.1512E-10, -1.301E-11, 27405.2093, 2.0285E-10, -1.0253E-11, 25724.5467,
    1.9059E-10, -7.496E-12, 24056.0293, 1.7833E-10, -4.739E-12, 22401.1752,
    1.6606E-10, -1.9819E-12, 20761.5025, 1.538E-10, 7.7515E-13, 19138.5294,
    1.4153E-10, 3.5322E-12, 17533.7741, 1.2927E-10, 6.2893E-12, 15948.7546,
    1.1701E-10, 9.0463E-12, 14384.9892, 1.0474E-10, 1.1803E-11, 12843.9959,
    9.2479E-11, 1.456E-11, 11327.2929, 8.0215E-11, 1.7317E-11, 9836.3984,
    6.7951E-11, 2.0075E-11, 8372.8305, 37659.217, 6.7092E-11, 2.887E-10,
    56543.5106, 6.555E-11, 2.7671E-10, 34885.9892, 6.4008E-11, 2.6471E-10,
    33229.986, 6.2466E-11, 2.5272E-10, 31577.0191, 6.0924E-11, 2.4073E-10,
    29928.6066, 5.9382E-11, 2.2874E-10, 28286.2667, 5.784E-11, 2.1674E-10,
    26651.5175, 5.6298E-11, 2.0475E-10, 25025.8771, 5.4756E-11, 1.9276E-10,
    23410.8638, 5.3214E-11, 1.8077E-10, 21807.9957, 5.1672E-11, 1.6877E-10,
    20218.7908, 5.013E-11, 1.5678E-10, 18644.7674, 4.8587E-11, 1.4479E-10,
    17087.4436, 4.7045E-11, 1.3279E-10, 15548.3375, 4.5503E-11, 1.208E-10,
    14028.9673, 4.3961E-11, 1.0881E-10, 12530.8511, 4.2419E-11, 9.6816E-11,
    11055.5071, 4.0877E-11, 8.4823E-11, 9604.4534, 3.9335E-11, 7.2831E-11,
    8179.2081, 3.7793E-11, 6.0838E-11, 6.7726E-11, 37659.217, -2.9552E-11,
    6.555E-11, 56543.5106, -2.941E-11, 6.3375E-11, 34885.9892, -2.9267E-11,
    6.1199E-11, 33229.986, -2.9124E-11, 5.9024E-11, 31577.0191, -2.8981E-11,
    5.6849E-11, 29928.6066, -2.8838E-11, 5.4673E-11, 28286.2667, -2.8695E-11,
    5.2498E-11, 26651.5175, -2.8552E-11, 5.0322E-11, 25025.8771, -2.8409E-11,
    4.8147E-11, 23410.8638, -2.8267E-11, 4.5972E-11, 21807.9957, -2.8124E-11,
    4.3796E-11, 20218.7908, -2.7981E-11, 4.1621E-11, 18644.7674, -2.7838E-11,
    3.9445E-11, 17087.4436, -2.7695E-11, 3.727E-11, 15548.3375, -2.7552E-11,
    3.5095E-11, 14028.9673, -2.7409E-11, 3.2919E-11, 12530.8511, -2.7266E-11,
    3.0744E-11, 11055.5071, -2.7124E-11, 2.8569E-11, 9604.4534, -2.6981E-11,
    2.6393E-11, 8179.2081, -2.6838E-11, 2.8846E-10, -3.2107E-11, 37659.217,
    2.7671E-10, -2.941E-11, 56543.5106, 2.6496E-10, -2.6712E-11, 34885.9892,
    2.5321E-10, -2.4014E-11, 33229.986, 2.4146E-10, -2.1316E-11, 31577.0191,
    2.2971E-10, -1.8619E-11, 29928.6066, 2.1796E-10, -1.5921E-11, 28286.2667,
    2.0621E-10, -1.3223E-11, 26651.5175, 1.9446E-10, -1.0526E-11, 25025.8771,
    1.827E-10, -7.8279E-12, 23410.8638, 1.7095E-10, -5.1302E-12, 21807.9957,
    1.592E-10, -2.4325E-12, 20218.7908, 1.4745E-10, 2.6521E-13, 18644.7674,
    1.357E-10, 2.9629E-12, 17087.4436, 1.2395E-10, 5.6606E-12, 15548.3375,
    1.122E-10, 8.3583E-12, 14028.9673, 1.0045E-10, 1.1056E-11, 12530.8511,
    8.8703E-11, 1.3754E-11, 11055.5071, 7.6953E-11, 1.6451E-11, 9604.4534,
    6.5202E-11, 1.9149E-11, 8179.2081, 35937.5645, 6.4878E-11, 2.7644E-10,
    34885.9892, 6.3375E-11, 2.6496E-10, 53834.4139, 6.1872E-11, 2.5348E-10,
    32241.0237, 6.0369E-11, 2.42E-10, 30649.1516, 5.8866E-11, 2.3052E-10,
    29060.3158, 5.7363E-11, 2.1904E-10, 27476.0344, 5.586E-11, 2.0756E-10,
    25897.8256, 5.4357E-11, 1.9608E-10, 24327.2075, 5.2854E-11, 1.846E-10,
    22765.6983, 5.1351E-11, 1.7312E-10, 21214.8161, 4.9848E-11, 1.6164E-10,
    19676.0791, 4.8345E-11, 1.5016E-10, 18151.0053, 4.6841E-11, 1.3869E-10,
    16641.113, 4.5338E-11, 1.2721E-10, 15147.9203, 4.3835E-11, 1.1573E-10,
    13672.9454, 4.2332E-11, 1.0425E-10, 12217.7063, 4.0829E-11, 9.2769E-11,
    10783.7212, 3.9326E-11, 8.129E-11, 9372.5083, 3.7823E-11, 6.9811E-11,
    7985.5857, 3.632E-11, 5.8331E-11, 6.6144E-11, 35937.5645, -2.6795E-11,
    6.4008E-11, 34885.9892, -2.6712E-11, 6.1872E-11, 53834.4139, -2.6628E-11,
    5.9735E-11, 32241.0237, -2.6545E-11, 5.7599E-11, 30649.1516, -2.6461E-11,
    5.5463E-11, 29060.3158, -2.6378E-11, 5.3326E-11, 27476.0344, -2.6294E-11,
    5.119E-11, 25897.8256, -2.6211E-11, 4.9054E-11, 24327.2075, -2.6127E-11,
    4.6917E-11, 22765.6983, -2.6044E-11, 4.4781E-11, 21214.8161, -2.596E-11,
    4.2645E-11, 19676.0791, -2.5877E-11, 4.0508E-11, 18151.0053, -2.5793E-11,
    3.8372E-11, 16641.113, -2.571E-11, 3.6236E-11, 15147.9203, -2.5626E-11,
    3.4099E-11, 13672.9454, -2.5542E-11, 3.1963E-11, 12217.7063, -2.5459E-11,
    2.9826E-11, 10783.7212, -2.5375E-11, 2.769E-11, 9372.5083, -2.5292E-11,
    2.5554E-11, 7985.5857, -2.5208E-11, 2.7595E-10, -3.1905E-11, 35937.5645,
    2.6471E-10, -2.9267E-11, 34885.9892, 2.5348E-10, -2.6628E-11, 53834.4139,
    2.4224E-10, -2.399E-11, 32241.0237, 2.31E-10, -2.1352E-11, 30649.1516,
    2.1977E-10, -1.8713E-11, 29060.3158, 2.0853E-10, -1.6075E-11, 27476.0344,
    1.9729E-10, -1.3437E-11, 25897.8256, 1.8606E-10, -1.0798E-11, 24327.2075,
    1.7482E-10, -8.1598E-12, 22765.6983, 1.6358E-10, -5.5214E-12, 21214.8161,
    1.5235E-10, -2.8831E-12, 19676.0791, 1.4111E-10, -2.4472E-13, 18151.0053,
    1.2987E-10, 2.3936E-12, 16641.113, 1.1864E-10, 5.032E-12, 15147.9203,
    1.074E-10, 7.6704E-12, 13672.9454, 9.6163E-11, 1.0309E-11, 12217.7063,
    8.4927E-11, 1.2947E-11, 10783.7212, 7.369E-11, 1.5585E-11, 9372.5083,
    6.2453E-11, 1.8224E-11, 7985.5857, 34218.9483, 6.2663E-11, 2.6417E-10,
    33229.986, 6.1199E-11, 2.5321E-10, 32241.0237, 5.9735E-11, 2.4224E-10,
    51252.0614, 5.8271E-11, 2.3128E-10, 29721.2841, 5.6807E-11, 2.2031E-10,
    28192.025, 5.5343E-11, 2.0934E-10, 26665.8021, 5.3879E-11, 1.9838E-10,
    25144.1337, 5.2415E-11, 1.8741E-10, 23628.5379, 5.0951E-11, 1.7645E-10,
    22120.5328, 4.9487E-11, 1.6548E-10, 20621.6366, 4.8023E-11, 1.5452E-10,
    19133.3673, 4.6559E-11, 1.4355E-10, 17657.2433, 4.5096E-11, 1.3258E-10,
    16194.7825, 4.3632E-11, 1.2162E-10, 14747.5032, 4.2168E-11, 1.1065E-10,
    13316.9235, 4.0704E-11, 9.9687E-11, 11904.5615, 3.924E-11, 8.8722E-11,
    10511.9354, 3.7776E-11, 7.7756E-11, 9140.5633, 3.6312E-11, 6.6791E-11,
    7791.9633, 3.4848E-11, 5.5825E-11, 6.4563E-11, 34218.9483, -2.4038E-11,
    6.2466E-11, 33229.986, -2.4014E-11, 6.0369E-11, 32241.0237, -2.399E-11,
    5.8271E-11, 51252.0614, -2.3966E-11, 5.6174E-11, 29721.2841, -2.3942E-11,
    5.4077E-11, 28192.025, -2.3917E-11, 5.1979E-11, 26665.8021, -2.3893E-11,
    4.9882E-11, 25144.1337, -2.3869E-11, 4.7785E-11, 23628.5379, -2.3845E-11,
    4.5688E-11, 22120.5328, -2.3821E-11, 4.359E-11, 20621.6366, -2.3797E-11,
    4.1493E-11, 19133.3673, -2.3772E-11, 3.9396E-11, 17657.2433, -2.3748E-11,
    3.7298E-11, 16194.7825, -2.3724E-11, 3.5201E-11, 14747.5032, -2.37E-11,
    3.3104E-11, 13316.9235, -2.3676E-11, 3.1006E-11, 11904.5615, -2.3651E-11,
    2.8909E-11, 10511.9354, -2.3627E-11, 2.6812E-11, 9140.5633, -2.3603E-11,
    2.4714E-11, 7791.9633, -2.3579E-11, 2.6344E-10, -3.1703E-11, 34218.9483,
    2.5272E-10, -2.9124E-11, 33229.986, 2.42E-10, -2.6545E-11, 32241.0237,
    2.3128E-10, -2.3966E-11, 51252.0614, 2.2055E-10, -2.1387E-11, 29721.2841,
    2.0983E-10, -1.8808E-11, 28192.025, 1.9911E-10, -1.6229E-11, 26665.8021,
    1.8838E-10, -1.365E-11, 25144.1337, 1.7766E-10, -1.1071E-11, 23628.5379,
    1.6694E-10, -8.4917E-12, 22120.5328, 1.5621E-10, -5.9127E-12, 20621.6366,
    1.4549E-10, -3.3337E-12, 19133.3673, 1.3477E-10, -7.5465E-13, 17657.2433,
    1.2404E-10, 1.8244E-12, 16194.7825, 1.1332E-10, 4.4034E-12, 14747.5032,
    1.026E-10, 6.9824E-12, 13316.9235, 9.1874E-11, 9.5614E-12, 11904.5615,
    8.1151E-11, 1.214E-11, 10511.9354, 7.0428E-11, 1.4719E-11, 9140.5633,
    5.9704E-11, 1.7298E-11, 7791.9633, 32504.8866, 6.0449E-11, 2.5191E-10,
    31577.0191, 5.9024E-11, 2.4146E-10, 30649.1516, 5.7599E-11, 2.31E-10,
    29721.2841, 5.6174E-11, 2.2055E-10, 48793.4166, 5.4749E-11, 2.101E-10,
    27323.7342, 5.3324E-11, 1.9965E-10, 25855.5699, 5.1899E-11, 1.892E-10,
    24390.4419, 5.0474E-11, 1.7874E-10, 22929.8683, 4.9049E-11, 1.6829E-10,
    21475.3673, 4.7624E-11, 1.5784E-10, 20028.457, 4.6199E-11, 1.4739E-10,
    18590.6556, 4.4774E-11, 1.3694E-10, 17163.4812, 4.335E-11, 1.2648E-10,
    15748.452, 4.1925E-11, 1.1603E-10, 14347.0861, 4.05E-11, 1.0558E-10,
    12960.9016, 3.9075E-11, 9.5127E-11, 11591.4167, 3.765E-11, 8.4675E-11,
    10240.1495, 3.6225E-11, 7.4223E-11, 8908.6182, 3.48E-11, 6.3771E-11,
    7598.341, 3.3375E-11, 5.3318E-11, 6.2982E-11, 32504.8866, -2.1281E-11,
    6.0924E-11, 31577.0191, -2.1316E-11, 5.8866E-11, 30649.1516, -2.1352E-11,
    5.6807E-11, 29721.2841, -2.1387E-11, 5.4749E-11, 48793.4166, -2.1422E-11,
    5.2691E-11, 27323.7342, -2.1457E-11, 5.0633E-11, 25855.5699, -2.1492E-11,
    4.8574E-11, 24390.4419, -2.1527E-11, 4.6516E-11, 22929.8683, -2.1563E-11,
    4.4458E-11, 21475.3673, -2.1598E-11, 4.2399E-11, 20028.457, -2.1633E-11,
    4.0341E-11, 18590.6556, -2.1668E-11, 3.8283E-11, 17163.4812, -2.1703E-11,
    3.6225E-11, 15748.452, -2.1738E-11, 3.4166E-11, 14347.0861, -2.1774E-11,
    3.2108E-11, 12960.9016, -2.1809E-11, 3.005E-11, 11591.4167, -2.1844E-11,
    2.7992E-11, 10240.1495, -2.1879E-11, 2.5933E-11, 8908.6182, -2.1914E-11,
    2.3875E-11, 7598.341, -2.1949E-11, 2.5094E-10, -3.1501E-11, 32504.8866,
    2.4073E-10, -2.8981E-11, 31577.0191, 2.3052E-10, -2.6461E-11, 30649.1516,
    2.2031E-10, -2.3942E-11, 29721.2841, 2.101E-10, -2.1422E-11, 48793.4166,
    1.9989E-10, -1.8902E-11, 27323.7342, 1.8968E-10, -1.6383E-11, 25855.5699,
    1.7947E-10, -1.3863E-11, 24390.4419, 1.6926E-10, -1.1343E-11, 22929.8683,
    1.5905E-10, -8.8236E-12, 21475.3673, 1.4884E-10, -6.3039E-12, 20028.457,
    1.3863E-10, -3.7843E-12, 18590.6556, 1.2842E-10, -1.2646E-12, 17163.4812,
    1.1821E-10, 1.2551E-12, 15748.452, 1.08E-10, 3.7748E-12, 14347.0861,
    9.7794E-11, 6.2944E-12, 12960.9016, 8.7584E-11, 8.8141E-12, 11591.4167,
    7.7375E-11, 1.1334E-11, 10240.1495, 6.7165E-11, 1.3853E-11, 8908.6182,
    5.6956E-11, 1.6373E-11, 7598.341, 30796.8974, 5.8235E-11, 2.3964E-10,
    29928.6066, 5.6849E-11, 2.2971E-10, 29060.3158, 5.5463E-11, 2.1977E-10,
    28192.025, 5.4077E-11, 2.0983E-10, 27323.7342, 5.2691E-11, 1.9989E-10,
    46455.4434, 5.1305E-11, 1.8995E-10, 25045.3376, 4.9919E-11, 1.8001E-10,
    23636.75, 4.8533E-11, 1.7007E-10, 22231.1987, 4.7147E-11, 1.6014E-10,
    20830.2018, 4.5761E-11, 1.502E-10, 19435.2775, 4.4375E-11, 1.4026E-10,
    18047.9439, 4.2989E-11, 1.3032E-10, 16669.7192, 4.1604E-11, 1.2038E-10,
    15302.1215, 4.0218E-11, 1.1044E-10, 13946.6689, 3.8832E-11, 1.005E-10,
    12604.8797, 3.7446E-11, 9.0566E-11, 11278.2719, 3.606E-11, 8.0628E-11,
    9968.3636, 3.4674E-11, 7.0689E-11, 8676.6732, 3.3288E-11, 6.0751E-11,
    7404.7186, 3.1902E-11, 5.0812E-11, 6.1401E-11, 30796.8974, -1.8524E-11,
    5.9382E-11, 29928.6066, -1.8619E-11, 5.7363E-11, 29060.3158, -1.8713E-11,
    5.5343E-11, 28192.025, -1.8808E-11, 5.3324E-11, 27323.7342, -1.8902E-11,
    5.1305E-11, 46455.4434, -1.8997E-11, 4.9286E-11, 25045.3376, -1.9091E-11,
    4.7266E-11, 23636.75, -1.9186E-11, 4.5247E-11, 22231.1987, -1.928E-11,
    4.3228E-11, 20830.2018, -1.9375E-11, 4.1209E-11, 19435.2775, -1.9469E-11,
    3.919E-11, 18047.9439, -1.9564E-11, 3.717E-11, 16669.7192, -1.9658E-11,
    3.5151E-11, 15302.1215, -1.9753E-11, 3.3132E-11, 13946.6689, -1.9847E-11,
    3.1113E-11, 12604.8797, -1.9942E-11, 2.9093E-11, 11278.2719, -2.0036E-11,
    2.7074E-11, 9968.3636, -2.0131E-11, 2.5055E-11, 8676.6732, -2.0225E-11,
    2.3036E-11, 7404.7186, -2.032E-11, 2.3843E-10, -3.1298E-11, 30796.8974,
    2.2874E-10, -2.8838E-11, 29928.6066, 2.1904E-10, -2.6378E-11, 29060.3158,
    2.0934E-10, -2.3917E-11, 28192.025, 1.9965E-10, -2.1457E-11, 27323.7342,
    1.8995E-10, -1.8997E-11, 46455.4434, 1.8026E-10, -1.6536E-11, 25045.3376,
    1.7056E-10, -1.4076E-11, 23636.75, 1.6086E-10, -1.1616E-11, 22231.1987,
    1.5117E-10, -9.1555E-12, 20830.2018, 1.4147E-10, -6.6952E-12, 19435.2775,
    1.3178E-10, -4.2348E-12, 18047.9439, 1.2208E-10, -1.7745E-12, 16669.7192,
    1.1238E-10, 6.858E-13, 15302.1215, 1.0269E-10, 3.1461E-12, 13946.6689,
    9.2991E-11, 5.6064E-12, 12604.8797, 8.3295E-11, 8.0668E-12, 11278.2719,
    7.3599E-11, 1.0527E-11, 9968.3636, 6.3903E-11, 1.2987E-11, 8676.6732,
    5.4207E-11, 1.5448E-11, 7404.7186, 29096.4989, 5.602E-11, 2.2738E-10,
    28286.2667, 5.4673E-11, 2.1796E-10, 27476.0344, 5.3326E-11, 2.0853E-10,
    26665.8021, 5.1979E-11, 1.9911E-10, 25855.5699, 5.0633E-11, 1.8968E-10,
    25045.3376, 4.9286E-11, 1.8026E-10, 44235.1054, 4.7939E-11, 1.7083E-10,
    22883.0582, 4.6592E-11, 1.6141E-10, 21532.5291, 4.5245E-11, 1.5198E-10,
    20185.0363, 4.3898E-11, 1.4256E-10, 18842.0979, 4.2551E-11, 1.3313E-10,
    17505.2322, 4.1204E-11, 1.2371E-10, 16175.9571, 3.9858E-11, 1.1428E-10,
    14855.7909, 3.8511E-11, 1.0486E-10, 13546.2518, 3.7164E-11, 9.5431E-11,
    12248.8578, 3.5817E-11, 8.6006E-11, 10965.1271, 3.447E-11, 7.6581E-11,
    9696.5778, 3.3123E-11, 6.7156E-11, 8444.7281, 3.1776E-11, 5.7731E-11,
    7211.0962, 3.0429E-11, 4.8306E-11, 5.982E-11, 29096.4989, -1.5767E-11,
    5.784E-11, 28286.2667, -1.5921E-11, 5.586E-11, 27476.0344, -1.6075E-11,
    5.3879E-11, 26665.8021, -1.6229E-11, 5.1899E-11, 25855.5699, -1.6383E-11,
    4.9919E-11, 25045.3376, -1.6536E-11, 4.7939E-11, 44235.1054, -1.669E-11,
    4.5959E-11, 22883.0582, -1.6844E-11, 4.3978E-11, 21532.5291, -1.6998E-11,
    4.1998E-11, 20185.0363, -1.7152E-11, 4.0018E-11, 18842.0979, -1.7306E-11,
    3.8038E-11, 17505.2322, -1.746E-11, 3.6058E-11, 16175.9571, -1.7613E-11,
    3.4077E-11, 14855.7909, -1.7767E-11, 3.2097E-11, 13546.2518, -1.7921E-11,
    3.0117E-11, 12248.8578, -1.8075E-11, 2.8137E-11, 10965.1271, -1.8229E-11,
    2.6157E-11, 9696.5778, -1.8383E-11, 2.4176E-11, 8444.7281, -1.8537E-11,
    2.2196E-11, 7211.0962, -1.869E-11, 2.2593E-10, -3.1096E-11, 29096.4989,
    2.1674E-10, -2.8695E-11, 28286.2667, 2.0756E-10, -2.6294E-11, 27476.0344,
    1.9838E-10, -2.3893E-11, 26665.8021, 1.892E-10, -2.1492E-11, 25855.5699,
    1.8001E-10, -1.9091E-11, 25045.3376, 1.7083E-10, -1.669E-11, 44235.1054,
    1.6165E-10, -1.4289E-11, 22883.0582, 1.5247E-10, -1.1888E-11, 21532.5291,
    1.4328E-10, -9.4874E-12, 20185.0363, 1.341E-10, -7.0864E-12, 18842.0979,
    1.2492E-10, -4.6854E-12, 17505.2322, 1.1574E-10, -2.2845E-12, 16175.9571,
    1.0655E-10, 1.1652E-13, 14855.7909, 9.737E-11, 2.5175E-12, 13546.2518,
    8.8188E-11, 4.9185E-12, 12248.8578, 7.9005E-11, 7.3194E-12, 10965.1271,
    6.9823E-11, 9.7204E-12, 9696.5778, 6.064E-11, 1.2121E-11, 8444.7281,
    5.1458E-11, 1.4522E-11, 7211.0962, 27405.2093, 5.3806E-11, 2.1512E-10,
    26651.5175, 5.2498E-11, 2.0621E-10, 25897.8256, 5.119E-11, 1.9729E-10,
    25144.1337, 4.9882E-11, 1.8838E-10, 24390.4419, 4.8574E-11, 1.7947E-10,
    23636.75, 4.7266E-11, 1.7056E-10, 22883.0582, 4.5959E-11, 1.6165E-10,
    42129.3663, 4.4651E-11, 1.5274E-10, 20833.8595, 4.3343E-11, 1.4383E-10,
    19539.8708, 4.2035E-11, 1.3491E-10, 18248.9184, 4.0727E-11, 1.26E-10,
    16962.5205, 3.9419E-11, 1.1709E-10, 15682.1951, 3.8112E-11, 1.0818E-10,
    14409.4604, 3.6804E-11, 9.9268E-11, 13145.8346, 3.5496E-11, 9.0356E-11,
    11892.8359, 3.4188E-11, 8.1445E-11, 10651.9823, 3.288E-11, 7.2534E-11,
    9424.7919, 3.1572E-11, 6.3622E-11, 8212.7831, 3.0265E-11, 5.4711E-11,
    7017.4738, 2.8957E-11, 4.5799E-11, 5.8239E-11, 27405.2093, -1.301E-11,
    5.6298E-11, 26651.5175, -1.3223E-11, 5.4357E-11, 25897.8256, -1.3437E-11,
    5.2415E-11, 25144.1337, -1.365E-11, 5.0474E-11, 24390.4419, -1.3863E-11,
    4.8533E-11, 23636.75, -1.4076E-11, 4.6592E-11, 22883.0582, -1.4289E-11,
    4.4651E-11, 42129.3663, -1.4503E-11, 4.271E-11, 20833.8595, -1.4716E-11,
    4.0768E-11, 19539.8708, -1.4929E-11, 3.8827E-11, 18248.9184, -1.5142E-11,
    3.6886E-11, 16962.5205, -1.5355E-11, 3.4945E-11, 15682.1951, -1.5569E-11,
    3.3004E-11, 14409.4604, -1.5782E-11, 3.1063E-11, 13145.8346, -1.5995E-11,
    2.9121E-11, 11892.8359, -1.6208E-11, 2.718E-11, 10651.9823, -1.6421E-11,
    2.5239E-11, 9424.7919, -1.6635E-11, 2.3298E-11, 8212.7831, -1.6848E-11,
    2.1357E-11, 7017.4738, -1.7061E-11, 2.1342E-10, -3.0894E-11, 27405.2093,
    2.0475E-10, -2.8552E-11, 26651.5175, 1.9608E-10, -2.6211E-11, 25897.8256,
    1.8741E-10, -2.3869E-11, 25144.1337, 1.7874E-10, -2.1527E-11, 24390.4419,
    1.7007E-10, -1.9186E-11, 23636.75, 1.6141E-10, -1.6844E-11, 22883.0582,
    1.5274E-10, -1.4503E-11, 42129.3663, 1.4407E-10, -1.2161E-11, 20833.8595,
    1.354E-10, -9.8193E-12, 19539.8708, 1.2673E-10, -7.4776E-12, 18248.9184,
    1.1806E-10, -5.136E-12, 16962.5205, 1.0939E-10, -2.7944E-12, 15682.1951,
    1.0072E-10, -4.5276E-13, 14409.4604, 9.2054E-11, 1.8889E-12, 13145.8346,
    8.3385E-11, 4.2305E-12, 11892.8359, 7.4716E-11, 6.5721E-12, 10651.9823,
    6.6047E-11, 8.9138E-12, 9424.7919, 5.7378E-11, 1.1255E-11, 8212.7831,
    4.8709E-11, 1.3597E-11, 7017.4738, 25724.5467, 5.1591E-11, 2.0285E-10,
    25025.8771, 5.0322E-11, 1.9446E-10, 24327.2075, 4.9054E-11, 1.8606E-10,
    23628.5379, 4.7785E-11, 1.7766E-10, 22929.8683, 4.6516E-11, 1.6926E-10,
    22231.1987, 4.5247E-11, 1.6086E-10, 21532.5291, 4.3978E-11, 1.5247E-10,
    20833.8595, 4.271E-11, 1.4407E-10, 40135.1899, 4.1441E-11, 1.3567E-10,
    18894.7053, 4.0172E-11, 1.2727E-10, 17655.7389, 3.8903E-11, 1.1887E-10,
    16419.8087, 3.7634E-11, 1.1048E-10, 15188.433, 3.6366E-11, 1.0208E-10,
    13963.1299, 3.5097E-11, 9.368E-11, 12745.4175, 3.3828E-11, 8.5282E-11,
    11536.814, 3.2559E-11, 7.6884E-11, 10338.8374, 3.129E-11, 6.8486E-11,
    9153.0061, 3.0022E-11, 6.0089E-11, 7980.838, 2.8753E-11, 5.1691E-11,
    6823.8514, 2.7484E-11, 4.3293E-11, 5.6658E-11, 25724.5467, -1.0253E-11,
    5.4756E-11, 25025.8771, -1.0526E-11, 5.2854E-11, 24327.2075, -1.0798E-11,
    5.0951E-11, 23628.5379, -1.1071E-11, 4.9049E-11, 22929.8683, -1.1343E-11,
    4.7147E-11, 22231.1987, -1.1616E-11, 4.5245E-11, 21532.5291, -1.1888E-11,
    4.3343E-11, 20833.8595, -1.2161E-11, 4.1441E-11, 40135.1899, -1.2433E-11,
    3.9539E-11, 18894.7053, -1.2706E-11, 3.7637E-11, 17655.7389, -1.2979E-11,
    3.5734E-11, 16419.8087, -1.3251E-11, 3.3832E-11, 15188.433, -1.3524E-11,
    3.193E-11, 13963.1299, -1.3796E-11, 3.0028E-11, 12745.4175, -1.4069E-11,
    2.8126E-11, 11536.814, -1.4341E-11, 2.6224E-11, 10338.8374, -1.4614E-11,
    2.4322E-11, 9153.0061, -1.4886E-11, 2.242E-11, 7980.838, -1.5159E-11,
    2.0517E-11, 6823.8514, -1.5431E-11, 2.0091E-10, -3.0692E-11, 25724.5467,
    1.9276E-10, -2.8409E-11, 25025.8771, 1.846E-10, -2.6127E-11, 24327.2075,
    1.7645E-10, -2.3845E-11, 23628.5379, 1.6829E-10, -2.1563E-11, 22929.8683,
    1.6014E-10, -1.928E-11, 22231.1987, 1.5198E-10, -1.6998E-11, 21532.5291,
    1.4383E-10, -1.4716E-11, 20833.8595, 1.3567E-10, -1.2433E-11, 40135.1899,
    1.2751E-10, -1.0151E-11, 18894.7053, 1.1936E-10, -7.8689E-12, 17655.7389,
    1.112E-10, -5.5866E-12, 16419.8087, 1.0305E-10, -3.3043E-12, 15188.433,
    9.4893E-11, -1.022E-12, 13963.1299, 8.6737E-11, 1.2602E-12, 12745.4175,
    7.8582E-11, 3.5425E-12, 11536.814, 7.0426E-11, 5.8248E-12, 10338.8374,
    6.2271E-11, 8.1071E-12, 9153.0061, 5.4115E-11, 1.0389E-11, 7980.838,
    4.596E-11, 1.2672E-11, 6823.8514, 24056.0293, 4.9377E-11, 1.9059E-10,
    23410.8638, 4.8147E-11, 1.827E-10, 22765.6983, 4.6917E-11, 1.7482E-10,
    22120.5328, 4.5688E-11, 1.6694E-10, 21475.3673, 4.4458E-11, 1.5905E-10,
    20830.2018, 4.3228E-11, 1.5117E-10, 20185.0363, 4.1998E-11, 1.4328E-10,
    19539.8708, 4.0768E-11, 1.354E-10, 18894.7053, 3.9539E-11, 1.2751E-10,
    38249.5398, 3.8309E-11, 1.1963E-10, 17062.5593, 3.7079E-11, 1.1175E-10,
    15877.097, 3.5849E-11, 1.0386E-10, 14694.671, 3.462E-11, 9.5977E-11,
    13516.7994, 3.339E-11, 8.8092E-11, 12345.0004, 3.216E-11, 8.0208E-11,
    11180.7921, 3.093E-11, 7.2324E-11, 10025.6926, 2.9701E-11, 6.4439E-11,
    8881.2202, 2.8471E-11, 5.6555E-11, 7748.893, 2.7241E-11, 4.8671E-11,
    6630.229, 2.6011E-11, 4.0786E-11, 5.5077E-11, 24056.0293, -7.496E-12,
    5.3214E-11, 23410.8638, -7.8279E-12, 5.1351E-11, 22765.6983, -8.1598E-12,
    4.9487E-11, 22120.5328, -8.4917E-12, 4.7624E-11, 21475.3673, -8.8236E-12,
    4.5761E-11, 20830.2018, -9.1555E-12, 4.3898E-11, 20185.0363, -9.4874E-12,
    4.2035E-11, 19539.8708, -9.8193E-12, 4.0172E-11, 18894.7053, -1.0151E-11,
    3.8309E-11, 38249.5398, -1.0483E-11, 3.6446E-11, 17062.5593, -1.0815E-11,
    3.4583E-11, 15877.097, -1.1147E-11, 3.272E-11, 14694.671, -1.1479E-11,
    3.0857E-11, 13516.7994, -1.1811E-11, 2.8994E-11, 12345.0004, -1.2143E-11,
    2.713E-11, 11180.7921, -1.2474E-11, 2.5267E-11, 10025.6926, -1.2806E-11,
    2.3404E-11, 8881.2202, -1.3138E-11, 2.1541E-11, 7748.893, -1.347E-11,
    1.9678E-11, 6630.229, -1.3802E-11, 1.8841E-10, -3.0489E-11, 24056.0293,
    1.8077E-10, -2.8267E-11, 23410.8638, 1.7312E-10, -2.6044E-11, 22765.6983,
    1.6548E-10, -2.3821E-11, 22120.5328, 1.5784E-10, -2.1598E-11, 21475.3673,
    1.502E-10, -1.9375E-11, 20830.2018, 1.4256E-10, -1.7152E-11, 20185.0363,
    1.3491E-10, -1.4929E-11, 19539.8708, 1.2727E-10, -1.2706E-11, 18894.7053,
    1.1963E-10, -1.0483E-11, 38249.5398, 1.1199E-10, -8.2601E-12, 17062.5593,
    1.0435E-10, -6.0372E-12, 15877.097, 9.6704E-11, -3.8143E-12, 14694.671,
    8.9062E-11, -1.5913E-12, 13516.7994, 8.142E-11, 6.3162E-13, 12345.0004,
    7.3779E-11, 2.8546E-12, 11180.7921, 6.6137E-11, 5.0775E-12, 10025.6926,
    5.8495E-11, 7.3004E-12, 8881.2202, 5.0853E-11, 9.5234E-12, 7748.893,
    4.3211E-11, 1.1746E-11, 6630.229, 22401.1752, 4.7162E-11, 1.7833E-10,
    21807.9957, 4.5972E-11, 1.7095E-10, 21214.8161, 4.4781E-11, 1.6358E-10,
    20621.6366, 4.359E-11, 1.5621E-10, 20028.457, 4.2399E-11, 1.4884E-10,
    19435.2775, 4.1209E-11, 1.4147E-10, 18842.0979, 4.0018E-11, 1.341E-10,
    18248.9184, 3.8827E-11, 1.2673E-10, 17655.7389, 3.7637E-11, 1.1936E-10,
    17062.5593, 3.6446E-11, 1.1199E-10, 36469.3798, 3.5255E-11, 1.0462E-10,
    15334.3853, 3.4064E-11, 9.7246E-11, 14200.9089, 3.2874E-11, 8.9876E-11,
    13070.4689, 3.1683E-11, 8.2505E-11, 11944.5832, 3.0492E-11, 7.5134E-11,
    10824.7702, 2.9302E-11, 6.7763E-11, 9712.5478, 2.8111E-11, 6.0392E-11,
    8609.4344, 2.692E-11, 5.3021E-11, 7516.9479, 2.5729E-11, 4.5651E-11,
    6436.6066, 2.4539E-11, 3.828E-11, 5.3496E-11, 22401.1752, -4.739E-12,
    5.1672E-11, 21807.9957, -5.1302E-12, 4.9848E-11, 21214.8161, -5.5214E-12,
    4.8023E-11, 20621.6366, -5.9127E-12, 4.6199E-11, 20028.457, -6.3039E-12,
    4.4375E-11, 19435.2775, -6.6952E-12, 4.2551E-11, 18842.0979, -7.0864E-12,
    4.0727E-11, 18248.9184, -7.4776E-12, 3.8903E-11, 17655.7389, -7.8689E-12,
    3.7079E-11, 17062.5593, -8.2601E-12, 3.5255E-11, 36469.3798, -8.6514E-12,
    3.3431E-11, 15334.3853, -9.0426E-12, 3.1607E-11, 14200.9089, -9.4339E-12,
    2.9783E-11, 13070.4689, -9.8251E-12, 2.7959E-11, 11944.5832, -1.0216E-11,
    2.6135E-11, 10824.7702, -1.0608E-11, 2.4311E-11, 9712.5478, -1.0999E-11,
    2.2487E-11, 8609.4344, -1.139E-11, 2.0663E-11, 7516.9479, -1.1781E-11,
    1.8839E-11, 6436.6066, -1.2173E-11, 1.759E-10, -3.0287E-11, 22401.1752,
    1.6877E-10, -2.8124E-11, 21807.9957, 1.6164E-10, -2.596E-11, 21214.8161,
    1.5452E-10, -2.3797E-11, 20621.6366, 1.4739E-10, -2.1633E-11, 20028.457,
    1.4026E-10, -1.9469E-11, 19435.2775, 1.3313E-10, -1.7306E-11, 18842.0979,
    1.26E-10, -1.5142E-11, 18248.9184, 1.1887E-10, -1.2979E-11, 17655.7389,
    1.1175E-10, -1.0815E-11, 17062.5593, 1.0462E-10, -8.6514E-12, 36469.3798,
    9.7489E-11, -6.4878E-12, 15334.3853, 9.036E-11, -4.3242E-12, 14200.9089,
    8.3232E-11, -2.1606E-12, 13070.4689, 7.6104E-11, 2.9897E-15, 11944.5832,
    6.8975E-11, 2.1666E-12, 10824.7702, 6.1847E-11, 4.3302E-12, 9712.5478,
    5.4719E-11, 6.4938E-12, 8609.4344, 4.759E-11, 8.6574E-12, 7516.9479,
    4.0462E-11, 1.0821E-11, 6436.6066, 20761.5025, 4.4948E-11, 1.6606E-10,
    20218.7908, 4.3796E-11, 1.592E-10, 19676.0791, 4.2645E-11, 1.5235E-10,
    19133.3673, 4.1493E-11, 1.4549E-10, 18590.6556, 4.0341E-11, 1.3863E-10,
    18047.9439, 3.919E-11, 1.3178E-10, 17505.2322, 3.8038E-11, 1.2492E-10,
    16962.5205, 3.6886E-11, 1.1806E-10, 16419.8087, 3.5734E-11, 1.112E-10,
    15877.097, 3.4583E-11, 1.0435E-10, 15334.3853, 3.3431E-11, 9.7489E-11,
    34791.6736, 3.2279E-11, 9.0632E-11, 13707.1469, 3.1128E-11, 8.3774E-11,
    12624.1383, 2.9976E-11, 7.6917E-11, 11544.1661, 2.8824E-11, 7.006E-11,
    10468.7483, 2.7673E-11, 6.3202E-11, 9399.403, 2.6521E-11, 5.6345E-11,
    8337.6485, 2.5369E-11, 4.9488E-11, 7285.0029, 2.4218E-11, 4.2631E-11,
    6242.9842, 2.3066E-11, 3.5773E-11, 5.1915E-11, 20761.5025, -1.9819E-12,
    5.013E-11, 20218.7908, -2.4325E-12, 4.8345E-11, 19676.0791, -2.8831E-12,
    4.6559E-11, 19133.3673, -3.3337E-12, 4.4774E-11, 18590.6556, -3.7843E-12,
    4.2989E-11, 18047.9439, -4.2348E-12, 4.1204E-11, 17505.2322, -4.6854E-12,
    3.9419E-11, 16962.5205, -5.136E-12, 3.7634E-11, 16419.8087, -5.5866E-12,
    3.5849E-11, 15877.097, -6.0372E-12, 3.4064E-11, 15334.3853, -6.4878E-12,
    3.2279E-11, 34791.6736, -6.9384E-12, 3.0494E-11, 13707.1469, -7.389E-12,
    2.8709E-11, 12624.1383, -7.8395E-12, 2.6924E-11, 11544.1661, -8.2901E-12,
    2.5139E-11, 10468.7483, -8.7407E-12, 2.3354E-11, 9399.403, -9.1913E-12,
    2.1569E-11, 8337.6485, -9.6419E-12, 1.9784E-11, 7285.0029, -1.0092E-11,
    1.7999E-11, 6242.9842, -1.0543E-11, 1.6339E-10, -3.0085E-11, 20761.5025,
    1.5678E-10, -2.7981E-11, 20218.7908, 1.5016E-10, -2.5877E-11, 19676.0791,
    1.4355E-10, -2.3772E-11, 19133.3673, 1.3694E-10, -2.1668E-11, 18590.6556,
    1.3032E-10, -1.9564E-11, 18047.9439, 1.2371E-10, -1.746E-11, 17505.2322,
    1.1709E-10, -1.5355E-11, 16962.5205, 1.1048E-10, -1.3251E-11, 16419.8087,
    1.0386E-10, -1.1147E-11, 15877.097, 9.7246E-11, -9.0426E-12, 15334.3853,
    9.0632E-11, -6.9384E-12, 34791.6736, 8.4017E-11, -4.8341E-12, 13707.1469,
    7.7402E-11, -2.7299E-12, 12624.1383, 7.0787E-11, -6.2564E-13, 11544.1661,
    6.4172E-11, 1.4786E-12, 10468.7483, 5.7558E-11, 3.5829E-12, 9399.403,
    5.0943E-11, 5.6871E-12, 8337.6485, 4.4328E-11, 7.7913E-12, 7285.0029,
    3.7713E-11, 9.8956E-12, 6242.9842, 19138.5294, 4.2734E-11, 1.538E-10,
    18644.7674, 4.1621E-11, 1.4745E-10, 18151.0053, 4.0508E-11, 1.4111E-10,
    17657.2433, 3.9396E-11, 1.3477E-10, 17163.4812, 3.8283E-11, 1.2842E-10,
    16669.7192, 3.717E-11, 1.2208E-10, 16175.9571, 3.6058E-11, 1.1574E-10,
    15682.1951, 3.4945E-11, 1.0939E-10, 15188.433, 3.3832E-11, 1.0305E-10,
    14694.671, 3.272E-11, 9.6704E-11, 14200.9089, 3.1607E-11, 9.036E-11,
    13707.1469, 3.0494E-11, 8.4017E-11, 33213.3848, 2.9382E-11, 7.7673E-11,
    12177.8078, 2.8269E-11, 7.1329E-11, 11143.7489, 2.7156E-11, 6.4986E-11,
    10112.7264, 2.6044E-11, 5.8642E-11, 9086.2582, 2.4931E-11, 5.2298E-11,
    8065.8627, 2.3818E-11, 4.5954E-11, 7053.0578, 2.2706E-11, 3.9611E-11,
    6049.3618, 2.1593E-11, 3.3267E-11, 5.0333E-11, 19138.5294, 7.7515E-13,
    4.8587E-11, 18644.7674, 2.6521E-13, 4.6841E-11, 18151.0053, -2.4472E-13,
    4.5096E-11, 17657.2433, -7.5465E-13, 4.335E-11, 17163.4812, -1.2646E-12,
    4.1604E-11, 16669.7192, -1.7745E-12, 3.9858E-11, 16175.9571, -2.2845E-12,
    3.8112E-11, 15682.1951, -2.7944E-12, 3.6366E-11, 15188.433, -3.3043E-12,
    3.462E-11, 14694.671, -3.8143E-12, 3.2874E-11, 14200.9089, -4.3242E-12,
    3.1128E-11, 13707.1469, -4.8341E-12, 2.9382E-11, 33213.3848, -5.3441E-12,
    2.7636E-11, 12177.8078, -5.854E-12, 2.589E-11, 11143.7489, -6.3639E-12,
    2.4144E-11, 10112.7264, -6.8739E-12, 2.2398E-11, 9086.2582, -7.3838E-12,
    2.0652E-11, 8065.8627, -7.8937E-12, 1.8906E-11, 7053.0578, -8.4037E-12,
    1.716E-11, 6049.3618, -8.9136E-12, 1.5089E-10, -2.9883E-11, 19138.5294,
    1.4479E-10, -2.7838E-11, 18644.7674, 1.3869E-10, -2.5793E-11, 18151.0053,
    1.3258E-10, -2.3748E-11, 17657.2433, 1.2648E-10, -2.1703E-11, 17163.4812,
    1.2038E-10, -1.9658E-11, 16669.7192, 1.1428E-10, -1.7613E-11, 16175.9571,
    1.0818E-10, -1.5569E-11, 15682.1951, 1.0208E-10, -1.3524E-11, 15188.433,
    9.5977E-11, -1.1479E-11, 14694.671, 8.9876E-11, -9.4339E-12, 14200.9089,
    8.3774E-11, -7.389E-12, 13707.1469, 7.7673E-11, -5.3441E-12, 33213.3848,
    7.1572E-11, -3.2992E-12, 12177.8078, 6.5471E-11, -1.2543E-12, 11143.7489,
    5.9369E-11, 7.9063E-13, 10112.7264, 5.3268E-11, 2.8355E-12, 9086.2582,
    4.7167E-11, 4.8804E-12, 8065.8627, 4.1065E-11, 6.9253E-12, 7053.0578,
    3.4964E-11, 8.9702E-12, 6049.3618, 17533.7741, 4.0519E-11, 1.4153E-10,
    17087.4436, 3.9445E-11, 1.357E-10, 16641.113, 3.8372E-11, 1.2987E-10,
    16194.7825, 3.7298E-11, 1.2404E-10, 15748.452, 3.6225E-11, 1.1821E-10,
    15302.1215, 3.5151E-11, 1.1238E-10, 14855.7909, 3.4077E-11, 1.0655E-10,
    14409.4604, 3.3004E-11, 1.0072E-10, 13963.1299, 3.193E-11, 9.4893E-11,
    13516.7994, 3.0857E-11, 8.9062E-11, 13070.4689, 2.9783E-11, 8.3232E-11,
    12624.1383, 2.8709E-11, 7.7402E-11, 12177.8078, 2.7636E-11, 7.1572E-11,
    31731.4773, 2.6562E-11, 6.5742E-11, 10743.3318, 2.5489E-11, 5.9911E-11,
    9756.7045, 2.4415E-11, 5.4081E-11, 8773.1134, 2.3341E-11, 4.8251E-11,
    7794.0768, 2.2268E-11, 4.2421E-11, 6821.1128, 2.1194E-11, 3.6591E-11,
    5855.7394, 2.012E-11, 3.076E-11, 4.8752E-11, 17533.7741, 3.5322E-12,
    4.7045E-11, 17087.4436, 2.9629E-12, 4.5338E-11, 16641.113, 2.3936E-12,
    4.3632E-11, 16194.7825, 1.8244E-12, 4.1925E-11, 15748.452, 1.2551E-12,
    4.0218E-11, 15302.1215, 6.858E-13, 3.8511E-11, 14855.7909, 1.1652E-13,
    3.6804E-11, 14409.4604, -4.5276E-13, 3.5097E-11, 13963.1299, -1.022E-12,
    3.339E-11, 13516.7994, -1.5913E-12, 3.1683E-11, 13070.4689, -2.1606E-12,
    2.9976E-11, 12624.1383, -2.7299E-12, 2.8269E-11, 12177.8078, -3.2992E-12,
    2.6562E-11, 31731.4773, -3.8684E-12, 2.4855E-11, 10743.3318, -4.4377E-12,
    2.3148E-11, 9756.7045, -5.007E-12, 2.1441E-11, 8773.1134, -5.5763E-12,
    1.9734E-11, 7794.0768, -6.1456E-12, 1.8027E-11, 6821.1128, -6.7148E-12,
    1.6321E-11, 5855.7394, -7.2841E-12, 1.3838E-10, -2.9681E-11, 17533.7741,
    1.3279E-10, -2.7695E-11, 17087.4436, 1.2721E-10, -2.571E-11, 16641.113,
    1.2162E-10, -2.3724E-11, 16194.7825, 1.1603E-10, -2.1738E-11, 15748.452,
    1.1044E-10, -1.9753E-11, 15302.1215, 1.0486E-10, -1.7767E-11, 14855.7909,
    9.9268E-11, -1.5782E-11, 14409.4604, 9.368E-11, -1.3796E-11, 13963.1299,
    8.8092E-11, -1.1811E-11, 13516.7994, 8.2505E-11, -9.8251E-12, 13070.4689,
    7.6917E-11, -7.8395E-12, 12624.1383, 7.1329E-11, -5.854E-12, 12177.8078,
    6.5742E-11, -3.8684E-12, 31731.4773, 6.0154E-11, -1.8829E-12, 10743.3318,
    5.4566E-11, 1.0266E-13, 9756.7045, 4.8978E-11, 2.0882E-12, 8773.1134,
    4.3391E-11, 4.0738E-12, 7794.0768, 3.7803E-11, 6.0593E-12, 6821.1128,
    3.2215E-11, 8.0449E-12, 5855.7394, 15948.7546, 3.8305E-11, 1.2927E-10,
    15548.3375, 3.727E-11, 1.2395E-10, 15147.9203, 3.6236E-11, 1.1864E-10,
    14747.5032, 3.5201E-11, 1.1332E-10, 14347.0861, 3.4166E-11, 1.08E-10,
    13946.6689, 3.3132E-11, 1.0269E-10, 13546.2518, 3.2097E-11, 9.737E-11,
    13145.8346, 3.1063E-11, 9.2054E-11, 12745.4175, 3.0028E-11, 8.6737E-11,
    12345.0004, 2.8994E-11, 8.142E-11, 11944.5832, 2.7959E-11, 7.6104E-11,
    11544.1661, 2.6924E-11, 7.0787E-11, 11143.7489, 2.589E-11, 6.5471E-11,
    10743.3318, 2.4855E-11, 6.0154E-11, 30342.9147, 2.3821E-11, 5.4837E-11,
    9400.6826, 2.2786E-11, 4.9521E-11, 8459.9686, 2.1751E-11, 4.4204E-11,
    7522.2909, 2.0717E-11, 3.8887E-11, 6589.1677, 1.9682E-11, 3.3571E-11,
    5662.1171, 1.8648E-11, 2.8254E-11, 4.7171E-11, 15948.7546, 6.2893E-12,
    4.5503E-11, 15548.3375, 5.6606E-12, 4.3835E-11, 15147.9203, 5.032E-12,
    4.2168E-11, 14747.5032, 4.4034E-12, 4.05E-11, 14347.0861, 3.7748E-12,
    3.8832E-11, 13946.6689, 3.1461E-12, 3.7164E-11, 13546.2518, 2.5175E-12,
    3.5496E-11, 13145.8346, 1.8889E-12, 3.3828E-11, 12745.4175, 1.2602E-12,
    3.216E-11, 12345.0004, 6.3162E-13, 3.0492E-11, 11944.5832, 2.9897E-15,
    2.8824E-11, 11544.1661, -6.2564E-13, 2.7156E-11, 11143.7489, -1.2543E-12,
    2.5489E-11, 10743.3318, -1.8829E-12, 2.3821E-11, 30342.9147, -2.5115E-12,
    2.2153E-11, 9400.6826, -3.1401E-12, 2.0485E-11, 8459.9686, -3.7688E-12,
    1.8817E-11, 7522.2909, -4.3974E-12, 1.7149E-11, 6589.1677, -5.026E-12,
    1.5481E-11, 5662.1171, -5.6547E-12, 1.2588E-10, -2.9478E-11, 15948.7546,
    1.208E-10, -2.7552E-11, 15548.3375, 1.1573E-10, -2.5626E-11, 15147.9203,
    1.1065E-10, -2.37E-11, 14747.5032, 1.0558E-10, -2.1774E-11, 14347.0861,
    1.005E-10, -1.9847E-11, 13946.6689, 9.5431E-11, -1.7921E-11, 13546.2518,
    9.0356E-11, -1.5995E-11, 13145.8346, 8.5282E-11, -1.4069E-11, 12745.4175,
    8.0208E-11, -1.2143E-11, 12345.0004, 7.5134E-11, -1.0216E-11, 11944.5832,
    7.006E-11, -8.2901E-12, 11544.1661, 6.4986E-11, -6.3639E-12, 11143.7489,
    5.9911E-11, -4.4377E-12, 10743.3318, 5.4837E-11, -2.5115E-12, 30342.9147,
    4.9763E-11, -5.8531E-13, 9400.6826, 4.4689E-11, 1.3409E-12, 8459.9686,
    3.9615E-11, 3.2671E-12, 7522.2909, 3.4541E-11, 5.1933E-12, 6589.1677,
    2.9466E-11, 7.1195E-12, 5662.1171, 14384.9892, 3.609E-11, 1.1701E-10,
    14028.9673, 3.5095E-11, 1.122E-10, 13672.9454, 3.4099E-11, 1.074E-10,
    13316.9235, 3.3104E-11, 1.026E-10, 12960.9016, 3.2108E-11, 9.7794E-11,
    12604.8797, 3.1113E-11, 9.2991E-11, 12248.8578, 3.0117E-11, 8.8188E-11,
    11892.8359, 2.9121E-11, 8.3385E-11, 11536.814, 2.8126E-11, 7.8582E-11,
    11180.7921, 2.713E-11, 7.3779E-11, 10824.7702, 2.6135E-11, 6.8975E-11,
    10468.7483, 2.5139E-11, 6.4172E-11, 10112.7264, 2.4144E-11, 5.9369E-11,
    9756.7045, 2.3148E-11, 5.4566E-11, 9400.6826, 2.2153E-11, 4.9763E-11,
    29044.6607, 2.1157E-11, 4.496E-11, 8146.8238, 2.0162E-11, 4.0157E-11,
    7250.5051, 1.9166E-11, 3.5354E-11, 6357.2227, 1.8171E-11, 3.0551E-11,
    5468.4947, 1.7175E-11, 2.5748E-11, 4.559E-11, 14384.9892, 9.0463E-12,
    4.3961E-11, 14028.9673, 8.3583E-12, 4.2332E-11, 13672.9454, 7.6704E-12,
    4.0704E-11, 13316.9235, 6.9824E-12, 3.9075E-11, 12960.9016, 6.2944E-12,
    3.7446E-11, 12604.8797, 5.6064E-12, 3.5817E-11, 12248.8578, 4.9185E-12,
    3.4188E-11, 11892.8359, 4.2305E-12, 3.2559E-11, 11536.814, 3.5425E-12,
    3.093E-11, 11180.7921, 2.8546E-12, 2.9302E-11, 10824.7702, 2.1666E-12,
    2.7673E-11, 10468.7483, 1.4786E-12, 2.6044E-11, 10112.7264, 7.9063E-13,
    2.4415E-11, 9756.7045, 1.0266E-13, 2.2786E-11, 9400.6826, -5.8531E-13,
    2.1157E-11, 29044.6607, -1.2733E-12, 1.9528E-11, 8146.8238, -1.9613E-12,
    1.7899E-11, 7250.5051, -2.6492E-12, 1.6271E-11, 6357.2227, -3.3372E-12,
    1.4642E-11, 5468.4947, -4.0252E-12, 1.1337E-10, -2.9276E-11, 14384.9892,
    1.0881E-10, -2.7409E-11, 14028.9673, 1.0425E-10, -2.5542E-11, 13672.9454,
    9.9687E-11, -2.3676E-11, 13316.9235, 9.5127E-11, -2.1809E-11, 12960.9016,
    9.0566E-11, -1.9942E-11, 12604.8797, 8.6006E-11, -1.8075E-11, 12248.8578,
    8.1445E-11, -1.6208E-11, 11892.8359, 7.6884E-11, -1.4341E-11, 11536.814,
    7.2324E-11, -1.2474E-11, 11180.7921, 6.7763E-11, -1.0608E-11, 10824.7702,
    6.3202E-11, -8.7407E-12, 10468.7483, 5.8642E-11, -6.8739E-12, 10112.7264,
    5.4081E-11, -5.007E-12, 9756.7045, 4.9521E-11, -3.1401E-12, 9400.6826,
    4.496E-11, -1.2733E-12, 29044.6607, 4.0399E-11, 5.9357E-13, 8146.8238,
    3.5839E-11, 2.4604E-12, 7250.5051, 3.1278E-11, 4.3273E-12, 6357.2227,
    2.6717E-11, 6.1941E-12, 5468.4947, 12843.9959, 3.3876E-11, 1.0474E-10,
    12530.8511, 3.2919E-11, 1.0045E-10, 12217.7063, 3.1963E-11, 9.6163E-11,
    11904.5615, 3.1006E-11, 9.1874E-11, 11591.4167, 3.005E-11, 8.7584E-11,
    11278.2719, 2.9093E-11, 8.3295E-11, 10965.1271, 2.8137E-11, 7.9005E-11,
    10651.9823, 2.718E-11, 7.4716E-11, 10338.8374, 2.6224E-11, 7.0426E-11,
    10025.6926, 2.5267E-11, 6.6137E-11, 9712.5478, 2.4311E-11, 6.1847E-11,
    9399.403, 2.3354E-11, 5.7558E-11, 9086.2582, 2.2398E-11, 5.3268E-11,
    8773.1134, 2.1441E-11, 4.8978E-11, 8459.9686, 2.0485E-11, 4.4689E-11,
    8146.8238, 1.9528E-11, 4.0399E-11, 27833.679, 1.8572E-11, 3.611E-11,
    6978.7192, 1.7615E-11, 3.182E-11, 6125.2776, 1.6659E-11, 2.7531E-11,
    5274.8723, 1.5702E-11, 2.3241E-11, 4.4009E-11, 12843.9959, 1.1803E-11,
    4.2419E-11, 12530.8511, 1.1056E-11, 4.0829E-11, 12217.7063, 1.0309E-11,
    3.924E-11, 11904.5615, 9.5614E-12, 3.765E-11, 11591.4167, 8.8141E-12,
    3.606E-11, 11278.2719, 8.0668E-12, 3.447E-11, 10965.1271, 7.3194E-12,
    3.288E-11, 10651.9823, 6.5721E-12, 3.129E-11, 10338.8374, 5.8248E-12,
    2.9701E-11, 10025.6926, 5.0775E-12, 2.8111E-11, 9712.5478, 4.3302E-12,
    2.6521E-11, 9399.403, 3.5829E-12, 2.4931E-11, 9086.2582, 2.8355E-12,
    2.3341E-11, 8773.1134, 2.0882E-12, 2.1751E-11, 8459.9686, 1.3409E-12,
    2.0162E-11, 8146.8238, 5.9357E-13, 1.8572E-11, 27833.679, -1.5375E-13,
    1.6982E-11, 6978.7192, -9.0107E-13, 1.5392E-11, 6125.2776, -1.6484E-12,
    1.3802E-11, 5274.8723, -2.3957E-12, 1.0086E-10, -2.9074E-11, 12843.9959,
    9.6816E-11, -2.7266E-11, 12530.8511, 9.2769E-11, -2.5459E-11, 12217.7063,
    8.8722E-11, -2.3651E-11, 11904.5615, 8.4675E-11, -2.1844E-11, 11591.4167,
    8.0628E-11, -2.0036E-11, 11278.2719, 7.6581E-11, -1.8229E-11, 10965.1271,
    7.2534E-11, -1.6421E-11, 10651.9823, 6.8486E-11, -1.4614E-11, 10338.8374,
    6.4439E-11, -1.2806E-11, 10025.6926, 6.0392E-11, -1.0999E-11, 9712.5478,
    5.6345E-11, -9.1913E-12, 9399.403, 5.2298E-11, -7.3838E-12, 9086.2582,
    4.8251E-11, -5.5763E-12, 8773.1134, 4.4204E-11, -3.7688E-12, 8459.9686,
    4.0157E-11, -1.9613E-12, 8146.8238, 3.611E-11, -1.5375E-13, 27833.679,
    3.2063E-11, 1.6538E-12, 6978.7192, 2.8016E-11, 3.4613E-12, 6125.2776,
    2.3969E-11, 5.2688E-12, 5274.8723, 11327.2929, 3.1661E-11, 9.2479E-11,
    11055.5071, 3.0744E-11, 8.8703E-11, 10783.7212, 2.9826E-11, 8.4927E-11,
    10511.9354, 2.8909E-11, 8.1151E-11, 10240.1495, 2.7992E-11, 7.7375E-11,
    9968.3636, 2.7074E-11, 7.3599E-11, 9696.5778, 2.6157E-11, 6.9823E-11,
    9424.7919, 2.5239E-11, 6.6047E-11, 9153.0061, 2.4322E-11, 6.2271E-11,
    8881.2202, 2.3404E-11, 5.8495E-11, 8609.4344, 2.2487E-11, 5.4719E-11,
    8337.6485, 2.1569E-11, 5.0943E-11, 8065.8627, 2.0652E-11, 4.7167E-11,
    7794.0768, 1.9734E-11, 4.3391E-11, 7522.2909, 1.8817E-11, 3.9615E-11,
    7250.5051, 1.7899E-11, 3.5839E-11, 6978.7192, 1.6982E-11, 3.2063E-11,
    26706.9334, 1.6065E-11, 2.8287E-11, 5893.3326, 1.5147E-11, 2.4511E-11,
    5081.2499, 1.423E-11, 2.0735E-11, 4.2428E-11, 11327.2929, 1.456E-11,
    4.0877E-11, 11055.5071, 1.3754E-11, 3.9326E-11, 10783.7212, 1.2947E-11,
    3.7776E-11, 10511.9354, 1.214E-11, 3.6225E-11, 10240.1495, 1.1334E-11,
    3.4674E-11, 9968.3636, 1.0527E-11, 3.3123E-11, 9696.5778, 9.7204E-12,
    3.1572E-11, 9424.7919, 8.9138E-12, 3.0022E-11, 9153.0061, 8.1071E-12,
    2.8471E-11, 8881.2202, 7.3004E-12, 2.692E-11, 8609.4344, 6.4938E-12,
    2.5369E-11, 8337.6485, 5.6871E-12, 2.3818E-11, 8065.8627, 4.8804E-12,
    2.2268E-11, 7794.0768, 4.0738E-12, 2.0717E-11, 7522.2909, 3.2671E-12,
    1.9166E-11, 7250.5051, 2.4604E-12, 1.7615E-11, 6978.7192, 1.6538E-12,
    1.6065E-11, 26706.9334, 8.471E-13, 1.4514E-11, 5893.3326, 4.0431E-14,
    1.2963E-11, 5081.2499, -7.6623E-13, 8.8357E-11, -2.8872E-11, 11327.2929,
    8.4823E-11, -2.7124E-11, 11055.5071, 8.129E-11, -2.5375E-11, 10783.7212,
    7.7756E-11, -2.3627E-11, 10511.9354, 7.4223E-11, -2.1879E-11, 10240.1495,
    7.0689E-11, -2.0131E-11, 9968.3636, 6.7156E-11, -1.8383E-11, 9696.5778,
    6.3622E-11, -1.6635E-11, 9424.7919, 6.0089E-11, -1.4886E-11, 9153.0061,
    5.6555E-11, -1.3138E-11, 8881.2202, 5.3021E-11, -1.139E-11, 8609.4344,
    4.9488E-11, -9.6419E-12, 8337.6485, 4.5954E-11, -7.8937E-12, 8065.8627,
    4.2421E-11, -6.1456E-12, 7794.0768, 3.8887E-11, -4.3974E-12, 7522.2909,
    3.5354E-11, -2.6492E-12, 7250.5051, 3.182E-11, -9.0107E-13, 6978.7192,
    2.8287E-11, 8.471E-13, 26706.9334, 2.4753E-11, 2.5953E-12, 5893.3326,
    2.122E-11, 4.3434E-12, 5081.2499, 9836.3984, 2.9447E-11, 8.0215E-11,
    9604.4534, 2.8569E-11, 7.6953E-11, 9372.5083, 2.769E-11, 7.369E-11,
    9140.5633, 2.6812E-11, 7.0428E-11, 8908.6182, 2.5933E-11, 6.7165E-11,
    8676.6732, 2.5055E-11, 6.3903E-11, 8444.7281, 2.4176E-11, 6.064E-11,
    8212.7831, 2.3298E-11, 5.7378E-11, 7980.838, 2.242E-11, 5.4115E-11, 7748.893,
    2.1541E-11, 5.0853E-11, 7516.9479, 2.0663E-11, 4.759E-11, 7285.0029,
    1.9784E-11, 4.4328E-11, 7053.0578, 1.8906E-11, 4.1065E-11, 6821.1128,
    1.8027E-11, 3.7803E-11, 6589.1677, 1.7149E-11, 3.4541E-11, 6357.2227,
    1.6271E-11, 3.1278E-11, 6125.2776, 1.5392E-11, 2.8016E-11, 5893.3326,
    1.4514E-11, 2.4753E-11, 25661.3875, 1.3635E-11, 2.1491E-11, 4887.6275,
    1.2757E-11, 1.8228E-11, 4.0847E-11, 9836.3984, 1.7317E-11, 3.9335E-11,
    9604.4534, 1.6451E-11, 3.7823E-11, 9372.5083, 1.5585E-11, 3.6312E-11,
    9140.5633, 1.4719E-11, 3.48E-11, 8908.6182, 1.3853E-11, 3.3288E-11,
    8676.6732, 1.2987E-11, 3.1776E-11, 8444.7281, 1.2121E-11, 3.0265E-11,
    8212.7831, 1.1255E-11, 2.8753E-11, 7980.838, 1.0389E-11, 2.7241E-11,
    7748.893, 9.5234E-12, 2.5729E-11, 7516.9479, 8.6574E-12, 2.4218E-11,
    7285.0029, 7.7913E-12, 2.2706E-11, 7053.0578, 6.9253E-12, 2.1194E-11,
    6821.1128, 6.0593E-12, 1.9682E-11, 6589.1677, 5.1933E-12, 1.8171E-11,
    6357.2227, 4.3273E-12, 1.6659E-11, 6125.2776, 3.4613E-12, 1.5147E-11,
    5893.3326, 2.5953E-12, 1.3635E-11, 25661.3875, 1.7293E-12, 1.2124E-11,
    4887.6275, 8.6324E-13, 7.5851E-11, -2.8669E-11, 9836.3984, 7.2831E-11,
    -2.6981E-11, 9604.4534, 6.9811E-11, -2.5292E-11, 9372.5083, 6.6791E-11,
    -2.3603E-11, 9140.5633, 6.3771E-11, -2.1914E-11, 8908.6182, 6.0751E-11,
    -2.0225E-11, 8676.6732, 5.7731E-11, -1.8537E-11, 8444.7281, 5.4711E-11,
    -1.6848E-11, 8212.7831, 5.1691E-11, -1.5159E-11, 7980.838, 4.8671E-11,
    -1.347E-11, 7748.893, 4.5651E-11, -1.1781E-11, 7516.9479, 4.2631E-11,
    -1.0092E-11, 7285.0029, 3.9611E-11, -8.4037E-12, 7053.0578, 3.6591E-11,
    -6.7148E-12, 6821.1128, 3.3571E-11, -5.026E-12, 6589.1677, 3.0551E-11,
    -3.3372E-12, 6357.2227, 2.7531E-11, -1.6484E-12, 6125.2776, 2.4511E-11,
    4.0431E-14, 5893.3326, 2.1491E-11, 1.7293E-12, 25661.3875, 1.8471E-11,
    3.4181E-12, 4887.6275, 8372.8305, 2.7232E-11, 6.7951E-11, 8179.2081,
    2.6393E-11, 6.5202E-11, 7985.5857, 2.5554E-11, 6.2453E-11, 7791.9633,
    2.4714E-11, 5.9704E-11, 7598.341, 2.3875E-11, 5.6956E-11, 7404.7186,
    2.3036E-11, 5.4207E-11, 7211.0962, 2.2196E-11, 5.1458E-11, 7017.4738,
    2.1357E-11, 4.8709E-11, 6823.8514, 2.0517E-11, 4.596E-11, 6630.229,
    1.9678E-11, 4.3211E-11, 6436.6066, 1.8839E-11, 4.0462E-11, 6242.9842,
    1.7999E-11, 3.7713E-11, 6049.3618, 1.716E-11, 3.4964E-11, 5855.7394,
    1.6321E-11, 3.2215E-11, 5662.1171, 1.5481E-11, 2.9466E-11, 5468.4947,
    1.4642E-11, 2.6717E-11, 5274.8723, 1.3802E-11, 2.3969E-11, 5081.2499,
    1.2963E-11, 2.122E-11, 4887.6275, 1.2124E-11, 1.8471E-11, 24694.0051,
    1.1284E-11, 1.5722E-11, 3.9266E-11, 8372.8305, 2.0075E-11, 3.7793E-11,
    8179.2081, 1.9149E-11, 3.632E-11, 7985.5857, 1.8224E-11, 3.4848E-11,
    7791.9633, 1.7298E-11, 3.3375E-11, 7598.341, 1.6373E-11, 3.1902E-11,
    7404.7186, 1.5448E-11, 3.0429E-11, 7211.0962, 1.4522E-11, 2.8957E-11,
    7017.4738, 1.3597E-11, 2.7484E-11, 6823.8514, 1.2672E-11, 2.6011E-11,
    6630.229, 1.1746E-11, 2.4539E-11, 6436.6066, 1.0821E-11, 2.3066E-11,
    6242.9842, 9.8956E-12, 2.1593E-11, 6049.3618, 8.9702E-12, 2.012E-11,
    5855.7394, 8.0449E-12, 1.8648E-11, 5662.1171, 7.1195E-12, 1.7175E-11,
    5468.4947, 6.1941E-12, 1.5702E-11, 5274.8723, 5.2688E-12, 1.423E-11,
    5081.2499, 4.3434E-12, 1.2757E-11, 4887.6275, 3.4181E-12, 1.1284E-11,
    24694.0051, 2.4927E-12, 6.3344E-11, -2.8467E-11, 8372.8305, 6.0838E-11,
    -2.6838E-11, 8179.2081, 5.8331E-11, -2.5208E-11, 7985.5857, 5.5825E-11,
    -2.3579E-11, 7791.9633, 5.3318E-11, -2.1949E-11, 7598.341, 5.0812E-11,
    -2.032E-11, 7404.7186, 4.8306E-11, -1.869E-11, 7211.0962, 4.5799E-11,
    -1.7061E-11, 7017.4738, 4.3293E-11, -1.5431E-11, 6823.8514, 4.0786E-11,
    -1.3802E-11, 6630.229, 3.828E-11, -1.2173E-11, 6436.6066, 3.5773E-11,
    -1.0543E-11, 6242.9842, 3.3267E-11, -8.9136E-12, 6049.3618, 3.076E-11,
    -7.2841E-12, 5855.7394, 2.8254E-11, -5.6547E-12, 5662.1171, 2.5748E-11,
    -4.0252E-12, 5468.4947, 2.3241E-11, -2.3957E-12, 5274.8723, 2.0735E-11,
    -7.6623E-13, 5081.2499, 1.8228E-11, 8.6324E-13, 4887.6275, 1.5722E-11,
    2.4927E-12, 24694.0051 };

  static const real_T b_a[3600] = { -59382.3876, -6.9307E-11, -3.0096E-10,
    -37659.217, -6.7726E-11, -2.8846E-10, -35937.5645, -6.6144E-11, -2.7595E-10,
    -34218.9483, -6.4563E-11, -2.6344E-10, -32504.8866, -6.2982E-11, -2.5094E-10,
    -30796.8974, -6.1401E-11, -2.3843E-10, -29096.4989, -5.982E-11, -2.2593E-10,
    -27405.2093, -5.8239E-11, -2.1342E-10, -25724.5467, -5.6658E-11, -2.0091E-10,
    -24056.0293, -5.5077E-11, -1.8841E-10, -22401.1752, -5.3496E-11, -1.759E-10,
    -20761.5025, -5.1915E-11, -1.6339E-10, -19138.5294, -5.0333E-11, -1.5089E-10,
    -17533.7741, -4.8752E-11, -1.3838E-10, -15948.7546, -4.7171E-11, -1.2588E-10,
    -14384.9892, -4.559E-11, -1.1337E-10, -12843.9959, -4.4009E-11, -1.0086E-10,
    -11327.2929, -4.2428E-11, -8.8357E-11, -9836.3984, -4.0847E-11, -7.5851E-11,
    -8372.8305, -3.9266E-11, -6.3344E-11, -6.9307E-11, -59382.3876, 3.231E-11,
    -6.7092E-11, -37659.217, 3.2107E-11, -6.4878E-11, -35937.5645, 3.1905E-11,
    -6.2663E-11, -34218.9483, 3.1703E-11, -6.0449E-11, -32504.8866, 3.1501E-11,
    -5.8235E-11, -30796.8974, 3.1298E-11, -5.602E-11, -29096.4989, 3.1096E-11,
    -5.3806E-11, -27405.2093, 3.0894E-11, -5.1591E-11, -25724.5467, 3.0692E-11,
    -4.9377E-11, -24056.0293, 3.0489E-11, -4.7162E-11, -22401.1752, 3.0287E-11,
    -4.4948E-11, -20761.5025, 3.0085E-11, -4.2734E-11, -19138.5294, 2.9883E-11,
    -4.0519E-11, -17533.7741, 2.9681E-11, -3.8305E-11, -15948.7546, 2.9478E-11,
    -3.609E-11, -14384.9892, 2.9276E-11, -3.3876E-11, -12843.9959, 2.9074E-11,
    -3.1661E-11, -11327.2929, 2.8872E-11, -2.9447E-11, -9836.3984, 2.8669E-11,
    -2.7232E-11, -8372.8305, 2.8467E-11, -3.0096E-10, 3.231E-11, -59382.3876,
    -2.887E-10, 2.9552E-11, -37659.217, -2.7644E-10, 2.6795E-11, -35937.5645,
    -2.6417E-10, 2.4038E-11, -34218.9483, -2.5191E-10, 2.1281E-11, -32504.8866,
    -2.3964E-10, 1.8524E-11, -30796.8974, -2.2738E-10, 1.5767E-11, -29096.4989,
    -2.1512E-10, 1.301E-11, -27405.2093, -2.0285E-10, 1.0253E-11, -25724.5467,
    -1.9059E-10, 7.496E-12, -24056.0293, -1.7833E-10, 4.739E-12, -22401.1752,
    -1.6606E-10, 1.9819E-12, -20761.5025, -1.538E-10, -7.7515E-13, -19138.5294,
    -1.4153E-10, -3.5322E-12, -17533.7741, -1.2927E-10, -6.2893E-12, -15948.7546,
    -1.1701E-10, -9.0463E-12, -14384.9892, -1.0474E-10, -1.1803E-11, -12843.9959,
    -9.2479E-11, -1.456E-11, -11327.2929, -8.0215E-11, -1.7317E-11, -9836.3984,
    -6.7951E-11, -2.0075E-11, -8372.8305, -37659.217, -6.7092E-11, -2.887E-10,
    -56543.5106, -6.555E-11, -2.7671E-10, -34885.9892, -6.4008E-11, -2.6471E-10,
    -33229.986, -6.2466E-11, -2.5272E-10, -31577.0191, -6.0924E-11, -2.4073E-10,
    -29928.6066, -5.9382E-11, -2.2874E-10, -28286.2667, -5.784E-11, -2.1674E-10,
    -26651.5175, -5.6298E-11, -2.0475E-10, -25025.8771, -5.4756E-11, -1.9276E-10,
    -23410.8638, -5.3214E-11, -1.8077E-10, -21807.9957, -5.1672E-11, -1.6877E-10,
    -20218.7908, -5.013E-11, -1.5678E-10, -18644.7674, -4.8587E-11, -1.4479E-10,
    -17087.4436, -4.7045E-11, -1.3279E-10, -15548.3375, -4.5503E-11, -1.208E-10,
    -14028.9673, -4.3961E-11, -1.0881E-10, -12530.8511, -4.2419E-11, -9.6816E-11,
    -11055.5071, -4.0877E-11, -8.4823E-11, -9604.4534, -3.9335E-11, -7.2831E-11,
    -8179.2081, -3.7793E-11, -6.0838E-11, -6.7726E-11, -37659.217, 2.9552E-11,
    -6.555E-11, -56543.5106, 2.941E-11, -6.3375E-11, -34885.9892, 2.9267E-11,
    -6.1199E-11, -33229.986, 2.9124E-11, -5.9024E-11, -31577.0191, 2.8981E-11,
    -5.6849E-11, -29928.6066, 2.8838E-11, -5.4673E-11, -28286.2667, 2.8695E-11,
    -5.2498E-11, -26651.5175, 2.8552E-11, -5.0322E-11, -25025.8771, 2.8409E-11,
    -4.8147E-11, -23410.8638, 2.8267E-11, -4.5972E-11, -21807.9957, 2.8124E-11,
    -4.3796E-11, -20218.7908, 2.7981E-11, -4.1621E-11, -18644.7674, 2.7838E-11,
    -3.9445E-11, -17087.4436, 2.7695E-11, -3.727E-11, -15548.3375, 2.7552E-11,
    -3.5095E-11, -14028.9673, 2.7409E-11, -3.2919E-11, -12530.8511, 2.7266E-11,
    -3.0744E-11, -11055.5071, 2.7124E-11, -2.8569E-11, -9604.4534, 2.6981E-11,
    -2.6393E-11, -8179.2081, 2.6838E-11, -2.8846E-10, 3.2107E-11, -37659.217,
    -2.7671E-10, 2.941E-11, -56543.5106, -2.6496E-10, 2.6712E-11, -34885.9892,
    -2.5321E-10, 2.4014E-11, -33229.986, -2.4146E-10, 2.1316E-11, -31577.0191,
    -2.2971E-10, 1.8619E-11, -29928.6066, -2.1796E-10, 1.5921E-11, -28286.2667,
    -2.0621E-10, 1.3223E-11, -26651.5175, -1.9446E-10, 1.0526E-11, -25025.8771,
    -1.827E-10, 7.8279E-12, -23410.8638, -1.7095E-10, 5.1302E-12, -21807.9957,
    -1.592E-10, 2.4325E-12, -20218.7908, -1.4745E-10, -2.6521E-13, -18644.7674,
    -1.357E-10, -2.9629E-12, -17087.4436, -1.2395E-10, -5.6606E-12, -15548.3375,
    -1.122E-10, -8.3583E-12, -14028.9673, -1.0045E-10, -1.1056E-11, -12530.8511,
    -8.8703E-11, -1.3754E-11, -11055.5071, -7.6953E-11, -1.6451E-11, -9604.4534,
    -6.5202E-11, -1.9149E-11, -8179.2081, -35937.5645, -6.4878E-11, -2.7644E-10,
    -34885.9892, -6.3375E-11, -2.6496E-10, -53834.4139, -6.1872E-11, -2.5348E-10,
    -32241.0237, -6.0369E-11, -2.42E-10, -30649.1516, -5.8866E-11, -2.3052E-10,
    -29060.3158, -5.7363E-11, -2.1904E-10, -27476.0344, -5.586E-11, -2.0756E-10,
    -25897.8256, -5.4357E-11, -1.9608E-10, -24327.2075, -5.2854E-11, -1.846E-10,
    -22765.6983, -5.1351E-11, -1.7312E-10, -21214.8161, -4.9848E-11, -1.6164E-10,
    -19676.0791, -4.8345E-11, -1.5016E-10, -18151.0053, -4.6841E-11, -1.3869E-10,
    -16641.113, -4.5338E-11, -1.2721E-10, -15147.9203, -4.3835E-11, -1.1573E-10,
    -13672.9454, -4.2332E-11, -1.0425E-10, -12217.7063, -4.0829E-11, -9.2769E-11,
    -10783.7212, -3.9326E-11, -8.129E-11, -9372.5083, -3.7823E-11, -6.9811E-11,
    -7985.5857, -3.632E-11, -5.8331E-11, -6.6144E-11, -35937.5645, 2.6795E-11,
    -6.4008E-11, -34885.9892, 2.6712E-11, -6.1872E-11, -53834.4139, 2.6628E-11,
    -5.9735E-11, -32241.0237, 2.6545E-11, -5.7599E-11, -30649.1516, 2.6461E-11,
    -5.5463E-11, -29060.3158, 2.6378E-11, -5.3326E-11, -27476.0344, 2.6294E-11,
    -5.119E-11, -25897.8256, 2.6211E-11, -4.9054E-11, -24327.2075, 2.6127E-11,
    -4.6917E-11, -22765.6983, 2.6044E-11, -4.4781E-11, -21214.8161, 2.596E-11,
    -4.2645E-11, -19676.0791, 2.5877E-11, -4.0508E-11, -18151.0053, 2.5793E-11,
    -3.8372E-11, -16641.113, 2.571E-11, -3.6236E-11, -15147.9203, 2.5626E-11,
    -3.4099E-11, -13672.9454, 2.5542E-11, -3.1963E-11, -12217.7063, 2.5459E-11,
    -2.9826E-11, -10783.7212, 2.5375E-11, -2.769E-11, -9372.5083, 2.5292E-11,
    -2.5554E-11, -7985.5857, 2.5208E-11, -2.7595E-10, 3.1905E-11, -35937.5645,
    -2.6471E-10, 2.9267E-11, -34885.9892, -2.5348E-10, 2.6628E-11, -53834.4139,
    -2.4224E-10, 2.399E-11, -32241.0237, -2.31E-10, 2.1352E-11, -30649.1516,
    -2.1977E-10, 1.8713E-11, -29060.3158, -2.0853E-10, 1.6075E-11, -27476.0344,
    -1.9729E-10, 1.3437E-11, -25897.8256, -1.8606E-10, 1.0798E-11, -24327.2075,
    -1.7482E-10, 8.1598E-12, -22765.6983, -1.6358E-10, 5.5214E-12, -21214.8161,
    -1.5235E-10, 2.8831E-12, -19676.0791, -1.4111E-10, 2.4472E-13, -18151.0053,
    -1.2987E-10, -2.3936E-12, -16641.113, -1.1864E-10, -5.032E-12, -15147.9203,
    -1.074E-10, -7.6704E-12, -13672.9454, -9.6163E-11, -1.0309E-11, -12217.7063,
    -8.4927E-11, -1.2947E-11, -10783.7212, -7.369E-11, -1.5585E-11, -9372.5083,
    -6.2453E-11, -1.8224E-11, -7985.5857, -34218.9483, -6.2663E-11, -2.6417E-10,
    -33229.986, -6.1199E-11, -2.5321E-10, -32241.0237, -5.9735E-11, -2.4224E-10,
    -51252.0614, -5.8271E-11, -2.3128E-10, -29721.2841, -5.6807E-11, -2.2031E-10,
    -28192.025, -5.5343E-11, -2.0934E-10, -26665.8021, -5.3879E-11, -1.9838E-10,
    -25144.1337, -5.2415E-11, -1.8741E-10, -23628.5379, -5.0951E-11, -1.7645E-10,
    -22120.5328, -4.9487E-11, -1.6548E-10, -20621.6366, -4.8023E-11, -1.5452E-10,
    -19133.3673, -4.6559E-11, -1.4355E-10, -17657.2433, -4.5096E-11, -1.3258E-10,
    -16194.7825, -4.3632E-11, -1.2162E-10, -14747.5032, -4.2168E-11, -1.1065E-10,
    -13316.9235, -4.0704E-11, -9.9687E-11, -11904.5615, -3.924E-11, -8.8722E-11,
    -10511.9354, -3.7776E-11, -7.7756E-11, -9140.5633, -3.6312E-11, -6.6791E-11,
    -7791.9633, -3.4848E-11, -5.5825E-11, -6.4563E-11, -34218.9483, 2.4038E-11,
    -6.2466E-11, -33229.986, 2.4014E-11, -6.0369E-11, -32241.0237, 2.399E-11,
    -5.8271E-11, -51252.0614, 2.3966E-11, -5.6174E-11, -29721.2841, 2.3942E-11,
    -5.4077E-11, -28192.025, 2.3917E-11, -5.1979E-11, -26665.8021, 2.3893E-11,
    -4.9882E-11, -25144.1337, 2.3869E-11, -4.7785E-11, -23628.5379, 2.3845E-11,
    -4.5688E-11, -22120.5328, 2.3821E-11, -4.359E-11, -20621.6366, 2.3797E-11,
    -4.1493E-11, -19133.3673, 2.3772E-11, -3.9396E-11, -17657.2433, 2.3748E-11,
    -3.7298E-11, -16194.7825, 2.3724E-11, -3.5201E-11, -14747.5032, 2.37E-11,
    -3.3104E-11, -13316.9235, 2.3676E-11, -3.1006E-11, -11904.5615, 2.3651E-11,
    -2.8909E-11, -10511.9354, 2.3627E-11, -2.6812E-11, -9140.5633, 2.3603E-11,
    -2.4714E-11, -7791.9633, 2.3579E-11, -2.6344E-10, 3.1703E-11, -34218.9483,
    -2.5272E-10, 2.9124E-11, -33229.986, -2.42E-10, 2.6545E-11, -32241.0237,
    -2.3128E-10, 2.3966E-11, -51252.0614, -2.2055E-10, 2.1387E-11, -29721.2841,
    -2.0983E-10, 1.8808E-11, -28192.025, -1.9911E-10, 1.6229E-11, -26665.8021,
    -1.8838E-10, 1.365E-11, -25144.1337, -1.7766E-10, 1.1071E-11, -23628.5379,
    -1.6694E-10, 8.4917E-12, -22120.5328, -1.5621E-10, 5.9127E-12, -20621.6366,
    -1.4549E-10, 3.3337E-12, -19133.3673, -1.3477E-10, 7.5465E-13, -17657.2433,
    -1.2404E-10, -1.8244E-12, -16194.7825, -1.1332E-10, -4.4034E-12, -14747.5032,
    -1.026E-10, -6.9824E-12, -13316.9235, -9.1874E-11, -9.5614E-12, -11904.5615,
    -8.1151E-11, -1.214E-11, -10511.9354, -7.0428E-11, -1.4719E-11, -9140.5633,
    -5.9704E-11, -1.7298E-11, -7791.9633, -32504.8866, -6.0449E-11, -2.5191E-10,
    -31577.0191, -5.9024E-11, -2.4146E-10, -30649.1516, -5.7599E-11, -2.31E-10,
    -29721.2841, -5.6174E-11, -2.2055E-10, -48793.4166, -5.4749E-11, -2.101E-10,
    -27323.7342, -5.3324E-11, -1.9965E-10, -25855.5699, -5.1899E-11, -1.892E-10,
    -24390.4419, -5.0474E-11, -1.7874E-10, -22929.8683, -4.9049E-11, -1.6829E-10,
    -21475.3673, -4.7624E-11, -1.5784E-10, -20028.457, -4.6199E-11, -1.4739E-10,
    -18590.6556, -4.4774E-11, -1.3694E-10, -17163.4812, -4.335E-11, -1.2648E-10,
    -15748.452, -4.1925E-11, -1.1603E-10, -14347.0861, -4.05E-11, -1.0558E-10,
    -12960.9016, -3.9075E-11, -9.5127E-11, -11591.4167, -3.765E-11, -8.4675E-11,
    -10240.1495, -3.6225E-11, -7.4223E-11, -8908.6182, -3.48E-11, -6.3771E-11,
    -7598.341, -3.3375E-11, -5.3318E-11, -6.2982E-11, -32504.8866, 2.1281E-11,
    -6.0924E-11, -31577.0191, 2.1316E-11, -5.8866E-11, -30649.1516, 2.1352E-11,
    -5.6807E-11, -29721.2841, 2.1387E-11, -5.4749E-11, -48793.4166, 2.1422E-11,
    -5.2691E-11, -27323.7342, 2.1457E-11, -5.0633E-11, -25855.5699, 2.1492E-11,
    -4.8574E-11, -24390.4419, 2.1527E-11, -4.6516E-11, -22929.8683, 2.1563E-11,
    -4.4458E-11, -21475.3673, 2.1598E-11, -4.2399E-11, -20028.457, 2.1633E-11,
    -4.0341E-11, -18590.6556, 2.1668E-11, -3.8283E-11, -17163.4812, 2.1703E-11,
    -3.6225E-11, -15748.452, 2.1738E-11, -3.4166E-11, -14347.0861, 2.1774E-11,
    -3.2108E-11, -12960.9016, 2.1809E-11, -3.005E-11, -11591.4167, 2.1844E-11,
    -2.7992E-11, -10240.1495, 2.1879E-11, -2.5933E-11, -8908.6182, 2.1914E-11,
    -2.3875E-11, -7598.341, 2.1949E-11, -2.5094E-10, 3.1501E-11, -32504.8866,
    -2.4073E-10, 2.8981E-11, -31577.0191, -2.3052E-10, 2.6461E-11, -30649.1516,
    -2.2031E-10, 2.3942E-11, -29721.2841, -2.101E-10, 2.1422E-11, -48793.4166,
    -1.9989E-10, 1.8902E-11, -27323.7342, -1.8968E-10, 1.6383E-11, -25855.5699,
    -1.7947E-10, 1.3863E-11, -24390.4419, -1.6926E-10, 1.1343E-11, -22929.8683,
    -1.5905E-10, 8.8236E-12, -21475.3673, -1.4884E-10, 6.3039E-12, -20028.457,
    -1.3863E-10, 3.7843E-12, -18590.6556, -1.2842E-10, 1.2646E-12, -17163.4812,
    -1.1821E-10, -1.2551E-12, -15748.452, -1.08E-10, -3.7748E-12, -14347.0861,
    -9.7794E-11, -6.2944E-12, -12960.9016, -8.7584E-11, -8.8141E-12, -11591.4167,
    -7.7375E-11, -1.1334E-11, -10240.1495, -6.7165E-11, -1.3853E-11, -8908.6182,
    -5.6956E-11, -1.6373E-11, -7598.341, -30796.8974, -5.8235E-11, -2.3964E-10,
    -29928.6066, -5.6849E-11, -2.2971E-10, -29060.3158, -5.5463E-11, -2.1977E-10,
    -28192.025, -5.4077E-11, -2.0983E-10, -27323.7342, -5.2691E-11, -1.9989E-10,
    -46455.4434, -5.1305E-11, -1.8995E-10, -25045.3376, -4.9919E-11, -1.8001E-10,
    -23636.75, -4.8533E-11, -1.7007E-10, -22231.1987, -4.7147E-11, -1.6014E-10,
    -20830.2018, -4.5761E-11, -1.502E-10, -19435.2775, -4.4375E-11, -1.4026E-10,
    -18047.9439, -4.2989E-11, -1.3032E-10, -16669.7192, -4.1604E-11, -1.2038E-10,
    -15302.1215, -4.0218E-11, -1.1044E-10, -13946.6689, -3.8832E-11, -1.005E-10,
    -12604.8797, -3.7446E-11, -9.0566E-11, -11278.2719, -3.606E-11, -8.0628E-11,
    -9968.3636, -3.4674E-11, -7.0689E-11, -8676.6732, -3.3288E-11, -6.0751E-11,
    -7404.7186, -3.1902E-11, -5.0812E-11, -6.1401E-11, -30796.8974, 1.8524E-11,
    -5.9382E-11, -29928.6066, 1.8619E-11, -5.7363E-11, -29060.3158, 1.8713E-11,
    -5.5343E-11, -28192.025, 1.8808E-11, -5.3324E-11, -27323.7342, 1.8902E-11,
    -5.1305E-11, -46455.4434, 1.8997E-11, -4.9286E-11, -25045.3376, 1.9091E-11,
    -4.7266E-11, -23636.75, 1.9186E-11, -4.5247E-11, -22231.1987, 1.928E-11,
    -4.3228E-11, -20830.2018, 1.9375E-11, -4.1209E-11, -19435.2775, 1.9469E-11,
    -3.919E-11, -18047.9439, 1.9564E-11, -3.717E-11, -16669.7192, 1.9658E-11,
    -3.5151E-11, -15302.1215, 1.9753E-11, -3.3132E-11, -13946.6689, 1.9847E-11,
    -3.1113E-11, -12604.8797, 1.9942E-11, -2.9093E-11, -11278.2719, 2.0036E-11,
    -2.7074E-11, -9968.3636, 2.0131E-11, -2.5055E-11, -8676.6732, 2.0225E-11,
    -2.3036E-11, -7404.7186, 2.032E-11, -2.3843E-10, 3.1298E-11, -30796.8974,
    -2.2874E-10, 2.8838E-11, -29928.6066, -2.1904E-10, 2.6378E-11, -29060.3158,
    -2.0934E-10, 2.3917E-11, -28192.025, -1.9965E-10, 2.1457E-11, -27323.7342,
    -1.8995E-10, 1.8997E-11, -46455.4434, -1.8026E-10, 1.6536E-11, -25045.3376,
    -1.7056E-10, 1.4076E-11, -23636.75, -1.6086E-10, 1.1616E-11, -22231.1987,
    -1.5117E-10, 9.1555E-12, -20830.2018, -1.4147E-10, 6.6952E-12, -19435.2775,
    -1.3178E-10, 4.2348E-12, -18047.9439, -1.2208E-10, 1.7745E-12, -16669.7192,
    -1.1238E-10, -6.858E-13, -15302.1215, -1.0269E-10, -3.1461E-12, -13946.6689,
    -9.2991E-11, -5.6064E-12, -12604.8797, -8.3295E-11, -8.0668E-12, -11278.2719,
    -7.3599E-11, -1.0527E-11, -9968.3636, -6.3903E-11, -1.2987E-11, -8676.6732,
    -5.4207E-11, -1.5448E-11, -7404.7186, -29096.4989, -5.602E-11, -2.2738E-10,
    -28286.2667, -5.4673E-11, -2.1796E-10, -27476.0344, -5.3326E-11, -2.0853E-10,
    -26665.8021, -5.1979E-11, -1.9911E-10, -25855.5699, -5.0633E-11, -1.8968E-10,
    -25045.3376, -4.9286E-11, -1.8026E-10, -44235.1054, -4.7939E-11, -1.7083E-10,
    -22883.0582, -4.6592E-11, -1.6141E-10, -21532.5291, -4.5245E-11, -1.5198E-10,
    -20185.0363, -4.3898E-11, -1.4256E-10, -18842.0979, -4.2551E-11, -1.3313E-10,
    -17505.2322, -4.1204E-11, -1.2371E-10, -16175.9571, -3.9858E-11, -1.1428E-10,
    -14855.7909, -3.8511E-11, -1.0486E-10, -13546.2518, -3.7164E-11, -9.5431E-11,
    -12248.8578, -3.5817E-11, -8.6006E-11, -10965.1271, -3.447E-11, -7.6581E-11,
    -9696.5778, -3.3123E-11, -6.7156E-11, -8444.7281, -3.1776E-11, -5.7731E-11,
    -7211.0962, -3.0429E-11, -4.8306E-11, -5.982E-11, -29096.4989, 1.5767E-11,
    -5.784E-11, -28286.2667, 1.5921E-11, -5.586E-11, -27476.0344, 1.6075E-11,
    -5.3879E-11, -26665.8021, 1.6229E-11, -5.1899E-11, -25855.5699, 1.6383E-11,
    -4.9919E-11, -25045.3376, 1.6536E-11, -4.7939E-11, -44235.1054, 1.669E-11,
    -4.5959E-11, -22883.0582, 1.6844E-11, -4.3978E-11, -21532.5291, 1.6998E-11,
    -4.1998E-11, -20185.0363, 1.7152E-11, -4.0018E-11, -18842.0979, 1.7306E-11,
    -3.8038E-11, -17505.2322, 1.746E-11, -3.6058E-11, -16175.9571, 1.7613E-11,
    -3.4077E-11, -14855.7909, 1.7767E-11, -3.2097E-11, -13546.2518, 1.7921E-11,
    -3.0117E-11, -12248.8578, 1.8075E-11, -2.8137E-11, -10965.1271, 1.8229E-11,
    -2.6157E-11, -9696.5778, 1.8383E-11, -2.4176E-11, -8444.7281, 1.8537E-11,
    -2.2196E-11, -7211.0962, 1.869E-11, -2.2593E-10, 3.1096E-11, -29096.4989,
    -2.1674E-10, 2.8695E-11, -28286.2667, -2.0756E-10, 2.6294E-11, -27476.0344,
    -1.9838E-10, 2.3893E-11, -26665.8021, -1.892E-10, 2.1492E-11, -25855.5699,
    -1.8001E-10, 1.9091E-11, -25045.3376, -1.7083E-10, 1.669E-11, -44235.1054,
    -1.6165E-10, 1.4289E-11, -22883.0582, -1.5247E-10, 1.1888E-11, -21532.5291,
    -1.4328E-10, 9.4874E-12, -20185.0363, -1.341E-10, 7.0864E-12, -18842.0979,
    -1.2492E-10, 4.6854E-12, -17505.2322, -1.1574E-10, 2.2845E-12, -16175.9571,
    -1.0655E-10, -1.1652E-13, -14855.7909, -9.737E-11, -2.5175E-12, -13546.2518,
    -8.8188E-11, -4.9185E-12, -12248.8578, -7.9005E-11, -7.3194E-12, -10965.1271,
    -6.9823E-11, -9.7204E-12, -9696.5778, -6.064E-11, -1.2121E-11, -8444.7281,
    -5.1458E-11, -1.4522E-11, -7211.0962, -27405.2093, -5.3806E-11, -2.1512E-10,
    -26651.5175, -5.2498E-11, -2.0621E-10, -25897.8256, -5.119E-11, -1.9729E-10,
    -25144.1337, -4.9882E-11, -1.8838E-10, -24390.4419, -4.8574E-11, -1.7947E-10,
    -23636.75, -4.7266E-11, -1.7056E-10, -22883.0582, -4.5959E-11, -1.6165E-10,
    -42129.3663, -4.4651E-11, -1.5274E-10, -20833.8595, -4.3343E-11, -1.4383E-10,
    -19539.8708, -4.2035E-11, -1.3491E-10, -18248.9184, -4.0727E-11, -1.26E-10,
    -16962.5205, -3.9419E-11, -1.1709E-10, -15682.1951, -3.8112E-11, -1.0818E-10,
    -14409.4604, -3.6804E-11, -9.9268E-11, -13145.8346, -3.5496E-11, -9.0356E-11,
    -11892.8359, -3.4188E-11, -8.1445E-11, -10651.9823, -3.288E-11, -7.2534E-11,
    -9424.7919, -3.1572E-11, -6.3622E-11, -8212.7831, -3.0265E-11, -5.4711E-11,
    -7017.4738, -2.8957E-11, -4.5799E-11, -5.8239E-11, -27405.2093, 1.301E-11,
    -5.6298E-11, -26651.5175, 1.3223E-11, -5.4357E-11, -25897.8256, 1.3437E-11,
    -5.2415E-11, -25144.1337, 1.365E-11, -5.0474E-11, -24390.4419, 1.3863E-11,
    -4.8533E-11, -23636.75, 1.4076E-11, -4.6592E-11, -22883.0582, 1.4289E-11,
    -4.4651E-11, -42129.3663, 1.4503E-11, -4.271E-11, -20833.8595, 1.4716E-11,
    -4.0768E-11, -19539.8708, 1.4929E-11, -3.8827E-11, -18248.9184, 1.5142E-11,
    -3.6886E-11, -16962.5205, 1.5355E-11, -3.4945E-11, -15682.1951, 1.5569E-11,
    -3.3004E-11, -14409.4604, 1.5782E-11, -3.1063E-11, -13145.8346, 1.5995E-11,
    -2.9121E-11, -11892.8359, 1.6208E-11, -2.718E-11, -10651.9823, 1.6421E-11,
    -2.5239E-11, -9424.7919, 1.6635E-11, -2.3298E-11, -8212.7831, 1.6848E-11,
    -2.1357E-11, -7017.4738, 1.7061E-11, -2.1342E-10, 3.0894E-11, -27405.2093,
    -2.0475E-10, 2.8552E-11, -26651.5175, -1.9608E-10, 2.6211E-11, -25897.8256,
    -1.8741E-10, 2.3869E-11, -25144.1337, -1.7874E-10, 2.1527E-11, -24390.4419,
    -1.7007E-10, 1.9186E-11, -23636.75, -1.6141E-10, 1.6844E-11, -22883.0582,
    -1.5274E-10, 1.4503E-11, -42129.3663, -1.4407E-10, 1.2161E-11, -20833.8595,
    -1.354E-10, 9.8193E-12, -19539.8708, -1.2673E-10, 7.4776E-12, -18248.9184,
    -1.1806E-10, 5.136E-12, -16962.5205, -1.0939E-10, 2.7944E-12, -15682.1951,
    -1.0072E-10, 4.5276E-13, -14409.4604, -9.2054E-11, -1.8889E-12, -13145.8346,
    -8.3385E-11, -4.2305E-12, -11892.8359, -7.4716E-11, -6.5721E-12, -10651.9823,
    -6.6047E-11, -8.9138E-12, -9424.7919, -5.7378E-11, -1.1255E-11, -8212.7831,
    -4.8709E-11, -1.3597E-11, -7017.4738, -25724.5467, -5.1591E-11, -2.0285E-10,
    -25025.8771, -5.0322E-11, -1.9446E-10, -24327.2075, -4.9054E-11, -1.8606E-10,
    -23628.5379, -4.7785E-11, -1.7766E-10, -22929.8683, -4.6516E-11, -1.6926E-10,
    -22231.1987, -4.5247E-11, -1.6086E-10, -21532.5291, -4.3978E-11, -1.5247E-10,
    -20833.8595, -4.271E-11, -1.4407E-10, -40135.1899, -4.1441E-11, -1.3567E-10,
    -18894.7053, -4.0172E-11, -1.2727E-10, -17655.7389, -3.8903E-11, -1.1887E-10,
    -16419.8087, -3.7634E-11, -1.1048E-10, -15188.433, -3.6366E-11, -1.0208E-10,
    -13963.1299, -3.5097E-11, -9.368E-11, -12745.4175, -3.3828E-11, -8.5282E-11,
    -11536.814, -3.2559E-11, -7.6884E-11, -10338.8374, -3.129E-11, -6.8486E-11,
    -9153.0061, -3.0022E-11, -6.0089E-11, -7980.838, -2.8753E-11, -5.1691E-11,
    -6823.8514, -2.7484E-11, -4.3293E-11, -5.6658E-11, -25724.5467, 1.0253E-11,
    -5.4756E-11, -25025.8771, 1.0526E-11, -5.2854E-11, -24327.2075, 1.0798E-11,
    -5.0951E-11, -23628.5379, 1.1071E-11, -4.9049E-11, -22929.8683, 1.1343E-11,
    -4.7147E-11, -22231.1987, 1.1616E-11, -4.5245E-11, -21532.5291, 1.1888E-11,
    -4.3343E-11, -20833.8595, 1.2161E-11, -4.1441E-11, -40135.1899, 1.2433E-11,
    -3.9539E-11, -18894.7053, 1.2706E-11, -3.7637E-11, -17655.7389, 1.2979E-11,
    -3.5734E-11, -16419.8087, 1.3251E-11, -3.3832E-11, -15188.433, 1.3524E-11,
    -3.193E-11, -13963.1299, 1.3796E-11, -3.0028E-11, -12745.4175, 1.4069E-11,
    -2.8126E-11, -11536.814, 1.4341E-11, -2.6224E-11, -10338.8374, 1.4614E-11,
    -2.4322E-11, -9153.0061, 1.4886E-11, -2.242E-11, -7980.838, 1.5159E-11,
    -2.0517E-11, -6823.8514, 1.5431E-11, -2.0091E-10, 3.0692E-11, -25724.5467,
    -1.9276E-10, 2.8409E-11, -25025.8771, -1.846E-10, 2.6127E-11, -24327.2075,
    -1.7645E-10, 2.3845E-11, -23628.5379, -1.6829E-10, 2.1563E-11, -22929.8683,
    -1.6014E-10, 1.928E-11, -22231.1987, -1.5198E-10, 1.6998E-11, -21532.5291,
    -1.4383E-10, 1.4716E-11, -20833.8595, -1.3567E-10, 1.2433E-11, -40135.1899,
    -1.2751E-10, 1.0151E-11, -18894.7053, -1.1936E-10, 7.8689E-12, -17655.7389,
    -1.112E-10, 5.5866E-12, -16419.8087, -1.0305E-10, 3.3043E-12, -15188.433,
    -9.4893E-11, 1.022E-12, -13963.1299, -8.6737E-11, -1.2602E-12, -12745.4175,
    -7.8582E-11, -3.5425E-12, -11536.814, -7.0426E-11, -5.8248E-12, -10338.8374,
    -6.2271E-11, -8.1071E-12, -9153.0061, -5.4115E-11, -1.0389E-11, -7980.838,
    -4.596E-11, -1.2672E-11, -6823.8514, -24056.0293, -4.9377E-11, -1.9059E-10,
    -23410.8638, -4.8147E-11, -1.827E-10, -22765.6983, -4.6917E-11, -1.7482E-10,
    -22120.5328, -4.5688E-11, -1.6694E-10, -21475.3673, -4.4458E-11, -1.5905E-10,
    -20830.2018, -4.3228E-11, -1.5117E-10, -20185.0363, -4.1998E-11, -1.4328E-10,
    -19539.8708, -4.0768E-11, -1.354E-10, -18894.7053, -3.9539E-11, -1.2751E-10,
    -38249.5398, -3.8309E-11, -1.1963E-10, -17062.5593, -3.7079E-11, -1.1175E-10,
    -15877.097, -3.5849E-11, -1.0386E-10, -14694.671, -3.462E-11, -9.5977E-11,
    -13516.7994, -3.339E-11, -8.8092E-11, -12345.0004, -3.216E-11, -8.0208E-11,
    -11180.7921, -3.093E-11, -7.2324E-11, -10025.6926, -2.9701E-11, -6.4439E-11,
    -8881.2202, -2.8471E-11, -5.6555E-11, -7748.893, -2.7241E-11, -4.8671E-11,
    -6630.229, -2.6011E-11, -4.0786E-11, -5.5077E-11, -24056.0293, 7.496E-12,
    -5.3214E-11, -23410.8638, 7.8279E-12, -5.1351E-11, -22765.6983, 8.1598E-12,
    -4.9487E-11, -22120.5328, 8.4917E-12, -4.7624E-11, -21475.3673, 8.8236E-12,
    -4.5761E-11, -20830.2018, 9.1555E-12, -4.3898E-11, -20185.0363, 9.4874E-12,
    -4.2035E-11, -19539.8708, 9.8193E-12, -4.0172E-11, -18894.7053, 1.0151E-11,
    -3.8309E-11, -38249.5398, 1.0483E-11, -3.6446E-11, -17062.5593, 1.0815E-11,
    -3.4583E-11, -15877.097, 1.1147E-11, -3.272E-11, -14694.671, 1.1479E-11,
    -3.0857E-11, -13516.7994, 1.1811E-11, -2.8994E-11, -12345.0004, 1.2143E-11,
    -2.713E-11, -11180.7921, 1.2474E-11, -2.5267E-11, -10025.6926, 1.2806E-11,
    -2.3404E-11, -8881.2202, 1.3138E-11, -2.1541E-11, -7748.893, 1.347E-11,
    -1.9678E-11, -6630.229, 1.3802E-11, -1.8841E-10, 3.0489E-11, -24056.0293,
    -1.8077E-10, 2.8267E-11, -23410.8638, -1.7312E-10, 2.6044E-11, -22765.6983,
    -1.6548E-10, 2.3821E-11, -22120.5328, -1.5784E-10, 2.1598E-11, -21475.3673,
    -1.502E-10, 1.9375E-11, -20830.2018, -1.4256E-10, 1.7152E-11, -20185.0363,
    -1.3491E-10, 1.4929E-11, -19539.8708, -1.2727E-10, 1.2706E-11, -18894.7053,
    -1.1963E-10, 1.0483E-11, -38249.5398, -1.1199E-10, 8.2601E-12, -17062.5593,
    -1.0435E-10, 6.0372E-12, -15877.097, -9.6704E-11, 3.8143E-12, -14694.671,
    -8.9062E-11, 1.5913E-12, -13516.7994, -8.142E-11, -6.3162E-13, -12345.0004,
    -7.3779E-11, -2.8546E-12, -11180.7921, -6.6137E-11, -5.0775E-12, -10025.6926,
    -5.8495E-11, -7.3004E-12, -8881.2202, -5.0853E-11, -9.5234E-12, -7748.893,
    -4.3211E-11, -1.1746E-11, -6630.229, -22401.1752, -4.7162E-11, -1.7833E-10,
    -21807.9957, -4.5972E-11, -1.7095E-10, -21214.8161, -4.4781E-11, -1.6358E-10,
    -20621.6366, -4.359E-11, -1.5621E-10, -20028.457, -4.2399E-11, -1.4884E-10,
    -19435.2775, -4.1209E-11, -1.4147E-10, -18842.0979, -4.0018E-11, -1.341E-10,
    -18248.9184, -3.8827E-11, -1.2673E-10, -17655.7389, -3.7637E-11, -1.1936E-10,
    -17062.5593, -3.6446E-11, -1.1199E-10, -36469.3798, -3.5255E-11, -1.0462E-10,
    -15334.3853, -3.4064E-11, -9.7246E-11, -14200.9089, -3.2874E-11, -8.9876E-11,
    -13070.4689, -3.1683E-11, -8.2505E-11, -11944.5832, -3.0492E-11, -7.5134E-11,
    -10824.7702, -2.9302E-11, -6.7763E-11, -9712.5478, -2.8111E-11, -6.0392E-11,
    -8609.4344, -2.692E-11, -5.3021E-11, -7516.9479, -2.5729E-11, -4.5651E-11,
    -6436.6066, -2.4539E-11, -3.828E-11, -5.3496E-11, -22401.1752, 4.739E-12,
    -5.1672E-11, -21807.9957, 5.1302E-12, -4.9848E-11, -21214.8161, 5.5214E-12,
    -4.8023E-11, -20621.6366, 5.9127E-12, -4.6199E-11, -20028.457, 6.3039E-12,
    -4.4375E-11, -19435.2775, 6.6952E-12, -4.2551E-11, -18842.0979, 7.0864E-12,
    -4.0727E-11, -18248.9184, 7.4776E-12, -3.8903E-11, -17655.7389, 7.8689E-12,
    -3.7079E-11, -17062.5593, 8.2601E-12, -3.5255E-11, -36469.3798, 8.6514E-12,
    -3.3431E-11, -15334.3853, 9.0426E-12, -3.1607E-11, -14200.9089, 9.4339E-12,
    -2.9783E-11, -13070.4689, 9.8251E-12, -2.7959E-11, -11944.5832, 1.0216E-11,
    -2.6135E-11, -10824.7702, 1.0608E-11, -2.4311E-11, -9712.5478, 1.0999E-11,
    -2.2487E-11, -8609.4344, 1.139E-11, -2.0663E-11, -7516.9479, 1.1781E-11,
    -1.8839E-11, -6436.6066, 1.2173E-11, -1.759E-10, 3.0287E-11, -22401.1752,
    -1.6877E-10, 2.8124E-11, -21807.9957, -1.6164E-10, 2.596E-11, -21214.8161,
    -1.5452E-10, 2.3797E-11, -20621.6366, -1.4739E-10, 2.1633E-11, -20028.457,
    -1.4026E-10, 1.9469E-11, -19435.2775, -1.3313E-10, 1.7306E-11, -18842.0979,
    -1.26E-10, 1.5142E-11, -18248.9184, -1.1887E-10, 1.2979E-11, -17655.7389,
    -1.1175E-10, 1.0815E-11, -17062.5593, -1.0462E-10, 8.6514E-12, -36469.3798,
    -9.7489E-11, 6.4878E-12, -15334.3853, -9.036E-11, 4.3242E-12, -14200.9089,
    -8.3232E-11, 2.1606E-12, -13070.4689, -7.6104E-11, -2.9897E-15, -11944.5832,
    -6.8975E-11, -2.1666E-12, -10824.7702, -6.1847E-11, -4.3302E-12, -9712.5478,
    -5.4719E-11, -6.4938E-12, -8609.4344, -4.759E-11, -8.6574E-12, -7516.9479,
    -4.0462E-11, -1.0821E-11, -6436.6066, -20761.5025, -4.4948E-11, -1.6606E-10,
    -20218.7908, -4.3796E-11, -1.592E-10, -19676.0791, -4.2645E-11, -1.5235E-10,
    -19133.3673, -4.1493E-11, -1.4549E-10, -18590.6556, -4.0341E-11, -1.3863E-10,
    -18047.9439, -3.919E-11, -1.3178E-10, -17505.2322, -3.8038E-11, -1.2492E-10,
    -16962.5205, -3.6886E-11, -1.1806E-10, -16419.8087, -3.5734E-11, -1.112E-10,
    -15877.097, -3.4583E-11, -1.0435E-10, -15334.3853, -3.3431E-11, -9.7489E-11,
    -34791.6736, -3.2279E-11, -9.0632E-11, -13707.1469, -3.1128E-11, -8.3774E-11,
    -12624.1383, -2.9976E-11, -7.6917E-11, -11544.1661, -2.8824E-11, -7.006E-11,
    -10468.7483, -2.7673E-11, -6.3202E-11, -9399.403, -2.6521E-11, -5.6345E-11,
    -8337.6485, -2.5369E-11, -4.9488E-11, -7285.0029, -2.4218E-11, -4.2631E-11,
    -6242.9842, -2.3066E-11, -3.5773E-11, -5.1915E-11, -20761.5025, 1.9819E-12,
    -5.013E-11, -20218.7908, 2.4325E-12, -4.8345E-11, -19676.0791, 2.8831E-12,
    -4.6559E-11, -19133.3673, 3.3337E-12, -4.4774E-11, -18590.6556, 3.7843E-12,
    -4.2989E-11, -18047.9439, 4.2348E-12, -4.1204E-11, -17505.2322, 4.6854E-12,
    -3.9419E-11, -16962.5205, 5.136E-12, -3.7634E-11, -16419.8087, 5.5866E-12,
    -3.5849E-11, -15877.097, 6.0372E-12, -3.4064E-11, -15334.3853, 6.4878E-12,
    -3.2279E-11, -34791.6736, 6.9384E-12, -3.0494E-11, -13707.1469, 7.389E-12,
    -2.8709E-11, -12624.1383, 7.8395E-12, -2.6924E-11, -11544.1661, 8.2901E-12,
    -2.5139E-11, -10468.7483, 8.7407E-12, -2.3354E-11, -9399.403, 9.1913E-12,
    -2.1569E-11, -8337.6485, 9.6419E-12, -1.9784E-11, -7285.0029, 1.0092E-11,
    -1.7999E-11, -6242.9842, 1.0543E-11, -1.6339E-10, 3.0085E-11, -20761.5025,
    -1.5678E-10, 2.7981E-11, -20218.7908, -1.5016E-10, 2.5877E-11, -19676.0791,
    -1.4355E-10, 2.3772E-11, -19133.3673, -1.3694E-10, 2.1668E-11, -18590.6556,
    -1.3032E-10, 1.9564E-11, -18047.9439, -1.2371E-10, 1.746E-11, -17505.2322,
    -1.1709E-10, 1.5355E-11, -16962.5205, -1.1048E-10, 1.3251E-11, -16419.8087,
    -1.0386E-10, 1.1147E-11, -15877.097, -9.7246E-11, 9.0426E-12, -15334.3853,
    -9.0632E-11, 6.9384E-12, -34791.6736, -8.4017E-11, 4.8341E-12, -13707.1469,
    -7.7402E-11, 2.7299E-12, -12624.1383, -7.0787E-11, 6.2564E-13, -11544.1661,
    -6.4172E-11, -1.4786E-12, -10468.7483, -5.7558E-11, -3.5829E-12, -9399.403,
    -5.0943E-11, -5.6871E-12, -8337.6485, -4.4328E-11, -7.7913E-12, -7285.0029,
    -3.7713E-11, -9.8956E-12, -6242.9842, -19138.5294, -4.2734E-11, -1.538E-10,
    -18644.7674, -4.1621E-11, -1.4745E-10, -18151.0053, -4.0508E-11, -1.4111E-10,
    -17657.2433, -3.9396E-11, -1.3477E-10, -17163.4812, -3.8283E-11, -1.2842E-10,
    -16669.7192, -3.717E-11, -1.2208E-10, -16175.9571, -3.6058E-11, -1.1574E-10,
    -15682.1951, -3.4945E-11, -1.0939E-10, -15188.433, -3.3832E-11, -1.0305E-10,
    -14694.671, -3.272E-11, -9.6704E-11, -14200.9089, -3.1607E-11, -9.036E-11,
    -13707.1469, -3.0494E-11, -8.4017E-11, -33213.3848, -2.9382E-11, -7.7673E-11,
    -12177.8078, -2.8269E-11, -7.1329E-11, -11143.7489, -2.7156E-11, -6.4986E-11,
    -10112.7264, -2.6044E-11, -5.8642E-11, -9086.2582, -2.4931E-11, -5.2298E-11,
    -8065.8627, -2.3818E-11, -4.5954E-11, -7053.0578, -2.2706E-11, -3.9611E-11,
    -6049.3618, -2.1593E-11, -3.3267E-11, -5.0333E-11, -19138.5294, -7.7515E-13,
    -4.8587E-11, -18644.7674, -2.6521E-13, -4.6841E-11, -18151.0053, 2.4472E-13,
    -4.5096E-11, -17657.2433, 7.5465E-13, -4.335E-11, -17163.4812, 1.2646E-12,
    -4.1604E-11, -16669.7192, 1.7745E-12, -3.9858E-11, -16175.9571, 2.2845E-12,
    -3.8112E-11, -15682.1951, 2.7944E-12, -3.6366E-11, -15188.433, 3.3043E-12,
    -3.462E-11, -14694.671, 3.8143E-12, -3.2874E-11, -14200.9089, 4.3242E-12,
    -3.1128E-11, -13707.1469, 4.8341E-12, -2.9382E-11, -33213.3848, 5.3441E-12,
    -2.7636E-11, -12177.8078, 5.854E-12, -2.589E-11, -11143.7489, 6.3639E-12,
    -2.4144E-11, -10112.7264, 6.8739E-12, -2.2398E-11, -9086.2582, 7.3838E-12,
    -2.0652E-11, -8065.8627, 7.8937E-12, -1.8906E-11, -7053.0578, 8.4037E-12,
    -1.716E-11, -6049.3618, 8.9136E-12, -1.5089E-10, 2.9883E-11, -19138.5294,
    -1.4479E-10, 2.7838E-11, -18644.7674, -1.3869E-10, 2.5793E-11, -18151.0053,
    -1.3258E-10, 2.3748E-11, -17657.2433, -1.2648E-10, 2.1703E-11, -17163.4812,
    -1.2038E-10, 1.9658E-11, -16669.7192, -1.1428E-10, 1.7613E-11, -16175.9571,
    -1.0818E-10, 1.5569E-11, -15682.1951, -1.0208E-10, 1.3524E-11, -15188.433,
    -9.5977E-11, 1.1479E-11, -14694.671, -8.9876E-11, 9.4339E-12, -14200.9089,
    -8.3774E-11, 7.389E-12, -13707.1469, -7.7673E-11, 5.3441E-12, -33213.3848,
    -7.1572E-11, 3.2992E-12, -12177.8078, -6.5471E-11, 1.2543E-12, -11143.7489,
    -5.9369E-11, -7.9063E-13, -10112.7264, -5.3268E-11, -2.8355E-12, -9086.2582,
    -4.7167E-11, -4.8804E-12, -8065.8627, -4.1065E-11, -6.9253E-12, -7053.0578,
    -3.4964E-11, -8.9702E-12, -6049.3618, -17533.7741, -4.0519E-11, -1.4153E-10,
    -17087.4436, -3.9445E-11, -1.357E-10, -16641.113, -3.8372E-11, -1.2987E-10,
    -16194.7825, -3.7298E-11, -1.2404E-10, -15748.452, -3.6225E-11, -1.1821E-10,
    -15302.1215, -3.5151E-11, -1.1238E-10, -14855.7909, -3.4077E-11, -1.0655E-10,
    -14409.4604, -3.3004E-11, -1.0072E-10, -13963.1299, -3.193E-11, -9.4893E-11,
    -13516.7994, -3.0857E-11, -8.9062E-11, -13070.4689, -2.9783E-11, -8.3232E-11,
    -12624.1383, -2.8709E-11, -7.7402E-11, -12177.8078, -2.7636E-11, -7.1572E-11,
    -31731.4773, -2.6562E-11, -6.5742E-11, -10743.3318, -2.5489E-11, -5.9911E-11,
    -9756.7045, -2.4415E-11, -5.4081E-11, -8773.1134, -2.3341E-11, -4.8251E-11,
    -7794.0768, -2.2268E-11, -4.2421E-11, -6821.1128, -2.1194E-11, -3.6591E-11,
    -5855.7394, -2.012E-11, -3.076E-11, -4.8752E-11, -17533.7741, -3.5322E-12,
    -4.7045E-11, -17087.4436, -2.9629E-12, -4.5338E-11, -16641.113, -2.3936E-12,
    -4.3632E-11, -16194.7825, -1.8244E-12, -4.1925E-11, -15748.452, -1.2551E-12,
    -4.0218E-11, -15302.1215, -6.858E-13, -3.8511E-11, -14855.7909, -1.1652E-13,
    -3.6804E-11, -14409.4604, 4.5276E-13, -3.5097E-11, -13963.1299, 1.022E-12,
    -3.339E-11, -13516.7994, 1.5913E-12, -3.1683E-11, -13070.4689, 2.1606E-12,
    -2.9976E-11, -12624.1383, 2.7299E-12, -2.8269E-11, -12177.8078, 3.2992E-12,
    -2.6562E-11, -31731.4773, 3.8684E-12, -2.4855E-11, -10743.3318, 4.4377E-12,
    -2.3148E-11, -9756.7045, 5.007E-12, -2.1441E-11, -8773.1134, 5.5763E-12,
    -1.9734E-11, -7794.0768, 6.1456E-12, -1.8027E-11, -6821.1128, 6.7148E-12,
    -1.6321E-11, -5855.7394, 7.2841E-12, -1.3838E-10, 2.9681E-11, -17533.7741,
    -1.3279E-10, 2.7695E-11, -17087.4436, -1.2721E-10, 2.571E-11, -16641.113,
    -1.2162E-10, 2.3724E-11, -16194.7825, -1.1603E-10, 2.1738E-11, -15748.452,
    -1.1044E-10, 1.9753E-11, -15302.1215, -1.0486E-10, 1.7767E-11, -14855.7909,
    -9.9268E-11, 1.5782E-11, -14409.4604, -9.368E-11, 1.3796E-11, -13963.1299,
    -8.8092E-11, 1.1811E-11, -13516.7994, -8.2505E-11, 9.8251E-12, -13070.4689,
    -7.6917E-11, 7.8395E-12, -12624.1383, -7.1329E-11, 5.854E-12, -12177.8078,
    -6.5742E-11, 3.8684E-12, -31731.4773, -6.0154E-11, 1.8829E-12, -10743.3318,
    -5.4566E-11, -1.0266E-13, -9756.7045, -4.8978E-11, -2.0882E-12, -8773.1134,
    -4.3391E-11, -4.0738E-12, -7794.0768, -3.7803E-11, -6.0593E-12, -6821.1128,
    -3.2215E-11, -8.0449E-12, -5855.7394, -15948.7546, -3.8305E-11, -1.2927E-10,
    -15548.3375, -3.727E-11, -1.2395E-10, -15147.9203, -3.6236E-11, -1.1864E-10,
    -14747.5032, -3.5201E-11, -1.1332E-10, -14347.0861, -3.4166E-11, -1.08E-10,
    -13946.6689, -3.3132E-11, -1.0269E-10, -13546.2518, -3.2097E-11, -9.737E-11,
    -13145.8346, -3.1063E-11, -9.2054E-11, -12745.4175, -3.0028E-11, -8.6737E-11,
    -12345.0004, -2.8994E-11, -8.142E-11, -11944.5832, -2.7959E-11, -7.6104E-11,
    -11544.1661, -2.6924E-11, -7.0787E-11, -11143.7489, -2.589E-11, -6.5471E-11,
    -10743.3318, -2.4855E-11, -6.0154E-11, -30342.9147, -2.3821E-11, -5.4837E-11,
    -9400.6826, -2.2786E-11, -4.9521E-11, -8459.9686, -2.1751E-11, -4.4204E-11,
    -7522.2909, -2.0717E-11, -3.8887E-11, -6589.1677, -1.9682E-11, -3.3571E-11,
    -5662.1171, -1.8648E-11, -2.8254E-11, -4.7171E-11, -15948.7546, -6.2893E-12,
    -4.5503E-11, -15548.3375, -5.6606E-12, -4.3835E-11, -15147.9203, -5.032E-12,
    -4.2168E-11, -14747.5032, -4.4034E-12, -4.05E-11, -14347.0861, -3.7748E-12,
    -3.8832E-11, -13946.6689, -3.1461E-12, -3.7164E-11, -13546.2518, -2.5175E-12,
    -3.5496E-11, -13145.8346, -1.8889E-12, -3.3828E-11, -12745.4175, -1.2602E-12,
    -3.216E-11, -12345.0004, -6.3162E-13, -3.0492E-11, -11944.5832, -2.9897E-15,
    -2.8824E-11, -11544.1661, 6.2564E-13, -2.7156E-11, -11143.7489, 1.2543E-12,
    -2.5489E-11, -10743.3318, 1.8829E-12, -2.3821E-11, -30342.9147, 2.5115E-12,
    -2.2153E-11, -9400.6826, 3.1401E-12, -2.0485E-11, -8459.9686, 3.7688E-12,
    -1.8817E-11, -7522.2909, 4.3974E-12, -1.7149E-11, -6589.1677, 5.026E-12,
    -1.5481E-11, -5662.1171, 5.6547E-12, -1.2588E-10, 2.9478E-11, -15948.7546,
    -1.208E-10, 2.7552E-11, -15548.3375, -1.1573E-10, 2.5626E-11, -15147.9203,
    -1.1065E-10, 2.37E-11, -14747.5032, -1.0558E-10, 2.1774E-11, -14347.0861,
    -1.005E-10, 1.9847E-11, -13946.6689, -9.5431E-11, 1.7921E-11, -13546.2518,
    -9.0356E-11, 1.5995E-11, -13145.8346, -8.5282E-11, 1.4069E-11, -12745.4175,
    -8.0208E-11, 1.2143E-11, -12345.0004, -7.5134E-11, 1.0216E-11, -11944.5832,
    -7.006E-11, 8.2901E-12, -11544.1661, -6.4986E-11, 6.3639E-12, -11143.7489,
    -5.9911E-11, 4.4377E-12, -10743.3318, -5.4837E-11, 2.5115E-12, -30342.9147,
    -4.9763E-11, 5.8531E-13, -9400.6826, -4.4689E-11, -1.3409E-12, -8459.9686,
    -3.9615E-11, -3.2671E-12, -7522.2909, -3.4541E-11, -5.1933E-12, -6589.1677,
    -2.9466E-11, -7.1195E-12, -5662.1171, -14384.9892, -3.609E-11, -1.1701E-10,
    -14028.9673, -3.5095E-11, -1.122E-10, -13672.9454, -3.4099E-11, -1.074E-10,
    -13316.9235, -3.3104E-11, -1.026E-10, -12960.9016, -3.2108E-11, -9.7794E-11,
    -12604.8797, -3.1113E-11, -9.2991E-11, -12248.8578, -3.0117E-11, -8.8188E-11,
    -11892.8359, -2.9121E-11, -8.3385E-11, -11536.814, -2.8126E-11, -7.8582E-11,
    -11180.7921, -2.713E-11, -7.3779E-11, -10824.7702, -2.6135E-11, -6.8975E-11,
    -10468.7483, -2.5139E-11, -6.4172E-11, -10112.7264, -2.4144E-11, -5.9369E-11,
    -9756.7045, -2.3148E-11, -5.4566E-11, -9400.6826, -2.2153E-11, -4.9763E-11,
    -29044.6607, -2.1157E-11, -4.496E-11, -8146.8238, -2.0162E-11, -4.0157E-11,
    -7250.5051, -1.9166E-11, -3.5354E-11, -6357.2227, -1.8171E-11, -3.0551E-11,
    -5468.4947, -1.7175E-11, -2.5748E-11, -4.559E-11, -14384.9892, -9.0463E-12,
    -4.3961E-11, -14028.9673, -8.3583E-12, -4.2332E-11, -13672.9454, -7.6704E-12,
    -4.0704E-11, -13316.9235, -6.9824E-12, -3.9075E-11, -12960.9016, -6.2944E-12,
    -3.7446E-11, -12604.8797, -5.6064E-12, -3.5817E-11, -12248.8578, -4.9185E-12,
    -3.4188E-11, -11892.8359, -4.2305E-12, -3.2559E-11, -11536.814, -3.5425E-12,
    -3.093E-11, -11180.7921, -2.8546E-12, -2.9302E-11, -10824.7702, -2.1666E-12,
    -2.7673E-11, -10468.7483, -1.4786E-12, -2.6044E-11, -10112.7264, -7.9063E-13,
    -2.4415E-11, -9756.7045, -1.0266E-13, -2.2786E-11, -9400.6826, 5.8531E-13,
    -2.1157E-11, -29044.6607, 1.2733E-12, -1.9528E-11, -8146.8238, 1.9613E-12,
    -1.7899E-11, -7250.5051, 2.6492E-12, -1.6271E-11, -6357.2227, 3.3372E-12,
    -1.4642E-11, -5468.4947, 4.0252E-12, -1.1337E-10, 2.9276E-11, -14384.9892,
    -1.0881E-10, 2.7409E-11, -14028.9673, -1.0425E-10, 2.5542E-11, -13672.9454,
    -9.9687E-11, 2.3676E-11, -13316.9235, -9.5127E-11, 2.1809E-11, -12960.9016,
    -9.0566E-11, 1.9942E-11, -12604.8797, -8.6006E-11, 1.8075E-11, -12248.8578,
    -8.1445E-11, 1.6208E-11, -11892.8359, -7.6884E-11, 1.4341E-11, -11536.814,
    -7.2324E-11, 1.2474E-11, -11180.7921, -6.7763E-11, 1.0608E-11, -10824.7702,
    -6.3202E-11, 8.7407E-12, -10468.7483, -5.8642E-11, 6.8739E-12, -10112.7264,
    -5.4081E-11, 5.007E-12, -9756.7045, -4.9521E-11, 3.1401E-12, -9400.6826,
    -4.496E-11, 1.2733E-12, -29044.6607, -4.0399E-11, -5.9357E-13, -8146.8238,
    -3.5839E-11, -2.4604E-12, -7250.5051, -3.1278E-11, -4.3273E-12, -6357.2227,
    -2.6717E-11, -6.1941E-12, -5468.4947, -12843.9959, -3.3876E-11, -1.0474E-10,
    -12530.8511, -3.2919E-11, -1.0045E-10, -12217.7063, -3.1963E-11, -9.6163E-11,
    -11904.5615, -3.1006E-11, -9.1874E-11, -11591.4167, -3.005E-11, -8.7584E-11,
    -11278.2719, -2.9093E-11, -8.3295E-11, -10965.1271, -2.8137E-11, -7.9005E-11,
    -10651.9823, -2.718E-11, -7.4716E-11, -10338.8374, -2.6224E-11, -7.0426E-11,
    -10025.6926, -2.5267E-11, -6.6137E-11, -9712.5478, -2.4311E-11, -6.1847E-11,
    -9399.403, -2.3354E-11, -5.7558E-11, -9086.2582, -2.2398E-11, -5.3268E-11,
    -8773.1134, -2.1441E-11, -4.8978E-11, -8459.9686, -2.0485E-11, -4.4689E-11,
    -8146.8238, -1.9528E-11, -4.0399E-11, -27833.679, -1.8572E-11, -3.611E-11,
    -6978.7192, -1.7615E-11, -3.182E-11, -6125.2776, -1.6659E-11, -2.7531E-11,
    -5274.8723, -1.5702E-11, -2.3241E-11, -4.4009E-11, -12843.9959, -1.1803E-11,
    -4.2419E-11, -12530.8511, -1.1056E-11, -4.0829E-11, -12217.7063, -1.0309E-11,
    -3.924E-11, -11904.5615, -9.5614E-12, -3.765E-11, -11591.4167, -8.8141E-12,
    -3.606E-11, -11278.2719, -8.0668E-12, -3.447E-11, -10965.1271, -7.3194E-12,
    -3.288E-11, -10651.9823, -6.5721E-12, -3.129E-11, -10338.8374, -5.8248E-12,
    -2.9701E-11, -10025.6926, -5.0775E-12, -2.8111E-11, -9712.5478, -4.3302E-12,
    -2.6521E-11, -9399.403, -3.5829E-12, -2.4931E-11, -9086.2582, -2.8355E-12,
    -2.3341E-11, -8773.1134, -2.0882E-12, -2.1751E-11, -8459.9686, -1.3409E-12,
    -2.0162E-11, -8146.8238, -5.9357E-13, -1.8572E-11, -27833.679, 1.5375E-13,
    -1.6982E-11, -6978.7192, 9.0107E-13, -1.5392E-11, -6125.2776, 1.6484E-12,
    -1.3802E-11, -5274.8723, 2.3957E-12, -1.0086E-10, 2.9074E-11, -12843.9959,
    -9.6816E-11, 2.7266E-11, -12530.8511, -9.2769E-11, 2.5459E-11, -12217.7063,
    -8.8722E-11, 2.3651E-11, -11904.5615, -8.4675E-11, 2.1844E-11, -11591.4167,
    -8.0628E-11, 2.0036E-11, -11278.2719, -7.6581E-11, 1.8229E-11, -10965.1271,
    -7.2534E-11, 1.6421E-11, -10651.9823, -6.8486E-11, 1.4614E-11, -10338.8374,
    -6.4439E-11, 1.2806E-11, -10025.6926, -6.0392E-11, 1.0999E-11, -9712.5478,
    -5.6345E-11, 9.1913E-12, -9399.403, -5.2298E-11, 7.3838E-12, -9086.2582,
    -4.8251E-11, 5.5763E-12, -8773.1134, -4.4204E-11, 3.7688E-12, -8459.9686,
    -4.0157E-11, 1.9613E-12, -8146.8238, -3.611E-11, 1.5375E-13, -27833.679,
    -3.2063E-11, -1.6538E-12, -6978.7192, -2.8016E-11, -3.4613E-12, -6125.2776,
    -2.3969E-11, -5.2688E-12, -5274.8723, -11327.2929, -3.1661E-11, -9.2479E-11,
    -11055.5071, -3.0744E-11, -8.8703E-11, -10783.7212, -2.9826E-11, -8.4927E-11,
    -10511.9354, -2.8909E-11, -8.1151E-11, -10240.1495, -2.7992E-11, -7.7375E-11,
    -9968.3636, -2.7074E-11, -7.3599E-11, -9696.5778, -2.6157E-11, -6.9823E-11,
    -9424.7919, -2.5239E-11, -6.6047E-11, -9153.0061, -2.4322E-11, -6.2271E-11,
    -8881.2202, -2.3404E-11, -5.8495E-11, -8609.4344, -2.2487E-11, -5.4719E-11,
    -8337.6485, -2.1569E-11, -5.0943E-11, -8065.8627, -2.0652E-11, -4.7167E-11,
    -7794.0768, -1.9734E-11, -4.3391E-11, -7522.2909, -1.8817E-11, -3.9615E-11,
    -7250.5051, -1.7899E-11, -3.5839E-11, -6978.7192, -1.6982E-11, -3.2063E-11,
    -26706.9334, -1.6065E-11, -2.8287E-11, -5893.3326, -1.5147E-11, -2.4511E-11,
    -5081.2499, -1.423E-11, -2.0735E-11, -4.2428E-11, -11327.2929, -1.456E-11,
    -4.0877E-11, -11055.5071, -1.3754E-11, -3.9326E-11, -10783.7212, -1.2947E-11,
    -3.7776E-11, -10511.9354, -1.214E-11, -3.6225E-11, -10240.1495, -1.1334E-11,
    -3.4674E-11, -9968.3636, -1.0527E-11, -3.3123E-11, -9696.5778, -9.7204E-12,
    -3.1572E-11, -9424.7919, -8.9138E-12, -3.0022E-11, -9153.0061, -8.1071E-12,
    -2.8471E-11, -8881.2202, -7.3004E-12, -2.692E-11, -8609.4344, -6.4938E-12,
    -2.5369E-11, -8337.6485, -5.6871E-12, -2.3818E-11, -8065.8627, -4.8804E-12,
    -2.2268E-11, -7794.0768, -4.0738E-12, -2.0717E-11, -7522.2909, -3.2671E-12,
    -1.9166E-11, -7250.5051, -2.4604E-12, -1.7615E-11, -6978.7192, -1.6538E-12,
    -1.6065E-11, -26706.9334, -8.471E-13, -1.4514E-11, -5893.3326, -4.0431E-14,
    -1.2963E-11, -5081.2499, 7.6623E-13, -8.8357E-11, 2.8872E-11, -11327.2929,
    -8.4823E-11, 2.7124E-11, -11055.5071, -8.129E-11, 2.5375E-11, -10783.7212,
    -7.7756E-11, 2.3627E-11, -10511.9354, -7.4223E-11, 2.1879E-11, -10240.1495,
    -7.0689E-11, 2.0131E-11, -9968.3636, -6.7156E-11, 1.8383E-11, -9696.5778,
    -6.3622E-11, 1.6635E-11, -9424.7919, -6.0089E-11, 1.4886E-11, -9153.0061,
    -5.6555E-11, 1.3138E-11, -8881.2202, -5.3021E-11, 1.139E-11, -8609.4344,
    -4.9488E-11, 9.6419E-12, -8337.6485, -4.5954E-11, 7.8937E-12, -8065.8627,
    -4.2421E-11, 6.1456E-12, -7794.0768, -3.8887E-11, 4.3974E-12, -7522.2909,
    -3.5354E-11, 2.6492E-12, -7250.5051, -3.182E-11, 9.0107E-13, -6978.7192,
    -2.8287E-11, -8.471E-13, -26706.9334, -2.4753E-11, -2.5953E-12, -5893.3326,
    -2.122E-11, -4.3434E-12, -5081.2499, -9836.3984, -2.9447E-11, -8.0215E-11,
    -9604.4534, -2.8569E-11, -7.6953E-11, -9372.5083, -2.769E-11, -7.369E-11,
    -9140.5633, -2.6812E-11, -7.0428E-11, -8908.6182, -2.5933E-11, -6.7165E-11,
    -8676.6732, -2.5055E-11, -6.3903E-11, -8444.7281, -2.4176E-11, -6.064E-11,
    -8212.7831, -2.3298E-11, -5.7378E-11, -7980.838, -2.242E-11, -5.4115E-11,
    -7748.893, -2.1541E-11, -5.0853E-11, -7516.9479, -2.0663E-11, -4.759E-11,
    -7285.0029, -1.9784E-11, -4.4328E-11, -7053.0578, -1.8906E-11, -4.1065E-11,
    -6821.1128, -1.8027E-11, -3.7803E-11, -6589.1677, -1.7149E-11, -3.4541E-11,
    -6357.2227, -1.6271E-11, -3.1278E-11, -6125.2776, -1.5392E-11, -2.8016E-11,
    -5893.3326, -1.4514E-11, -2.4753E-11, -25661.3875, -1.3635E-11, -2.1491E-11,
    -4887.6275, -1.2757E-11, -1.8228E-11, -4.0847E-11, -9836.3984, -1.7317E-11,
    -3.9335E-11, -9604.4534, -1.6451E-11, -3.7823E-11, -9372.5083, -1.5585E-11,
    -3.6312E-11, -9140.5633, -1.4719E-11, -3.48E-11, -8908.6182, -1.3853E-11,
    -3.3288E-11, -8676.6732, -1.2987E-11, -3.1776E-11, -8444.7281, -1.2121E-11,
    -3.0265E-11, -8212.7831, -1.1255E-11, -2.8753E-11, -7980.838, -1.0389E-11,
    -2.7241E-11, -7748.893, -9.5234E-12, -2.5729E-11, -7516.9479, -8.6574E-12,
    -2.4218E-11, -7285.0029, -7.7913E-12, -2.2706E-11, -7053.0578, -6.9253E-12,
    -2.1194E-11, -6821.1128, -6.0593E-12, -1.9682E-11, -6589.1677, -5.1933E-12,
    -1.8171E-11, -6357.2227, -4.3273E-12, -1.6659E-11, -6125.2776, -3.4613E-12,
    -1.5147E-11, -5893.3326, -2.5953E-12, -1.3635E-11, -25661.3875, -1.7293E-12,
    -1.2124E-11, -4887.6275, -8.6324E-13, -7.5851E-11, 2.8669E-11, -9836.3984,
    -7.2831E-11, 2.6981E-11, -9604.4534, -6.9811E-11, 2.5292E-11, -9372.5083,
    -6.6791E-11, 2.3603E-11, -9140.5633, -6.3771E-11, 2.1914E-11, -8908.6182,
    -6.0751E-11, 2.0225E-11, -8676.6732, -5.7731E-11, 1.8537E-11, -8444.7281,
    -5.4711E-11, 1.6848E-11, -8212.7831, -5.1691E-11, 1.5159E-11, -7980.838,
    -4.8671E-11, 1.347E-11, -7748.893, -4.5651E-11, 1.1781E-11, -7516.9479,
    -4.2631E-11, 1.0092E-11, -7285.0029, -3.9611E-11, 8.4037E-12, -7053.0578,
    -3.6591E-11, 6.7148E-12, -6821.1128, -3.3571E-11, 5.026E-12, -6589.1677,
    -3.0551E-11, 3.3372E-12, -6357.2227, -2.7531E-11, 1.6484E-12, -6125.2776,
    -2.4511E-11, -4.0431E-14, -5893.3326, -2.1491E-11, -1.7293E-12, -25661.3875,
    -1.8471E-11, -3.4181E-12, -4887.6275, -8372.8305, -2.7232E-11, -6.7951E-11,
    -8179.2081, -2.6393E-11, -6.5202E-11, -7985.5857, -2.5554E-11, -6.2453E-11,
    -7791.9633, -2.4714E-11, -5.9704E-11, -7598.341, -2.3875E-11, -5.6956E-11,
    -7404.7186, -2.3036E-11, -5.4207E-11, -7211.0962, -2.2196E-11, -5.1458E-11,
    -7017.4738, -2.1357E-11, -4.8709E-11, -6823.8514, -2.0517E-11, -4.596E-11,
    -6630.229, -1.9678E-11, -4.3211E-11, -6436.6066, -1.8839E-11, -4.0462E-11,
    -6242.9842, -1.7999E-11, -3.7713E-11, -6049.3618, -1.716E-11, -3.4964E-11,
    -5855.7394, -1.6321E-11, -3.2215E-11, -5662.1171, -1.5481E-11, -2.9466E-11,
    -5468.4947, -1.4642E-11, -2.6717E-11, -5274.8723, -1.3802E-11, -2.3969E-11,
    -5081.2499, -1.2963E-11, -2.122E-11, -4887.6275, -1.2124E-11, -1.8471E-11,
    -24694.0051, -1.1284E-11, -1.5722E-11, -3.9266E-11, -8372.8305, -2.0075E-11,
    -3.7793E-11, -8179.2081, -1.9149E-11, -3.632E-11, -7985.5857, -1.8224E-11,
    -3.4848E-11, -7791.9633, -1.7298E-11, -3.3375E-11, -7598.341, -1.6373E-11,
    -3.1902E-11, -7404.7186, -1.5448E-11, -3.0429E-11, -7211.0962, -1.4522E-11,
    -2.8957E-11, -7017.4738, -1.3597E-11, -2.7484E-11, -6823.8514, -1.2672E-11,
    -2.6011E-11, -6630.229, -1.1746E-11, -2.4539E-11, -6436.6066, -1.0821E-11,
    -2.3066E-11, -6242.9842, -9.8956E-12, -2.1593E-11, -6049.3618, -8.9702E-12,
    -2.012E-11, -5855.7394, -8.0449E-12, -1.8648E-11, -5662.1171, -7.1195E-12,
    -1.7175E-11, -5468.4947, -6.1941E-12, -1.5702E-11, -5274.8723, -5.2688E-12,
    -1.423E-11, -5081.2499, -4.3434E-12, -1.2757E-11, -4887.6275, -3.4181E-12,
    -1.1284E-11, -24694.0051, -2.4927E-12, -6.3344E-11, 2.8467E-11, -8372.8305,
    -6.0838E-11, 2.6838E-11, -8179.2081, -5.8331E-11, 2.5208E-11, -7985.5857,
    -5.5825E-11, 2.3579E-11, -7791.9633, -5.3318E-11, 2.1949E-11, -7598.341,
    -5.0812E-11, 2.032E-11, -7404.7186, -4.8306E-11, 1.869E-11, -7211.0962,
    -4.5799E-11, 1.7061E-11, -7017.4738, -4.3293E-11, 1.5431E-11, -6823.8514,
    -4.0786E-11, 1.3802E-11, -6630.229, -3.828E-11, 1.2173E-11, -6436.6066,
    -3.5773E-11, 1.0543E-11, -6242.9842, -3.3267E-11, 8.9136E-12, -6049.3618,
    -3.076E-11, 7.2841E-12, -5855.7394, -2.8254E-11, 5.6547E-12, -5662.1171,
    -2.5748E-11, 4.0252E-12, -5468.4947, -2.3241E-11, 2.3957E-12, -5274.8723,
    -2.0735E-11, 7.6623E-13, -5081.2499, -1.8228E-11, -8.6324E-13, -4887.6275,
    -1.5722E-11, -2.4927E-12, -24694.0051 };

  static const real_T dv[60] = { 0.075295, -1.8179E-17, -0.075295, 0.060393,
    -1.5125E-17, -0.060393, 0.047849, -1.2479E-17, -0.047849, 0.037314,
    -1.0168E-17, -0.037314, 0.02849, -8.125E-18, -0.02849, 0.021119, -6.2903E-18,
    -0.021119, 0.014984, -4.6095E-18, -0.014984, 0.0098978, -3.0317E-18,
    -0.0098978, 0.0057007, -1.5086E-18, -0.0057007, 0.0022566, 7.1654E-21,
    -0.0022566, -0.00055045, 1.5626E-18, 0.00055045, -0.0028196, 3.2058E-18,
    0.0028196, -0.0046354, 4.9874E-18, 0.0046354, -0.0060696, 6.9611E-18,
    0.0060696, -0.0071836, 9.1859E-18, 0.0071836, -0.0080296, 1.1727E-17,
    0.0080296, -0.0086518, 1.4658E-17, 0.0086518, -0.0090878, 1.8057E-17,
    0.0090878, -0.0093698, 2.2024E-17, 0.0093698, -0.0095247, 2.6665E-17,
    0.0095247 };

  static const int8_T A[7200] = { -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 };

  static const int8_T a[7200] = { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 };

  static const int8_T c_a[7200] = { -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 };

  static const int8_T iv[720] = { 1, 0, 0, 2, 0, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 0,
    0, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 4, 0, 0,
    0, 1, 0, 0, 4, 0, 0, 0, 1, 0, 0, 4, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 1, 1, 0, 0, 6, 0, 0, 0, 1, 0, 0, 6, 0, 0, 0, 1, 0, 0, 6, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 8, 0, 0, 0, 1, 0, 0, 8, 0,
    0, 0, 1, 0, 0, 8, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0,
    0, 10, 0, 0, 0, 1, 0, 0, 10, 0, 0, 0, 1, 0, 0, 10, 0, 0, 0, 1, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 12, 0, 0, 0, 1, 0, 0, 12, 0, 0, 0, 1, 0,
    0, 12, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 14, 0,
    0, 0, 1, 0, 0, 14, 0, 0, 0, 1, 0, 0, 14, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 0, 1, 1, 0, 0, 16, 0, 0, 0, 1, 0, 0, 16, 0, 0, 0, 1, 0, 0, 16, 0,
    0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 18, 0, 0, 0, 1,
    0, 0, 18, 0, 0, 0, 1, 0, 0, 18, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 1, 1, 0, 0, 20, 0, 0, 0, 1, 0, 0, 20, 0, 0, 0, 1, 0, 0, 20, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 22, 0, 0, 0, 1, 0, 0, 22,
    0, 0, 0, 1, 0, 0, 22, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1,
    1, 0, 0, 24, 0, 0, 0, 1, 0, 0, 24, 0, 0, 0, 1, 0, 0, 24, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 26, 0, 0, 0, 1, 0, 0, 26, 0, 0, 0,
    1, 0, 0, 26, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0,
    28, 0, 0, 0, 1, 0, 0, 28, 0, 0, 0, 1, 0, 0, 28, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 30, 0, 0, 0, 1, 0, 0, 30, 0, 0, 0, 1, 0, 0,
    30, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 32, 0, 0,
    0, 1, 0, 0, 32, 0, 0, 0, 1, 0, 0, 32, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 1, 1, 0, 0, 34, 0, 0, 0, 1, 0, 0, 34, 0, 0, 0, 1, 0, 0, 34, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 36, 0, 0, 0, 1, 0,
    0, 36, 0, 0, 0, 1, 0, 0, 36, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
    0, 1, 1, 0, 0, 38, 0, 0, 0, 1, 0, 0, 38, 0, 0, 0, 1, 0, 0, 38, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 40, 0, 0, 0, 1, 0, 0, 40, 0,
    0, 0, 1, 0, 0, 40, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1 };

  real_T IGA[7200];
  real_T A_data[3600];
  double y_tmp[3600];
  real_T b_del_lam[240];
  real_T result[240];
  real_T varargin_1_data[240];
  real_T del_lam[120];
  real_T esig[120];
  real_T ftol[120];
  real_T igr2[120];
  real_T ilam[120];
  real_T lam[120];
  real_T mesil[120];
  real_T nlamt[120];
  real_T H[60];
  double del_z[60];
  real_T r1[60];
  real_T d;
  real_T mu;
  real_T mu_old;
  real_T ssq;
  int32_T b_i;
  int32_T exitflag;
  int32_T i;
  int32_T i1;
  int32_T idx;
  int32_T iy;
  int32_T unusedU0;
  uint8_T ii_data[240];
  boolean_T x[240];

  //  A=[1 0 0 2 0 0 ;0 1 0 0 2 0 ;0 0 1 0 0 2 ;0 0 0 1 0 0 ;0 0 0 0 1 0 ;0 0 0 0 0 1 ]; 
  //  B=[0.20877 0 0 ;0 0.20877 0 ;0 0 0.20877 ;0.20877 0 0 ;0 0.20877 0 ;0 0 0.20877]; 
  //  Q_MPC=[100 0 0 0 0 0 ;0 100 0 0 0 0 ;0 0 100 0 0 0 ;0 0 0 100000 0 0 ;0 0 0 0 100000 0 ;0 0 0 0 0 100000]; 
  //  % P_MPC=[1567.5765 1.5888e-12 2.7019e-12 8536.6551 1.0896e-11 3.5199e-11 ;1.5888e-12 1567.5765 9.8891e-12 6.7324e-11 8536.6551 6.1982e-11 ;2.7019e-12 9.8891e-12 1567.5765 1.4677e-11 9.7557e-11 8536.6551 ;8536.6551 6.7324e-11 1.4677e-11 423802.2092 1.1408e-09 6.5756e-10 ;1.0896e-11 8536.6551 9.7557e-11 1.1408e-09 423802.2092 6.8455e-10 ;3.5199e-11 6.1982e-11 8536.6551 6.5756e-10 6.8455e-10 423802.2092 ;]; 
  //  R_MPC=[100 0 0 ;0 100 0 ;0 0 100];
  //  V_TRMPC=[0.25752 ;0.25752 ;0.25752 ;0.25752 ;0.25752 ;0.25752];
  //  Hp = 20;
  //  %  X_QP=[0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ]; 
  //  X_QP = zeros(3*Hp,1);
  //  solution matrices x_traj = A_tilde x_0 + B_tilde u_traj
  //  n=6;
  //  n = size(A,1);
  //  m=3;
  //  m= size(B,2);
  //  A_tilde = zeros((Hp)*n,n);
  //  % A_tilde = [1 0 0 2 0 0 ;0 1 0 0 2 0 ;0 0 1 0 0 2 ;0 0 0 1 0 0 ;0 0 0 0 1 0 ;0 0 0 0 0 1 ;1 0 0 4 0 0 ;0 1 0 0 4 0 ;0 0 1 0 0 4 ;0 0 0 1 0 0 ;0 0 0 0 1 0 ;0 0 0 0 0 1 ;1 0 0 6 0 0 ;0 1 0 0 6 0 ;0 0 1 0 0 6 ;0 0 0 1 0 0 ;0 0 0 0 1 0 ;0 0 0 0 0 1 ;1 0 0 8 0 0 ;0 1 0 0 8 0 ;0 0 1 0 0 8 ;0 0 0 1 0 0 ;0 0 0 0 1 0 ;0 0 0 0 0 1 ;1 0 0 10 0 0 ;0 1 0 0 10 0 ;0 0 1 0 0 10 ;0 0 0 1 0 0 ;0 0 0 0 1 0 ;0 0 0 0 0 1 ;1 0 0 12 0 0 ;0 1 0 0 12 0 ;0 0 1 0 0 12 ;0 0 0 1 0 0 ;0 0 0 0 1 0 ;0 0 0 0 0 1 ;1 0 0 14 0 0 ;0 1 0 0 14 0 ;0 0 1 0 0 14 ;0 0 0 1 0 0 ;0 0 0 0 1 0 ;0 0 0 0 0 1 ;1 0 0 16 0 0 ;0 1 0 0 16 0 ;0 0 1 0 0 16 ;0 0 0 1 0 0 ;0 0 0 0 1 0 ;0 0 0 0 0 1 ;1 0 0 18 0 0 ;0 1 0 0 18 0 ;0 0 1 0 0 18 ;0 0 0 1 0 0 ;0 0 0 0 1 0 ;0 0 0 0 0 1 ;1 0 0 20 0 0 ;0 1 0 0 20 0 ;0 0 1 0 0 20 ;0 0 0 1 0 0 ;0 0 0 0 1 0 ;0 0 0 0 0 1 ;1 0 0 22 0 0 ;0 1 0 0 22 0 ;0 0 1 0 0 22 ;0 0 0 1 0 0 ;0 0 0 0 1 0 ;0 0 0 0 0 1 ;1 0 0 24 0 0 ;0 1 0 0 24 0 ;0 0 1 0 0 24 ;0 0 0 1 0 0 ;0 0 0 0 1 0 ;0 0 0 0 0 1 ;1 0 0 26 0 0 ;0 1 0 0 26 0 ;0 0 1 0 0 26 ;0 0 0 1 0 0 ;0 0 0 0 1 0 ;0 0 0 0 0 1 ;1 0 0 28 0 0 ;0 1 0 0 28 0 ;0 0 1 0 0 28 ;0 0 0 1 0 0 ;0 0 0 0 1 0 ;0 0 0 0 0 1 ;1 0 0 30 0 0 ;0 1 0 0 30 0 ;0 0 1 0 0 30 ;0 0 0 1 0 0 ;0 0 0 0 1 0 ;0 0 0 0 0 1 ;1 0 0 32 0 0 ;0 1 0 0 32 0 ;0 0 1 0 0 32 ;0 0 0 1 0 0 ;0 0 0 0 1 0 ;0 0 0 0 0 1 ;1 0 0 34 0 0 ;0 1 0 0 34 0 ;0 0 1 0 0 34 ;0 0 0 1 0 0 ;0 0 0 0 1 0 ;0 0 0 0 0 1 ;1 0 0 36 0 0 ;0 1 0 0 36 0 ;0 0 1 0 0 36 ;0 0 0 1 0 0 ;0 0 0 0 1 0 ;0 0 0 0 0 1 ;1 0 0 38 0 0 ;0 1 0 0 38 0 ;0 0 1 0 0 38 ;0 0 0 1 0 0 ;0 0 0 0 1 0 ;0 0 0 0 0 1 ;1 0 0 40 0 0 ;0 1 0 0 40 0 ;0 0 1 0 0 40 ;0 0 0 1 0 0 ;0 0 0 0 1 0 ;0 0 0 0 0 1 ;]; 
  //  A_tilde(1:n,:) = A;
  //  for tt=1:Hp-1
  //    A_tilde(n*tt+1:n*(tt+1),:) = A * A_tilde(n*(tt-1)+1:n*tt,:);
  //  end
  //  B_tilde = zeros(length(A)*Hp,size(B,2)*Hp);
  //
  //  for ii = 1:Hp
  //      for jj = 1:ii
  //              B_tilde(n*(ii-1)+1:n*ii,m*(jj-1)+1:m*(jj)) = (A^(ii-jj))*B;
  //      end
  //  end
  //  % B_tilde=[0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;2.4802 0 0 2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 2.4802 0 0 2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 2.4802 0 0 2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;2.6786 0 0 2.4802 0 0 2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 2.6786 0 0 2.4802 0 0 2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 2.6786 0 0 2.4802 0 0 2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;2.877 0 0 2.6786 0 0 2.4802 0 0 2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 2.877 0 0 2.6786 0 0 2.4802 0 0 2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 2.877 0 0 2.6786 0 0 2.4802 0 0 2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;3.0754 0 0 2.877 0 0 2.6786 0 0 2.4802 0 0 2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 3.0754 0 0 2.877 0 0 2.6786 0 0 2.4802 0 0 2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 3.0754 0 0 2.877 0 0 2.6786 0 0 2.4802 0 0 2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 ;0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 ;3.2738 0 0 3.0754 0 0 2.877 0 0 2.6786 0 0 2.4802 0 0 2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 ;0 3.2738 0 0 3.0754 0 0 2.877 0 0 2.6786 0 0 2.4802 0 0 2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 ;0 0 3.2738 0 0 3.0754 0 0 2.877 0 0 2.6786 0 0 2.4802 0 0 2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 ;0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 ;0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 ;0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 ;3.4722 0 0 3.2738 0 0 3.0754 0 0 2.877 0 0 2.6786 0 0 2.4802 0 0 2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 ;0 3.4722 0 0 3.2738 0 0 3.0754 0 0 2.877 0 0 2.6786 0 0 2.4802 0 0 2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 ;0 0 3.4722 0 0 3.2738 0 0 3.0754 0 0 2.877 0 0 2.6786 0 0 2.4802 0 0 2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 ;0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 ;0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 ;0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 ;3.6706 0 0 3.4722 0 0 3.2738 0 0 3.0754 0 0 2.877 0 0 2.6786 0 0 2.4802 0 0 2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 ;0 3.6706 0 0 3.4722 0 0 3.2738 0 0 3.0754 0 0 2.877 0 0 2.6786 0 0 2.4802 0 0 2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 ;0 0 3.6706 0 0 3.4722 0 0 3.2738 0 0 3.0754 0 0 2.877 0 0 2.6786 0 0 2.4802 0 0 2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 ;0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 ;0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 ;0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 ;3.869 0 0 3.6706 0 0 3.4722 0 0 3.2738 0 0 3.0754 0 0 2.877 0 0 2.6786 0 0 2.4802 0 0 2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 ;0 3.869 0 0 3.6706 0 0 3.4722 0 0 3.2738 0 0 3.0754 0 0 2.877 0 0 2.6786 0 0 2.4802 0 0 2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 ;0 0 3.869 0 0 3.6706 0 0 3.4722 0 0 3.2738 0 0 3.0754 0 0 2.877 0 0 2.6786 0 0 2.4802 0 0 2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 ;0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 ;0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 ;0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 ;]; 
  //  cost
  //  Q_tilde = blkdiag(kron(eye(Hp-1),Q_MPC),P_MPC);
  //  % Q_tilde = [70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1567.5765 1.5888e-12 2.7019e-12 8536.6551 1.0896e-11 3.5199e-11 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.5888e-12 1567.5765 9.8891e-12 6.7324e-11 8536.6551 6.1982e-11 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2.7019e-12 9.8891e-12 1567.5765 1.4677e-11 9.7557e-11 8536.6551 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8536.6551 6.7324e-11 1.4677e-11 423802.2092 1.1408e-09 6.5756e-10 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.0896e-11 8536.6551 9.7557e-11 1.1408e-09 423802.2092 6.8455e-10 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3.5199e-11 6.1982e-11 8536.6551 6.5756e-10 6.8455e-10 423802.2092 ;]; 
  //  R_tilde = kron(eye(Hp),R_MPC);
  //  R_tilde=[100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 ]; 
  //  F_tilde = kron([ones(m-1,1);zeros(m*(Hp-1)+1,1)],F_MPC);
  //  L3 = blkdiag(kron(eye(Hp-1),Q_MPC),zeros(n,n));
  //  % L3=[70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;]; 
  //  L4 = zeros(n,Hp*n);
  //  L4(:,end-(n-1):end) = eye(n);
  //  % L4 = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 ;]; 
  //  PHI = blkdiag(kron(eye(Hp),eye(6,6)));
  //  % PHI = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 ;]; 
  //  constraints
  //  CONTROL INPUTS
  //  Define and Solve QP
  //  GAMMA = zeros(n*Hp,1);%Dimension of GAMMA is n*horizon x n
  //  for i = 1:Hp
  //      GAMMA(i*n-5:i*n,1) = xfinal;
  //  end
  //  % GAMMA = [0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;]; 
  //  QQ = 2*R_tilde + 2*OMEGA'*Q_tilde*OMEGA;
  //  QQ = [488965.6068 0 0 462030.3957 0 0 435330.0526 0 0 408699.4453 0 0 382173.4419 0 0 355786.9101 0 0 329574.718 0 0 303571.7335 0 0 277812.8244 0 0 252332.8589 0 0 227166.7047 0 0 202349.2298 0 0 177915.3021 0 0 153899.7896 0 0 130337.5602 0 0 107263.4818 0 0 84712.4224 0 0 62719.2499 0 0 41318.8321 0 0 20546.0372 0 0 ;0 488965.6068 0 0 462030.3957 0 0 435330.0526 0 0 408699.4453 0 0 382173.4419 0 0 355786.9101 0 0 329574.718 0 0 303571.7335 0 0 277812.8244 0 0 252332.8589 0 0 227166.7047 0 0 202349.2298 0 0 177915.3021 0 0 153899.7896 0 0 130337.5602 0 0 107263.4818 0 0 84712.4224 0 0 62719.2499 0 0 41318.8321 0 0 20546.0372 0 ;0 0 488965.6068 0 0 462030.3957 0 0 435330.0526 0 0 408699.4453 0 0 382173.4419 0 0 355786.9101 0 0 329574.718 0 0 303571.7335 0 0 277812.8244 0 0 252332.8589 0 0 227166.7047 0 0 202349.2298 0 0 177915.3021 0 0 153899.7896 0 0 130337.5602 0 0 107263.4818 0 0 84712.4224 0 0 62719.2499 0 0 41318.8321 0 0 20546.0372 ;462030.3957 0 0 445435.4433 0 0 419732.2252 0 0 394263.875 0 0 368865.2607 0 0 343571.2502 0 0 318416.7114 0 0 293436.5123 0 0 268665.5207 0 0 244138.6047 0 0 219890.6321 0 0 195956.4709 0 0 172370.989 0 0 149169.0543 0 0 126385.5347 0 0 104055.2983 0 0 82213.2129 0 0 60894.1465 0 0 40132.9669 0 0 19964.5421 0 0 ;0 462030.3957 0 0 445435.4433 0 0 419732.2252 0 0 394263.875 0 0 368865.2607 0 0 343571.2502 0 0 318416.7114 0 0 293436.5123 0 0 268665.5207 0 0 244138.6047 0 0 219890.6321 0 0 195956.4709 0 0 172370.989 0 0 149169.0543 0 0 126385.5347 0 0 104055.2983 0 0 82213.2129 0 0 60894.1465 0 0 40132.9669 0 0 19964.5421 0 ;0 0 462030.3957 0 0 445435.4433 0 0 419732.2252 0 0 394263.875 0 0 368865.2607 0 0 343571.2502 0 0 318416.7114 0 0 293436.5123 0 0 268665.5207 0 0 244138.6047 0 0 219890.6321 0 0 195956.4709 0 0 172370.989 0 0 149169.0543 0 0 126385.5347 0 0 104055.2983 0 0 82213.2129 0 0 60894.1465 0 0 40132.9669 0 0 19964.5421 ;435330.0526 0 0 419732.2252 0 0 404334.3977 0 0 379828.3047 0 0 355557.0796 0 0 331355.5903 0 0 307258.7049 0 0 283301.2911 0 0 259518.217 0 0 235944.3505 0 0 212614.5595 0 0 189563.712 0 0 166826.6758 0 0 144438.3189 0 0 122433.5093 0 0 100847.1148 0 0 79714.0034 0 0 59069.0431 0 0 38947.1016 0 0 19383.0471 0 0 ;0 435330.0526 0 0 419732.2252 0 0 404334.3977 0 0 379828.3047 0 0 355557.0796 0 0 331355.5903 0 0 307258.7049 0 0 283301.2911 0 0 259518.217 0 0 235944.3505 0 0 212614.5595 0 0 189563.712 0 0 166826.6758 0 0 144438.3189 0 0 122433.5093 0 0 100847.1148 0 0 79714.0034 0 0 59069.0431 0 0 38947.1016 0 0 19383.0471 0 ;0 0 435330.0526 0 0 419732.2252 0 0 404334.3977 0 0 379828.3047 0 0 355557.0796 0 0 331355.5903 0 0 307258.7049 0 0 283301.2911 0 0 259518.217 0 0 235944.3505 0 0 212614.5595 0 0 189563.712 0 0 166826.6758 0 0 144438.3189 0 0 122433.5093 0 0 100847.1148 0 0 79714.0034 0 0 59069.0431 0 0 38947.1016 0 0 19383.0471 ;408699.4453 0 0 394263.875 0 0 379828.3047 0 0 365592.7343 0 0 342248.8984 0 0 319139.9304 0 0 296100.6983 0 0 273166.0699 0 0 250370.9133 0 0 227750.0963 0 0 205338.4869 0 0 183170.9531 0 0 161282.3626 0 0 139707.5836 0 0 118481.4838 0 0 97638.9313 0 0 77214.7939 0 0 57243.9396 0 0 37761.2364 0 0 18801.5521 0 0 ;0 408699.4453 0 0 394263.875 0 0 379828.3047 0 0 365592.7343 0 0 342248.8984 0 0 319139.9304 0 0 296100.6983 0 0 273166.0699 0 0 250370.9133 0 0 227750.0963 0 0 205338.4869 0 0 183170.9531 0 0 161282.3626 0 0 139707.5836 0 0 118481.4838 0 0 97638.9313 0 0 77214.7939 0 0 57243.9396 0 0 37761.2364 0 0 18801.5521 0 ;0 0 408699.4453 0 0 394263.875 0 0 379828.3047 0 0 365592.7343 0 0 342248.8984 0 0 319139.9304 0 0 296100.6983 0 0 273166.0699 0 0 250370.9133 0 0 227750.0963 0 0 205338.4869 0 0 183170.9531 0 0 161282.3626 0 0 139707.5836 0 0 118481.4838 0 0 97638.9313 0 0 77214.7939 0 0 57243.9396 0 0 37761.2364 0 0 18801.5521 ;382173.4419 0 0 368865.2607 0 0 355557.0796 0 0 342248.8984 0 0 329140.7173 0 0 306924.2705 0 0 284942.6917 0 0 263030.8488 0 0 241223.6096 0 0 219555.8422 0 0 198062.4144 0 0 176778.1942 0 0 155738.0495 0 0 134976.8482 0 0 114529.4583 0 0 94430.7478 0 0 74715.5844 0 0 55418.8362 0 0 36575.3711 0 0 18220.0571 0 0 ;0 382173.4419 0 0 368865.2607 0 0 355557.0796 0 0 342248.8984 0 0 329140.7173 0 0 306924.2705 0 0 284942.6917 0 0 263030.8488 0 0 241223.6096 0 0 219555.8422 0 0 198062.4144 0 0 176778.1942 0 0 155738.0495 0 0 134976.8482 0 0 114529.4583 0 0 94430.7478 0 0 74715.5844 0 0 55418.8362 0 0 36575.3711 0 0 18220.0571 0 ;0 0 382173.4419 0 0 368865.2607 0 0 355557.0796 0 0 342248.8984 0 0 329140.7173 0 0 306924.2705 0 0 284942.6917 0 0 263030.8488 0 0 241223.6096 0 0 219555.8422 0 0 198062.4144 0 0 176778.1942 0 0 155738.0495 0 0 134976.8482 0 0 114529.4583 0 0 94430.7478 0 0 74715.5844 0 0 55418.8362 0 0 36575.3711 0 0 18220.0571 ;355786.9101 0 0 343571.2502 0 0 331355.5903 0 0 319139.9304 0 0 306924.2705 0 0 294908.6106 0 0 273784.6851 0 0 252895.6276 0 0 232076.3059 0 0 211361.588 0 0 190786.3418 0 0 170385.4353 0 0 150193.7363 0 0 130246.1129 0 0 110577.4329 0 0 91222.5642 0 0 72216.3749 0 0 53593.7328 0 0 35389.5059 0 0 17638.5621 0 0 ;0 355786.9101 0 0 343571.2502 0 0 331355.5903 0 0 319139.9304 0 0 306924.2705 0 0 294908.6106 0 0 273784.6851 0 0 252895.6276 0 0 232076.3059 0 0 211361.588 0 0 190786.3418 0 0 170385.4353 0 0 150193.7363 0 0 130246.1129 0 0 110577.4329 0 0 91222.5642 0 0 72216.3749 0 0 53593.7328 0 0 35389.5059 0 0 17638.5621 0 ;0 0 355786.9101 0 0 343571.2502 0 0 331355.5903 0 0 319139.9304 0 0 306924.2705 0 0 294908.6106 0 0 273784.6851 0 0 252895.6276 0 0 232076.3059 0 0 211361.588 0 0 190786.3418 0 0 170385.4353 0 0 150193.7363 0 0 130246.1129 0 0 110577.4329 0 0 91222.5642 0 0 72216.3749 0 0 53593.7328 0 0 35389.5059 0 0 17638.5621 ;329574.718 0 0 318416.7114 0 0 307258.7049 0 0 296100.6983 0 0 284942.6917 0 0 273784.6851 0 0 262826.6786 0 0 242760.4064 0 0 222929.0022 0 0 203167.3338 0 0 183510.2692 0 0 163992.6764 0 0 144649.4231 0 0 125515.3775 0 0 106625.4074 0 0 88014.3807 0 0 69717.1654 0 0 51768.6294 0 0 34203.6407 0 0 17057.067 0 0 ;0 329574.718 0 0 318416.7114 0 0 307258.7049 0 0 296100.6983 0 0 284942.6917 0 0 273784.6851 0 0 262826.6786 0 0 242760.4064 0 0 222929.0022 0 0 203167.3338 0 0 183510.2692 0 0 163992.6764 0 0 144649.4231 0 0 125515.3775 0 0 106625.4074 0 0 88014.3807 0 0 69717.1654 0 0 51768.6294 0 0 34203.6407 0 0 17057.067 0 ;0 0 329574.718 0 0 318416.7114 0 0 307258.7049 0 0 296100.6983 0 0 284942.6917 0 0 273784.6851 0 0 262826.6786 0 0 242760.4064 0 0 222929.0022 0 0 203167.3338 0 0 183510.2692 0 0 163992.6764 0 0 144649.4231 0 0 125515.3775 0 0 106625.4074 0 0 88014.3807 0 0 69717.1654 0 0 51768.6294 0 0 34203.6407 0 0 17057.067 ;303571.7335 0 0 293436.5123 0 0 283301.2911 0 0 273166.0699 0 0 263030.8488 0 0 252895.6276 0 0 242760.4064 0 0 232825.1852 0 0 213781.6985 0 0 194973.0796 0 0 176234.1966 0 0 157599.9174 0 0 139105.11 0 0 120784.6422 0 0 102673.3819 0 0 84806.1972 0 0 67217.9559 0 0 49943.526 0 0 33017.7754 0 0 16475.572 0 0 ;0 303571.7335 0 0 293436.5123 0 0 283301.2911 0 0 273166.0699 0 0 263030.8488 0 0 252895.6276 0 0 242760.4064 0 0 232825.1852 0 0 213781.6985 0 0 194973.0796 0 0 176234.1966 0 0 157599.9174 0 0 139105.11 0 0 120784.6422 0 0 102673.3819 0 0 84806.1972 0 0 67217.9559 0 0 49943.526 0 0 33017.7754 0 0 16475.572 0 ;0 0 303571.7335 0 0 293436.5123 0 0 283301.2911 0 0 273166.0699 0 0 263030.8488 0 0 252895.6276 0 0 242760.4064 0 0 232825.1852 0 0 213781.6985 0 0 194973.0796 0 0 176234.1966 0 0 157599.9174 0 0 139105.11 0 0 120784.6422 0 0 102673.3819 0 0 84806.1972 0 0 67217.9559 0 0 49943.526 0 0 33017.7754 0 0 16475.572 ;277812.8244 0 0 268665.5207 0 0 259518.217 0 0 250370.9133 0 0 241223.6096 0 0 232076.3059 0 0 222929.0022 0 0 213781.6985 0 0 204834.3947 0 0 186778.8254 0 0 168958.1241 0 0 151207.1585 0 0 133560.7968 0 0 116053.9068 0 0 98721.3565 0 0 81598.0137 0 0 64718.7464 0 0 48118.4226 0 0 31831.9102 0 0 15894.077 0 0 ;0 277812.8244 0 0 268665.5207 0 0 259518.217 0 0 250370.9133 0 0 241223.6096 0 0 232076.3059 0 0 222929.0022 0 0 213781.6985 0 0 204834.3947 0 0 186778.8254 0 0 168958.1241 0 0 151207.1585 0 0 133560.7968 0 0 116053.9068 0 0 98721.3565 0 0 81598.0137 0 0 64718.7464 0 0 48118.4226 0 0 31831.9102 0 0 15894.077 0 ;0 0 277812.8244 0 0 268665.5207 0 0 259518.217 0 0 250370.9133 0 0 241223.6096 0 0 232076.3059 0 0 222929.0022 0 0 213781.6985 0 0 204834.3947 0 0 186778.8254 0 0 168958.1241 0 0 151207.1585 0 0 133560.7968 0 0 116053.9068 0 0 98721.3565 0 0 81598.0137 0 0 64718.7464 0 0 48118.4226 0 0 31831.9102 0 0 15894.077 ;252332.8589 0 0 244138.6047 0 0 235944.3505 0 0 227750.0963 0 0 219555.8422 0 0 211361.588 0 0 203167.3338 0 0 194973.0796 0 0 186778.8254 0 0 178784.5713 0 0 161682.0515 0 0 144814.3996 0 0 128016.4836 0 0 111323.1715 0 0 94769.331 0 0 78389.8302 0 0 62219.5369 0 0 46293.3192 0 0 30646.0449 0 0 15312.582 0 0 ;0 252332.8589 0 0 244138.6047 0 0 235944.3505 0 0 227750.0963 0 0 219555.8422 0 0 211361.588 0 0 203167.3338 0 0 194973.0796 0 0 186778.8254 0 0 178784.5713 0 0 161682.0515 0 0 144814.3996 0 0 128016.4836 0 0 111323.1715 0 0 94769.331 0 0 78389.8302 0 0 62219.5369 0 0 46293.3192 0 0 30646.0449 0 0 15312.582 0 ;0 0 252332.8589 0 0 244138.6047 0 0 235944.3505 0 0 227750.0963 0 0 219555.8422 0 0 211361.588 0 0 203167.3338 0 0 194973.0796 0 0 186778.8254 0 0 178784.5713 0 0 161682.0515 0 0 144814.3996 0 0 128016.4836 0 0 111323.1715 0 0 94769.331 0 0 78389.8302 0 0 62219.5369 0 0 46293.3192 0 0 30646.0449 0 0 15312.582 ;227166.7047 0 0 219890.6321 0 0 212614.5595 0 0 205338.4869 0 0 198062.4144 0 0 190786.3418 0 0 183510.2692 0 0 176234.1966 0 0 168958.1241 0 0 161682.0515 0 0 154605.9789 0 0 138421.6407 0 0 122472.1705 0 0 106592.4361 0 0 90817.3055 0 0 75181.6466 0 0 59720.3274 0 0 44468.2158 0 0 29460.1797 0 0 14731.087 0 0 ;0 227166.7047 0 0 219890.6321 0 0 212614.5595 0 0 205338.4869 0 0 198062.4144 0 0 190786.3418 0 0 183510.2692 0 0 176234.1966 0 0 168958.1241 0 0 161682.0515 0 0 154605.9789 0 0 138421.6407 0 0 122472.1705 0 0 106592.4361 0 0 90817.3055 0 0 75181.6466 0 0 59720.3274 0 0 44468.2158 0 0 29460.1797 0 0 14731.087 0 ;0 0 227166.7047 0 0 219890.6321 0 0 212614.5595 0 0 205338.4869 0 0 198062.4144 0 0 190786.3418 0 0 183510.2692 0 0 176234.1966 0 0 168958.1241 0 0 161682.0515 0 0 154605.9789 0 0 138421.6407 0 0 122472.1705 0 0 106592.4361 0 0 90817.3055 0 0 75181.6466 0 0 59720.3274 0 0 44468.2158 0 0 29460.1797 0 0 14731.087 ;202349.2298 0 0 195956.4709 0 0 189563.712 0 0 183170.9531 0 0 176778.1942 0 0 170385.4353 0 0 163992.6764 0 0 157599.9174 0 0 151207.1585 0 0 144814.3996 0 0 138421.6407 0 0 132228.8818 0 0 116927.8573 0 0 101861.7008 0 0 86865.28 0 0 71973.4631 0 0 57221.1179 0 0 42643.1124 0 0 28274.3144 0 0 14149.592 0 0 ;0 202349.2298 0 0 195956.4709 0 0 189563.712 0 0 183170.9531 0 0 176778.1942 0 0 170385.4353 0 0 163992.6764 0 0 157599.9174 0 0 151207.1585 0 0 144814.3996 0 0 138421.6407 0 0 132228.8818 0 0 116927.8573 0 0 101861.7008 0 0 86865.28 0 0 71973.4631 0 0 57221.1179 0 0 42643.1124 0 0 28274.3144 0 0 14149.592 0 ;0 0 202349.2298 0 0 195956.4709 0 0 189563.712 0 0 183170.9531 0 0 176778.1942 0 0 170385.4353 0 0 163992.6764 0 0 157599.9174 0 0 151207.1585 0 0 144814.3996 0 0 138421.6407 0 0 132228.8818 0 0 116927.8573 0 0 101861.7008 0 0 86865.28 0 0 71973.4631 0 0 57221.1179 0 0 42643.1124 0 0 28274.3144 0 0 14149.592 ;177915.3021 0 0 172370.989 0 0 166826.6758 0 0 161282.3626 0 0 155738.0495 0 0 150193.7363 0 0 144649.4231 0 0 139105.11 0 0 133560.7968 0 0 128016.4836 0 0 122472.1705 0 0 116927.8573 0 0 111583.5442 0 0 97130.9654 0 0 82913.2546 0 0 68765.2796 0 0 54721.9084 0 0 40818.009 0 0 27088.4492 0 0 13568.0969 0 0 ;0 177915.3021 0 0 172370.989 0 0 166826.6758 0 0 161282.3626 0 0 155738.0495 0 0 150193.7363 0 0 144649.4231 0 0 139105.11 0 0 133560.7968 0 0 128016.4836 0 0 122472.1705 0 0 116927.8573 0 0 111583.5442 0 0 97130.9654 0 0 82913.2546 0 0 68765.2796 0 0 54721.9084 0 0 40818.009 0 0 27088.4492 0 0 13568.0969 0 ;0 0 177915.3021 0 0 172370.989 0 0 166826.6758 0 0 161282.3626 0 0 155738.0495 0 0 150193.7363 0 0 144649.4231 0 0 139105.11 0 0 133560.7968 0 0 128016.4836 0 0 122472.1705 0 0 116927.8573 0 0 111583.5442 0 0 97130.9654 0 0 82913.2546 0 0 68765.2796 0 0 54721.9084 0 0 40818.009 0 0 27088.4492 0 0 13568.0969 ;153899.7896 0 0 149169.0543 0 0 144438.3189 0 0 139707.5836 0 0 134976.8482 0 0 130246.1129 0 0 125515.3775 0 0 120784.6422 0 0 116053.9068 0 0 111323.1715 0 0 106592.4361 0 0 101861.7008 0 0 97130.9654 0 0 92600.23 0 0 78961.2291 0 0 65557.0961 0 0 52222.6989 0 0 38992.9056 0 0 25902.5839 0 0 12986.6019 0 0 ;0 153899.7896 0 0 149169.0543 0 0 144438.3189 0 0 139707.5836 0 0 134976.8482 0 0 130246.1129 0 0 125515.3775 0 0 120784.6422 0 0 116053.9068 0 0 111323.1715 0 0 106592.4361 0 0 101861.7008 0 0 97130.9654 0 0 92600.23 0 0 78961.2291 0 0 65557.0961 0 0 52222.6989 0 0 38992.9056 0 0 25902.5839 0 0 12986.6019 0 ;0 0 153899.7896 0 0 149169.0543 0 0 144438.3189 0 0 139707.5836 0 0 134976.8482 0 0 130246.1129 0 0 125515.3775 0 0 120784.6422 0 0 116053.9068 0 0 111323.1715 0 0 106592.4361 0 0 101861.7008 0 0 97130.9654 0 0 92600.23 0 0 78961.2291 0 0 65557.0961 0 0 52222.6989 0 0 38992.9056 0 0 25902.5839 0 0 12986.6019 ;130337.5602 0 0 126385.5347 0 0 122433.5093 0 0 118481.4838 0 0 114529.4583 0 0 110577.4329 0 0 106625.4074 0 0 102673.3819 0 0 98721.3565 0 0 94769.331 0 0 90817.3055 0 0 86865.28 0 0 82913.2546 0 0 78961.2291 0 0 75209.2036 0 0 62348.9126 0 0 49723.4894 0 0 37167.8021 0 0 24716.7187 0 0 12405.1069 0 0 ;0 130337.5602 0 0 126385.5347 0 0 122433.5093 0 0 118481.4838 0 0 114529.4583 0 0 110577.4329 0 0 106625.4074 0 0 102673.3819 0 0 98721.3565 0 0 94769.331 0 0 90817.3055 0 0 86865.28 0 0 82913.2546 0 0 78961.2291 0 0 75209.2036 0 0 62348.9126 0 0 49723.4894 0 0 37167.8021 0 0 24716.7187 0 0 12405.1069 0 ;0 0 130337.5602 0 0 126385.5347 0 0 122433.5093 0 0 118481.4838 0 0 114529.4583 0 0 110577.4329 0 0 106625.4074 0 0 102673.3819 0 0 98721.3565 0 0 94769.331 0 0 90817.3055 0 0 86865.28 0 0 82913.2546 0 0 78961.2291 0 0 75209.2036 0 0 62348.9126 0 0 49723.4894 0 0 37167.8021 0 0 24716.7187 0 0 12405.1069 ;107263.4818 0 0 104055.2983 0 0 100847.1148 0 0 97638.9313 0 0 94430.7478 0 0 91222.5642 0 0 88014.3807 0 0 84806.1972 0 0 81598.0137 0 0 78389.8302 0 0 75181.6466 0 0 71973.4631 0 0 68765.2796 0 0 65557.0961 0 0 62348.9126 0 0 59340.729 0 0 47224.2799 0 0 35342.6987 0 0 23530.8534 0 0 11823.6119 0 0 ;0 107263.4818 0 0 104055.2983 0 0 100847.1148 0 0 97638.9313 0 0 94430.7478 0 0 91222.5642 0 0 88014.3807 0 0 84806.1972 0 0 81598.0137 0 0 78389.8302 0 0 75181.6466 0 0 71973.4631 0 0 68765.2796 0 0 65557.0961 0 0 62348.9126 0 0 59340.729 0 0 47224.2799 0 0 35342.6987 0 0 23530.8534 0 0 11823.6119 0 ;0 0 107263.4818 0 0 104055.2983 0 0 100847.1148 0 0 97638.9313 0 0 94430.7478 0 0 91222.5642 0 0 88014.3807 0 0 84806.1972 0 0 81598.0137 0 0 78389.8302 0 0 75181.6466 0 0 71973.4631 0 0 68765.2796 0 0 65557.0961 0 0 62348.9126 0 0 59340.729 0 0 47224.2799 0 0 35342.6987 0 0 23530.8534 0 0 11823.6119 ;84712.4224 0 0 82213.2129 0 0 79714.0034 0 0 77214.7939 0 0 74715.5844 0 0 72216.3749 0 0 69717.1654 0 0 67217.9559 0 0 64718.7464 0 0 62219.5369 0 0 59720.3274 0 0 57221.1179 0 0 54721.9084 0 0 52222.6989 0 0 49723.4894 0 0 47224.2799 0 0 44925.0704 0 0 33517.5953 0 0 22344.9882 0 0 11242.1169 0 0 ;0 84712.4224 0 0 82213.2129 0 0 79714.0034 0 0 77214.7939 0 0 74715.5844 0 0 72216.3749 0 0 69717.1654 0 0 67217.9559 0 0 64718.7464 0 0 62219.5369 0 0 59720.3274 0 0 57221.1179 0 0 54721.9084 0 0 52222.6989 0 0 49723.4894 0 0 47224.2799 0 0 44925.0704 0 0 33517.5953 0 0 22344.9882 0 0 11242.1169 0 ;0 0 84712.4224 0 0 82213.2129 0 0 79714.0034 0 0 77214.7939 0 0 74715.5844 0 0 72216.3749 0 0 69717.1654 0 0 67217.9559 0 0 64718.7464 0 0 62219.5369 0 0 59720.3274 0 0 57221.1179 0 0 54721.9084 0 0 52222.6989 0 0 49723.4894 0 0 47224.2799 0 0 44925.0704 0 0 33517.5953 0 0 22344.9882 0 0 11242.1169 ;62719.2499 0 0 60894.1465 0 0 59069.0431 0 0 57243.9396 0 0 55418.8362 0 0 53593.7328 0 0 51768.6294 0 0 49943.526 0 0 48118.4226 0 0 46293.3192 0 0 44468.2158 0 0 42643.1124 0 0 40818.009 0 0 38992.9056 0 0 37167.8021 0 0 35342.6987 0 0 33517.5953 0 0 31892.4919 0 0 21159.1229 0 0 10660.6218 0 0 ;0 62719.2499 0 0 60894.1465 0 0 59069.0431 0 0 57243.9396 0 0 55418.8362 0 0 53593.7328 0 0 51768.6294 0 0 49943.526 0 0 48118.4226 0 0 46293.3192 0 0 44468.2158 0 0 42643.1124 0 0 40818.009 0 0 38992.9056 0 0 37167.8021 0 0 35342.6987 0 0 33517.5953 0 0 31892.4919 0 0 21159.1229 0 0 10660.6218 0 ;0 0 62719.2499 0 0 60894.1465 0 0 59069.0431 0 0 57243.9396 0 0 55418.8362 0 0 53593.7328 0 0 51768.6294 0 0 49943.526 0 0 48118.4226 0 0 46293.3192 0 0 44468.2158 0 0 42643.1124 0 0 40818.009 0 0 38992.9056 0 0 37167.8021 0 0 35342.6987 0 0 33517.5953 0 0 31892.4919 0 0 21159.1229 0 0 10660.6218 ;41318.8321 0 0 40132.9669 0 0 38947.1016 0 0 37761.2364 0 0 36575.3711 0 0 35389.5059 0 0 34203.6407 0 0 33017.7754 0 0 31831.9102 0 0 30646.0449 0 0 29460.1797 0 0 28274.3144 0 0 27088.4492 0 0 25902.5839 0 0 24716.7187 0 0 23530.8534 0 0 22344.9882 0 0 21159.1229 0 0 20173.2577 0 0 10079.1268 0 0 ;0 41318.8321 0 0 40132.9669 0 0 38947.1016 0 0 37761.2364 0 0 36575.3711 0 0 35389.5059 0 0 34203.6407 0 0 33017.7754 0 0 31831.9102 0 0 30646.0449 0 0 29460.1797 0 0 28274.3144 0 0 27088.4492 0 0 25902.5839 0 0 24716.7187 0 0 23530.8534 0 0 22344.9882 0 0 21159.1229 0 0 20173.2577 0 0 10079.1268 0 ;0 0 41318.8321 0 0 40132.9669 0 0 38947.1016 0 0 37761.2364 0 0 36575.3711 0 0 35389.5059 0 0 34203.6407 0 0 33017.7754 0 0 31831.9102 0 0 30646.0449 0 0 29460.1797 0 0 28274.3144 0 0 27088.4492 0 0 25902.5839 0 0 24716.7187 0 0 23530.8534 0 0 22344.9882 0 0 21159.1229 0 0 20173.2577 0 0 10079.1268 ;20546.0372 0 0 19964.5421 0 0 19383.0471 0 0 18801.5521 0 0 18220.0571 0 0 17638.5621 0 0 17057.067 0 0 16475.572 0 0 15894.077 0 0 15312.582 0 0 14731.087 0 0 14149.592 0 0 13568.0969 0 0 12986.6019 0 0 12405.1069 0 0 11823.6119 0 0 11242.1169 0 0 10660.6218 0 0 10079.1268 0 0 9697.6318 0 0 ;0 20546.0372 0 0 19964.5421 0 0 19383.0471 0 0 18801.5521 0 0 18220.0571 0 0 17638.5621 0 0 17057.067 0 0 16475.572 0 0 15894.077 0 0 15312.582 0 0 14731.087 0 0 14149.592 0 0 13568.0969 0 0 12986.6019 0 0 12405.1069 0 0 11823.6119 0 0 11242.1169 0 0 10660.6218 0 0 10079.1268 0 0 9697.6318 0 ;0 0 20546.0372 0 0 19964.5421 0 0 19383.0471 0 0 18801.5521 0 0 18220.0571 0 0 17638.5621 0 0 17057.067 0 0 16475.572 0 0 15894.077 0 0 15312.582 0 0 14731.087 0 0 14149.592 0 0 13568.0969 0 0 12986.6019 0 0 12405.1069 0 0 11823.6119 0 0 11242.1169 0 0 10660.6218 0 0 10079.1268 0 0 9697.6318]; 
  for (i = 0; i < 120; i++) {
    d = 0.0;
    for (i1 = 0; i1 < 6; i1++) {
      d += (2.0 * x0[i1]) * (static_cast<real_T>(iv[i1 + (6 * i)]));
    }

    lam[i] = d;
  }

  for (i = 0; i < 120; i++) {
    d = 0.0;
    for (i1 = 0; i1 < 120; i1++) {
      d += lam[i1] * b[i1 + (120 * i)];
    }

    ftol[i] = d;
  }

  //  + F_tilde';
  //  L = [-eye(m*Hp);eye(m*Hp)];
  //  % L = [-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 ;1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 ;]; 
  //  b = [kron(ones(Hp,1), V_TRMPC(1:m,1)); kron(ones(Hp,1), V_TRMPC(m+1:end,1))]; 
  //  % b=[0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;]; 
  //  Call QP Solver
  //  Call QP Solver
  //  QQ = (QQ+QQ')/2;
  //  % QQ = [138371.1425 1.3007e-10 1.3212e-10 112666.5366 1.2498e-10 1.274e-10 106967.442 1.1989e-10 1.2267e-10 101279.3704 1.148e-10 1.1794e-10 95607.833 1.0971e-10 1.1322e-10 89958.3413 1.0462e-10 1.0849e-10 84336.4067 9.953e-11 1.0377e-10 78747.5407 9.4441e-11 9.904e-11 73197.2546 8.9351e-11 9.4314e-11 67691.06 8.4261e-11 8.9588e-11 62234.4682 7.9172e-11 8.4862e-11 56832.9906 7.4082e-11 8.0136e-11 51492.1387 6.8992e-11 7.541e-11 46217.4239 6.3903e-11 7.0684e-11 41014.3576 5.8813e-11 6.5958e-11 35888.4513 5.3723e-11 6.1232e-11 30845.2164 4.8633e-11 5.6505e-11 25890.1642 4.3544e-11 5.1779e-11 21028.8063 3.8454e-11 4.7053e-11 16266.654 3.3364e-11 4.2327e-11 ;1.3007e-10 138371.1425 4.3201e-10 1.272e-10 112666.5366 4.1299e-10 1.2433e-10 106967.442 3.9397e-10 1.2146e-10 101279.3704 3.7494e-10 1.186e-10 95607.833 3.5592e-10 1.1573e-10 89958.3413 3.369e-10 1.1286e-10 84336.4067 3.1787e-10 1.0999e-10 78747.5407 2.9885e-10 1.0712e-10 73197.2546 2.7983e-10 1.0425e-10 67691.06 2.608e-10 1.0139e-10 62234.4682 2.4178e-10 9.8518e-11 56832.9906 2.2276e-10 9.5649e-11 51492.1387 2.0373e-10 9.2781e-11 46217.4239 1.8471e-10 8.9913e-11 41014.3576 1.6569e-10 8.7044e-11 35888.4513 1.4666e-10 8.4176e-11 30845.2164 1.2764e-10 8.1308e-11 25890.1642 1.0861e-10 7.844e-11 21028.8063 8.9591e-11 7.5571e-11 16266.654 7.0567e-11 ;1.3212e-10 4.3201e-10 138371.1425 1.2659e-10 4.1439e-10 112666.5366 1.2105e-10 3.9677e-10 106967.442 1.1552e-10 3.7915e-10 101279.3704 1.0999e-10 3.6152e-10 95607.833 1.0445e-10 3.439e-10 89958.3413 9.8919e-11 3.2628e-10 84336.4067 9.3385e-11 3.0865e-10 78747.5407 8.7851e-11 2.9103e-10 73197.2546 8.2317e-11 2.7341e-10 67691.06 7.6783e-11 2.5578e-10 62234.4682 7.1249e-11 2.3816e-10 56832.9906 6.5715e-11 2.2054e-10 51492.1387 6.0181e-11 2.0292e-10 46217.4239 5.4647e-11 1.8529e-10 41014.3576 4.9113e-11 1.6767e-10 35888.4513 4.3579e-11 1.5005e-10 30845.2164 3.8045e-11 1.3242e-10 25890.1642 3.2511e-11 1.148e-10 21028.8063 2.6977e-11 9.7177e-11 16266.654 ;112666.5366 1.272e-10 1.2659e-10 129154.3028 1.2224e-10 1.2208e-10 103675.0808 1.1727e-10 1.1756e-10 98201.3702 1.1231e-10 1.1305e-10 92738.6825 1.0734e-10 1.0854e-10 87292.529 1.0238e-10 1.0402e-10 81868.4212 9.7413e-11 9.9509e-11 76471.8706 9.2448e-11 9.4995e-11 71108.3885 8.7483e-11 9.0482e-11 65783.4864 8.2519e-11 8.5969e-11 60502.6757 7.7554e-11 8.1455e-11 55271.4678 7.259e-11 7.6942e-11 50095.3742 6.7625e-11 7.2429e-11 44979.9063 6.2661e-11 6.7915e-11 39930.5754 5.7696e-11 6.3402e-11 34952.8931 5.2731e-11 5.8889e-11 30052.3707 4.7767e-11 5.4375e-11 25234.5197 4.2802e-11 4.9862e-11 20504.8515 3.7838e-11 4.5348e-11 15868.8775 3.2873e-11 4.0835e-11 ;1.2498e-10 112666.5366 4.1439e-10 1.2224e-10 129154.3028 3.9615e-10 1.1949e-10 103675.0808 3.779e-10 1.1675e-10 98201.3702 3.5966e-10 1.1401e-10 92738.6825 3.4141e-10 1.1126e-10 87292.529 3.2317e-10 1.0852e-10 81868.4212 3.0492e-10 1.0578e-10 76471.8706 2.8668e-10 1.0303e-10 71108.3885 2.6843e-10 1.0029e-10 65783.4864 2.5019e-10 9.7547e-11 60502.6757 2.3194e-10 9.4804e-11 55271.4678 2.137e-10 9.2061e-11 50095.3742 1.9545e-10 8.9318e-11 44979.9063 1.7721e-10 8.6574e-11 39930.5754 1.5896e-10 8.3831e-11 34952.8931 1.4072e-10 8.1088e-11 30052.3707 1.2247e-10 7.8345e-11 25234.5197 1.0423e-10 7.5602e-11 20504.8515 8.5983e-11 7.2858e-11 15868.8775 6.7738e-11 ;1.274e-10 4.1299e-10 112666.5366 1.2208e-10 3.9615e-10 129154.3028 1.1675e-10 3.793e-10 103675.0808 1.1143e-10 3.6246e-10 98201.3702 1.0611e-10 3.4561e-10 92738.6825 1.0079e-10 3.2877e-10 87292.529 9.5469e-11 3.1192e-10 81868.4212 9.0148e-11 2.9508e-10 76471.8706 8.4827e-11 2.7824e-10 71108.3885 7.9506e-11 2.6139e-10 65783.4864 7.4184e-11 2.4455e-10 60502.6757 6.8863e-11 2.277e-10 55271.4678 6.3542e-11 2.1086e-10 50095.3742 5.8221e-11 1.9401e-10 44979.9063 5.2899e-11 1.7717e-10 39930.5754 4.7578e-11 1.6032e-10 34952.8931 4.2257e-11 1.4348e-10 30052.3707 3.6936e-11 1.2664e-10 25234.5197 3.1614e-11 1.0979e-10 20504.8515 2.6293e-11 9.2947e-11 15868.8775 ;106967.442 1.2433e-10 1.2105e-10 103675.0808 1.1949e-10 1.1675e-10 120382.7195 1.1465e-10 1.1245e-10 95123.37 1.0981e-10 1.0815e-10 89869.532 1.0497e-10 1.0385e-10 84626.7167 1.0013e-10 9.9552e-11 79400.4358 9.5295e-11 9.5251e-11 74196.2006 9.0455e-11 9.0951e-11 69019.5224 8.5616e-11 8.665e-11 63875.9129 8.0776e-11 8.2349e-11 58770.8833 7.5937e-11 7.8049e-11 53709.9451 7.1097e-11 7.3748e-11 48698.6097 6.6258e-11 6.9447e-11 43742.3886 6.1418e-11 6.5147e-11 38846.7932 5.6579e-11 6.0846e-11 34017.3348 5.1739e-11 5.6546e-11 29259.525 4.69e-11 5.2245e-11 24578.8751 4.206e-11 4.7944e-11 19980.8967 3.7221e-11 4.3644e-11 15471.101 3.2381e-11 3.9343e-11 ;1.1989e-10 106967.442 3.9677e-10 1.1727e-10 103675.0808 3.793e-10 1.1465e-10 120382.7195 3.6184e-10 1.1203e-10 95123.37 3.4437e-10 1.0942e-10 89869.532 3.269e-10 1.068e-10 84626.7167 3.0944e-10 1.0418e-10 79400.4358 2.9197e-10 1.0156e-10 74196.2006 2.745e-10 9.8944e-11 69019.5224 2.5704e-10 9.6326e-11 63875.9129 2.3957e-10 9.3708e-11 58770.8833 2.2211e-10 9.109e-11 53709.9451 2.0464e-10 8.8472e-11 48698.6097 1.8717e-10 8.5854e-11 43742.3886 1.6971e-10 8.3236e-11 38846.7932 1.5224e-10 8.0618e-11 34017.3348 1.3477e-10 7.8e-11 29259.525 1.1731e-10 7.5382e-11 24578.8751 9.9841e-11 7.2764e-11 19980.8967 8.2375e-11 7.0146e-11 15471.101 6.4909e-11 ;1.2267e-10 3.9397e-10 106967.442 1.1756e-10 3.779e-10 103675.0808 1.1245e-10 3.6184e-10 120382.7195 1.0735e-10 3.4577e-10 95123.37 1.0224e-10 3.297e-10 89869.532 9.7128e-11 3.1364e-10 84626.7167 9.202e-11 2.9757e-10 79400.4358 8.6911e-11 2.8151e-10 74196.2006 8.1803e-11 2.6544e-10 69019.5224 7.6694e-11 2.4938e-10 63875.9129 7.1586e-11 2.3331e-10 58770.8833 6.6477e-11 2.1724e-10 53709.9451 6.1368e-11 2.0118e-10 48698.6097 5.626e-11 1.8511e-10 43742.3886 5.1151e-11 1.6905e-10 38846.7932 4.6043e-11 1.5298e-10 34017.3348 4.0934e-11 1.3691e-10 29259.525 3.5826e-11 1.2085e-10 24578.8751 3.0717e-11 1.0478e-10 19980.8967 2.5609e-11 8.8717e-11 15471.101 ;101279.3704 1.2146e-10 1.1552e-10 98201.3702 1.1675e-10 1.1143e-10 95123.37 1.1203e-10 1.0735e-10 112045.3699 1.0732e-10 1.0326e-10 87000.3814 1.0261e-10 9.9169e-11 81960.9045 9.7892e-11 9.5081e-11 76932.4503 9.3177e-11 9.0994e-11 71920.5305 8.8463e-11 8.6906e-11 66930.6563 8.3748e-11 8.2818e-11 61968.3393 7.9034e-11 7.873e-11 57039.0908 7.432e-11 7.4642e-11 52148.4223 6.9605e-11 7.0554e-11 47301.8452 6.4891e-11 6.6466e-11 42504.871 6.0176e-11 6.2378e-11 37763.0109 5.5462e-11 5.829e-11 33081.7766 5.0748e-11 5.4203e-11 28466.6793 4.6033e-11 5.0115e-11 23923.2306 4.1319e-11 4.6027e-11 19456.9418 3.6604e-11 4.1939e-11 15073.3244 3.189e-11 3.7851e-11 ;1.148e-10 101279.3704 3.7915e-10 1.1231e-10 98201.3702 3.6246e-10 1.0981e-10 95123.37 3.4577e-10 1.0732e-10 112045.3699 3.2908e-10 1.0483e-10 87000.3814 3.1239e-10 1.0233e-10 81960.9045 2.9571e-10 9.9842e-11 76932.4503 2.7902e-10 9.7349e-11 71920.5305 2.6233e-10 9.4856e-11 66930.6563 2.4564e-10 9.2363e-11 61968.3393 2.2896e-10 8.987e-11 57039.0908 2.1227e-10 8.7377e-11 52148.4223 1.9558e-10 8.4884e-11 47301.8452 1.7889e-10 8.2391e-11 42504.871 1.6221e-10 7.9898e-11 37763.0109 1.4552e-10 7.7405e-11 33081.7766 1.2883e-10 7.4912e-11 28466.6793 1.1214e-10 7.2419e-11 23923.2306 9.5455e-11 6.9926e-11 19456.9418 7.8767e-11 6.7433e-11 15073.3244 6.2079e-11 ;1.1794e-10 3.7494e-10 101279.3704 1.1305e-10 3.5966e-10 98201.3702 1.0815e-10 3.4437e-10 95123.37 1.0326e-10 3.2908e-10 112045.3699 9.8361e-11 3.138e-10 87000.3814 9.3466e-11 2.9851e-10 81960.9045 8.857e-11 2.8322e-10 76932.4503 8.3674e-11 2.6793e-10 71920.5305 7.8778e-11 2.5265e-10 66930.6563 7.3883e-11 2.3736e-10 61968.3393 6.8987e-11 2.2207e-10 57039.0908 6.4091e-11 2.0678e-10 52148.4223 5.9195e-11 1.915e-10 47301.8452 5.4299e-11 1.7621e-10 42504.871 4.9404e-11 1.6092e-10 37763.0109 4.4508e-11 1.4564e-10 33081.7766 3.9612e-11 1.3035e-10 28466.6793 3.4716e-11 1.1506e-10 23923.2306 2.982e-11 9.9774e-11 19456.9418 2.4925e-11 8.4487e-11 15073.3244 ;95607.833 1.186e-10 1.0999e-10 92738.6825 1.1401e-10 1.0611e-10 89869.532 1.0942e-10 1.0224e-10 87000.3814 1.0483e-10 9.8361e-11 104131.2309 1.0024e-10 9.4486e-11 79295.0922 9.5649e-11 9.0611e-11 74464.4649 9.106e-11 8.6736e-11 69644.8604 8.647e-11 8.2861e-11 64841.7902 8.1881e-11 7.8986e-11 60060.7658 7.7292e-11 7.5111e-11 55307.2984 7.2702e-11 7.1235e-11 50586.8996 6.8113e-11 6.736e-11 45905.0807 6.3524e-11 6.3485e-11 41267.3533 5.8934e-11 5.961e-11 36679.2287 5.4345e-11 5.5735e-11 32146.2183 4.9756e-11 5.186e-11 27673.8337 4.5166e-11 4.7984e-11 23267.5861 4.0577e-11 4.4109e-11 18932.987 3.5988e-11 4.0234e-11 14675.5479 3.1398e-11 3.6359e-11 ;1.0971e-10 95607.833 3.6152e-10 1.0734e-10 92738.6825 3.4561e-10 1.0497e-10 89869.532 3.297e-10 1.0261e-10 87000.3814 3.138e-10 1.0024e-10 104131.2309 2.9789e-10 9.787e-11 79295.0922 2.8198e-10 9.5502e-11 74464.4649 2.6607e-10 9.3135e-11 69644.8604 2.5016e-10 9.0767e-11 64841.7902 2.3425e-10 8.8399e-11 60060.7658 2.1834e-10 8.6031e-11 55307.2984 2.0243e-10 8.3663e-11 50586.8996 1.8652e-10 8.1295e-11 45905.0807 1.7061e-10 7.8927e-11 41267.3533 1.547e-10 7.6559e-11 36679.2287 1.388e-10 7.4191e-11 32146.2183 1.2289e-10 7.1823e-11 27673.8337 1.0698e-10 6.9456e-11 23267.5861 9.1068e-11 6.7088e-11 18932.987 7.5159e-11 6.472e-11 14675.5479 5.925e-11 ;1.1322e-10 3.5592e-10 95607.833 1.0854e-10 3.4141e-10 92738.6825 1.0385e-10 3.269e-10 89869.532 9.9169e-11 3.1239e-10 87000.3814 9.4486e-11 2.9789e-10 104131.2309 8.9803e-11 2.8338e-10 79295.0922 8.512e-11 2.6887e-10 74464.4649 8.0437e-11 2.5436e-10 69644.8604 7.5754e-11 2.3985e-10 64841.7902 7.1071e-11 2.2534e-10 60060.7658 6.6388e-11 2.1083e-10 55307.2984 6.1705e-11 1.9633e-10 50586.8996 5.7022e-11 1.8182e-10 45905.0807 5.2339e-11 1.6731e-10 41267.3533 4.7656e-11 1.528e-10 36679.2287 4.2973e-11 1.3829e-10 32146.2183 3.829e-11 1.2378e-10 27673.8337 3.3607e-11 1.0927e-10 23267.5861 2.8924e-11 9.4766e-11 18932.987 2.4241e-11 8.0257e-11 14675.5479 ;89958.3413 1.1573e-10 1.0445e-10 87292.529 1.1126e-10 1.0079e-10 84626.7167 1.068e-10 9.7128e-11 81960.9045 1.0233e-10 9.3466e-11 79295.0922 9.787e-11 8.9803e-11 96629.2799 9.3406e-11 8.6141e-11 71996.4794 8.8942e-11 8.2478e-11 67369.1904 8.4478e-11 7.8816e-11 62752.9241 8.0013e-11 7.5154e-11 58153.1922 7.5549e-11 7.1491e-11 53575.5059 7.1085e-11 6.7829e-11 49025.3768 6.6621e-11 6.4166e-11 44508.3163 6.2157e-11 6.0504e-11 40029.8357 5.7692e-11 5.6841e-11 35595.4465 5.3228e-11 5.3179e-11 31210.6601 4.8764e-11 4.9517e-11 26880.988 4.43e-11 4.5854e-11 22611.9415 3.9835e-11 4.2192e-11 18409.0322 3.5371e-11 3.8529e-11 14277.7714 3.0907e-11 3.4867e-11 ;1.0462e-10 89958.3413 3.439e-10 1.0238e-10 87292.529 3.2877e-10 1.0013e-10 84626.7167 3.1364e-10 9.7892e-11 81960.9045 2.9851e-10 9.5649e-11 79295.0922 2.8338e-10 9.3406e-11 96629.2799 2.6825e-10 9.1163e-11 71996.4794 2.5312e-10 8.892e-11 67369.1904 2.3799e-10 8.6678e-11 62752.9241 2.2286e-10 8.4435e-11 58153.1922 2.0773e-10 8.2192e-11 53575.5059 1.9259e-10 7.9949e-11 49025.3768 1.7746e-10 7.7706e-11 44508.3163 1.6233e-10 7.5464e-11 40029.8357 1.472e-10 7.3221e-11 35595.4465 1.3207e-10 7.0978e-11 31210.6601 1.1694e-10 6.8735e-11 26880.988 1.0181e-10 6.6492e-11 22611.9415 8.6682e-11 6.425e-11 18409.0322 7.1551e-11 6.2007e-11 14277.7714 5.6421e-11 ;1.0849e-10 3.369e-10 89958.3413 1.0402e-10 3.2317e-10 87292.529 9.9552e-11 3.0944e-10 84626.7167 9.5081e-11 2.9571e-10 81960.9045 9.0611e-11 2.8198e-10 79295.0922 8.6141e-11 2.6825e-10 96629.2799 8.1671e-11 2.5452e-10 71996.4794 7.72e-11 2.4079e-10 67369.1904 7.273e-11 2.2706e-10 62752.9241 6.826e-11 2.1333e-10 58153.1922 6.3789e-11 1.996e-10 53575.5059 5.9319e-11 1.8587e-10 49025.3768 5.4849e-11 1.7214e-10 44508.3163 5.0378e-11 1.5841e-10 40029.8357 4.5908e-11 1.4468e-10 35595.4465 4.1438e-11 1.3095e-10 31210.6601 3.6967e-11 1.1722e-10 26880.988 3.2497e-11 1.0349e-10 22611.9415 2.8027e-11 8.9758e-11 18409.0322 2.3556e-11 7.6028e-11 14277.7714 ;84336.4067 1.1286e-10 9.8919e-11 81868.4212 1.0852e-10 9.5469e-11 79400.4358 1.0418e-10 9.202e-11 76932.4503 9.9842e-11 8.857e-11 74464.4649 9.5502e-11 8.512e-11 71996.4794 9.1163e-11 8.1671e-11 89528.494 8.6824e-11 7.8221e-11 65093.5203 8.2485e-11 7.4771e-11 60664.058 7.8146e-11 7.1321e-11 56245.6186 7.3807e-11 6.7872e-11 51843.7135 6.9468e-11 6.4422e-11 47463.8541 6.5128e-11 6.0972e-11 43111.5518 6.0789e-11 5.7523e-11 38792.318 5.645e-11 5.4073e-11 34511.6642 5.2111e-11 5.0623e-11 30275.1018 4.7772e-11 4.7174e-11 26088.1423 4.3433e-11 4.3724e-11 21956.297 3.9094e-11 4.0274e-11 17885.0774 3.4755e-11 3.6825e-11 13879.9948 3.0415e-11 3.3375e-11 ;9.953e-11 84336.4067 3.2628e-10 9.7413e-11 81868.4212 3.1192e-10 9.5295e-11 79400.4358 2.9757e-10 9.3177e-11 76932.4503 2.8322e-10 9.106e-11 74464.4649 2.6887e-10 8.8942e-11 71996.4794 2.5452e-10 8.6824e-11 89528.494 2.4017e-10 8.4706e-11 65093.5203 2.2581e-10 8.2589e-11 60664.058 2.1146e-10 8.0471e-11 56245.6186 1.9711e-10 7.8353e-11 51843.7135 1.8276e-10 7.6236e-11 47463.8541 1.6841e-10 7.4118e-11 43111.5518 1.5405e-10 7.2e-11 38792.318 1.397e-10 6.9882e-11 34511.6642 1.2535e-10 6.7765e-11 30275.1018 1.11e-10 6.5647e-11 26088.1423 9.6647e-11 6.3529e-11 21956.297 8.2295e-11 6.1412e-11 17885.0774 6.7943e-11 5.9294e-11 13879.9948 5.3591e-11 ;1.0377e-10 3.1787e-10 84336.4067 9.9509e-11 3.0492e-10 81868.4212 9.5251e-11 2.9197e-10 79400.4358 9.0994e-11 2.7902e-10 76932.4503 8.6736e-11 2.6607e-10 74464.4649 8.2478e-11 2.5312e-10 71996.4794 7.8221e-11 2.4017e-10 89528.494 7.3963e-11 2.2721e-10 65093.5203 6.9706e-11 2.1426e-10 60664.058 6.5448e-11 2.0131e-10 56245.6186 6.119e-11 1.8836e-10 51843.7135 5.6933e-11 1.7541e-10 47463.8541 5.2675e-11 1.6246e-10 43111.5518 4.8418e-11 1.4951e-10 38792.318 4.416e-11 1.3655e-10 34511.6642 3.9903e-11 1.236e-10 30275.1018 3.5645e-11 1.1065e-10 26088.1423 3.1387e-11 9.77e-11 21956.297 2.713e-11 8.4749e-11 17885.0774 2.2872e-11 7.1798e-11 13879.9948 ;78747.5407 1.0999e-10 9.3385e-11 76471.8706 1.0578e-10 9.0148e-11 74196.2006 1.0156e-10 8.6911e-11 71920.5305 9.7349e-11 8.3674e-11 69644.8604 9.3135e-11 8.0437e-11 67369.1904 8.892e-11 7.72e-11 65093.5203 8.4706e-11 7.3963e-11 82817.8502 8.0492e-11 7.0726e-11 58575.192 7.6278e-11 6.7489e-11 54338.0451 7.2064e-11 6.4252e-11 50111.9211 6.785e-11 6.1015e-11 45902.3313 6.3636e-11 5.7778e-11 41714.7873 5.9422e-11 5.4542e-11 37554.8004 5.5208e-11 5.1305e-11 33427.882 5.0994e-11 4.8068e-11 29339.5436 4.678e-11 4.4831e-11 25295.2966 4.2566e-11 4.1594e-11 21300.6524 3.8352e-11 3.8357e-11 17361.1225 3.4138e-11 3.512e-11 13482.2183 2.9924e-11 3.1883e-11 ;9.4441e-11 78747.5407 3.0865e-10 9.2448e-11 76471.8706 2.9508e-10 9.0455e-11 74196.2006 2.8151e-10 8.8463e-11 71920.5305 2.6793e-10 8.647e-11 69644.8604 2.5436e-10 8.4478e-11 67369.1904 2.4079e-10 8.2485e-11 65093.5203 2.2721e-10 8.0492e-11 82817.8502 2.1364e-10 7.85e-11 58575.192 2.0007e-10 7.6507e-11 54338.0451 1.8649e-10 7.4515e-11 50111.9211 1.7292e-10 7.2522e-11 45902.3313 1.5935e-10 7.0529e-11 41714.7873 1.4577e-10 6.8537e-11 37554.8004 1.322e-10 6.6544e-11 33427.882 1.1863e-10 6.4551e-11 29339.5436 1.0505e-10 6.2559e-11 25295.2966 9.1482e-11 6.0566e-11 21300.6524 7.7908e-11 5.8574e-11 17361.1225 6.4335e-11 5.6581e-11 13482.2183 5.0762e-11 ;9.904e-11 2.9885e-10 78747.5407 9.4995e-11 2.8668e-10 76471.8706 9.0951e-11 2.745e-10 74196.2006 8.6906e-11 2.6233e-10 71920.5305 8.2861e-11 2.5016e-10 69644.8604 7.8816e-11 2.3799e-10 67369.1904 7.4771e-11 2.2581e-10 65093.5203 7.0726e-11 2.1364e-10 82817.8502 6.6681e-11 2.0147e-10 58575.192 6.2637e-11 1.893e-10 54338.0451 5.8592e-11 1.7712e-10 50111.9211 5.4547e-11 1.6495e-10 45902.3313 5.0502e-11 1.5278e-10 41714.7873 4.6457e-11 1.406e-10 37554.8004 4.2412e-11 1.2843e-10 33427.882 3.8367e-11 1.1626e-10 29339.5436 3.4323e-11 1.0409e-10 25295.2966 3.0278e-11 9.1913e-11 21300.6524 2.6233e-11 7.9741e-11 17361.1225 2.2188e-11 6.7568e-11 13482.2183 ;73197.2546 1.0712e-10 8.7851e-11 71108.3885 1.0303e-10 8.4827e-11 69019.5224 9.8944e-11 8.1803e-11 66930.6563 9.4856e-11 7.8778e-11 64841.7902 9.0767e-11 7.5754e-11 62752.9241 8.6678e-11 7.273e-11 60664.058 8.2589e-11 6.9706e-11 58575.192 7.85e-11 6.6681e-11 76486.3259 7.4411e-11 6.3657e-11 52430.4715 7.0322e-11 6.0633e-11 48380.1286 6.6233e-11 5.7609e-11 44340.8086 6.2144e-11 5.4585e-11 40318.0228 5.8055e-11 5.156e-11 36317.2827 5.3966e-11 4.8536e-11 32344.0998 4.9877e-11 4.5512e-11 28403.9854 4.5788e-11 4.2488e-11 24502.4509 4.1699e-11 3.9463e-11 20645.0079 3.761e-11 3.6439e-11 16837.1677 3.3521e-11 3.3415e-11 13084.4418 2.9432e-11 3.0391e-11 ;8.9351e-11 73197.2546 2.9103e-10 8.7483e-11 71108.3885 2.7824e-10 8.5616e-11 69019.5224 2.6544e-10 8.3748e-11 66930.6563 2.5265e-10 8.1881e-11 64841.7902 2.3985e-10 8.0013e-11 62752.9241 2.2706e-10 7.8146e-11 60664.058 2.1426e-10 7.6278e-11 58575.192 2.0147e-10 7.4411e-11 76486.3259 1.8867e-10 7.2543e-11 52430.4715 1.7588e-10 7.0676e-11 48380.1286 1.6308e-10 6.8808e-11 44340.8086 1.5029e-10 6.6941e-11 40318.0228 1.3749e-10 6.5073e-11 36317.2827 1.247e-10 6.3206e-11 32344.0998 1.1191e-10 6.1338e-11 28403.9854 9.9111e-11 5.9471e-11 24502.4509 8.6316e-11 5.7603e-11 20645.0079 7.3522e-11 5.5736e-11 16837.1677 6.0727e-11 5.3868e-11 13084.4418 4.7933e-11 ;9.4314e-11 2.7983e-10 73197.2546 9.0482e-11 2.6843e-10 71108.3885 8.665e-11 2.5704e-10 69019.5224 8.2818e-11 2.4564e-10 66930.6563 7.8986e-11 2.3425e-10 64841.7902 7.5154e-11 2.2286e-10 62752.9241 7.1321e-11 2.1146e-10 60664.058 6.7489e-11 2.0007e-10 58575.192 6.3657e-11 1.8867e-10 76486.3259 5.9825e-11 1.7728e-10 52430.4715 5.5993e-11 1.6589e-10 48380.1286 5.2161e-11 1.5449e-10 44340.8086 4.8329e-11 1.431e-10 40318.0228 4.4497e-11 1.317e-10 36317.2827 4.0664e-11 1.2031e-10 32344.0998 3.6832e-11 1.0891e-10 28403.9854 3.3e-11 9.752e-11 24502.4509 2.9168e-11 8.6126e-11 20645.0079 2.5336e-11 7.4732e-11 16837.1677 2.1504e-11 6.3338e-11 13084.4418 ;67691.06 1.0425e-10 8.2317e-11 65783.4864 1.0029e-10 7.9506e-11 63875.9129 9.6326e-11 7.6694e-11 61968.3393 9.2363e-11 7.3883e-11 60060.7658 8.8399e-11 7.1071e-11 58153.1922 8.4435e-11 6.826e-11 56245.6186 8.0471e-11 6.5448e-11 54338.0451 7.6507e-11 6.2637e-11 52430.4715 7.2543e-11 5.9825e-11 70522.898 6.8579e-11 5.7014e-11 46648.3362 6.4616e-11 5.4202e-11 42779.2858 6.0652e-11 5.1391e-11 38921.2583 5.6688e-11 4.8579e-11 35079.7651 5.2724e-11 4.5768e-11 31260.3175 4.876e-11 4.2956e-11 27468.4271 4.4796e-11 4.0145e-11 23709.6053 4.0833e-11 3.7333e-11 19989.3634 3.6869e-11 3.4522e-11 16313.2129 3.2905e-11 3.171e-11 12686.6652 2.8941e-11 2.8899e-11 ;8.4261e-11 67691.06 2.7341e-10 8.2519e-11 65783.4864 2.6139e-10 8.0776e-11 63875.9129 2.4938e-10 7.9034e-11 61968.3393 2.3736e-10 7.7292e-11 60060.7658 2.2534e-10 7.5549e-11 58153.1922 2.1333e-10 7.3807e-11 56245.6186 2.0131e-10 7.2064e-11 54338.0451 1.893e-10 7.0322e-11 52430.4715 1.7728e-10 6.8579e-11 70522.898 1.6526e-10 6.6837e-11 46648.3362 1.5325e-10 6.5095e-11 42779.2858 1.4123e-10 6.3352e-11 38921.2583 1.2922e-10 6.161e-11 35079.7651 1.172e-10 5.9867e-11 31260.3175 1.0518e-10 5.8125e-11 27468.4271 9.3167e-11 5.6382e-11 23709.6053 8.1151e-11 5.464e-11 19989.3634 6.9135e-11 5.2898e-11 16313.2129 5.7119e-11 5.1155e-11 12686.6652 4.5103e-11 ;8.9588e-11 2.608e-10 67691.06 8.5969e-11 2.5019e-10 65783.4864 8.2349e-11 2.3957e-10 63875.9129 7.873e-11 2.2896e-10 61968.3393 7.5111e-11 2.1834e-10 60060.7658 7.1491e-11 2.0773e-10 58153.1922 6.7872e-11 1.9711e-10 56245.6186 6.4252e-11 1.8649e-10 54338.0451 6.0633e-11 1.7588e-10 52430.4715 5.7014e-11 1.6526e-10 70522.898 5.3394e-11 1.5465e-10 46648.3362 4.9775e-11 1.4403e-10 42779.2858 4.6155e-11 1.3342e-10 38921.2583 4.2536e-11 1.228e-10 35079.7651 3.8917e-11 1.1219e-10 31260.3175 3.5297e-11 1.0157e-10 27468.4271 3.1678e-11 9.0955e-11 23709.6053 2.8058e-11 8.0339e-11 19989.3634 2.4439e-11 6.9724e-11 16313.2129 2.082e-11 5.9108e-11 12686.6652 ;62234.4682 1.0139e-10 7.6783e-11 60502.6757 9.7547e-11 7.4184e-11 58770.8833 9.3708e-11 7.1586e-11 57039.0908 8.987e-11 6.8987e-11 55307.2984 8.6031e-11 6.6388e-11 53575.5059 8.2192e-11 6.3789e-11 51843.7135 7.8353e-11 6.119e-11 50111.9211 7.4515e-11 5.8592e-11 48380.1286 7.0676e-11 5.5993e-11 46648.3362 6.6837e-11 5.3394e-11 64916.5437 6.2998e-11 5.0795e-11 41217.7631 5.916e-11 4.8197e-11 37524.4938 5.5321e-11 4.5598e-11 33842.2474 5.1482e-11 4.2999e-11 30176.5353 4.7643e-11 4.04e-11 26532.8689 4.3805e-11 3.7802e-11 22916.7596 3.9966e-11 3.5203e-11 19333.7188 3.6127e-11 3.2604e-11 15789.2581 3.2288e-11 3.0005e-11 12288.8887 2.8449e-11 2.7407e-11 ;7.9172e-11 62234.4682 2.5578e-10 7.7554e-11 60502.6757 2.4455e-10 7.5937e-11 58770.8833 2.3331e-10 7.432e-11 57039.0908 2.2207e-10 7.2702e-11 55307.2984 2.1083e-10 7.1085e-11 53575.5059 1.996e-10 6.9468e-11 51843.7135 1.8836e-10 6.785e-11 50111.9211 1.7712e-10 6.6233e-11 48380.1286 1.6589e-10 6.4616e-11 46648.3362 1.5465e-10 6.2998e-11 64916.5437 1.4341e-10 6.1381e-11 41217.7631 1.3217e-10 5.9764e-11 37524.4938 1.2094e-10 5.8146e-11 33842.2474 1.097e-10 5.6529e-11 30176.5353 9.8461e-11 5.4912e-11 26532.8689 8.7223e-11 5.3294e-11 22916.7596 7.5986e-11 5.1677e-11 19333.7188 6.4749e-11 5.006e-11 15789.2581 5.3511e-11 4.8442e-11 12288.8887 4.2274e-11 ;8.4862e-11 2.4178e-10 62234.4682 8.1455e-11 2.3194e-10 60502.6757 7.8049e-11 2.2211e-10 58770.8833 7.4642e-11 2.1227e-10 57039.0908 7.1235e-11 2.0243e-10 55307.2984 6.7829e-11 1.9259e-10 53575.5059 6.4422e-11 1.8276e-10 51843.7135 6.1015e-11 1.7292e-10 50111.9211 5.7609e-11 1.6308e-10 48380.1286 5.4202e-11 1.5325e-10 46648.3362 5.0795e-11 1.4341e-10 64916.5437 4.7389e-11 1.3357e-10 41217.7631 4.3982e-11 1.2374e-10 37524.4938 4.0575e-11 1.139e-10 33842.2474 3.7169e-11 1.0406e-10 30176.5353 3.3762e-11 9.4226e-11 26532.8689 3.0356e-11 8.4389e-11 22916.7596 2.6949e-11 7.4552e-11 19333.7188 2.3542e-11 6.4715e-11 15789.2581 2.0136e-11 5.4878e-11 12288.8887 ;56832.9906 9.8518e-11 7.1249e-11 55271.4678 9.4804e-11 6.8863e-11 53709.9451 9.109e-11 6.6477e-11 52148.4223 8.7377e-11 6.4091e-11 50586.8996 8.3663e-11 6.1705e-11 49025.3768 7.9949e-11 5.9319e-11 47463.8541 7.6236e-11 5.6933e-11 45902.3313 7.2522e-11 5.4547e-11 44340.8086 6.8808e-11 5.2161e-11 42779.2858 6.5095e-11 4.9775e-11 41217.7631 6.1381e-11 4.7389e-11 59656.2403 5.7667e-11 4.5003e-11 36127.7293 5.3954e-11 4.2617e-11 32604.7298 5.024e-11 4.0231e-11 29092.7531 4.6526e-11 3.7845e-11 25597.3106 4.2813e-11 3.5459e-11 22123.9139 3.9099e-11 3.3073e-11 18678.0743 3.5385e-11 3.0687e-11 15265.3032 3.1672e-11 2.8301e-11 11891.1121 2.7958e-11 2.5915e-11 ;7.4082e-11 56832.9906 2.3816e-10 7.259e-11 55271.4678 2.277e-10 7.1097e-11 53709.9451 2.1724e-10 6.9605e-11 52148.4223 2.0678e-10 6.8113e-11 50586.8996 1.9633e-10 6.6621e-11 49025.3768 1.8587e-10 6.5128e-11 47463.8541 1.7541e-10 6.3636e-11 45902.3313 1.6495e-10 6.2144e-11 44340.8086 1.5449e-10 6.0652e-11 42779.2858 1.4403e-10 5.916e-11 41217.7631 1.3357e-10 5.7667e-11 59656.2403 1.2311e-10 5.6175e-11 36127.7293 1.1266e-10 5.4683e-11 32604.7298 1.022e-10 5.3191e-11 29092.7531 9.1738e-11 5.1698e-11 25597.3106 8.128e-11 5.0206e-11 22123.9139 7.0821e-11 4.8714e-11 18678.0743 6.0362e-11 4.7222e-11 15265.3032 4.9903e-11 4.5729e-11 11891.1121 3.9444e-11 ;8.0136e-11 2.2276e-10 56832.9906 7.6942e-11 2.137e-10 55271.4678 7.3748e-11 2.0464e-10 53709.9451 7.0554e-11 1.9558e-10 52148.4223 6.736e-11 1.8652e-10 50586.8996 6.4166e-11 1.7746e-10 49025.3768 6.0972e-11 1.6841e-10 47463.8541 5.7778e-11 1.5935e-10 45902.3313 5.4585e-11 1.5029e-10 44340.8086 5.1391e-11 1.4123e-10 42779.2858 4.8197e-11 1.3217e-10 41217.7631 4.5003e-11 1.2311e-10 59656.2403 4.1809e-11 1.1406e-10 36127.7293 3.8615e-11 1.05e-10 32604.7298 3.5421e-11 9.594e-11 29092.7531 3.2227e-11 8.6882e-11 25597.3106 2.9033e-11 7.7823e-11 22123.9139 2.5839e-11 6.8765e-11 18678.0743 2.2645e-11 5.9707e-11 15265.3032 1.9451e-11 5.0648e-11 11891.1121 ;51492.1387 9.5649e-11 6.5715e-11 50095.3742 9.2061e-11 6.3542e-11 48698.6097 8.8472e-11 6.1368e-11 47301.8452 8.4884e-11 5.9195e-11 45905.0807 8.1295e-11 5.7022e-11 44508.3163 7.7706e-11 5.4849e-11 43111.5518 7.4118e-11 5.2675e-11 41714.7873 7.0529e-11 5.0502e-11 40318.0228 6.6941e-11 4.8329e-11 38921.2583 6.3352e-11 4.6155e-11 37524.4938 5.9764e-11 4.3982e-11 36127.7293 5.6175e-11 4.1809e-11 54730.9648 5.2586e-11 3.9636e-11 31367.2121 4.8998e-11 3.7462e-11 28008.9708 4.5409e-11 3.5289e-11 24661.7524 4.1821e-11 3.3116e-11 21331.0682 3.8232e-11 3.0942e-11 18022.4297 3.4644e-11 2.8769e-11 14741.3484 3.1055e-11 2.6596e-11 11493.3356 2.7467e-11 2.4423e-11 ;6.8992e-11 51492.1387 2.2054e-10 6.7625e-11 50095.3742 2.1086e-10 6.6258e-11 48698.6097 2.0118e-10 6.4891e-11 47301.8452 1.915e-10 6.3524e-11 45905.0807 1.8182e-10 6.2157e-11 44508.3163 1.7214e-10 6.0789e-11 43111.5518 1.6246e-10 5.9422e-11 41714.7873 1.5278e-10 5.8055e-11 40318.0228 1.431e-10 5.6688e-11 38921.2583 1.3342e-10 5.5321e-11 37524.4938 1.2374e-10 5.3954e-11 36127.7293 1.1406e-10 5.2586e-11 54730.9648 1.0438e-10 5.1219e-11 31367.2121 9.4696e-11 4.9852e-11 28008.9708 8.5016e-11 4.8485e-11 24661.7524 7.5336e-11 4.7118e-11 21331.0682 6.5656e-11 4.5751e-11 18022.4297 5.5975e-11 4.4384e-11 14741.3484 4.6295e-11 4.3016e-11 11493.3356 3.6615e-11 ;7.541e-11 2.0373e-10 51492.1387 7.2429e-11 1.9545e-10 50095.3742 6.9447e-11 1.8717e-10 48698.6097 6.6466e-11 1.7889e-10 47301.8452 6.3485e-11 1.7061e-10 45905.0807 6.0504e-11 1.6233e-10 44508.3163 5.7523e-11 1.5405e-10 43111.5518 5.4542e-11 1.4577e-10 41714.7873 5.156e-11 1.3749e-10 40318.0228 4.8579e-11 1.2922e-10 38921.2583 4.5598e-11 1.2094e-10 37524.4938 4.2617e-11 1.1266e-10 36127.7293 3.9636e-11 1.0438e-10 54730.9648 3.6654e-11 9.6097e-11 31367.2121 3.3673e-11 8.7817e-11 28008.9708 3.0692e-11 7.9537e-11 24661.7524 2.7711e-11 7.1257e-11 21331.0682 2.473e-11 6.2978e-11 18022.4297 2.1748e-11 5.4698e-11 14741.3484 1.8767e-11 4.6418e-11 11493.3356 ;46217.4239 9.2781e-11 6.0181e-11 44979.9063 8.9318e-11 5.8221e-11 43742.3886 8.5854e-11 5.626e-11 42504.871 8.2391e-11 5.4299e-11 41267.3533 7.8927e-11 5.2339e-11 40029.8357 7.5464e-11 5.0378e-11 38792.318 7.2e-11 4.8418e-11 37554.8004 6.8537e-11 4.6457e-11 36317.2827 6.5073e-11 4.4497e-11 35079.7651 6.161e-11 4.2536e-11 33842.2474 5.8146e-11 4.0575e-11 32604.7298 5.4683e-11 3.8615e-11 31367.2121 5.1219e-11 3.6654e-11 50129.6945 4.7756e-11 3.4694e-11 26925.1886 4.4292e-11 3.2733e-11 23726.1941 4.0829e-11 3.0773e-11 20538.2225 3.7365e-11 2.8812e-11 17366.7852 3.3902e-11 2.6852e-11 14217.3936 3.0438e-11 2.4891e-11 11095.5591 2.6975e-11 2.293e-11 ;6.3903e-11 46217.4239 2.0292e-10 6.2661e-11 44979.9063 1.9401e-10 6.1418e-11 43742.3886 1.8511e-10 6.0176e-11 42504.871 1.7621e-10 5.8934e-11 41267.3533 1.6731e-10 5.7692e-11 40029.8357 1.5841e-10 5.645e-11 38792.318 1.4951e-10 5.5208e-11 37554.8004 1.406e-10 5.3966e-11 36317.2827 1.317e-10 5.2724e-11 35079.7651 1.228e-10 5.1482e-11 33842.2474 1.139e-10 5.024e-11 32604.7298 1.05e-10 4.8998e-11 31367.2121 9.6097e-11 4.7756e-11 50129.6945 8.7195e-11 4.6514e-11 26925.1886 7.8293e-11 4.5272e-11 23726.1941 6.9392e-11 4.403e-11 20538.2225 6.049e-11 4.2788e-11 17366.7852 5.1589e-11 4.1546e-11 14217.3936 4.2687e-11 4.0304e-11 11095.5591 3.3786e-11 ;7.0684e-11 1.8471e-10 46217.4239 6.7915e-11 1.7721e-10 44979.9063 6.5147e-11 1.6971e-10 43742.3886 6.2378e-11 1.6221e-10 42504.871 5.961e-11 1.547e-10 41267.3533 5.6841e-11 1.472e-10 40029.8357 5.4073e-11 1.397e-10 38792.318 5.1305e-11 1.322e-10 37554.8004 4.8536e-11 1.247e-10 36317.2827 4.5768e-11 1.172e-10 35079.7651 4.2999e-11 1.097e-10 33842.2474 4.0231e-11 1.022e-10 32604.7298 3.7462e-11 9.4696e-11 31367.2121 3.4694e-11 8.7195e-11 50129.6945 3.1925e-11 7.9694e-11 26925.1886 2.9157e-11 7.2193e-11 23726.1941 2.6388e-11 6.4692e-11 20538.2225 2.362e-11 5.7191e-11 17366.7852 2.0852e-11 4.969e-11 14217.3936 1.8083e-11 4.2189e-11 11095.5591 ;41014.3576 8.9913e-11 5.4647e-11 39930.5754 8.6574e-11 5.2899e-11 38846.7932 8.3236e-11 5.1151e-11 37763.0109 7.9898e-11 4.9404e-11 36679.2287 7.6559e-11 4.7656e-11 35595.4465 7.3221e-11 4.5908e-11 34511.6642 6.9882e-11 4.416e-11 33427.882 6.6544e-11 4.2412e-11 32344.0998 6.3206e-11 4.0664e-11 31260.3175 5.9867e-11 3.8917e-11 30176.5353 5.6529e-11 3.7169e-11 29092.7531 5.3191e-11 3.5421e-11 28008.9708 4.9852e-11 3.3673e-11 26925.1886 4.6514e-11 3.1925e-11 45841.4064 4.3175e-11 3.0178e-11 22790.6359 3.9837e-11 2.843e-11 19745.3769 3.6499e-11 2.6682e-11 16711.1407 3.316e-11 2.4934e-11 13693.4388 2.9822e-11 2.3186e-11 10697.7825 2.6484e-11 2.1438e-11 ;5.8813e-11 41014.3576 1.8529e-10 5.7696e-11 39930.5754 1.7717e-10 5.6579e-11 38846.7932 1.6905e-10 5.5462e-11 37763.0109 1.6092e-10 5.4345e-11 36679.2287 1.528e-10 5.3228e-11 35595.4465 1.4468e-10 5.2111e-11 34511.6642 1.3655e-10 5.0994e-11 33427.882 1.2843e-10 4.9877e-11 32344.0998 1.2031e-10 4.876e-11 31260.3175 1.1219e-10 4.7643e-11 30176.5353 1.0406e-10 4.6526e-11 29092.7531 9.594e-11 4.5409e-11 28008.9708 8.7817e-11 4.4292e-11 26925.1886 7.9694e-11 4.3175e-11 45841.4064 7.1571e-11 4.2058e-11 22790.6359 6.3448e-11 4.0941e-11 19745.3769 5.5325e-11 3.9825e-11 16711.1407 4.7202e-11 3.8708e-11 13693.4388 3.9079e-11 3.7591e-11 10697.7825 3.0956e-11 ;6.5958e-11 1.6569e-10 41014.3576 6.3402e-11 1.5896e-10 39930.5754 6.0846e-11 1.5224e-10 38846.7932 5.829e-11 1.4552e-10 37763.0109 5.5735e-11 1.388e-10 36679.2287 5.3179e-11 1.3207e-10 35595.4465 5.0623e-11 1.2535e-10 34511.6642 4.8068e-11 1.1863e-10 33427.882 4.5512e-11 1.1191e-10 32344.0998 4.2956e-11 1.0518e-10 31260.3175 4.04e-11 9.8461e-11 30176.5353 3.7845e-11 9.1738e-11 29092.7531 3.5289e-11 8.5016e-11 28008.9708 3.2733e-11 7.8293e-11 26925.1886 3.0178e-11 7.1571e-11 45841.4064 2.7622e-11 6.4849e-11 22790.6359 2.5066e-11 5.8126e-11 19745.3769 2.251e-11 5.1404e-11 16711.1407 1.9955e-11 4.4681e-11 13693.4388 1.7399e-11 3.7959e-11 10697.7825 ;35888.4513 8.7044e-11 4.9113e-11 34952.8931 8.3831e-11 4.7578e-11 34017.3348 8.0618e-11 4.6043e-11 33081.7766 7.7405e-11 4.4508e-11 32146.2183 7.4191e-11 4.2973e-11 31210.6601 7.0978e-11 4.1438e-11 30275.1018 6.7765e-11 3.9903e-11 29339.5436 6.4551e-11 3.8367e-11 28403.9854 6.1338e-11 3.6832e-11 27468.4271 5.8125e-11 3.5297e-11 26532.8689 5.4912e-11 3.3762e-11 25597.3106 5.1698e-11 3.2227e-11 24661.7524 4.8485e-11 3.0692e-11 23726.1941 4.5272e-11 2.9157e-11 22790.6359 4.2058e-11 2.7622e-11 41855.0777 3.8845e-11 2.6087e-11 18952.5312 3.5632e-11 2.4552e-11 16055.4961 3.2419e-11 2.3017e-11 13169.4839 2.9205e-11 2.1481e-11 10300.006 2.5992e-11 1.9946e-11 ;5.3723e-11 35888.4513 1.6767e-10 5.2731e-11 34952.8931 1.6032e-10 5.1739e-11 34017.3348 1.5298e-10 5.0748e-11 33081.7766 1.4564e-10 4.9756e-11 32146.2183 1.3829e-10 4.8764e-11 31210.6601 1.3095e-10 4.7772e-11 30275.1018 1.236e-10 4.678e-11 29339.5436 1.1626e-10 4.5788e-11 28403.9854 1.0891e-10 4.4796e-11 27468.4271 1.0157e-10 4.3805e-11 26532.8689 9.4226e-11 4.2813e-11 25597.3106 8.6882e-11 4.1821e-11 24661.7524 7.9537e-11 4.0829e-11 23726.1941 7.2193e-11 3.9837e-11 22790.6359 6.4849e-11 3.8845e-11 41855.0777 5.7504e-11 3.7853e-11 18952.5312 5.016e-11 3.6861e-11 16055.4961 4.2816e-11 3.587e-11 13169.4839 3.5471e-11 3.4878e-11 10300.006 2.8127e-11 ;6.1232e-11 1.4666e-10 35888.4513 5.8889e-11 1.4072e-10 34952.8931 5.6546e-11 1.3477e-10 34017.3348 5.4203e-11 1.2883e-10 33081.7766 5.186e-11 1.2289e-10 32146.2183 4.9517e-11 1.1694e-10 31210.6601 4.7174e-11 1.11e-10 30275.1018 4.4831e-11 1.0505e-10 29339.5436 4.2488e-11 9.9111e-11 28403.9854 4.0145e-11 9.3167e-11 27468.4271 3.7802e-11 8.7223e-11 26532.8689 3.5459e-11 8.128e-11 25597.3106 3.3116e-11 7.5336e-11 24661.7524 3.0773e-11 6.9392e-11 23726.1941 2.843e-11 6.3448e-11 22790.6359 2.6087e-11 5.7504e-11 41855.0777 2.3744e-11 5.156e-11 18952.5312 2.1401e-11 4.5617e-11 16055.4961 1.9058e-11 3.9673e-11 13169.4839 1.6715e-11 3.3729e-11 10300.006 ;30845.2164 8.4176e-11 4.3579e-11 30052.3707 8.1088e-11 4.2257e-11 29259.525 7.8e-11 4.0934e-11 28466.6793 7.4912e-11 3.9612e-11 27673.8337 7.1823e-11 3.829e-11 26880.988 6.8735e-11 3.6967e-11 26088.1423 6.5647e-11 3.5645e-11 25295.2966 6.2559e-11 3.4323e-11 24502.4509 5.9471e-11 3.3e-11 23709.6053 5.6382e-11 3.1678e-11 22916.7596 5.3294e-11 3.0356e-11 22123.9139 5.0206e-11 2.9033e-11 21331.0682 4.7118e-11 2.7711e-11 20538.2225 4.403e-11 2.6388e-11 19745.3769 4.0941e-11 2.5066e-11 18952.5312 3.7853e-11 2.3744e-11 38159.6855 3.4765e-11 2.2421e-11 15399.8516 3.1677e-11 2.1099e-11 12645.5291 2.8589e-11 1.9777e-11 9902.2295 2.5501e-11 1.8454e-11 ;4.8633e-11 30845.2164 1.5005e-10 4.7767e-11 30052.3707 1.4348e-10 4.69e-11 29259.525 1.3691e-10 4.6033e-11 28466.6793 1.3035e-10 4.5166e-11 27673.8337 1.2378e-10 4.43e-11 26880.988 1.1722e-10 4.3433e-11 26088.1423 1.1065e-10 4.2566e-11 25295.2966 1.0409e-10 4.1699e-11 24502.4509 9.752e-11 4.0833e-11 23709.6053 9.0955e-11 3.9966e-11 22916.7596 8.4389e-11 3.9099e-11 22123.9139 7.7823e-11 3.8232e-11 21331.0682 7.1257e-11 3.7365e-11 20538.2225 6.4692e-11 3.6499e-11 19745.3769 5.8126e-11 3.5632e-11 18952.5312 5.156e-11 3.4765e-11 38159.6855 4.4995e-11 3.3898e-11 15399.8516 3.8429e-11 3.3032e-11 12645.5291 3.1863e-11 3.2165e-11 9902.2295 2.5298e-11 ;5.6505e-11 1.2764e-10 30845.2164 5.4375e-11 1.2247e-10 30052.3707 5.2245e-11 1.1731e-10 29259.525 5.0115e-11 1.1214e-10 28466.6793 4.7984e-11 1.0698e-10 27673.8337 4.5854e-11 1.0181e-10 26880.988 4.3724e-11 9.6647e-11 26088.1423 4.1594e-11 9.1482e-11 25295.2966 3.9463e-11 8.6316e-11 24502.4509 3.7333e-11 8.1151e-11 23709.6053 3.5203e-11 7.5986e-11 22916.7596 3.3073e-11 7.0821e-11 22123.9139 3.0942e-11 6.5656e-11 21331.0682 2.8812e-11 6.049e-11 20538.2225 2.6682e-11 5.5325e-11 19745.3769 2.4552e-11 5.016e-11 18952.5312 2.2421e-11 4.4995e-11 38159.6855 2.0291e-11 3.9829e-11 15399.8516 1.8161e-11 3.4664e-11 12645.5291 1.6031e-11 2.9499e-11 9902.2295 ;25890.1642 8.1308e-11 3.8045e-11 25234.5197 7.8345e-11 3.6936e-11 24578.8751 7.5382e-11 3.5826e-11 23923.2306 7.2419e-11 3.4716e-11 23267.5861 6.9456e-11 3.3607e-11 22611.9415 6.6492e-11 3.2497e-11 21956.297 6.3529e-11 3.1387e-11 21300.6524 6.0566e-11 3.0278e-11 20645.0079 5.7603e-11 2.9168e-11 19989.3634 5.464e-11 2.8058e-11 19333.7188 5.1677e-11 2.6949e-11 18678.0743 4.8714e-11 2.5839e-11 18022.4297 4.5751e-11 2.473e-11 17366.7852 4.2788e-11 2.362e-11 16711.1407 3.9825e-11 2.251e-11 16055.4961 3.6861e-11 2.1401e-11 15399.8516 3.3898e-11 2.0291e-11 34744.207 3.0935e-11 1.9181e-11 12121.5743 2.7972e-11 1.8072e-11 9504.4529 2.5009e-11 1.6962e-11 ;4.3544e-11 25890.1642 1.3242e-10 4.2802e-11 25234.5197 1.2664e-10 4.206e-11 24578.8751 1.2085e-10 4.1319e-11 23923.2306 1.1506e-10 4.0577e-11 23267.5861 1.0927e-10 3.9835e-11 22611.9415 1.0349e-10 3.9094e-11 21956.297 9.77e-11 3.8352e-11 21300.6524 9.1913e-11 3.761e-11 20645.0079 8.6126e-11 3.6869e-11 19989.3634 8.0339e-11 3.6127e-11 19333.7188 7.4552e-11 3.5385e-11 18678.0743 6.8765e-11 3.4644e-11 18022.4297 6.2978e-11 3.3902e-11 17366.7852 5.7191e-11 3.316e-11 16711.1407 5.1404e-11 3.2419e-11 16055.4961 4.5617e-11 3.1677e-11 15399.8516 3.9829e-11 3.0935e-11 34744.207 3.4042e-11 3.0194e-11 12121.5743 2.8255e-11 2.9452e-11 9504.4529 2.2468e-11 ;5.1779e-11 1.0861e-10 25890.1642 4.9862e-11 1.0423e-10 25234.5197 4.7944e-11 9.9841e-11 24578.8751 4.6027e-11 9.5455e-11 23923.2306 4.4109e-11 9.1068e-11 23267.5861 4.2192e-11 8.6682e-11 22611.9415 4.0274e-11 8.2295e-11 21956.297 3.8357e-11 7.7908e-11 21300.6524 3.6439e-11 7.3522e-11 20645.0079 3.4522e-11 6.9135e-11 19989.3634 3.2604e-11 6.4749e-11 19333.7188 3.0687e-11 6.0362e-11 18678.0743 2.8769e-11 5.5975e-11 18022.4297 2.6852e-11 5.1589e-11 17366.7852 2.4934e-11 4.7202e-11 16711.1407 2.3017e-11 4.2816e-11 16055.4961 2.1099e-11 3.8429e-11 15399.8516 1.9181e-11 3.4042e-11 34744.207 1.7264e-11 2.9656e-11 12121.5743 1.5346e-11 2.5269e-11 9504.4529 ;21028.8063 7.844e-11 3.2511e-11 20504.8515 7.5602e-11 3.1614e-11 19980.8967 7.2764e-11 3.0717e-11 19456.9418 6.9926e-11 2.982e-11 18932.987 6.7088e-11 2.8924e-11 18409.0322 6.425e-11 2.8027e-11 17885.0774 6.1412e-11 2.713e-11 17361.1225 5.8574e-11 2.6233e-11 16837.1677 5.5736e-11 2.5336e-11 16313.2129 5.2898e-11 2.4439e-11 15789.2581 5.006e-11 2.3542e-11 15265.3032 4.7222e-11 2.2645e-11 14741.3484 4.4384e-11 2.1748e-11 14217.3936 4.1546e-11 2.0852e-11 13693.4388 3.8708e-11 1.9955e-11 13169.4839 3.587e-11 1.9058e-11 12645.5291 3.3032e-11 1.8161e-11 12121.5743 3.0194e-11 1.7264e-11 31597.6194 2.7356e-11 1.6367e-11 9106.6764 2.4518e-11 1.547e-11 ;3.8454e-11 21028.8063 1.148e-10 3.7838e-11 20504.8515 1.0979e-10 3.7221e-11 19980.8967 1.0478e-10 3.6604e-11 19456.9418 9.9774e-11 3.5988e-11 18932.987 9.4766e-11 3.5371e-11 18409.0322 8.9758e-11 3.4755e-11 17885.0774 8.4749e-11 3.4138e-11 17361.1225 7.9741e-11 3.3521e-11 16837.1677 7.4732e-11 3.2905e-11 16313.2129 6.9724e-11 3.2288e-11 15789.2581 6.4715e-11 3.1672e-11 15265.3032 5.9707e-11 3.1055e-11 14741.3484 5.4698e-11 3.0438e-11 14217.3936 4.969e-11 2.9822e-11 13693.4388 4.4681e-11 2.9205e-11 13169.4839 3.9673e-11 2.8589e-11 12645.5291 3.4664e-11 2.7972e-11 12121.5743 2.9656e-11 2.7356e-11 31597.6194 2.4647e-11 2.6739e-11 9106.6764 1.9639e-11 ;4.7053e-11 8.9591e-11 21028.8063 4.5348e-11 8.5983e-11 20504.8515 4.3644e-11 8.2375e-11 19980.8967 4.1939e-11 7.8767e-11 19456.9418 4.0234e-11 7.5159e-11 18932.987 3.8529e-11 7.1551e-11 18409.0322 3.6825e-11 6.7943e-11 17885.0774 3.512e-11 6.4335e-11 17361.1225 3.3415e-11 6.0727e-11 16837.1677 3.171e-11 5.7119e-11 16313.2129 3.0005e-11 5.3511e-11 15789.2581 2.8301e-11 4.9903e-11 15265.3032 2.6596e-11 4.6295e-11 14741.3484 2.4891e-11 4.2687e-11 14217.3936 2.3186e-11 3.9079e-11 13693.4388 2.1481e-11 3.5471e-11 13169.4839 1.9777e-11 3.1863e-11 12645.5291 1.8072e-11 2.8255e-11 12121.5743 1.6367e-11 2.4647e-11 31597.6194 1.4662e-11 2.1039e-11 9106.6764 ;16266.654 7.5571e-11 2.6977e-11 15868.8775 7.2858e-11 2.6293e-11 15471.101 7.0146e-11 2.5609e-11 15073.3244 6.7433e-11 2.4925e-11 14675.5479 6.472e-11 2.4241e-11 14277.7714 6.2007e-11 2.3556e-11 13879.9948 5.9294e-11 2.2872e-11 13482.2183 5.6581e-11 2.2188e-11 13084.4418 5.3868e-11 2.1504e-11 12686.6652 5.1155e-11 2.082e-11 12288.8887 4.8442e-11 2.0136e-11 11891.1121 4.5729e-11 1.9451e-11 11493.3356 4.3016e-11 1.8767e-11 11095.5591 4.0304e-11 1.8083e-11 10697.7825 3.7591e-11 1.7399e-11 10300.006 3.4878e-11 1.6715e-11 9902.2295 3.2165e-11 1.6031e-11 9504.4529 2.9452e-11 1.5346e-11 9106.6764 2.6739e-11 1.4662e-11 28708.8999 2.4026e-11 1.3978e-11 ;3.3364e-11 16266.654 9.7177e-11 3.2873e-11 15868.8775 9.2947e-11 3.2381e-11 15471.101 8.8717e-11 3.189e-11 15073.3244 8.4487e-11 3.1398e-11 14675.5479 8.0257e-11 3.0907e-11 14277.7714 7.6028e-11 3.0415e-11 13879.9948 7.1798e-11 2.9924e-11 13482.2183 6.7568e-11 2.9432e-11 13084.4418 6.3338e-11 2.8941e-11 12686.6652 5.9108e-11 2.8449e-11 12288.8887 5.4878e-11 2.7958e-11 11891.1121 5.0648e-11 2.7467e-11 11493.3356 4.6418e-11 2.6975e-11 11095.5591 4.2189e-11 2.6484e-11 10697.7825 3.7959e-11 2.5992e-11 10300.006 3.3729e-11 2.5501e-11 9902.2295 2.9499e-11 2.5009e-11 9504.4529 2.5269e-11 2.4518e-11 9106.6764 2.1039e-11 2.4026e-11 28708.8999 1.6809e-11 ;4.2327e-11 7.0567e-11 16266.654 4.0835e-11 6.7738e-11 15868.8775 3.9343e-11 6.4909e-11 15471.101 3.7851e-11 6.2079e-11 15073.3244 3.6359e-11 5.925e-11 14675.5479 3.4867e-11 5.6421e-11 14277.7714 3.3375e-11 5.3591e-11 13879.9948 3.1883e-11 5.0762e-11 13482.2183 3.0391e-11 4.7933e-11 13084.4418 2.8899e-11 4.5103e-11 12686.6652 2.7407e-11 4.2274e-11 12288.8887 2.5915e-11 3.9444e-11 11891.1121 2.4423e-11 3.6615e-11 11493.3356 2.293e-11 3.3786e-11 11095.5591 2.1438e-11 3.0956e-11 10697.7825 1.9946e-11 2.8127e-11 10300.006 1.8454e-11 2.5298e-11 9902.2295 1.6962e-11 2.2468e-11 9504.4529 1.547e-11 1.9639e-11 9106.6764 1.3978e-11 1.6809e-11 28708.8999 ;]; 
  for (b_i = 0; b_i < 60; b_i++) {
    d = 0.0;
    for (i = 0; i < 120; i++) {
      d += ftol[i] * b_b[i + (120 * b_i)];
    }

    H[b_i] = d;
    X_QP[b_i] = dv[b_i];
  }

  //  Solve quadratic programming problem using Wright's (1997) Method
  //  Minimise J(x) = 1/2x'Hx + f'x
  //  Subject to: Ax <= b
  //  Reference: S. J. Wright, "Applying New Optimization Algorithms to Model
  //  Predictive Control," in Chemical Process Control-V, CACHE, AIChE
  //  Symposium, 1997, pp. 147-155.
  // Number of decision variables
  // Test for Cold Start
  // Warm Start
  // to tune
  // to tune
  // Default Values
  for (b_i = 0; b_i < 120; b_i++) {
    lam[b_i] = 100.0;
    ftol[b_i] = 100.0;
    esig[b_i] = 0.001;
  }

  mu = 10000.0;

  //  %Linsolve options
  //  opU.UT = true;
  //  opUT.UT = true;
  //  opUT.TRANSA = true;
  // Begin Searching
  //  for iter = 1:maxiter
  unusedU0 = 0;
  exitflag = 0;
  while ((unusedU0 <= 100) && (exitflag != 1)) {
    int32_T ii_size_idx_0;
    int32_T info;
    int32_T j;
    boolean_T exitg1;
    boolean_T y;

    // Create common matrices
    for (b_i = 0; b_i < 120; b_i++) {
      d = lam[b_i];
      ssq = 1.0 / d;
      ilam[b_i] = ssq;
      nlamt[b_i] = (-d) / ftol[b_i];
      mesil[b_i] = (mu * esig[b_i]) * ssq;
    }

    // RHS
    for (b_i = 0; b_i < 60; b_i++) {
      for (i = 0; i < 120; i++) {
        idx = i + (120 * b_i);
        IGA[idx] = nlamt[i] * (static_cast<real_T>(A[idx]));
      }

      d = 0.0;
      for (i = 0; i < 60; i++) {
        d += b_a[b_i + (60 * i)] * X_QP[i];
      }

      ssq = 0.0;
      for (i = 0; i < 120; i++) {
        ssq += (static_cast<real_T>(c_a[b_i + (60 * i)])) * lam[i];
      }

      r1[b_i] = (d - ssq) - H[b_i];
    }

    for (i = 0; i < 120; i++) {
      d = 0.0;
      for (i1 = 0; i1 < 60; i1++) {
        d += (static_cast<real_T>(a[i + (120 * i1)])) * X_QP[i1];
      }

      igr2[i] = nlamt[i] * ((d + 0.17982) - mesil[i]);
    }

    // Solve
    for (i = 0; i < 60; i++) {
      for (i1 = 0; i1 < 60; i1++) {
        d = 0.0;
        for (b_i = 0; b_i < 120; b_i++) {
          d += (static_cast<real_T>(c_a[i + (60 * b_i)])) * IGA[b_i + (120 * i1)];
        }

        y_tmp[i + (60 * i1)] = d;
      }
    }

    for (i = 0; i < 3600; i++) {
      A_data[i] = b_H[i] - y_tmp[i];
    }

    info = -1;
    j = 0;
    exitg1 = false;
    while ((!exitg1) && (j < 60)) {
      int32_T idxA1j;
      int32_T idxAjj;
      int32_T ix;
      idxA1j = j * 60;
      idxAjj = idxA1j + j;
      ssq = 0.0;
      if (j >= 1) {
        ix = idxA1j;
        iy = idxA1j;
        for (b_i = 0; b_i < j; b_i++) {
          ssq += A_data[ix] * A_data[iy];
          ix++;
          iy++;
        }
      }

      ssq = A_data[idxAjj] - ssq;
      if (ssq > 0.0) {
        ssq = std::sqrt(ssq);
        A_data[idxAjj] = ssq;
        if ((j + 1) < 60) {
          int32_T idxAjjp1;
          b_i = idxA1j + 61;
          idxAjjp1 = idxAjj + 61;
          if (j != 0) {
            iy = idxAjj + 60;
            i = (idxA1j + (60 * (58 - j))) + 61;
            for (idx = b_i; idx <= i; idx += 60) {
              ix = idxA1j;
              mu_old = 0.0;
              i1 = (idx + j) - 1;
              for (ii_size_idx_0 = idx; ii_size_idx_0 <= i1; ii_size_idx_0++) {
                mu_old += A_data[ii_size_idx_0 - 1] * A_data[ix];
                ix++;
              }

              A_data[iy] += -mu_old;
              iy += 60;
            }
          }

          ssq = 1.0 / ssq;
          i = (idxAjj + (60 * (58 - j))) + 61;
          for (b_i = idxAjjp1; b_i <= i; b_i += 60) {
            A_data[b_i - 1] *= ssq;
          }
        }

        j++;
      } else {
        info = j;
        exitg1 = true;
      }
    }

    // [R] = chol(H-At*IGA);
    if ((info + 1) == 0) {
      for (i = 0; i < 60; i++) {
        d = 0.0;
        for (i1 = 0; i1 < 120; i1++) {
          d += (static_cast<real_T>(c_a[i + (60 * i1)])) * igr2[i1];
        }

        del_z[i] = r1[i] - d;
      }

      for (i = 0; i < 3600; i++) {
        y_tmp[i] = b_H[i] - y_tmp[i];
      }

      mldivide(y_tmp, del_z);

      // old method (LU?)
      //   del_z = linsolve (R, linsolve (R, (r1-At*igr2), opUT), opU); %exploit matrix properties for solving 
    } else {
      // Not Positive Definite (problem? eg infeasible)
      for (i = 0; i < 60; i++) {
        d = 0.0;
        for (i1 = 0; i1 < 120; i1++) {
          d += (static_cast<real_T>(c_a[i + (60 * i1)])) * igr2[i1];
        }

        del_z[i] = r1[i] - d;
      }

      for (i = 0; i < 3600; i++) {
        y_tmp[i] = b_H[i] - y_tmp[i];
      }

      mldivide(y_tmp, del_z);

      // old method (LU?)
    }

    // Decide on suitable alpha (from Wright's paper)
    // Try Max Increment (alpha = 1)
    // Check lam and ftol > 0
    for (b_i = 0; b_i < 120; b_i++) {
      d = 0.0;
      for (i = 0; i < 60; i++) {
        d += IGA[b_i + (120 * i)] * del_z[i];
      }

      d = igr2[b_i] - d;
      del_lam[b_i] = d;
      ssq = ftol[b_i];
      mu_old = ((-ssq) + mesil[b_i]) - ((ilam[b_i] * ssq) * d);
      mesil[b_i] = mu_old;
      d += lam[b_i];
      nlamt[b_i] = d;
      ssq += mu_old;
      ilam[b_i] = ssq;
      x[b_i] = (d < 2.2204460492503131E-16);
      x[b_i + 120] = (ssq < 2.2204460492503131E-16);
    }

    y = false;
    b_i = 0;
    exitg1 = false;
    while ((!exitg1) && (b_i < 240)) {
      if (!x[b_i]) {
        b_i++;
      } else {
        y = true;
        exitg1 = true;
      }
    }

    if (!y) {
      // KKT met
      std::memcpy(&lam[0], &nlamt[0], 120U * (sizeof(real_T)));
      std::memcpy(&ftol[0], &ilam[0], 120U * (sizeof(real_T)));
      for (i = 0; i < 60; i++) {
        X_QP[i] += del_z[i];
      }
    } else {
      // KKT failed - solve by finding minimum ratio
      for (i = 0; i < 120; i++) {
        result[i] = nlamt[i];
        result[i + 120] = ilam[i];
      }

      idx = 0;
      b_i = 0;
      exitg1 = false;
      while ((!exitg1) && (b_i < 240)) {
        if (result[b_i] < 2.2204460492503131E-16) {
          idx++;
          ii_data[idx - 1] = static_cast<uint8_T>(static_cast<int32_T>(b_i + 1));
          if (idx >= 240) {
            exitg1 = true;
          } else {
            b_i++;
          }
        } else {
          b_i++;
        }
      }

      if (1 > idx) {
        ii_size_idx_0 = 0;
      } else {
        ii_size_idx_0 = idx;
      }

      // detects elements breaking KKT condition
      for (i = 0; i < 120; i++) {
        b_del_lam[i] = del_lam[i];
        b_del_lam[i + 120] = mesil[i];
      }

      for (i = 0; i < ii_size_idx_0; i++) {
        b_i = (static_cast<int32_T>(ii_data[i])) - 1;
        varargin_1_data[i] = 1.0 - (result[b_i] / b_del_lam[b_i]);
      }

      if (ii_size_idx_0 <= 2) {
        if (ii_size_idx_0 == 1) {
          ssq = varargin_1_data[0];
        } else if (varargin_1_data[0] > varargin_1_data[1]) {
          ssq = varargin_1_data[1];
        } else if (rtIsNaN(varargin_1_data[0])) {
          if (!rtIsNaN(varargin_1_data[1])) {
            ssq = varargin_1_data[1];
          } else {
            ssq = varargin_1_data[0];
          }
        } else {
          ssq = varargin_1_data[0];
        }
      } else {
        if (!rtIsNaN(varargin_1_data[0])) {
          idx = 1;
        } else {
          idx = 0;
          b_i = 2;
          exitg1 = false;
          while ((!exitg1) && (b_i <= ii_size_idx_0)) {
            if (!rtIsNaN(varargin_1_data[b_i - 1])) {
              idx = b_i;
              exitg1 = true;
            } else {
              b_i++;
            }
          }
        }

        if (idx == 0) {
          ssq = varargin_1_data[0];
        } else {
          ssq = varargin_1_data[idx - 1];
          i = idx + 1;
          for (b_i = i; b_i <= ii_size_idx_0; b_i++) {
            d = varargin_1_data[b_i - 1];
            if (ssq > d) {
              ssq = d;
            }
          }
        }
      }

      ssq *= 0.995;

      // solves for min ratio (max value of alpha allowed)
      // Increment
      for (i = 0; i < 120; i++) {
        lam[i] += ssq * del_lam[i];
        ftol[i] += ssq * mesil[i];
      }

      for (i = 0; i < 60; i++) {
        X_QP[i] += ssq * del_z[i];
      }
    }

    // Complimentary Gap
    mu_old = mu;
    ssq = 0.0;
    for (i = 0; i < 120; i++) {
      ssq += ftol[i] * lam[i];
    }

    mu = ssq / 120.0;
    if (mu < 0.001) {
      exitflag = 1;
    } else {
      // Solve for new Sigma
      ssq = mu / mu_old;
      if (ssq > 0.1) {
        // to tune
        ssq = 0.1;
      }

      for (b_i = 0; b_i < 120; b_i++) {
        esig[b_i] = ssq;
      }
    }

    unusedU0++;
  }

  // Check for failure
  //  num_iter = iter;
  //  solution to warm-start next iteration
  //  Extract first control
  Fx = X_QP[0];
  Fy = X_QP[1];
  Fz = X_QP[2];
}

template<typename T>
void CoordinatorBase<T>::mldivide(double A[3600], double B[60])
  {
    double b_A[3600];
    int i;
    int ix;
    int iy;
    int jA;
    int k;
    signed char ipiv[60];
    std::memcpy(&b_A[0], &A[0], 3600U * sizeof(double));
    for (i = 0; i < 60; i++) {
      ipiv[i] = static_cast<signed char>(i + 1);
    }

    for (int j = 0; j < 59; j++) {
      double smax;
      int b_tmp;
      int jp1j;
      int mmj_tmp;
      signed char i1;
      mmj_tmp = 58 - j;
      b_tmp = j * 61;
      jp1j = b_tmp + 2;
      iy = 60 - j;
      jA = 0;
      ix = b_tmp;
      smax = std::abs(b_A[b_tmp]);
      for (k = 2; k <= iy; k++) {
        double s;
        ix++;
        s = std::abs(b_A[ix]);
        if (s > smax) {
          jA = k - 1;
          smax = s;
        }
      }

      if (b_A[b_tmp + jA] != 0.0) {
        if (jA != 0) {
          iy = j + jA;
          ipiv[j] = static_cast<signed char>(iy + 1);
          ix = j;
          for (k = 0; k < 60; k++) {
            smax = b_A[ix];
            b_A[ix] = b_A[iy];
            b_A[iy] = smax;
            ix += 60;
            iy += 60;
          }
        }

        i = (b_tmp - j) + 60;
        for (jA = jp1j; jA <= i; jA++) {
          b_A[jA - 1] /= b_A[b_tmp];
        }
      }

      iy = b_tmp + 60;
      jA = b_tmp;
      for (k = 0; k <= mmj_tmp; k++) {
        smax = b_A[iy];
        if (b_A[iy] != 0.0) {
          ix = b_tmp + 1;
          i = jA + 62;
          jp1j = (jA - j) + 120;
          for (int ijA = i; ijA <= jp1j; ijA++) {
            b_A[ijA - 1] += b_A[ix] * -smax;
            ix++;
          }
        }

        iy += 60;
        jA += 60;
      }

      i1 = ipiv[j];
      if (i1 != j + 1) {
        smax = B[j];
        B[j] = B[i1 - 1];
        B[i1 - 1] = smax;
      }
    }

    for (k = 0; k < 60; k++) {
      iy = 60 * k;
      if (B[k] != 0.0) {
        i = k + 2;
        for (jA = i; jA < 61; jA++) {
          B[jA - 1] -= B[k] * b_A[(jA + iy) - 1];
        }
      }
    }

    for (k = 59; k >= 0; k--) {
      iy = 60 * k;
      if (B[k] != 0.0) {
        B[k] /= b_A[k + iy];
        for (jA = 0; jA < k; jA++) {
          B[jA] -= B[k] * b_A[jA + iy];
        }
      }
    }
  }

