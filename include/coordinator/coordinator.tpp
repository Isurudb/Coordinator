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
    double r=0, p=0, y=225;  // Rotate the previous pose by 45* about Z

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
    if (sqrt(q_e.getX()*q_e.getX()+q_e.getY()*q_e.getY()+q_e.getZ()*q_e.getZ())>0.05){
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

   double a[18] = { 0.069074, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.0,
    0.069074, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.0, 0.069074, 0.0, 0.0,
    0.069074 };

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
 double Omega2= arg_omegay;
 double Omega3= arg_omegaz;

 static const real_T a1[18] = { -0.1171782276523663, -0.0, -0.0, -0.0,
    -0.42978465701183144, -0.0, -0.0, -0.0, -0.35586060535946518,
    -1.171782276523663, -0.0, -0.0, -0.0, -4.2978465701183142, -0.0, -0.0, -0.0,
    -3.5586060535946515 };

  static const real_T b_a[9] = { 0.0, -0.048539486885480869, 0.0,
    0.048539486885480869, 0.0, 0.097078973770961738, -0.0, -0.097078973770961738,
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
    dv1[i] = ((dv[i] * fx) + (dv[i + 3] * fy)) + (dv[i + 6] * fz);
    d = 0.0;
    for (int32_T i1 = 0; i1 < 6; i1++) {
      d += a1[i + (3 * i1)] * b_qx[i1];
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
    0.0, 0.0, 0.0, 0.0, 1648.7366, -5.9399E-14, -5.506E-12, 12184.5822,
    -3.3197E-11, -1.4564E-10, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -5.9399E-14, 1648.7366, -2.7043E-11, 4.0889E-12, 12184.5822, -1.9781E-10,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -5.506E-12, -2.7043E-11,
    1648.7366, -4.7997E-11, -1.9479E-10, 12184.5822, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 12184.5822, 4.0889E-12, -4.7997E-11, 611791.602,
    -9.1704E-10, -1.9541E-9, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -3.3197E-11, 12184.5822, -1.9479E-10, -9.1704E-10, 611791.602, -7.7746E-10,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.4564E-10, -1.9781E-10,
    12184.5822, -1.9541E-9, -7.7746E-10, 611791.602 };

  static const real_T b_b[7200] = { 0.069074, 0.0, 0.0, 0.069074, 0.0, 0.0,
    0.20722, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.34537, 0.0, 0.0, 0.069074, 0.0, 0.0,
    0.48352, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.62167, 0.0, 0.0, 0.069074, 0.0, 0.0,
    0.75981, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.89796, 0.0, 0.0, 0.069074, 0.0, 0.0,
    1.0361, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.1743, 0.0, 0.0, 0.069074, 0.0, 0.0,
    1.3124, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.4506, 0.0, 0.0, 0.069074, 0.0, 0.0,
    1.5887, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.7268, 0.0, 0.0, 0.069074, 0.0, 0.0,
    1.865, 0.0, 0.0, 0.069074, 0.0, 0.0, 2.0031, 0.0, 0.0, 0.069074, 0.0, 0.0,
    2.1413, 0.0, 0.0, 0.069074, 0.0, 0.0, 2.2794, 0.0, 0.0, 0.069074, 0.0, 0.0,
    2.4176, 0.0, 0.0, 0.069074, 0.0, 0.0, 2.5557, 0.0, 0.0, 0.069074, 0.0, 0.0,
    2.6939, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.20722, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.34537, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.48352, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.62167, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.75981, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.89796, 0.0, 0.0, 0.069074, 0.0,
    0.0, 1.0361, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.1743, 0.0, 0.0, 0.069074, 0.0,
    0.0, 1.3124, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.4506, 0.0, 0.0, 0.069074, 0.0,
    0.0, 1.5887, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.7268, 0.0, 0.0, 0.069074, 0.0,
    0.0, 1.865, 0.0, 0.0, 0.069074, 0.0, 0.0, 2.0031, 0.0, 0.0, 0.069074, 0.0,
    0.0, 2.1413, 0.0, 0.0, 0.069074, 0.0, 0.0, 2.2794, 0.0, 0.0, 0.069074, 0.0,
    0.0, 2.4176, 0.0, 0.0, 0.069074, 0.0, 0.0, 2.5557, 0.0, 0.0, 0.069074, 0.0,
    0.0, 2.6939, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.069074,
    0.0, 0.0, 0.20722, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.34537, 0.0, 0.0, 0.069074,
    0.0, 0.0, 0.48352, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.62167, 0.0, 0.0, 0.069074,
    0.0, 0.0, 0.75981, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.89796, 0.0, 0.0, 0.069074,
    0.0, 0.0, 1.0361, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.1743, 0.0, 0.0, 0.069074,
    0.0, 0.0, 1.3124, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.4506, 0.0, 0.0, 0.069074,
    0.0, 0.0, 1.5887, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.7268, 0.0, 0.0, 0.069074,
    0.0, 0.0, 1.865, 0.0, 0.0, 0.069074, 0.0, 0.0, 2.0031, 0.0, 0.0, 0.069074,
    0.0, 0.0, 2.1413, 0.0, 0.0, 0.069074, 0.0, 0.0, 2.2794, 0.0, 0.0, 0.069074,
    0.0, 0.0, 2.4176, 0.0, 0.0, 0.069074, 0.0, 0.0, 2.5557, 0.0, 0.0, 0.069074,
    0.0, 0.0, 2.6939, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.069074,
    0.0, 0.0, 0.069074, 0.0, 0.0, 0.20722, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.34537,
    0.0, 0.0, 0.069074, 0.0, 0.0, 0.48352, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.62167,
    0.0, 0.0, 0.069074, 0.0, 0.0, 0.75981, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.89796,
    0.0, 0.0, 0.069074, 0.0, 0.0, 1.0361, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.1743,
    0.0, 0.0, 0.069074, 0.0, 0.0, 1.3124, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.4506,
    0.0, 0.0, 0.069074, 0.0, 0.0, 1.5887, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.7268,
    0.0, 0.0, 0.069074, 0.0, 0.0, 1.865, 0.0, 0.0, 0.069074, 0.0, 0.0, 2.0031,
    0.0, 0.0, 0.069074, 0.0, 0.0, 2.1413, 0.0, 0.0, 0.069074, 0.0, 0.0, 2.2794,
    0.0, 0.0, 0.069074, 0.0, 0.0, 2.4176, 0.0, 0.0, 0.069074, 0.0, 0.0, 2.5557,
    0.0, 0.0, 0.069074, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.069074,
    0.0, 0.0, 0.069074, 0.0, 0.0, 0.20722, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.34537,
    0.0, 0.0, 0.069074, 0.0, 0.0, 0.48352, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.62167,
    0.0, 0.0, 0.069074, 0.0, 0.0, 0.75981, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.89796,
    0.0, 0.0, 0.069074, 0.0, 0.0, 1.0361, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.1743,
    0.0, 0.0, 0.069074, 0.0, 0.0, 1.3124, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.4506,
    0.0, 0.0, 0.069074, 0.0, 0.0, 1.5887, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.7268,
    0.0, 0.0, 0.069074, 0.0, 0.0, 1.865, 0.0, 0.0, 0.069074, 0.0, 0.0, 2.0031,
    0.0, 0.0, 0.069074, 0.0, 0.0, 2.1413, 0.0, 0.0, 0.069074, 0.0, 0.0, 2.2794,
    0.0, 0.0, 0.069074, 0.0, 0.0, 2.4176, 0.0, 0.0, 0.069074, 0.0, 0.0, 2.5557,
    0.0, 0.0, 0.069074, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.069074,
    0.0, 0.0, 0.069074, 0.0, 0.0, 0.20722, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.34537,
    0.0, 0.0, 0.069074, 0.0, 0.0, 0.48352, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.62167,
    0.0, 0.0, 0.069074, 0.0, 0.0, 0.75981, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.89796,
    0.0, 0.0, 0.069074, 0.0, 0.0, 1.0361, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.1743,
    0.0, 0.0, 0.069074, 0.0, 0.0, 1.3124, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.4506,
    0.0, 0.0, 0.069074, 0.0, 0.0, 1.5887, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.7268,
    0.0, 0.0, 0.069074, 0.0, 0.0, 1.865, 0.0, 0.0, 0.069074, 0.0, 0.0, 2.0031,
    0.0, 0.0, 0.069074, 0.0, 0.0, 2.1413, 0.0, 0.0, 0.069074, 0.0, 0.0, 2.2794,
    0.0, 0.0, 0.069074, 0.0, 0.0, 2.4176, 0.0, 0.0, 0.069074, 0.0, 0.0, 2.5557,
    0.0, 0.0, 0.069074, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.069074, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.20722, 0.0, 0.0, 0.069074,
    0.0, 0.0, 0.34537, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.48352, 0.0, 0.0, 0.069074,
    0.0, 0.0, 0.62167, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.75981, 0.0, 0.0, 0.069074,
    0.0, 0.0, 0.89796, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.0361, 0.0, 0.0, 0.069074,
    0.0, 0.0, 1.1743, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.3124, 0.0, 0.0, 0.069074,
    0.0, 0.0, 1.4506, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.5887, 0.0, 0.0, 0.069074,
    0.0, 0.0, 1.7268, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.865, 0.0, 0.0, 0.069074,
    0.0, 0.0, 2.0031, 0.0, 0.0, 0.069074, 0.0, 0.0, 2.1413, 0.0, 0.0, 0.069074,
    0.0, 0.0, 2.2794, 0.0, 0.0, 0.069074, 0.0, 0.0, 2.4176, 0.0, 0.0, 0.069074,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.069074, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.20722, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.34537, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.48352, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.62167, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.75981, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.89796, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.0361, 0.0, 0.0, 0.069074, 0.0,
    0.0, 1.1743, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.3124, 0.0, 0.0, 0.069074, 0.0,
    0.0, 1.4506, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.5887, 0.0, 0.0, 0.069074, 0.0,
    0.0, 1.7268, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.865, 0.0, 0.0, 0.069074, 0.0,
    0.0, 2.0031, 0.0, 0.0, 0.069074, 0.0, 0.0, 2.1413, 0.0, 0.0, 0.069074, 0.0,
    0.0, 2.2794, 0.0, 0.0, 0.069074, 0.0, 0.0, 2.4176, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.069074, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.20722, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.34537, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.48352, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.62167, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.75981, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.89796, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.0361, 0.0, 0.0, 0.069074, 0.0,
    0.0, 1.1743, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.3124, 0.0, 0.0, 0.069074, 0.0,
    0.0, 1.4506, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.5887, 0.0, 0.0, 0.069074, 0.0,
    0.0, 1.7268, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.865, 0.0, 0.0, 0.069074, 0.0,
    0.0, 2.0031, 0.0, 0.0, 0.069074, 0.0, 0.0, 2.1413, 0.0, 0.0, 0.069074, 0.0,
    0.0, 2.2794, 0.0, 0.0, 0.069074, 0.0, 0.0, 2.4176, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.069074, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.20722, 0.0, 0.0,
    0.069074, 0.0, 0.0, 0.34537, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.48352, 0.0, 0.0,
    0.069074, 0.0, 0.0, 0.62167, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.75981, 0.0, 0.0,
    0.069074, 0.0, 0.0, 0.89796, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.0361, 0.0, 0.0,
    0.069074, 0.0, 0.0, 1.1743, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.3124, 0.0, 0.0,
    0.069074, 0.0, 0.0, 1.4506, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.5887, 0.0, 0.0,
    0.069074, 0.0, 0.0, 1.7268, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.865, 0.0, 0.0,
    0.069074, 0.0, 0.0, 2.0031, 0.0, 0.0, 0.069074, 0.0, 0.0, 2.1413, 0.0, 0.0,
    0.069074, 0.0, 0.0, 2.2794, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.069074, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.20722, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.34537, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.48352, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.62167, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.75981, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.89796, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.0361, 0.0, 0.0, 0.069074, 0.0,
    0.0, 1.1743, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.3124, 0.0, 0.0, 0.069074, 0.0,
    0.0, 1.4506, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.5887, 0.0, 0.0, 0.069074, 0.0,
    0.0, 1.7268, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.865, 0.0, 0.0, 0.069074, 0.0,
    0.0, 2.0031, 0.0, 0.0, 0.069074, 0.0, 0.0, 2.1413, 0.0, 0.0, 0.069074, 0.0,
    0.0, 2.2794, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.069074, 0.0, 0.0, 0.20722, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.34537, 0.0,
    0.0, 0.069074, 0.0, 0.0, 0.48352, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.62167, 0.0,
    0.0, 0.069074, 0.0, 0.0, 0.75981, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.89796, 0.0,
    0.0, 0.069074, 0.0, 0.0, 1.0361, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.1743, 0.0,
    0.0, 0.069074, 0.0, 0.0, 1.3124, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.4506, 0.0,
    0.0, 0.069074, 0.0, 0.0, 1.5887, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.7268, 0.0,
    0.0, 0.069074, 0.0, 0.0, 1.865, 0.0, 0.0, 0.069074, 0.0, 0.0, 2.0031, 0.0,
    0.0, 0.069074, 0.0, 0.0, 2.1413, 0.0, 0.0, 0.069074, 0.0, 0.0, 2.2794, 0.0,
    0.0, 0.069074, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.069074, 0.0, 0.0, 0.20722, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.34537, 0.0,
    0.0, 0.069074, 0.0, 0.0, 0.48352, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.62167, 0.0,
    0.0, 0.069074, 0.0, 0.0, 0.75981, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.89796, 0.0,
    0.0, 0.069074, 0.0, 0.0, 1.0361, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.1743, 0.0,
    0.0, 0.069074, 0.0, 0.0, 1.3124, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.4506, 0.0,
    0.0, 0.069074, 0.0, 0.0, 1.5887, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.7268, 0.0,
    0.0, 0.069074, 0.0, 0.0, 1.865, 0.0, 0.0, 0.069074, 0.0, 0.0, 2.0031, 0.0,
    0.0, 0.069074, 0.0, 0.0, 2.1413, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.069074, 0.0, 0.0,
    0.20722, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.34537, 0.0, 0.0, 0.069074, 0.0, 0.0,
    0.48352, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.62167, 0.0, 0.0, 0.069074, 0.0, 0.0,
    0.75981, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.89796, 0.0, 0.0, 0.069074, 0.0, 0.0,
    1.0361, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.1743, 0.0, 0.0, 0.069074, 0.0, 0.0,
    1.3124, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.4506, 0.0, 0.0, 0.069074, 0.0, 0.0,
    1.5887, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.7268, 0.0, 0.0, 0.069074, 0.0, 0.0,
    1.865, 0.0, 0.0, 0.069074, 0.0, 0.0, 2.0031, 0.0, 0.0, 0.069074, 0.0, 0.0,
    2.1413, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.069074, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.20722, 0.0, 0.0,
    0.069074, 0.0, 0.0, 0.34537, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.48352, 0.0, 0.0,
    0.069074, 0.0, 0.0, 0.62167, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.75981, 0.0, 0.0,
    0.069074, 0.0, 0.0, 0.89796, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.0361, 0.0, 0.0,
    0.069074, 0.0, 0.0, 1.1743, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.3124, 0.0, 0.0,
    0.069074, 0.0, 0.0, 1.4506, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.5887, 0.0, 0.0,
    0.069074, 0.0, 0.0, 1.7268, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.865, 0.0, 0.0,
    0.069074, 0.0, 0.0, 2.0031, 0.0, 0.0, 0.069074, 0.0, 0.0, 2.1413, 0.0, 0.0,
    0.069074, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.069074, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.20722, 0.0, 0.0,
    0.069074, 0.0, 0.0, 0.34537, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.48352, 0.0, 0.0,
    0.069074, 0.0, 0.0, 0.62167, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.75981, 0.0, 0.0,
    0.069074, 0.0, 0.0, 0.89796, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.0361, 0.0, 0.0,
    0.069074, 0.0, 0.0, 1.1743, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.3124, 0.0, 0.0,
    0.069074, 0.0, 0.0, 1.4506, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.5887, 0.0, 0.0,
    0.069074, 0.0, 0.0, 1.7268, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.865, 0.0, 0.0,
    0.069074, 0.0, 0.0, 2.0031, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.069074, 0.0, 0.0, 0.20722, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.34537, 0.0,
    0.0, 0.069074, 0.0, 0.0, 0.48352, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.62167, 0.0,
    0.0, 0.069074, 0.0, 0.0, 0.75981, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.89796, 0.0,
    0.0, 0.069074, 0.0, 0.0, 1.0361, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.1743, 0.0,
    0.0, 0.069074, 0.0, 0.0, 1.3124, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.4506, 0.0,
    0.0, 0.069074, 0.0, 0.0, 1.5887, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.7268, 0.0,
    0.0, 0.069074, 0.0, 0.0, 1.865, 0.0, 0.0, 0.069074, 0.0, 0.0, 2.0031, 0.0,
    0.0, 0.069074, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.069074, 0.0, 0.0,
    0.20722, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.34537, 0.0, 0.0, 0.069074, 0.0, 0.0,
    0.48352, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.62167, 0.0, 0.0, 0.069074, 0.0, 0.0,
    0.75981, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.89796, 0.0, 0.0, 0.069074, 0.0, 0.0,
    1.0361, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.1743, 0.0, 0.0, 0.069074, 0.0, 0.0,
    1.3124, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.4506, 0.0, 0.0, 0.069074, 0.0, 0.0,
    1.5887, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.7268, 0.0, 0.0, 0.069074, 0.0, 0.0,
    1.865, 0.0, 0.0, 0.069074, 0.0, 0.0, 2.0031, 0.0, 0.0, 0.069074, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.20722, 0.0,
    0.0, 0.069074, 0.0, 0.0, 0.34537, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.48352, 0.0,
    0.0, 0.069074, 0.0, 0.0, 0.62167, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.75981, 0.0,
    0.0, 0.069074, 0.0, 0.0, 0.89796, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.0361, 0.0,
    0.0, 0.069074, 0.0, 0.0, 1.1743, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.3124, 0.0,
    0.0, 0.069074, 0.0, 0.0, 1.4506, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.5887, 0.0,
    0.0, 0.069074, 0.0, 0.0, 1.7268, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.865, 0.0,
    0.0, 0.069074, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.069074, 0.0, 0.0, 0.20722, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.34537, 0.0,
    0.0, 0.069074, 0.0, 0.0, 0.48352, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.62167, 0.0,
    0.0, 0.069074, 0.0, 0.0, 0.75981, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.89796, 0.0,
    0.0, 0.069074, 0.0, 0.0, 1.0361, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.1743, 0.0,
    0.0, 0.069074, 0.0, 0.0, 1.3124, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.4506, 0.0,
    0.0, 0.069074, 0.0, 0.0, 1.5887, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.7268, 0.0,
    0.0, 0.069074, 0.0, 0.0, 1.865, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.20722, 0.0,
    0.0, 0.069074, 0.0, 0.0, 0.34537, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.48352, 0.0,
    0.0, 0.069074, 0.0, 0.0, 0.62167, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.75981, 0.0,
    0.0, 0.069074, 0.0, 0.0, 0.89796, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.0361, 0.0,
    0.0, 0.069074, 0.0, 0.0, 1.1743, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.3124, 0.0,
    0.0, 0.069074, 0.0, 0.0, 1.4506, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.5887, 0.0,
    0.0, 0.069074, 0.0, 0.0, 1.7268, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.865, 0.0,
    0.0, 0.069074, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.069074, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.20722, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.34537, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.48352, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.62167, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.75981, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.89796, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.0361, 0.0, 0.0, 0.069074, 0.0,
    0.0, 1.1743, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.3124, 0.0, 0.0, 0.069074, 0.0,
    0.0, 1.4506, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.5887, 0.0, 0.0, 0.069074, 0.0,
    0.0, 1.7268, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.069074, 0.0, 0.0,
    0.20722, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.34537, 0.0, 0.0, 0.069074, 0.0, 0.0,
    0.48352, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.62167, 0.0, 0.0, 0.069074, 0.0, 0.0,
    0.75981, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.89796, 0.0, 0.0, 0.069074, 0.0, 0.0,
    1.0361, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.1743, 0.0, 0.0, 0.069074, 0.0, 0.0,
    1.3124, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.4506, 0.0, 0.0, 0.069074, 0.0, 0.0,
    1.5887, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.7268, 0.0, 0.0, 0.069074, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.069074,
    0.0, 0.0, 0.069074, 0.0, 0.0, 0.20722, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.34537,
    0.0, 0.0, 0.069074, 0.0, 0.0, 0.48352, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.62167,
    0.0, 0.0, 0.069074, 0.0, 0.0, 0.75981, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.89796,
    0.0, 0.0, 0.069074, 0.0, 0.0, 1.0361, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.1743,
    0.0, 0.0, 0.069074, 0.0, 0.0, 1.3124, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.4506,
    0.0, 0.0, 0.069074, 0.0, 0.0, 1.5887, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.7268,
    0.0, 0.0, 0.069074, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.069074, 0.0, 0.0,
    0.20722, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.34537, 0.0, 0.0, 0.069074, 0.0, 0.0,
    0.48352, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.62167, 0.0, 0.0, 0.069074, 0.0, 0.0,
    0.75981, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.89796, 0.0, 0.0, 0.069074, 0.0, 0.0,
    1.0361, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.1743, 0.0, 0.0, 0.069074, 0.0, 0.0,
    1.3124, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.4506, 0.0, 0.0, 0.069074, 0.0, 0.0,
    1.5887, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.069074, 0.0, 0.0,
    0.069074, 0.0, 0.0, 0.20722, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.34537, 0.0, 0.0,
    0.069074, 0.0, 0.0, 0.48352, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.62167, 0.0, 0.0,
    0.069074, 0.0, 0.0, 0.75981, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.89796, 0.0, 0.0,
    0.069074, 0.0, 0.0, 1.0361, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.1743, 0.0, 0.0,
    0.069074, 0.0, 0.0, 1.3124, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.4506, 0.0, 0.0,
    0.069074, 0.0, 0.0, 1.5887, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.069074, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.20722, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.34537, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.48352, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.62167, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.75981, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.89796, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.0361, 0.0, 0.0, 0.069074, 0.0,
    0.0, 1.1743, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.3124, 0.0, 0.0, 0.069074, 0.0,
    0.0, 1.4506, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.5887, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.20722, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.34537, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.48352, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.62167, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.75981, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.89796, 0.0, 0.0, 0.069074, 0.0,
    0.0, 1.0361, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.1743, 0.0, 0.0, 0.069074, 0.0,
    0.0, 1.3124, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.4506, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.069074, 0.0, 0.0,
    0.069074, 0.0, 0.0, 0.20722, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.34537, 0.0, 0.0,
    0.069074, 0.0, 0.0, 0.48352, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.62167, 0.0, 0.0,
    0.069074, 0.0, 0.0, 0.75981, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.89796, 0.0, 0.0,
    0.069074, 0.0, 0.0, 1.0361, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.1743, 0.0, 0.0,
    0.069074, 0.0, 0.0, 1.3124, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.4506, 0.0, 0.0,
    0.069074, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.069074, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.20722, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.34537, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.48352, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.62167, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.75981, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.89796, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.0361, 0.0, 0.0, 0.069074, 0.0,
    0.0, 1.1743, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.3124, 0.0, 0.0, 0.069074, 0.0,
    0.0, 1.4506, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.069074, 0.0, 0.0,
    0.20722, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.34537, 0.0, 0.0, 0.069074, 0.0, 0.0,
    0.48352, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.62167, 0.0, 0.0, 0.069074, 0.0, 0.0,
    0.75981, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.89796, 0.0, 0.0, 0.069074, 0.0, 0.0,
    1.0361, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.1743, 0.0, 0.0, 0.069074, 0.0, 0.0,
    1.3124, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.20722, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.34537, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.48352, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.62167, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.75981, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.89796, 0.0, 0.0, 0.069074, 0.0,
    0.0, 1.0361, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.1743, 0.0, 0.0, 0.069074, 0.0,
    0.0, 1.3124, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.069074,
    0.0, 0.0, 0.20722, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.34537, 0.0, 0.0, 0.069074,
    0.0, 0.0, 0.48352, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.62167, 0.0, 0.0, 0.069074,
    0.0, 0.0, 0.75981, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.89796, 0.0, 0.0, 0.069074,
    0.0, 0.0, 1.0361, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.1743, 0.0, 0.0, 0.069074,
    0.0, 0.0, 1.3124, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.069074,
    0.0, 0.0, 0.069074, 0.0, 0.0, 0.20722, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.34537,
    0.0, 0.0, 0.069074, 0.0, 0.0, 0.48352, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.62167,
    0.0, 0.0, 0.069074, 0.0, 0.0, 0.75981, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.89796,
    0.0, 0.0, 0.069074, 0.0, 0.0, 1.0361, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.1743,
    0.0, 0.0, 0.069074, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.069074,
    0.0, 0.0, 0.069074, 0.0, 0.0, 0.20722, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.34537,
    0.0, 0.0, 0.069074, 0.0, 0.0, 0.48352, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.62167,
    0.0, 0.0, 0.069074, 0.0, 0.0, 0.75981, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.89796,
    0.0, 0.0, 0.069074, 0.0, 0.0, 1.0361, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.1743,
    0.0, 0.0, 0.069074, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.069074,
    0.0, 0.0, 0.069074, 0.0, 0.0, 0.20722, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.34537,
    0.0, 0.0, 0.069074, 0.0, 0.0, 0.48352, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.62167,
    0.0, 0.0, 0.069074, 0.0, 0.0, 0.75981, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.89796,
    0.0, 0.0, 0.069074, 0.0, 0.0, 1.0361, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.1743,
    0.0, 0.0, 0.069074, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.069074, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.20722, 0.0, 0.0, 0.069074,
    0.0, 0.0, 0.34537, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.48352, 0.0, 0.0, 0.069074,
    0.0, 0.0, 0.62167, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.75981, 0.0, 0.0, 0.069074,
    0.0, 0.0, 0.89796, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.0361, 0.0, 0.0, 0.069074,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.069074, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.20722, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.34537, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.48352, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.62167, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.75981, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.89796, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.0361, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.069074, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.20722, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.34537, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.48352, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.62167, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.75981, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.89796, 0.0, 0.0, 0.069074, 0.0, 0.0, 1.0361, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.069074, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.20722, 0.0, 0.0,
    0.069074, 0.0, 0.0, 0.34537, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.48352, 0.0, 0.0,
    0.069074, 0.0, 0.0, 0.62167, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.75981, 0.0, 0.0,
    0.069074, 0.0, 0.0, 0.89796, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.069074, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.20722, 0.0, 0.0, 0.069074,
    0.0, 0.0, 0.34537, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.48352, 0.0, 0.0, 0.069074,
    0.0, 0.0, 0.62167, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.75981, 0.0, 0.0, 0.069074,
    0.0, 0.0, 0.89796, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.069074, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.20722, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.34537, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.48352, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.62167, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.75981, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.89796, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.069074, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.20722, 0.0, 0.0, 0.069074,
    0.0, 0.0, 0.34537, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.48352, 0.0, 0.0, 0.069074,
    0.0, 0.0, 0.62167, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.75981, 0.0, 0.0, 0.069074,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.069074, 0.0, 0.0, 0.20722, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.34537, 0.0,
    0.0, 0.069074, 0.0, 0.0, 0.48352, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.62167, 0.0,
    0.0, 0.069074, 0.0, 0.0, 0.75981, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.20722, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.34537, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.48352, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.62167, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.75981, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.069074, 0.0, 0.0,
    0.20722, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.34537, 0.0, 0.0, 0.069074, 0.0, 0.0,
    0.48352, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.62167, 0.0, 0.0, 0.069074, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.069074, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.20722, 0.0, 0.0, 0.069074,
    0.0, 0.0, 0.34537, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.48352, 0.0, 0.0, 0.069074,
    0.0, 0.0, 0.62167, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.069074, 0.0, 0.0,
    0.069074, 0.0, 0.0, 0.20722, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.34537, 0.0, 0.0,
    0.069074, 0.0, 0.0, 0.48352, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.62167, 0.0, 0.0,
    0.069074, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.20722, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.34537, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.48352, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.069074, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.20722, 0.0, 0.0, 0.069074,
    0.0, 0.0, 0.34537, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.48352, 0.0, 0.0, 0.069074,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.069074,
    0.0, 0.0, 0.20722, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.34537, 0.0, 0.0, 0.069074,
    0.0, 0.0, 0.48352, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.20722,
    0.0, 0.0, 0.069074, 0.0, 0.0, 0.34537, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.069074, 0.0, 0.0, 0.20722, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.34537, 0.0,
    0.0, 0.069074, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.20722, 0.0, 0.0,
    0.069074, 0.0, 0.0, 0.34537, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.069074,
    0.0, 0.0, 0.069074, 0.0, 0.0, 0.20722, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.20722, 0.0, 0.0,
    0.069074, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.069074, 0.0,
    0.0, 0.20722, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.069074, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.069074, 0.0, 0.0, 0.069074 };

  static const real_T b_H[3600] = { 83072.1302, -2.0446E-11, -1.7062E-10,
    60176.5678, -2.0479E-11, -1.6561E-10, 57283.6773, -2.0513E-11, -1.606E-10,
    54396.1306, -2.0547E-11, -1.5558E-10, 51516.5995, -2.0581E-11, -1.5057E-10,
    48647.756, -2.0615E-11, -1.4555E-10, 45792.2718, -2.0649E-11, -1.4054E-10,
    42952.819, -2.0682E-11, -1.3552E-10, 40132.0693, -2.0716E-11, -1.3051E-10,
    37332.6947, -2.075E-11, -1.255E-10, 34557.367, -2.0784E-11, -1.2048E-10,
    31808.7581, -2.0818E-11, -1.1547E-10, 29089.5399, -2.0852E-11, -1.1045E-10,
    26402.3844, -2.0885E-11, -1.0544E-10, 23749.9632, -2.0919E-11, -1.0043E-10,
    21134.9484, -2.0953E-11, -9.5411E-11, 18560.0119, -2.0987E-11, -9.0397E-11,
    16027.8254, -2.1021E-11, -8.5383E-11, 13541.061, -2.1054E-11, -8.0369E-11,
    11102.3904, -2.1088E-11, -7.5355E-11, -2.0446E-11, 83072.1302, -5.4603E-10,
    -1.9768E-11, 60176.5678, -5.2218E-10, -1.909E-11, 57283.6773, -4.9834E-10,
    -1.8412E-11, 54396.1306, -4.7449E-10, -1.7735E-11, 51516.5995, -4.5065E-10,
    -1.7057E-11, 48647.756, -4.268E-10, -1.6379E-11, 45792.2718, -4.0295E-10,
    -1.5701E-11, 42952.819, -3.7911E-10, -1.5023E-11, 40132.0693, -3.5526E-10,
    -1.4346E-11, 37332.6947, -3.3142E-10, -1.3668E-11, 34557.367, -3.0757E-10,
    -1.299E-11, 31808.7581, -2.8373E-10, -1.2312E-11, 29089.5399, -2.5988E-10,
    -1.1635E-11, 26402.3844, -2.3603E-10, -1.0957E-11, 23749.9632, -2.1219E-10,
    -1.0279E-11, 21134.9484, -1.8834E-10, -9.6013E-12, 18560.0119, -1.645E-10,
    -8.9235E-12, 16027.8254, -1.4065E-10, -8.2457E-12, 13541.061, -1.168E-10,
    -7.568E-12, 11102.3904, -9.2958E-11, -1.7062E-10, -5.4603E-10, 83072.1302,
    -1.6375E-10, -5.2213E-10, 60176.5678, -1.5687E-10, -4.9822E-10, 57283.6773,
    -1.4999E-10, -4.7432E-10, 54396.1306, -1.4311E-10, -4.5042E-10, 51516.5995,
    -1.3624E-10, -4.2651E-10, 48647.756, -1.2936E-10, -4.0261E-10, 45792.2718,
    -1.2248E-10, -3.7871E-10, 42952.819, -1.156E-10, -3.548E-10, 40132.0693,
    -1.0872E-10, -3.309E-10, 37332.6947, -1.0185E-10, -3.0699E-10, 34557.367,
    -9.497E-11, -2.8309E-10, 31808.7581, -8.8092E-11, -2.5919E-10, 29089.5399,
    -8.1214E-11, -2.3528E-10, 26402.3844, -7.4336E-11, -2.1138E-10, 23749.9632,
    -6.7459E-11, -1.8748E-10, 21134.9484, -6.0581E-11, -1.6357E-10, 18560.0119,
    -5.3703E-11, -1.3967E-10, 16027.8254, -4.6826E-11, -1.1577E-10, 13541.061,
    -3.9948E-11, -9.1863E-11, 11102.3904, 60176.5678, -1.9768E-11, -1.6375E-10,
    78346.9426, -1.9804E-11, -1.5894E-10, 55563.7419, -1.984E-11, -1.5414E-10,
    52783.2131, -1.9876E-11, -1.4933E-10, 50008.028, -1.9912E-11, -1.4453E-10,
    47240.8586, -1.9948E-11, -1.3973E-10, 44484.3767, -1.9984E-11, -1.3492E-10,
    41741.2542, -2.0021E-11, -1.3012E-10, 39014.1631, -2.0057E-11, -1.2531E-10,
    36305.775, -2.0093E-11, -1.2051E-10, 33618.7621, -2.0129E-11, -1.1571E-10,
    30955.796, -2.0165E-11, -1.109E-10, 28319.5488, -2.0201E-11, -1.061E-10,
    25712.6923, -2.0237E-11, -1.0129E-10, 23137.8984, -2.0273E-11, -9.649E-11,
    20597.8389, -2.0309E-11, -9.1686E-11, 18095.1858, -2.0345E-11, -8.6882E-11,
    15632.6109, -2.0381E-11, -8.2078E-11, 13212.7861, -2.0418E-11, -7.7274E-11,
    10838.3833, -2.0454E-11, -7.247E-11, -2.0479E-11, 60176.5678, -5.2213E-10,
    -1.9804E-11, 78346.9426, -4.9931E-10, -1.9128E-11, 55563.7419, -4.765E-10,
    -1.8453E-11, 52783.2131, -4.5369E-10, -1.7777E-11, 50008.028, -4.3087E-10,
    -1.7102E-11, 47240.8586, -4.0806E-10, -1.6426E-11, 44484.3767, -3.8524E-10,
    -1.5751E-11, 41741.2542, -3.6243E-10, -1.5075E-11, 39014.1631, -3.3962E-10,
    -1.44E-11, 36305.775, -3.168E-10, -1.3724E-11, 33618.7621, -2.9399E-10,
    -1.3049E-11, 30955.796, -2.7118E-10, -1.2373E-11, 28319.5488, -2.4836E-10,
    -1.1698E-11, 25712.6923, -2.2555E-10, -1.1022E-11, 23137.8984, -2.0273E-10,
    -1.0347E-11, 20597.8389, -1.7992E-10, -9.6714E-12, 18095.1858, -1.5711E-10,
    -8.9959E-12, 15632.6109, -1.3429E-10, -8.3204E-12, 13212.7861, -1.1148E-10,
    -7.6449E-12, 10838.3833, -8.8667E-11, -1.6561E-10, -5.2218E-10, 60176.5678,
    -1.5894E-10, -4.9931E-10, 78346.9426, -1.5227E-10, -4.7644E-10, 55563.7419,
    -1.4561E-10, -4.5357E-10, 52783.2131, -1.3894E-10, -4.307E-10, 50008.028,
    -1.3227E-10, -4.0783E-10, 47240.8586, -1.256E-10, -3.8496E-10, 44484.3767,
    -1.1894E-10, -3.6208E-10, 41741.2542, -1.1227E-10, -3.3921E-10, 39014.1631,
    -1.056E-10, -3.1634E-10, 36305.775, -9.8935E-11, -2.9347E-10, 33618.7621,
    -9.2267E-11, -2.706E-10, 30955.796, -8.56E-11, -2.4773E-10, 28319.5488,
    -7.8932E-11, -2.2486E-10, 25712.6923, -7.2265E-11, -2.0199E-10, 23137.8984,
    -6.5597E-11, -1.7911E-10, 20597.8389, -5.893E-11, -1.5624E-10, 18095.1858,
    -5.2262E-11, -1.3337E-10, 15632.6109, -4.5594E-11, -1.105E-10, 13212.7861,
    -3.8927E-11, -8.7629E-11, 10838.3833, 57283.6773, -1.909E-11, -1.5687E-10,
    55563.7419, -1.9128E-11, -1.5227E-10, 73843.8065, -1.9167E-11, -1.4768E-10,
    51170.2956, -1.9205E-11, -1.4309E-10, 48499.4565, -1.9244E-11, -1.3849E-10,
    45833.9613, -1.9282E-11, -1.339E-10, 43176.4816, -1.932E-11, -1.2931E-10,
    40529.6895, -1.9359E-11, -1.2471E-10, 37896.2568, -1.9397E-11, -1.2012E-10,
    35278.8554, -1.9435E-11, -1.1552E-10, 32680.1571, -1.9474E-11, -1.1093E-10,
    30102.834, -1.9512E-11, -1.0634E-10, 27549.5577, -1.955E-11, -1.0174E-10,
    25023.0003, -1.9589E-11, -9.7149E-11, 22525.8335, -1.9627E-11, -9.2555E-11,
    20060.7294, -1.9665E-11, -8.7961E-11, 17630.3597, -1.9704E-11, -8.3367E-11,
    15237.3963, -1.9742E-11, -7.8773E-11, 12884.5112, -1.9781E-11, -7.4179E-11,
    10574.3762, -1.9819E-11, -6.9585E-11, -2.0513E-11, 57283.6773, -4.9822E-10,
    -1.984E-11, 55563.7419, -4.7644E-10, -1.9167E-11, 73843.8065, -4.5466E-10,
    -1.8494E-11, 51170.2956, -4.3288E-10, -1.782E-11, 48499.4565, -4.111E-10,
    -1.7147E-11, 45833.9613, -3.8932E-10, -1.6474E-11, 43176.4816, -3.6753E-10,
    -1.5801E-11, 40529.6895, -3.4575E-10, -1.5127E-11, 37896.2568, -3.2397E-10,
    -1.4454E-11, 35278.8554, -3.0219E-10, -1.3781E-11, 32680.1571, -2.8041E-10,
    -1.3108E-11, 30102.834, -2.5863E-10, -1.2434E-11, 27549.5577, -2.3685E-10,
    -1.1761E-11, 25023.0003, -2.1506E-10, -1.1088E-11, 22525.8335, -1.9328E-10,
    -1.0415E-11, 20060.7294, -1.715E-10, -9.7415E-12, 17630.3597, -1.4972E-10,
    -9.0683E-12, 15237.3963, -1.2794E-10, -8.395E-12, 12884.5112, -1.0616E-10,
    -7.7218E-12, 10574.3762, -8.4375E-11, -1.606E-10, -4.9834E-10, 57283.6773,
    -1.5414E-10, -4.765E-10, 55563.7419, -1.4768E-10, -4.5466E-10, 73843.8065,
    -1.4122E-10, -4.3282E-10, 51170.2956, -1.3477E-10, -4.1098E-10, 48499.4565,
    -1.2831E-10, -3.8914E-10, 45833.9613, -1.2185E-10, -3.673E-10, 43176.4816,
    -1.1539E-10, -3.4546E-10, 40529.6895, -1.0894E-10, -3.2363E-10, 37896.2568,
    -1.0248E-10, -3.0179E-10, 35278.8554, -9.6022E-11, -2.7995E-10, 32680.1571,
    -8.9565E-11, -2.5811E-10, 30102.834, -8.3107E-11, -2.3627E-10, 27549.5577,
    -7.665E-11, -2.1443E-10, 25023.0003, -7.0193E-11, -1.9259E-10, 22525.8335,
    -6.3735E-11, -1.7075E-10, 20060.7294, -5.7278E-11, -1.4891E-10, 17630.3597,
    -5.0821E-11, -1.2707E-10, 15237.3963, -4.4363E-11, -1.0523E-10, 12884.5112,
    -3.7906E-11, -8.3395E-11, 10574.3762, 54396.1306, -1.8412E-11, -1.4999E-10,
    52783.2131, -1.8453E-11, -1.4561E-10, 51170.2956, -1.8494E-11, -1.4122E-10,
    69557.3781, -1.8534E-11, -1.3684E-10, 46990.885, -1.8575E-11, -1.3246E-10,
    44427.0639, -1.8615E-11, -1.2807E-10, 41868.5865, -1.8656E-11, -1.2369E-10,
    39318.1248, -1.8697E-11, -1.1931E-10, 36778.3505, -1.8737E-11, -1.1492E-10,
    34251.9357, -1.8778E-11, -1.1054E-10, 31741.5522, -1.8819E-11, -1.0615E-10,
    29249.8719, -1.8859E-11, -1.0177E-10, 26779.5666, -1.89E-11, -9.7387E-11,
    24333.3082, -1.894E-11, -9.3003E-11, 21913.7687, -1.8981E-11, -8.8619E-11,
    19523.6198, -1.9022E-11, -8.4236E-11, 17165.5336, -1.9062E-11, -7.9852E-11,
    14842.1818, -1.9103E-11, -7.5468E-11, 12556.2363, -1.9144E-11, -7.1085E-11,
    10310.3691, -1.9184E-11, -6.6701E-11, -2.0547E-11, 54396.1306, -4.7432E-10,
    -1.9876E-11, 52783.2131, -4.5357E-10, -1.9205E-11, 51170.2956, -4.3282E-10,
    -1.8534E-11, 69557.3781, -4.1207E-10, -1.7863E-11, 46990.885, -3.9132E-10,
    -1.7192E-11, 44427.0639, -3.7057E-10, -1.6521E-11, 41868.5865, -3.4982E-10,
    -1.585E-11, 39318.1248, -3.2907E-10, -1.5179E-11, 36778.3505, -3.0833E-10,
    -1.4508E-11, 34251.9357, -2.8758E-10, -1.3837E-11, 31741.5522, -2.6683E-10,
    -1.3166E-11, 29249.8719, -2.4608E-10, -1.2495E-11, 26779.5666, -2.2533E-10,
    -1.1825E-11, 24333.3082, -2.0458E-10, -1.1154E-11, 21913.7687, -1.8383E-10,
    -1.0483E-11, 19523.6198, -1.6308E-10, -9.8116E-12, 17165.5336, -1.4233E-10,
    -9.1406E-12, 14842.1818, -1.2158E-10, -8.4697E-12, 12556.2363, -1.0083E-10,
    -7.7987E-12, 10310.3691, -8.0084E-11, -1.5558E-10, -4.7449E-10, 54396.1306,
    -1.4933E-10, -4.5369E-10, 52783.2131, -1.4309E-10, -4.3288E-10, 51170.2956,
    -1.3684E-10, -4.1207E-10, 69557.3781, -1.3059E-10, -3.9126E-10, 46990.885,
    -1.2435E-10, -3.7046E-10, 44427.0639, -1.181E-10, -3.4965E-10, 41868.5865,
    -1.1185E-10, -3.2884E-10, 39318.1248, -1.056E-10, -3.0804E-10, 36778.3505,
    -9.9357E-11, -2.8723E-10, 34251.9357, -9.311E-11, -2.6642E-10, 31741.5522,
    -8.6862E-11, -2.4562E-10, 29249.8719, -8.0615E-11, -2.2481E-10, 26779.5666,
    -7.4368E-11, -2.04E-10, 24333.3082, -6.8121E-11, -1.832E-10, 21913.7687,
    -6.1874E-11, -1.6239E-10, 19523.6198, -5.5626E-11, -1.4158E-10, 17165.5336,
    -4.9379E-11, -1.2078E-10, 14842.1818, -4.3132E-11, -9.9969E-11, 12556.2363,
    -3.6885E-11, -7.9162E-11, 10310.3691, 51516.5995, -1.7735E-11, -1.4311E-10,
    50008.028, -1.7777E-11, -1.3894E-10, 48499.4565, -1.782E-11, -1.3477E-10,
    46990.885, -1.7863E-11, -1.3059E-10, 65482.3136, -1.7906E-11, -1.2642E-10,
    43020.1665, -1.7949E-11, -1.2225E-10, 40560.6914, -1.7992E-11, -1.1807E-10,
    38106.56, -1.8035E-11, -1.139E-10, 35660.4443, -1.8078E-11, -1.0973E-10,
    33225.0161, -1.8121E-11, -1.0555E-10, 30802.9473, -1.8163E-11, -1.0138E-10,
    28396.9098, -1.8206E-11, -9.7205E-11, 26009.5755, -1.8249E-11, -9.3031E-11,
    23643.6162, -1.8292E-11, -8.8857E-11, 21301.7038, -1.8335E-11, -8.4684E-11,
    18986.5103, -1.8378E-11, -8.051E-11, 16700.7075, -1.8421E-11, -7.6337E-11,
    14446.9673, -1.8464E-11, -7.2163E-11, 12227.9615, -1.8507E-11, -6.799E-11,
    10046.362, -1.855E-11, -6.3816E-11, -2.0581E-11, 51516.5995, -4.5042E-10,
    -1.9912E-11, 50008.028, -4.307E-10, -1.9244E-11, 48499.4565, -4.1098E-10,
    -1.8575E-11, 46990.885, -3.9126E-10, -1.7906E-11, 65482.3136, -3.7155E-10,
    -1.7237E-11, 43020.1665, -3.5183E-10, -1.6569E-11, 40560.6914, -3.3211E-10,
    -1.59E-11, 38106.56, -3.124E-10, -1.5231E-11, 35660.4443, -2.9268E-10,
    -1.4563E-11, 33225.0161, -2.7296E-10, -1.3894E-11, 30802.9473, -2.5325E-10,
    -1.3225E-11, 28396.9098, -2.3353E-10, -1.2557E-11, 26009.5755, -2.1381E-10,
    -1.1888E-11, 23643.6162, -1.9409E-10, -1.1219E-11, 21301.7038, -1.7438E-10,
    -1.055E-11, 18986.5103, -1.5466E-10, -9.8817E-12, 16700.7075, -1.3494E-10,
    -9.213E-12, 14446.9673, -1.1523E-10, -8.5443E-12, 12227.9615, -9.551E-11,
    -7.8756E-12, 10046.362, -7.5793E-11, -1.5057E-10, -4.5065E-10, 51516.5995,
    -1.4453E-10, -4.3087E-10, 50008.028, -1.3849E-10, -4.111E-10, 48499.4565,
    -1.3246E-10, -3.9132E-10, 46990.885, -1.2642E-10, -3.7155E-10, 65482.3136,
    -1.2038E-10, -3.5177E-10, 43020.1665, -1.1435E-10, -3.32E-10, 40560.6914,
    -1.0831E-10, -3.1222E-10, 38106.56, -1.0227E-10, -2.9245E-10, 35660.4443,
    -9.6234E-11, -2.7267E-10, 33225.0161, -9.0197E-11, -2.529E-10, 30802.9473,
    -8.416E-11, -2.3313E-10, 28396.9098, -7.8123E-11, -2.1335E-10, 26009.5755,
    -7.2086E-11, -1.9358E-10, 23643.6162, -6.6049E-11, -1.738E-10, 21301.7038,
    -6.0012E-11, -1.5403E-10, 18986.5103, -5.3975E-11, -1.3425E-10, 16700.7075,
    -4.7938E-11, -1.1448E-10, 14446.9673, -4.1901E-11, -9.4703E-11, 12227.9615,
    -3.5864E-11, -7.4928E-11, 10046.362, 48647.756, -1.7057E-11, -1.3624E-10,
    47240.8586, -1.7102E-11, -1.3227E-10, 45833.9613, -1.7147E-11, -1.2831E-10,
    44427.0639, -1.7192E-11, -1.2435E-10, 43020.1665, -1.7237E-11, -1.2038E-10,
    61613.2692, -1.7283E-11, -1.1642E-10, 39252.7963, -1.7328E-11, -1.1246E-10,
    36894.9953, -1.7373E-11, -1.0849E-10, 34542.538, -1.7418E-11, -1.0453E-10,
    32198.0964, -1.7463E-11, -1.0057E-10, 29864.3424, -1.7508E-11, -9.6602E-11,
    27543.9477, -1.7554E-11, -9.2639E-11, 25239.5843, -1.7599E-11, -8.8675E-11,
    22953.9241, -1.7644E-11, -8.4712E-11, 20689.639, -1.7689E-11, -8.0748E-11,
    18449.4008, -1.7734E-11, -7.6785E-11, 16235.8814, -1.7779E-11, -7.2822E-11,
    14051.7527, -1.7825E-11, -6.8858E-11, 11899.6866, -1.787E-11, -6.4895E-11,
    9782.3549, -1.7915E-11, -6.0932E-11, -2.0615E-11, 48647.756, -4.2651E-10,
    -1.9948E-11, 47240.8586, -4.0783E-10, -1.9282E-11, 45833.9613, -3.8914E-10,
    -1.8615E-11, 44427.0639, -3.7046E-10, -1.7949E-11, 43020.1665, -3.5177E-10,
    -1.7283E-11, 61613.2692, -3.3309E-10, -1.6616E-11, 39252.7963, -3.144E-10,
    -1.595E-11, 36894.9953, -2.9572E-10, -1.5283E-11, 34542.538, -2.7703E-10,
    -1.4617E-11, 32198.0964, -2.5835E-10, -1.395E-11, 29864.3424, -2.3966E-10,
    -1.3284E-11, 27543.9477, -2.2098E-10, -1.2618E-11, 25239.5843, -2.0229E-10,
    -1.1951E-11, 22953.9241, -1.8361E-10, -1.1285E-11, 20689.639, -1.6493E-10,
    -1.0618E-11, 18449.4008, -1.4624E-10, -9.9518E-12, 16235.8814, -1.2756E-10,
    -9.2854E-12, 14051.7527, -1.0887E-10, -8.6189E-12, 11899.6866, -9.0186E-11,
    -7.9525E-12, 9782.3549, -7.1501E-11, -1.4555E-10, -4.268E-10, 48647.756,
    -1.3973E-10, -4.0806E-10, 47240.8586, -1.339E-10, -3.8932E-10, 45833.9613,
    -1.2807E-10, -3.7057E-10, 44427.0639, -1.2225E-10, -3.5183E-10, 43020.1665,
    -1.1642E-10, -3.3309E-10, 61613.2692, -1.1059E-10, -3.1435E-10, 39252.7963,
    -1.0477E-10, -2.956E-10, 36894.9953, -9.8938E-11, -2.7686E-10, 34542.538,
    -9.3111E-11, -2.5812E-10, 32198.0964, -8.7284E-11, -2.3938E-10, 29864.3424,
    -8.1458E-11, -2.2063E-10, 27543.9477, -7.5631E-11, -2.0189E-10, 25239.5843,
    -6.9804E-11, -1.8315E-10, 22953.9241, -6.3977E-11, -1.6441E-10, 20689.639,
    -5.815E-11, -1.4566E-10, 18449.4008, -5.2323E-11, -1.2692E-10, 16235.8814,
    -4.6496E-11, -1.0818E-10, 14051.7527, -4.0669E-11, -8.9437E-11, 11899.6866,
    -3.4843E-11, -7.0694E-11, 9782.3549, 45792.2718, -1.6379E-11, -1.2936E-10,
    44484.3767, -1.6426E-11, -1.256E-10, 43176.4816, -1.6474E-11, -1.2185E-10,
    41868.5865, -1.6521E-11, -1.181E-10, 40560.6914, -1.6569E-11, -1.1435E-10,
    39252.7963, -1.6616E-11, -1.1059E-10, 57944.9012, -1.6664E-11, -1.0684E-10,
    35683.4305, -1.6711E-11, -1.0309E-10, 33424.6318, -1.6758E-11, -9.9332E-11,
    31171.1768, -1.6806E-11, -9.5579E-11, 28925.7374, -1.6853E-11, -9.1826E-11,
    26690.9856, -1.6901E-11, -8.8073E-11, 24469.5932, -1.6948E-11, -8.4319E-11,
    22264.2321, -1.6996E-11, -8.0566E-11, 20077.5741, -1.7043E-11, -7.6813E-11,
    17912.2913, -1.709E-11, -7.306E-11, 15771.0553, -1.7138E-11, -6.9307E-11,
    13656.5382, -1.7185E-11, -6.5553E-11, 11571.4117, -1.7233E-11, -6.18E-11,
    9518.3479, -1.728E-11, -5.8047E-11, -2.0649E-11, 45792.2718, -4.0261E-10,
    -1.9984E-11, 44484.3767, -3.8496E-10, -1.932E-11, 43176.4816, -3.673E-10,
    -1.8656E-11, 41868.5865, -3.4965E-10, -1.7992E-11, 40560.6914, -3.32E-10,
    -1.7328E-11, 39252.7963, -3.1435E-10, -1.6664E-11, 57944.9012, -2.9669E-10,
    -1.5999E-11, 35683.4305, -2.7904E-10, -1.5335E-11, 33424.6318, -2.6139E-10,
    -1.4671E-11, 31171.1768, -2.4374E-10, -1.4007E-11, 28925.7374, -2.2608E-10,
    -1.3343E-11, 26690.9856, -2.0843E-10, -1.2679E-11, 24469.5932, -1.9078E-10,
    -1.2014E-11, 22264.2321, -1.7313E-10, -1.135E-11, 20077.5741, -1.5547E-10,
    -1.0686E-11, 17912.2913, -1.3782E-10, -1.0022E-11, 15771.0553, -1.2017E-10,
    -9.3577E-12, 13656.5382, -1.0252E-10, -8.6936E-12, 11571.4117, -8.4863E-11,
    -8.0294E-12, 9518.3479, -6.721E-11, -1.4054E-10, -4.0295E-10, 45792.2718,
    -1.3492E-10, -3.8524E-10, 44484.3767, -1.2931E-10, -3.6753E-10, 43176.4816,
    -1.2369E-10, -3.4982E-10, 41868.5865, -1.1807E-10, -3.3211E-10, 40560.6914,
    -1.1246E-10, -3.144E-10, 39252.7963, -1.0684E-10, -2.9669E-10, 57944.9012,
    -1.0122E-10, -2.7898E-10, 35683.4305, -9.5605E-11, -2.6127E-10, 33424.6318,
    -8.9989E-11, -2.4356E-10, 31171.1768, -8.4372E-11, -2.2585E-10, 28925.7374,
    -7.8755E-11, -2.0814E-10, 26690.9856, -7.3138E-11, -1.9043E-10, 24469.5932,
    -6.7522E-11, -1.7272E-10, 22264.2321, -6.1905E-11, -1.5501E-10, 20077.5741,
    -5.6288E-11, -1.373E-10, 17912.2913, -5.0672E-11, -1.1959E-10, 15771.0553,
    -4.5055E-11, -1.0188E-10, 13656.5382, -3.9438E-11, -8.4171E-11, 11571.4117,
    -3.3821E-11, -6.6461E-11, 9518.3479, 42952.819, -1.5701E-11, -1.2248E-10,
    41741.2542, -1.5751E-11, -1.1894E-10, 40529.6895, -1.5801E-11, -1.1539E-10,
    39318.1248, -1.585E-11, -1.1185E-10, 38106.56, -1.59E-11, -1.0831E-10,
    36894.9953, -1.595E-11, -1.0477E-10, 35683.4305, -1.5999E-11, -1.0122E-10,
    54471.8658, -1.6049E-11, -9.7679E-11, 32306.7255, -1.6099E-11, -9.4136E-11,
    30144.2571, -1.6148E-11, -9.0593E-11, 27987.1325, -1.6198E-11, -8.705E-11,
    25838.0235, -1.6248E-11, -8.3507E-11, 23699.6021, -1.6298E-11, -7.9964E-11,
    21574.54, -1.6347E-11, -7.6421E-11, 19465.5093, -1.6397E-11, -7.2878E-11,
    17375.1817, -1.6447E-11, -6.9335E-11, 15306.2292, -1.6496E-11, -6.5792E-11,
    13261.3236, -1.6546E-11, -6.2248E-11, 11243.1368, -1.6596E-11, -5.8705E-11,
    9254.3408, -1.6645E-11, -5.5162E-11, -2.0682E-11, 42952.819, -3.7871E-10,
    -2.0021E-11, 41741.2542, -3.6208E-10, -1.9359E-11, 40529.6895, -3.4546E-10,
    -1.8697E-11, 39318.1248, -3.2884E-10, -1.8035E-11, 38106.56, -3.1222E-10,
    -1.7373E-11, 36894.9953, -2.956E-10, -1.6711E-11, 35683.4305, -2.7898E-10,
    -1.6049E-11, 54471.8658, -2.6236E-10, -1.5387E-11, 32306.7255, -2.4574E-10,
    -1.4725E-11, 30144.2571, -2.2912E-10, -1.4063E-11, 27987.1325, -2.125E-10,
    -1.3402E-11, 25838.0235, -1.9588E-10, -1.274E-11, 23699.6021, -1.7926E-10,
    -1.2078E-11, 21574.54, -1.6264E-10, -1.1416E-11, 19465.5093, -1.4602E-10,
    -1.0754E-11, 17375.1817, -1.294E-10, -1.0092E-11, 15306.2292, -1.1278E-10,
    -9.4301E-12, 13261.3236, -9.6159E-11, -8.7682E-12, 11243.1368, -7.9539E-11,
    -8.1063E-12, 9254.3408, -6.2919E-11, -1.3552E-10, -3.7911E-10, 42952.819,
    -1.3012E-10, -3.6243E-10, 41741.2542, -1.2471E-10, -3.4575E-10, 40529.6895,
    -1.1931E-10, -3.2907E-10, 39318.1248, -1.139E-10, -3.124E-10, 38106.56,
    -1.0849E-10, -2.9572E-10, 36894.9953, -1.0309E-10, -2.7904E-10, 35683.4305,
    -9.7679E-11, -2.6236E-10, 54471.8658, -9.2272E-11, -2.4568E-10, 32306.7255,
    -8.6866E-11, -2.2901E-10, 30144.2571, -8.1459E-11, -2.1233E-10, 27987.1325,
    -7.6053E-11, -1.9565E-10, 25838.0235, -7.0646E-11, -1.7897E-10, 23699.6021,
    -6.524E-11, -1.623E-10, 21574.54, -5.9833E-11, -1.4562E-10, 19465.5093,
    -5.4427E-11, -1.2894E-10, 17375.1817, -4.902E-11, -1.1226E-10, 15306.2292,
    -4.3613E-11, -9.5583E-11, 13261.3236, -3.8207E-11, -7.8905E-11, 11243.1368,
    -3.28E-11, -6.2227E-11, 9254.3408, 40132.0693, -1.5023E-11, -1.156E-10,
    39014.1631, -1.5075E-11, -1.1227E-10, 37896.2568, -1.5127E-11, -1.0894E-10,
    36778.3505, -1.5179E-11, -1.056E-10, 35660.4443, -1.5231E-11, -1.0227E-10,
    34542.538, -1.5283E-11, -9.8938E-11, 33424.6318, -1.5335E-11, -9.5605E-11,
    32306.7255, -1.5387E-11, -9.2272E-11, 51188.8193, -1.5439E-11, -8.894E-11,
    29117.3375, -1.5491E-11, -8.5607E-11, 27048.5276, -1.5543E-11, -8.2274E-11,
    24985.0615, -1.5595E-11, -7.8941E-11, 22929.611, -1.5647E-11, -7.5608E-11,
    20884.848, -1.5699E-11, -7.2275E-11, 18853.4444, -1.5751E-11, -6.8942E-11,
    16838.0722, -1.5803E-11, -6.5609E-11, 14841.4031, -1.5855E-11, -6.2276E-11,
    12866.1091, -1.5907E-11, -5.8944E-11, 10914.862, -1.5959E-11, -5.5611E-11,
    8990.3337, -1.6011E-11, -5.2278E-11, -2.0716E-11, 40132.0693, -3.548E-10,
    -2.0057E-11, 39014.1631, -3.3921E-10, -1.9397E-11, 37896.2568, -3.2363E-10,
    -1.8737E-11, 36778.3505, -3.0804E-10, -1.8078E-11, 35660.4443, -2.9245E-10,
    -1.7418E-11, 34542.538, -2.7686E-10, -1.6758E-11, 33424.6318, -2.6127E-10,
    -1.6099E-11, 32306.7255, -2.4568E-10, -1.5439E-11, 51188.8193, -2.301E-10,
    -1.478E-11, 29117.3375, -2.1451E-10, -1.412E-11, 27048.5276, -1.9892E-10,
    -1.346E-11, 24985.0615, -1.8333E-10, -1.2801E-11, 22929.611, -1.6774E-10,
    -1.2141E-11, 20884.848, -1.5216E-10, -1.1481E-11, 18853.4444, -1.3657E-10,
    -1.0822E-11, 16838.0722, -1.2098E-10, -1.0162E-11, 14841.4031, -1.0539E-10,
    -9.5025E-12, 12866.1091, -8.9804E-11, -8.8428E-12, 10914.862, -7.4216E-11,
    -8.1832E-12, 8990.3337, -5.8627E-11, -1.3051E-10, -3.5526E-10, 40132.0693,
    -1.2531E-10, -3.3962E-10, 39014.1631, -1.2012E-10, -3.2397E-10, 37896.2568,
    -1.1492E-10, -3.0833E-10, 36778.3505, -1.0973E-10, -2.9268E-10, 35660.4443,
    -1.0453E-10, -2.7703E-10, 34542.538, -9.9332E-11, -2.6139E-10, 33424.6318,
    -9.4136E-11, -2.4574E-10, 32306.7255, -8.894E-11, -2.301E-10, 51188.8193,
    -8.3743E-11, -2.1445E-10, 29117.3375, -7.8547E-11, -1.9881E-10, 27048.5276,
    -7.335E-11, -1.8316E-10, 24985.0615, -6.8154E-11, -1.6751E-10, 22929.611,
    -6.2958E-11, -1.5187E-10, 20884.848, -5.7761E-11, -1.3622E-10, 18853.4444,
    -5.2565E-11, -1.2058E-10, 16838.0722, -4.7368E-11, -1.0493E-10, 14841.4031,
    -4.2172E-11, -8.9285E-11, 12866.1091, -3.6976E-11, -7.3639E-11, 10914.862,
    -3.1779E-11, -5.7993E-11, 8990.3337, 37332.6947, -1.4346E-11, -1.0872E-10,
    36305.775, -1.44E-11, -1.056E-10, 35278.8554, -1.4454E-11, -1.0248E-10,
    34251.9357, -1.4508E-11, -9.9357E-11, 33225.0161, -1.4563E-11, -9.6234E-11,
    32198.0964, -1.4617E-11, -9.3111E-11, 31171.1768, -1.4671E-11, -8.9989E-11,
    30144.2571, -1.4725E-11, -8.6866E-11, 29117.3375, -1.478E-11, -8.3743E-11,
    48090.4178, -1.4834E-11, -8.062E-11, 26109.9227, -1.4888E-11, -7.7498E-11,
    24132.0994, -1.4942E-11, -7.4375E-11, 22159.6198, -1.4996E-11, -7.1252E-11,
    20195.156, -1.5051E-11, -6.813E-11, 18241.3796, -1.5105E-11, -6.5007E-11,
    16300.9627, -1.5159E-11, -6.1884E-11, 14376.577, -1.5213E-11, -5.8761E-11,
    12470.8945, -1.5268E-11, -5.5639E-11, 10586.5871, -1.5322E-11, -5.2516E-11,
    8726.3266, -1.5376E-11, -4.9393E-11, -2.075E-11, 37332.6947, -3.309E-10,
    -2.0093E-11, 36305.775, -3.1634E-10, -1.9435E-11, 35278.8554, -3.0179E-10,
    -1.8778E-11, 34251.9357, -2.8723E-10, -1.8121E-11, 33225.0161, -2.7267E-10,
    -1.7463E-11, 32198.0964, -2.5812E-10, -1.6806E-11, 31171.1768, -2.4356E-10,
    -1.6148E-11, 30144.2571, -2.2901E-10, -1.5491E-11, 29117.3375, -2.1445E-10,
    -1.4834E-11, 48090.4178, -1.999E-10, -1.4176E-11, 26109.9227, -1.8534E-10,
    -1.3519E-11, 24132.0994, -1.7078E-10, -1.2862E-11, 22159.6198, -1.5623E-10,
    -1.2204E-11, 20195.156, -1.4167E-10, -1.1547E-11, 18241.3796, -1.2712E-10,
    -1.089E-11, 16300.9627, -1.1256E-10, -1.0232E-11, 14376.577, -9.8004E-11,
    -9.5748E-12, 12470.8945, -8.3448E-11, -8.9175E-12, 10586.5871, -6.8892E-11,
    -8.2601E-12, 8726.3266, -5.4336E-11, -1.255E-10, -3.3142E-10, 37332.6947,
    -1.2051E-10, -3.168E-10, 36305.775, -1.1552E-10, -3.0219E-10, 35278.8554,
    -1.1054E-10, -2.8758E-10, 34251.9357, -1.0555E-10, -2.7296E-10, 33225.0161,
    -1.0057E-10, -2.5835E-10, 32198.0964, -9.5579E-11, -2.4374E-10, 31171.1768,
    -9.0593E-11, -2.2912E-10, 30144.2571, -8.5607E-11, -2.1451E-10, 29117.3375,
    -8.062E-11, -1.999E-10, 48090.4178, -7.5634E-11, -1.8528E-10, 26109.9227,
    -7.0648E-11, -1.7067E-10, 24132.0994, -6.5662E-11, -1.5605E-10, 22159.6198,
    -6.0676E-11, -1.4144E-10, 20195.156, -5.5689E-11, -1.2683E-10, 18241.3796,
    -5.0703E-11, -1.1221E-10, 16300.9627, -4.5717E-11, -9.76E-11, 14376.577,
    -4.0731E-11, -8.2987E-11, 12470.8945, -3.5744E-11, -6.8373E-11, 10586.5871,
    -3.0758E-11, -5.376E-11, 8726.3266, 34557.367, -1.3668E-11, -1.0185E-10,
    33618.7621, -1.3724E-11, -9.8935E-11, 32680.1571, -1.3781E-11, -9.6022E-11,
    31741.5522, -1.3837E-11, -9.311E-11, 30802.9473, -1.3894E-11, -9.0197E-11,
    29864.3424, -1.395E-11, -8.7284E-11, 28925.7374, -1.4007E-11, -8.4372E-11,
    27987.1325, -1.4063E-11, -8.1459E-11, 27048.5276, -1.412E-11, -7.8547E-11,
    26109.9227, -1.4176E-11, -7.5634E-11, 45171.3177, -1.4233E-11, -7.2722E-11,
    23279.1373, -1.4289E-11, -6.9809E-11, 21389.6287, -1.4346E-11, -6.6897E-11,
    19505.4639, -1.4402E-11, -6.3984E-11, 17629.3148, -1.4459E-11, -6.1071E-11,
    15763.8531, -1.4515E-11, -5.8159E-11, 13911.7509, -1.4572E-11, -5.5246E-11,
    12075.68, -1.4628E-11, -5.2334E-11, 10258.3122, -1.4685E-11, -4.9421E-11,
    8462.3195, -1.4741E-11, -4.6509E-11, -2.0784E-11, 34557.367, -3.0699E-10,
    -2.0129E-11, 33618.7621, -2.9347E-10, -1.9474E-11, 32680.1571, -2.7995E-10,
    -1.8819E-11, 31741.5522, -2.6642E-10, -1.8163E-11, 30802.9473, -2.529E-10,
    -1.7508E-11, 29864.3424, -2.3938E-10, -1.6853E-11, 28925.7374, -2.2585E-10,
    -1.6198E-11, 27987.1325, -2.1233E-10, -1.5543E-11, 27048.5276, -1.9881E-10,
    -1.4888E-11, 26109.9227, -1.8528E-10, -1.4233E-11, 45171.3177, -1.7176E-10,
    -1.3578E-11, 23279.1373, -1.5823E-10, -1.2923E-11, 21389.6287, -1.4471E-10,
    -1.2268E-11, 19505.4639, -1.3119E-10, -1.1612E-11, 17629.3148, -1.1766E-10,
    -1.0957E-11, 15763.8531, -1.0414E-10, -1.0302E-11, 13911.7509, -9.0616E-11,
    -9.6472E-12, 12075.68, -7.7092E-11, -8.9921E-12, 10258.3122, -6.3569E-11,
    -8.337E-12, 8462.3195, -5.0045E-11, -1.2048E-10, -3.0757E-10, 34557.367,
    -1.1571E-10, -2.9399E-10, 33618.7621, -1.1093E-10, -2.8041E-10, 32680.1571,
    -1.0615E-10, -2.6683E-10, 31741.5522, -1.0138E-10, -2.5325E-10, 30802.9473,
    -9.6602E-11, -2.3966E-10, 29864.3424, -9.1826E-11, -2.2608E-10, 28925.7374,
    -8.705E-11, -2.125E-10, 27987.1325, -8.2274E-11, -1.9892E-10, 27048.5276,
    -7.7498E-11, -1.8534E-10, 26109.9227, -7.2722E-11, -1.7176E-10, 45171.3177,
    -6.7946E-11, -1.5818E-10, 23279.1373, -6.317E-11, -1.446E-10, 21389.6287,
    -5.8393E-11, -1.3101E-10, 19505.4639, -5.3617E-11, -1.1743E-10, 17629.3148,
    -4.8841E-11, -1.0385E-10, 15763.8531, -4.4065E-11, -9.027E-11, 13911.7509,
    -3.9289E-11, -7.6689E-11, 12075.68, -3.4513E-11, -6.3107E-11, 10258.3122,
    -2.9737E-11, -4.9526E-11, 8462.3195, 31808.7581, -1.299E-11, -9.497E-11,
    30955.796, -1.3049E-11, -9.2267E-11, 30102.834, -1.3108E-11, -8.9565E-11,
    29249.8719, -1.3166E-11, -8.6862E-11, 28396.9098, -1.3225E-11, -8.416E-11,
    27543.9477, -1.3284E-11, -8.1458E-11, 26690.9856, -1.3343E-11, -7.8755E-11,
    25838.0235, -1.3402E-11, -7.6053E-11, 24985.0615, -1.346E-11, -7.335E-11,
    24132.0994, -1.3519E-11, -7.0648E-11, 23279.1373, -1.3578E-11, -6.7946E-11,
    42426.1752, -1.3637E-11, -6.5243E-11, 20619.6376, -1.3695E-11, -6.2541E-11,
    18815.7719, -1.3754E-11, -5.9838E-11, 17017.2499, -1.3813E-11, -5.7136E-11,
    15226.7436, -1.3872E-11, -5.4434E-11, 13446.9248, -1.393E-11, -5.1731E-11,
    11680.4654, -1.3989E-11, -4.9029E-11, 9930.0373, -1.4048E-11, -4.6326E-11,
    8198.3124, -1.4107E-11, -4.3624E-11, -2.0818E-11, 31808.7581, -2.8309E-10,
    -2.0165E-11, 30955.796, -2.706E-10, -1.9512E-11, 30102.834, -2.5811E-10,
    -1.8859E-11, 29249.8719, -2.4562E-10, -1.8206E-11, 28396.9098, -2.3313E-10,
    -1.7554E-11, 27543.9477, -2.2063E-10, -1.6901E-11, 26690.9856, -2.0814E-10,
    -1.6248E-11, 25838.0235, -1.9565E-10, -1.5595E-11, 24985.0615, -1.8316E-10,
    -1.4942E-11, 24132.0994, -1.7067E-10, -1.4289E-11, 23279.1373, -1.5818E-10,
    -1.3637E-11, 42426.1752, -1.4569E-10, -1.2984E-11, 20619.6376, -1.3319E-10,
    -1.2331E-11, 18815.7719, -1.207E-10, -1.1678E-11, 17017.2499, -1.0821E-10,
    -1.1025E-11, 15226.7436, -9.5719E-11, -1.0372E-11, 13446.9248, -8.3228E-11,
    -9.7196E-12, 11680.4654, -7.0736E-11, -9.0667E-12, 9930.0373, -5.8245E-11,
    -8.4139E-12, 8198.3124, -4.5754E-11, -1.1547E-10, -2.8373E-10, 31808.7581,
    -1.109E-10, -2.7118E-10, 30955.796, -1.0634E-10, -2.5863E-10, 30102.834,
    -1.0177E-10, -2.4608E-10, 29249.8719, -9.7205E-11, -2.3353E-10, 28396.9098,
    -9.2639E-11, -2.2098E-10, 27543.9477, -8.8073E-11, -2.0843E-10, 26690.9856,
    -8.3507E-11, -1.9588E-10, 25838.0235, -7.8941E-11, -1.8333E-10, 24985.0615,
    -7.4375E-11, -1.7078E-10, 24132.0994, -6.9809E-11, -1.5823E-10, 23279.1373,
    -6.5243E-11, -1.4569E-10, 42426.1752, -6.0677E-11, -1.3314E-10, 20619.6376,
    -5.6111E-11, -1.2059E-10, 18815.7719, -5.1545E-11, -1.0804E-10, 17017.2499,
    -4.698E-11, -9.5489E-11, 15226.7436, -4.2414E-11, -8.294E-11, 13446.9248,
    -3.7848E-11, -7.0391E-11, 11680.4654, -3.3282E-11, -5.7842E-11, 9930.0373,
    -2.8716E-11, -4.5293E-11, 8198.3124, 29089.5399, -1.2312E-11, -8.8092E-11,
    28319.5488, -1.2373E-11, -8.56E-11, 27549.5577, -1.2434E-11, -8.3107E-11,
    26779.5666, -1.2495E-11, -8.0615E-11, 26009.5755, -1.2557E-11, -7.8123E-11,
    25239.5843, -1.2618E-11, -7.5631E-11, 24469.5932, -1.2679E-11, -7.3138E-11,
    23699.6021, -1.274E-11, -7.0646E-11, 22929.611, -1.2801E-11, -6.8154E-11,
    22159.6198, -1.2862E-11, -6.5662E-11, 21389.6287, -1.2923E-11, -6.317E-11,
    20619.6376, -1.2984E-11, -6.0677E-11, 39849.6465, -1.3045E-11, -5.8185E-11,
    18126.0798, -1.3106E-11, -5.5693E-11, 16405.1851, -1.3167E-11, -5.3201E-11,
    14689.6341, -1.3228E-11, -5.0708E-11, 12982.0987, -1.3289E-11, -4.8216E-11,
    11285.2509, -1.335E-11, -4.5724E-11, 9601.7625, -1.3411E-11, -4.3232E-11,
    7934.3053, -1.3472E-11, -4.0739E-11, -2.0852E-11, 29089.5399, -2.5919E-10,
    -2.0201E-11, 28319.5488, -2.4773E-10, -1.955E-11, 27549.5577, -2.3627E-10,
    -1.89E-11, 26779.5666, -2.2481E-10, -1.8249E-11, 26009.5755, -2.1335E-10,
    -1.7599E-11, 25239.5843, -2.0189E-10, -1.6948E-11, 24469.5932, -1.9043E-10,
    -1.6298E-11, 23699.6021, -1.7897E-10, -1.5647E-11, 22929.611, -1.6751E-10,
    -1.4996E-11, 22159.6198, -1.5605E-10, -1.4346E-11, 21389.6287, -1.446E-10,
    -1.3695E-11, 20619.6376, -1.3314E-10, -1.3045E-11, 39849.6465, -1.2168E-10,
    -1.2394E-11, 18126.0798, -1.1022E-10, -1.1744E-11, 16405.1851, -9.8758E-11,
    -1.1093E-11, 14689.6341, -8.7299E-11, -1.0443E-11, 12982.0987, -7.584E-11,
    -9.7919E-12, 11285.2509, -6.4381E-11, -9.1414E-12, 9601.7625, -5.2922E-11,
    -8.4908E-12, 7934.3053, -4.1462E-11, -1.1045E-10, -2.5988E-10, 29089.5399,
    -1.061E-10, -2.4836E-10, 28319.5488, -1.0174E-10, -2.3685E-10, 27549.5577,
    -9.7387E-11, -2.2533E-10, 26779.5666, -9.3031E-11, -2.1381E-10, 26009.5755,
    -8.8675E-11, -2.0229E-10, 25239.5843, -8.4319E-11, -1.9078E-10, 24469.5932,
    -7.9964E-11, -1.7926E-10, 23699.6021, -7.5608E-11, -1.6774E-10, 22929.611,
    -7.1252E-11, -1.5623E-10, 22159.6198, -6.6897E-11, -1.4471E-10, 21389.6287,
    -6.2541E-11, -1.3319E-10, 20619.6376, -5.8185E-11, -1.2168E-10, 39849.6465,
    -5.3829E-11, -1.1016E-10, 18126.0798, -4.9474E-11, -9.8643E-11, 16405.1851,
    -4.5118E-11, -8.7126E-11, 14689.6341, -4.0762E-11, -7.5609E-11, 12982.0987,
    -3.6406E-11, -6.4093E-11, 11285.2509, -3.2051E-11, -5.2576E-11, 9601.7625,
    -2.7695E-11, -4.1059E-11, 7934.3053, 26402.3844, -1.1635E-11, -8.1214E-11,
    25712.6923, -1.1698E-11, -7.8932E-11, 25023.0003, -1.1761E-11, -7.665E-11,
    24333.3082, -1.1825E-11, -7.4368E-11, 23643.6162, -1.1888E-11, -7.2086E-11,
    22953.9241, -1.1951E-11, -6.9804E-11, 22264.2321, -1.2014E-11, -6.7522E-11,
    21574.54, -1.2078E-11, -6.524E-11, 20884.848, -1.2141E-11, -6.2958E-11,
    20195.156, -1.2204E-11, -6.0676E-11, 19505.4639, -1.2268E-11, -5.8393E-11,
    18815.7719, -1.2331E-11, -5.6111E-11, 18126.0798, -1.2394E-11, -5.3829E-11,
    37436.3878, -1.2457E-11, -5.1547E-11, 15793.1202, -1.2521E-11, -4.9265E-11,
    14152.5245, -1.2584E-11, -4.6983E-11, 12517.2726, -1.2647E-11, -4.4701E-11,
    10890.0363, -1.2711E-11, -4.2419E-11, 9273.4876, -1.2774E-11, -4.0137E-11,
    7670.2982, -1.2837E-11, -3.7855E-11, -2.0885E-11, 26402.3844, -2.3528E-10,
    -2.0237E-11, 25712.6923, -2.2486E-10, -1.9589E-11, 25023.0003, -2.1443E-10,
    -1.894E-11, 24333.3082, -2.04E-10, -1.8292E-11, 23643.6162, -1.9358E-10,
    -1.7644E-11, 22953.9241, -1.8315E-10, -1.6996E-11, 22264.2321, -1.7272E-10,
    -1.6347E-11, 21574.54, -1.623E-10, -1.5699E-11, 20884.848, -1.5187E-10,
    -1.5051E-11, 20195.156, -1.4144E-10, -1.4402E-11, 19505.4639, -1.3101E-10,
    -1.3754E-11, 18815.7719, -1.2059E-10, -1.3106E-11, 18126.0798, -1.1016E-10,
    -1.2457E-11, 37436.3878, -9.9733E-11, -1.1809E-11, 15793.1202, -8.9306E-11,
    -1.1161E-11, 14152.5245, -7.8879E-11, -1.0513E-11, 12517.2726, -6.8452E-11,
    -9.8643E-12, 10890.0363, -5.8025E-11, -9.216E-12, 9273.4876, -4.7598E-11,
    -8.5677E-12, 7670.2982, -3.7171E-11, -1.0544E-10, -2.3603E-10, 26402.3844,
    -1.0129E-10, -2.2555E-10, 25712.6923, -9.7149E-11, -2.1506E-10, 25023.0003,
    -9.3003E-11, -2.0458E-10, 24333.3082, -8.8857E-11, -1.9409E-10, 23643.6162,
    -8.4712E-11, -1.8361E-10, 22953.9241, -8.0566E-11, -1.7313E-10, 22264.2321,
    -7.6421E-11, -1.6264E-10, 21574.54, -7.2275E-11, -1.5216E-10, 20884.848,
    -6.813E-11, -1.4167E-10, 20195.156, -6.3984E-11, -1.3119E-10, 19505.4639,
    -5.9838E-11, -1.207E-10, 18815.7719, -5.5693E-11, -1.1022E-10, 18126.0798,
    -5.1547E-11, -9.9733E-11, 37436.3878, -4.7402E-11, -8.9248E-11, 15793.1202,
    -4.3256E-11, -7.8764E-11, 14152.5245, -3.911E-11, -6.8279E-11, 12517.2726,
    -3.4965E-11, -5.7794E-11, 10890.0363, -3.0819E-11, -4.731E-11, 9273.4876,
    -2.6674E-11, -3.6825E-11, 7670.2982, 23749.9632, -1.0957E-11, -7.4336E-11,
    23137.8984, -1.1022E-11, -7.2265E-11, 22525.8335, -1.1088E-11, -7.0193E-11,
    21913.7687, -1.1154E-11, -6.8121E-11, 21301.7038, -1.1219E-11, -6.6049E-11,
    20689.639, -1.1285E-11, -6.3977E-11, 20077.5741, -1.135E-11, -6.1905E-11,
    19465.5093, -1.1416E-11, -5.9833E-11, 18853.4444, -1.1481E-11, -5.7761E-11,
    18241.3796, -1.1547E-11, -5.5689E-11, 17629.3148, -1.1612E-11, -5.3617E-11,
    17017.2499, -1.1678E-11, -5.1545E-11, 16405.1851, -1.1744E-11, -4.9474E-11,
    15793.1202, -1.1809E-11, -4.7402E-11, 35181.0554, -1.1875E-11, -4.533E-11,
    13615.415, -1.194E-11, -4.3258E-11, 12052.4465, -1.2006E-11, -4.1186E-11,
    10494.8218, -1.2071E-11, -3.9114E-11, 8945.2127, -1.2137E-11, -3.7042E-11,
    7406.2911, -1.2203E-11, -3.497E-11, -2.0919E-11, 23749.9632, -2.1138E-10,
    -2.0273E-11, 23137.8984, -2.0199E-10, -1.9627E-11, 22525.8335, -1.9259E-10,
    -1.8981E-11, 21913.7687, -1.832E-10, -1.8335E-11, 21301.7038, -1.738E-10,
    -1.7689E-11, 20689.639, -1.6441E-10, -1.7043E-11, 20077.5741, -1.5501E-10,
    -1.6397E-11, 19465.5093, -1.4562E-10, -1.5751E-11, 18853.4444, -1.3622E-10,
    -1.5105E-11, 18241.3796, -1.2683E-10, -1.4459E-11, 17629.3148, -1.1743E-10,
    -1.3813E-11, 17017.2499, -1.0804E-10, -1.3167E-11, 16405.1851, -9.8643E-11,
    -1.2521E-11, 15793.1202, -8.9248E-11, -1.1875E-11, 35181.0554, -7.9854E-11,
    -1.1229E-11, 13615.415, -7.0459E-11, -1.0583E-11, 12052.4465, -6.1064E-11,
    -9.9367E-12, 10494.8218, -5.1669E-11, -9.2906E-12, 8945.2127, -4.2275E-11,
    -8.6446E-12, 7406.2911, -3.288E-11, -1.0043E-10, -2.1219E-10, 23749.9632,
    -9.649E-11, -2.0273E-10, 23137.8984, -9.2555E-11, -1.9328E-10, 22525.8335,
    -8.8619E-11, -1.8383E-10, 21913.7687, -8.4684E-11, -1.7438E-10, 21301.7038,
    -8.0748E-11, -1.6493E-10, 20689.639, -7.6813E-11, -1.5547E-10, 20077.5741,
    -7.2878E-11, -1.4602E-10, 19465.5093, -6.8942E-11, -1.3657E-10, 18853.4444,
    -6.5007E-11, -1.2712E-10, 18241.3796, -6.1071E-11, -1.1766E-10, 17629.3148,
    -5.7136E-11, -1.0821E-10, 17017.2499, -5.3201E-11, -9.8758E-11, 16405.1851,
    -4.9265E-11, -8.9306E-11, 15793.1202, -4.533E-11, -7.9854E-11, 35181.0554,
    -4.1394E-11, -7.0401E-11, 13615.415, -3.7459E-11, -6.0949E-11, 12052.4465,
    -3.3523E-11, -5.1496E-11, 10494.8218, -2.9588E-11, -4.2044E-11, 8945.2127,
    -2.5653E-11, -3.2592E-11, 7406.2911, 21134.9484, -1.0279E-11, -6.7459E-11,
    20597.8389, -1.0347E-11, -6.5597E-11, 20060.7294, -1.0415E-11, -6.3735E-11,
    19523.6198, -1.0483E-11, -6.1874E-11, 18986.5103, -1.055E-11, -6.0012E-11,
    18449.4008, -1.0618E-11, -5.815E-11, 17912.2913, -1.0686E-11, -5.6288E-11,
    17375.1817, -1.0754E-11, -5.4427E-11, 16838.0722, -1.0822E-11, -5.2565E-11,
    16300.9627, -1.089E-11, -5.0703E-11, 15763.8531, -1.0957E-11, -4.8841E-11,
    15226.7436, -1.1025E-11, -4.698E-11, 14689.6341, -1.1093E-11, -4.5118E-11,
    14152.5245, -1.1161E-11, -4.3256E-11, 13615.415, -1.1229E-11, -4.1394E-11,
    33078.3055, -1.1297E-11, -3.9533E-11, 11587.6204, -1.1364E-11, -3.7671E-11,
    10099.6072, -1.1432E-11, -3.5809E-11, 8616.9378, -1.15E-11, -3.3947E-11,
    7142.2841, -1.1568E-11, -3.2086E-11, -2.0953E-11, 21134.9484, -1.8748E-10,
    -2.0309E-11, 20597.8389, -1.7911E-10, -1.9665E-11, 20060.7294, -1.7075E-10,
    -1.9022E-11, 19523.6198, -1.6239E-10, -1.8378E-11, 18986.5103, -1.5403E-10,
    -1.7734E-11, 18449.4008, -1.4566E-10, -1.709E-11, 17912.2913, -1.373E-10,
    -1.6447E-11, 17375.1817, -1.2894E-10, -1.5803E-11, 16838.0722, -1.2058E-10,
    -1.5159E-11, 16300.9627, -1.1221E-10, -1.4515E-11, 15763.8531, -1.0385E-10,
    -1.3872E-11, 15226.7436, -9.5489E-11, -1.3228E-11, 14689.6341, -8.7126E-11,
    -1.2584E-11, 14152.5245, -7.8764E-11, -1.194E-11, 13615.415, -7.0401E-11,
    -1.1297E-11, 33078.3055, -6.2039E-11, -1.0653E-11, 11587.6204, -5.3676E-11,
    -1.0009E-11, 10099.6072, -4.5314E-11, -9.3653E-12, 8616.9378, -3.6951E-11,
    -8.7215E-12, 7142.2841, -2.8588E-11, -9.5411E-11, -1.8834E-10, 21134.9484,
    -9.1686E-11, -1.7992E-10, 20597.8389, -8.7961E-11, -1.715E-10, 20060.7294,
    -8.4236E-11, -1.6308E-10, 19523.6198, -8.051E-11, -1.5466E-10, 18986.5103,
    -7.6785E-11, -1.4624E-10, 18449.4008, -7.306E-11, -1.3782E-10, 17912.2913,
    -6.9335E-11, -1.294E-10, 17375.1817, -6.5609E-11, -1.2098E-10, 16838.0722,
    -6.1884E-11, -1.1256E-10, 16300.9627, -5.8159E-11, -1.0414E-10, 15763.8531,
    -5.4434E-11, -9.5719E-11, 15226.7436, -5.0708E-11, -8.7299E-11, 14689.6341,
    -4.6983E-11, -7.8879E-11, 14152.5245, -4.3258E-11, -7.0459E-11, 13615.415,
    -3.9533E-11, -6.2039E-11, 33078.3055, -3.5807E-11, -5.3618E-11, 11587.6204,
    -3.2082E-11, -4.5198E-11, 10099.6072, -2.8357E-11, -3.6778E-11, 8616.9378,
    -2.4632E-11, -2.8358E-11, 7142.2841, 18560.0119, -9.6013E-12, -6.0581E-11,
    18095.1858, -9.6714E-12, -5.893E-11, 17630.3597, -9.7415E-12, -5.7278E-11,
    17165.5336, -9.8116E-12, -5.5626E-11, 16700.7075, -9.8817E-12, -5.3975E-11,
    16235.8814, -9.9518E-12, -5.2323E-11, 15771.0553, -1.0022E-11, -5.0672E-11,
    15306.2292, -1.0092E-11, -4.902E-11, 14841.4031, -1.0162E-11, -4.7368E-11,
    14376.577, -1.0232E-11, -4.5717E-11, 13911.7509, -1.0302E-11, -4.4065E-11,
    13446.9248, -1.0372E-11, -4.2414E-11, 12982.0987, -1.0443E-11, -4.0762E-11,
    12517.2726, -1.0513E-11, -3.911E-11, 12052.4465, -1.0583E-11, -3.7459E-11,
    11587.6204, -1.0653E-11, -3.5807E-11, 31122.7943, -1.0723E-11, -3.4156E-11,
    9704.3927, -1.0793E-11, -3.2504E-11, 8288.6629, -1.0863E-11, -3.0853E-11,
    6878.277, -1.0933E-11, -2.9201E-11, -2.0987E-11, 18560.0119, -1.6357E-10,
    -2.0345E-11, 18095.1858, -1.5624E-10, -1.9704E-11, 17630.3597, -1.4891E-10,
    -1.9062E-11, 17165.5336, -1.4158E-10, -1.8421E-11, 16700.7075, -1.3425E-10,
    -1.7779E-11, 16235.8814, -1.2692E-10, -1.7138E-11, 15771.0553, -1.1959E-10,
    -1.6496E-11, 15306.2292, -1.1226E-10, -1.5855E-11, 14841.4031, -1.0493E-10,
    -1.5213E-11, 14376.577, -9.76E-11, -1.4572E-11, 13911.7509, -9.027E-11,
    -1.393E-11, 13446.9248, -8.294E-11, -1.3289E-11, 12982.0987, -7.5609E-11,
    -1.2647E-11, 12517.2726, -6.8279E-11, -1.2006E-11, 12052.4465, -6.0949E-11,
    -1.1364E-11, 11587.6204, -5.3618E-11, -1.0723E-11, 31122.7943, -4.6288E-11,
    -1.0081E-11, 9704.3927, -3.8958E-11, -9.4399E-12, 8288.6629, -3.1627E-11,
    -8.7984E-12, 6878.277, -2.4297E-11, -9.0397E-11, -1.645E-10, 18560.0119,
    -8.6882E-11, -1.5711E-10, 18095.1858, -8.3367E-11, -1.4972E-10, 17630.3597,
    -7.9852E-11, -1.4233E-10, 17165.5336, -7.6337E-11, -1.3494E-10, 16700.7075,
    -7.2822E-11, -1.2756E-10, 16235.8814, -6.9307E-11, -1.2017E-10, 15771.0553,
    -6.5792E-11, -1.1278E-10, 15306.2292, -6.2276E-11, -1.0539E-10, 14841.4031,
    -5.8761E-11, -9.8004E-11, 14376.577, -5.5246E-11, -9.0616E-11, 13911.7509,
    -5.1731E-11, -8.3228E-11, 13446.9248, -4.8216E-11, -7.584E-11, 12982.0987,
    -4.4701E-11, -6.8452E-11, 12517.2726, -4.1186E-11, -6.1064E-11, 12052.4465,
    -3.7671E-11, -5.3676E-11, 11587.6204, -3.4156E-11, -4.6288E-11, 31122.7943,
    -3.0641E-11, -3.89E-11, 9704.3927, -2.7126E-11, -3.1512E-11, 8288.6629,
    -2.361E-11, -2.4124E-11, 6878.277, 16027.8254, -8.9235E-12, -5.3703E-11,
    15632.6109, -8.9959E-12, -5.2262E-11, 15237.3963, -9.0683E-12, -5.0821E-11,
    14842.1818, -9.1406E-12, -4.9379E-11, 14446.9673, -9.213E-12, -4.7938E-11,
    14051.7527, -9.2854E-12, -4.6496E-11, 13656.5382, -9.3577E-12, -4.5055E-11,
    13261.3236, -9.4301E-12, -4.3613E-11, 12866.1091, -9.5025E-12, -4.2172E-11,
    12470.8945, -9.5748E-12, -4.0731E-11, 12075.68, -9.6472E-12, -3.9289E-11,
    11680.4654, -9.7196E-12, -3.7848E-11, 11285.2509, -9.7919E-12, -3.6406E-11,
    10890.0363, -9.8643E-12, -3.4965E-11, 10494.8218, -9.9367E-12, -3.3523E-11,
    10099.6072, -1.0009E-11, -3.2082E-11, 9704.3927, -1.0081E-11, -3.0641E-11,
    29309.1781, -1.0154E-11, -2.9199E-11, 7960.3881, -1.0226E-11, -2.7758E-11,
    6614.2699, -1.0299E-11, -2.6316E-11, -2.1021E-11, 16027.8254, -1.3967E-10,
    -2.0381E-11, 15632.6109, -1.3337E-10, -1.9742E-11, 15237.3963, -1.2707E-10,
    -1.9103E-11, 14842.1818, -1.2078E-10, -1.8464E-11, 14446.9673, -1.1448E-10,
    -1.7825E-11, 14051.7527, -1.0818E-10, -1.7185E-11, 13656.5382, -1.0188E-10,
    -1.6546E-11, 13261.3236, -9.5583E-11, -1.5907E-11, 12866.1091, -8.9285E-11,
    -1.5268E-11, 12470.8945, -8.2987E-11, -1.4628E-11, 12075.68, -7.6689E-11,
    -1.3989E-11, 11680.4654, -7.0391E-11, -1.335E-11, 11285.2509, -6.4093E-11,
    -1.2711E-11, 10890.0363, -5.7794E-11, -1.2071E-11, 10494.8218, -5.1496E-11,
    -1.1432E-11, 10099.6072, -4.5198E-11, -1.0793E-11, 9704.3927, -3.89E-11,
    -1.0154E-11, 29309.1781, -3.2602E-11, -9.5145E-12, 7960.3881, -2.6304E-11,
    -8.8753E-12, 6614.2699, -2.0006E-11, -8.5383E-11, -1.4065E-10, 16027.8254,
    -8.2078E-11, -1.3429E-10, 15632.6109, -7.8773E-11, -1.2794E-10, 15237.3963,
    -7.5468E-11, -1.2158E-10, 14842.1818, -7.2163E-11, -1.1523E-10, 14446.9673,
    -6.8858E-11, -1.0887E-10, 14051.7527, -6.5553E-11, -1.0252E-10, 13656.5382,
    -6.2248E-11, -9.6159E-11, 13261.3236, -5.8944E-11, -8.9804E-11, 12866.1091,
    -5.5639E-11, -8.3448E-11, 12470.8945, -5.2334E-11, -7.7092E-11, 12075.68,
    -4.9029E-11, -7.0736E-11, 11680.4654, -4.5724E-11, -6.4381E-11, 11285.2509,
    -4.2419E-11, -5.8025E-11, 10890.0363, -3.9114E-11, -5.1669E-11, 10494.8218,
    -3.5809E-11, -4.5314E-11, 10099.6072, -3.2504E-11, -3.8958E-11, 9704.3927,
    -2.9199E-11, -3.2602E-11, 29309.1781, -2.5894E-11, -2.6246E-11, 7960.3881,
    -2.2589E-11, -1.9891E-11, 6614.2699, 13541.061, -8.2457E-12, -4.6826E-11,
    13212.7861, -8.3204E-12, -4.5594E-11, 12884.5112, -8.395E-12, -4.4363E-11,
    12556.2363, -8.4697E-12, -4.3132E-11, 12227.9615, -8.5443E-12, -4.1901E-11,
    11899.6866, -8.6189E-12, -4.0669E-11, 11571.4117, -8.6936E-12, -3.9438E-11,
    11243.1368, -8.7682E-12, -3.8207E-11, 10914.862, -8.8428E-12, -3.6976E-11,
    10586.5871, -8.9175E-12, -3.5744E-11, 10258.3122, -8.9921E-12, -3.4513E-11,
    9930.0373, -9.0667E-12, -3.3282E-11, 9601.7625, -9.1414E-12, -3.2051E-11,
    9273.4876, -9.216E-12, -3.0819E-11, 8945.2127, -9.2906E-12, -2.9588E-11,
    8616.9378, -9.3653E-12, -2.8357E-11, 8288.6629, -9.4399E-12, -2.7126E-11,
    7960.3881, -9.5145E-12, -2.5894E-11, 27632.1132, -9.5892E-12, -2.4663E-11,
    6350.2628, -9.6638E-12, -2.3432E-11, -2.1054E-11, 13541.061, -1.1577E-10,
    -2.0418E-11, 13212.7861, -1.105E-10, -1.9781E-11, 12884.5112, -1.0523E-10,
    -1.9144E-11, 12556.2363, -9.9969E-11, -1.8507E-11, 12227.9615, -9.4703E-11,
    -1.787E-11, 11899.6866, -8.9437E-11, -1.7233E-11, 11571.4117, -8.4171E-11,
    -1.6596E-11, 11243.1368, -7.8905E-11, -1.5959E-11, 10914.862, -7.3639E-11,
    -1.5322E-11, 10586.5871, -6.8373E-11, -1.4685E-11, 10258.3122, -6.3107E-11,
    -1.4048E-11, 9930.0373, -5.7842E-11, -1.3411E-11, 9601.7625, -5.2576E-11,
    -1.2774E-11, 9273.4876, -4.731E-11, -1.2137E-11, 8945.2127, -4.2044E-11,
    -1.15E-11, 8616.9378, -3.6778E-11, -1.0863E-11, 8288.6629, -3.1512E-11,
    -1.0226E-11, 7960.3881, -2.6246E-11, -9.5892E-12, 27632.1132, -2.098E-11,
    -8.9522E-12, 6350.2628, -1.5715E-11, -8.0369E-11, -1.168E-10, 13541.061,
    -7.7274E-11, -1.1148E-10, 13212.7861, -7.4179E-11, -1.0616E-10, 12884.5112,
    -7.1085E-11, -1.0083E-10, 12556.2363, -6.799E-11, -9.551E-11, 12227.9615,
    -6.4895E-11, -9.0186E-11, 11899.6866, -6.18E-11, -8.4863E-11, 11571.4117,
    -5.8705E-11, -7.9539E-11, 11243.1368, -5.5611E-11, -7.4216E-11, 10914.862,
    -5.2516E-11, -6.8892E-11, 10586.5871, -4.9421E-11, -6.3569E-11, 10258.3122,
    -4.6326E-11, -5.8245E-11, 9930.0373, -4.3232E-11, -5.2922E-11, 9601.7625,
    -4.0137E-11, -4.7598E-11, 9273.4876, -3.7042E-11, -4.2275E-11, 8945.2127,
    -3.3947E-11, -3.6951E-11, 8616.9378, -3.0853E-11, -3.1627E-11, 8288.6629,
    -2.7758E-11, -2.6304E-11, 7960.3881, -2.4663E-11, -2.098E-11, 27632.1132,
    -2.1568E-11, -1.5657E-11, 6350.2628, 11102.3904, -7.568E-12, -3.9948E-11,
    10838.3833, -7.6449E-12, -3.8927E-11, 10574.3762, -7.7218E-12, -3.7906E-11,
    10310.3691, -7.7987E-12, -3.6885E-11, 10046.362, -7.8756E-12, -3.5864E-11,
    9782.3549, -7.9525E-12, -3.4843E-11, 9518.3479, -8.0294E-12, -3.3821E-11,
    9254.3408, -8.1063E-12, -3.28E-11, 8990.3337, -8.1832E-12, -3.1779E-11,
    8726.3266, -8.2601E-12, -3.0758E-11, 8462.3195, -8.337E-12, -2.9737E-11,
    8198.3124, -8.4139E-12, -2.8716E-11, 7934.3053, -8.4908E-12, -2.7695E-11,
    7670.2982, -8.5677E-12, -2.6674E-11, 7406.2911, -8.6446E-12, -2.5653E-11,
    7142.2841, -8.7215E-12, -2.4632E-11, 6878.277, -8.7984E-12, -2.361E-11,
    6614.2699, -8.8753E-12, -2.2589E-11, 6350.2628, -8.9522E-12, -2.1568E-11,
    26086.2557, -9.0291E-12, -2.0547E-11, -2.1088E-11, 11102.3904, -9.1863E-11,
    -2.0454E-11, 10838.3833, -8.7629E-11, -1.9819E-11, 10574.3762, -8.3395E-11,
    -1.9184E-11, 10310.3691, -7.9162E-11, -1.855E-11, 10046.362, -7.4928E-11,
    -1.7915E-11, 9782.3549, -7.0694E-11, -1.728E-11, 9518.3479, -6.6461E-11,
    -1.6645E-11, 9254.3408, -6.2227E-11, -1.6011E-11, 8990.3337, -5.7993E-11,
    -1.5376E-11, 8726.3266, -5.376E-11, -1.4741E-11, 8462.3195, -4.9526E-11,
    -1.4107E-11, 8198.3124, -4.5293E-11, -1.3472E-11, 7934.3053, -4.1059E-11,
    -1.2837E-11, 7670.2982, -3.6825E-11, -1.2203E-11, 7406.2911, -3.2592E-11,
    -1.1568E-11, 7142.2841, -2.8358E-11, -1.0933E-11, 6878.277, -2.4124E-11,
    -1.0299E-11, 6614.2699, -1.9891E-11, -9.6638E-12, 6350.2628, -1.5657E-11,
    -9.0291E-12, 26086.2557, -1.1423E-11, -7.5355E-11, -9.2958E-11, 11102.3904,
    -7.247E-11, -8.8667E-11, 10838.3833, -6.9585E-11, -8.4375E-11, 10574.3762,
    -6.6701E-11, -8.0084E-11, 10310.3691, -6.3816E-11, -7.5793E-11, 10046.362,
    -6.0932E-11, -7.1501E-11, 9782.3549, -5.8047E-11, -6.721E-11, 9518.3479,
    -5.5162E-11, -6.2919E-11, 9254.3408, -5.2278E-11, -5.8627E-11, 8990.3337,
    -4.9393E-11, -5.4336E-11, 8726.3266, -4.6509E-11, -5.0045E-11, 8462.3195,
    -4.3624E-11, -4.5754E-11, 8198.3124, -4.0739E-11, -4.1462E-11, 7934.3053,
    -3.7855E-11, -3.7171E-11, 7670.2982, -3.497E-11, -3.288E-11, 7406.2911,
    -3.2086E-11, -2.8588E-11, 7142.2841, -2.9201E-11, -2.4297E-11, 6878.277,
    -2.6316E-11, -2.0006E-11, 6614.2699, -2.3432E-11, -1.5715E-11, 6350.2628,
    -2.0547E-11, -1.1423E-11, 26086.2557 };

  static const real_T b_a[3600] = { -83072.1302, 2.0446E-11, 1.7062E-10,
    -60176.5678, 2.0479E-11, 1.6561E-10, -57283.6773, 2.0513E-11, 1.606E-10,
    -54396.1306, 2.0547E-11, 1.5558E-10, -51516.5995, 2.0581E-11, 1.5057E-10,
    -48647.756, 2.0615E-11, 1.4555E-10, -45792.2718, 2.0649E-11, 1.4054E-10,
    -42952.819, 2.0682E-11, 1.3552E-10, -40132.0693, 2.0716E-11, 1.3051E-10,
    -37332.6947, 2.075E-11, 1.255E-10, -34557.367, 2.0784E-11, 1.2048E-10,
    -31808.7581, 2.0818E-11, 1.1547E-10, -29089.5399, 2.0852E-11, 1.1045E-10,
    -26402.3844, 2.0885E-11, 1.0544E-10, -23749.9632, 2.0919E-11, 1.0043E-10,
    -21134.9484, 2.0953E-11, 9.5411E-11, -18560.0119, 2.0987E-11, 9.0397E-11,
    -16027.8254, 2.1021E-11, 8.5383E-11, -13541.061, 2.1054E-11, 8.0369E-11,
    -11102.3904, 2.1088E-11, 7.5355E-11, 2.0446E-11, -83072.1302, 5.4603E-10,
    1.9768E-11, -60176.5678, 5.2218E-10, 1.909E-11, -57283.6773, 4.9834E-10,
    1.8412E-11, -54396.1306, 4.7449E-10, 1.7735E-11, -51516.5995, 4.5065E-10,
    1.7057E-11, -48647.756, 4.268E-10, 1.6379E-11, -45792.2718, 4.0295E-10,
    1.5701E-11, -42952.819, 3.7911E-10, 1.5023E-11, -40132.0693, 3.5526E-10,
    1.4346E-11, -37332.6947, 3.3142E-10, 1.3668E-11, -34557.367, 3.0757E-10,
    1.299E-11, -31808.7581, 2.8373E-10, 1.2312E-11, -29089.5399, 2.5988E-10,
    1.1635E-11, -26402.3844, 2.3603E-10, 1.0957E-11, -23749.9632, 2.1219E-10,
    1.0279E-11, -21134.9484, 1.8834E-10, 9.6013E-12, -18560.0119, 1.645E-10,
    8.9235E-12, -16027.8254, 1.4065E-10, 8.2457E-12, -13541.061, 1.168E-10,
    7.568E-12, -11102.3904, 9.2958E-11, 1.7062E-10, 5.4603E-10, -83072.1302,
    1.6375E-10, 5.2213E-10, -60176.5678, 1.5687E-10, 4.9822E-10, -57283.6773,
    1.4999E-10, 4.7432E-10, -54396.1306, 1.4311E-10, 4.5042E-10, -51516.5995,
    1.3624E-10, 4.2651E-10, -48647.756, 1.2936E-10, 4.0261E-10, -45792.2718,
    1.2248E-10, 3.7871E-10, -42952.819, 1.156E-10, 3.548E-10, -40132.0693,
    1.0872E-10, 3.309E-10, -37332.6947, 1.0185E-10, 3.0699E-10, -34557.367,
    9.497E-11, 2.8309E-10, -31808.7581, 8.8092E-11, 2.5919E-10, -29089.5399,
    8.1214E-11, 2.3528E-10, -26402.3844, 7.4336E-11, 2.1138E-10, -23749.9632,
    6.7459E-11, 1.8748E-10, -21134.9484, 6.0581E-11, 1.6357E-10, -18560.0119,
    5.3703E-11, 1.3967E-10, -16027.8254, 4.6826E-11, 1.1577E-10, -13541.061,
    3.9948E-11, 9.1863E-11, -11102.3904, -60176.5678, 1.9768E-11, 1.6375E-10,
    -78346.9426, 1.9804E-11, 1.5894E-10, -55563.7419, 1.984E-11, 1.5414E-10,
    -52783.2131, 1.9876E-11, 1.4933E-10, -50008.028, 1.9912E-11, 1.4453E-10,
    -47240.8586, 1.9948E-11, 1.3973E-10, -44484.3767, 1.9984E-11, 1.3492E-10,
    -41741.2542, 2.0021E-11, 1.3012E-10, -39014.1631, 2.0057E-11, 1.2531E-10,
    -36305.775, 2.0093E-11, 1.2051E-10, -33618.7621, 2.0129E-11, 1.1571E-10,
    -30955.796, 2.0165E-11, 1.109E-10, -28319.5488, 2.0201E-11, 1.061E-10,
    -25712.6923, 2.0237E-11, 1.0129E-10, -23137.8984, 2.0273E-11, 9.649E-11,
    -20597.8389, 2.0309E-11, 9.1686E-11, -18095.1858, 2.0345E-11, 8.6882E-11,
    -15632.6109, 2.0381E-11, 8.2078E-11, -13212.7861, 2.0418E-11, 7.7274E-11,
    -10838.3833, 2.0454E-11, 7.247E-11, 2.0479E-11, -60176.5678, 5.2213E-10,
    1.9804E-11, -78346.9426, 4.9931E-10, 1.9128E-11, -55563.7419, 4.765E-10,
    1.8453E-11, -52783.2131, 4.5369E-10, 1.7777E-11, -50008.028, 4.3087E-10,
    1.7102E-11, -47240.8586, 4.0806E-10, 1.6426E-11, -44484.3767, 3.8524E-10,
    1.5751E-11, -41741.2542, 3.6243E-10, 1.5075E-11, -39014.1631, 3.3962E-10,
    1.44E-11, -36305.775, 3.168E-10, 1.3724E-11, -33618.7621, 2.9399E-10,
    1.3049E-11, -30955.796, 2.7118E-10, 1.2373E-11, -28319.5488, 2.4836E-10,
    1.1698E-11, -25712.6923, 2.2555E-10, 1.1022E-11, -23137.8984, 2.0273E-10,
    1.0347E-11, -20597.8389, 1.7992E-10, 9.6714E-12, -18095.1858, 1.5711E-10,
    8.9959E-12, -15632.6109, 1.3429E-10, 8.3204E-12, -13212.7861, 1.1148E-10,
    7.6449E-12, -10838.3833, 8.8667E-11, 1.6561E-10, 5.2218E-10, -60176.5678,
    1.5894E-10, 4.9931E-10, -78346.9426, 1.5227E-10, 4.7644E-10, -55563.7419,
    1.4561E-10, 4.5357E-10, -52783.2131, 1.3894E-10, 4.307E-10, -50008.028,
    1.3227E-10, 4.0783E-10, -47240.8586, 1.256E-10, 3.8496E-10, -44484.3767,
    1.1894E-10, 3.6208E-10, -41741.2542, 1.1227E-10, 3.3921E-10, -39014.1631,
    1.056E-10, 3.1634E-10, -36305.775, 9.8935E-11, 2.9347E-10, -33618.7621,
    9.2267E-11, 2.706E-10, -30955.796, 8.56E-11, 2.4773E-10, -28319.5488,
    7.8932E-11, 2.2486E-10, -25712.6923, 7.2265E-11, 2.0199E-10, -23137.8984,
    6.5597E-11, 1.7911E-10, -20597.8389, 5.893E-11, 1.5624E-10, -18095.1858,
    5.2262E-11, 1.3337E-10, -15632.6109, 4.5594E-11, 1.105E-10, -13212.7861,
    3.8927E-11, 8.7629E-11, -10838.3833, -57283.6773, 1.909E-11, 1.5687E-10,
    -55563.7419, 1.9128E-11, 1.5227E-10, -73843.8065, 1.9167E-11, 1.4768E-10,
    -51170.2956, 1.9205E-11, 1.4309E-10, -48499.4565, 1.9244E-11, 1.3849E-10,
    -45833.9613, 1.9282E-11, 1.339E-10, -43176.4816, 1.932E-11, 1.2931E-10,
    -40529.6895, 1.9359E-11, 1.2471E-10, -37896.2568, 1.9397E-11, 1.2012E-10,
    -35278.8554, 1.9435E-11, 1.1552E-10, -32680.1571, 1.9474E-11, 1.1093E-10,
    -30102.834, 1.9512E-11, 1.0634E-10, -27549.5577, 1.955E-11, 1.0174E-10,
    -25023.0003, 1.9589E-11, 9.7149E-11, -22525.8335, 1.9627E-11, 9.2555E-11,
    -20060.7294, 1.9665E-11, 8.7961E-11, -17630.3597, 1.9704E-11, 8.3367E-11,
    -15237.3963, 1.9742E-11, 7.8773E-11, -12884.5112, 1.9781E-11, 7.4179E-11,
    -10574.3762, 1.9819E-11, 6.9585E-11, 2.0513E-11, -57283.6773, 4.9822E-10,
    1.984E-11, -55563.7419, 4.7644E-10, 1.9167E-11, -73843.8065, 4.5466E-10,
    1.8494E-11, -51170.2956, 4.3288E-10, 1.782E-11, -48499.4565, 4.111E-10,
    1.7147E-11, -45833.9613, 3.8932E-10, 1.6474E-11, -43176.4816, 3.6753E-10,
    1.5801E-11, -40529.6895, 3.4575E-10, 1.5127E-11, -37896.2568, 3.2397E-10,
    1.4454E-11, -35278.8554, 3.0219E-10, 1.3781E-11, -32680.1571, 2.8041E-10,
    1.3108E-11, -30102.834, 2.5863E-10, 1.2434E-11, -27549.5577, 2.3685E-10,
    1.1761E-11, -25023.0003, 2.1506E-10, 1.1088E-11, -22525.8335, 1.9328E-10,
    1.0415E-11, -20060.7294, 1.715E-10, 9.7415E-12, -17630.3597, 1.4972E-10,
    9.0683E-12, -15237.3963, 1.2794E-10, 8.395E-12, -12884.5112, 1.0616E-10,
    7.7218E-12, -10574.3762, 8.4375E-11, 1.606E-10, 4.9834E-10, -57283.6773,
    1.5414E-10, 4.765E-10, -55563.7419, 1.4768E-10, 4.5466E-10, -73843.8065,
    1.4122E-10, 4.3282E-10, -51170.2956, 1.3477E-10, 4.1098E-10, -48499.4565,
    1.2831E-10, 3.8914E-10, -45833.9613, 1.2185E-10, 3.673E-10, -43176.4816,
    1.1539E-10, 3.4546E-10, -40529.6895, 1.0894E-10, 3.2363E-10, -37896.2568,
    1.0248E-10, 3.0179E-10, -35278.8554, 9.6022E-11, 2.7995E-10, -32680.1571,
    8.9565E-11, 2.5811E-10, -30102.834, 8.3107E-11, 2.3627E-10, -27549.5577,
    7.665E-11, 2.1443E-10, -25023.0003, 7.0193E-11, 1.9259E-10, -22525.8335,
    6.3735E-11, 1.7075E-10, -20060.7294, 5.7278E-11, 1.4891E-10, -17630.3597,
    5.0821E-11, 1.2707E-10, -15237.3963, 4.4363E-11, 1.0523E-10, -12884.5112,
    3.7906E-11, 8.3395E-11, -10574.3762, -54396.1306, 1.8412E-11, 1.4999E-10,
    -52783.2131, 1.8453E-11, 1.4561E-10, -51170.2956, 1.8494E-11, 1.4122E-10,
    -69557.3781, 1.8534E-11, 1.3684E-10, -46990.885, 1.8575E-11, 1.3246E-10,
    -44427.0639, 1.8615E-11, 1.2807E-10, -41868.5865, 1.8656E-11, 1.2369E-10,
    -39318.1248, 1.8697E-11, 1.1931E-10, -36778.3505, 1.8737E-11, 1.1492E-10,
    -34251.9357, 1.8778E-11, 1.1054E-10, -31741.5522, 1.8819E-11, 1.0615E-10,
    -29249.8719, 1.8859E-11, 1.0177E-10, -26779.5666, 1.89E-11, 9.7387E-11,
    -24333.3082, 1.894E-11, 9.3003E-11, -21913.7687, 1.8981E-11, 8.8619E-11,
    -19523.6198, 1.9022E-11, 8.4236E-11, -17165.5336, 1.9062E-11, 7.9852E-11,
    -14842.1818, 1.9103E-11, 7.5468E-11, -12556.2363, 1.9144E-11, 7.1085E-11,
    -10310.3691, 1.9184E-11, 6.6701E-11, 2.0547E-11, -54396.1306, 4.7432E-10,
    1.9876E-11, -52783.2131, 4.5357E-10, 1.9205E-11, -51170.2956, 4.3282E-10,
    1.8534E-11, -69557.3781, 4.1207E-10, 1.7863E-11, -46990.885, 3.9132E-10,
    1.7192E-11, -44427.0639, 3.7057E-10, 1.6521E-11, -41868.5865, 3.4982E-10,
    1.585E-11, -39318.1248, 3.2907E-10, 1.5179E-11, -36778.3505, 3.0833E-10,
    1.4508E-11, -34251.9357, 2.8758E-10, 1.3837E-11, -31741.5522, 2.6683E-10,
    1.3166E-11, -29249.8719, 2.4608E-10, 1.2495E-11, -26779.5666, 2.2533E-10,
    1.1825E-11, -24333.3082, 2.0458E-10, 1.1154E-11, -21913.7687, 1.8383E-10,
    1.0483E-11, -19523.6198, 1.6308E-10, 9.8116E-12, -17165.5336, 1.4233E-10,
    9.1406E-12, -14842.1818, 1.2158E-10, 8.4697E-12, -12556.2363, 1.0083E-10,
    7.7987E-12, -10310.3691, 8.0084E-11, 1.5558E-10, 4.7449E-10, -54396.1306,
    1.4933E-10, 4.5369E-10, -52783.2131, 1.4309E-10, 4.3288E-10, -51170.2956,
    1.3684E-10, 4.1207E-10, -69557.3781, 1.3059E-10, 3.9126E-10, -46990.885,
    1.2435E-10, 3.7046E-10, -44427.0639, 1.181E-10, 3.4965E-10, -41868.5865,
    1.1185E-10, 3.2884E-10, -39318.1248, 1.056E-10, 3.0804E-10, -36778.3505,
    9.9357E-11, 2.8723E-10, -34251.9357, 9.311E-11, 2.6642E-10, -31741.5522,
    8.6862E-11, 2.4562E-10, -29249.8719, 8.0615E-11, 2.2481E-10, -26779.5666,
    7.4368E-11, 2.04E-10, -24333.3082, 6.8121E-11, 1.832E-10, -21913.7687,
    6.1874E-11, 1.6239E-10, -19523.6198, 5.5626E-11, 1.4158E-10, -17165.5336,
    4.9379E-11, 1.2078E-10, -14842.1818, 4.3132E-11, 9.9969E-11, -12556.2363,
    3.6885E-11, 7.9162E-11, -10310.3691, -51516.5995, 1.7735E-11, 1.4311E-10,
    -50008.028, 1.7777E-11, 1.3894E-10, -48499.4565, 1.782E-11, 1.3477E-10,
    -46990.885, 1.7863E-11, 1.3059E-10, -65482.3136, 1.7906E-11, 1.2642E-10,
    -43020.1665, 1.7949E-11, 1.2225E-10, -40560.6914, 1.7992E-11, 1.1807E-10,
    -38106.56, 1.8035E-11, 1.139E-10, -35660.4443, 1.8078E-11, 1.0973E-10,
    -33225.0161, 1.8121E-11, 1.0555E-10, -30802.9473, 1.8163E-11, 1.0138E-10,
    -28396.9098, 1.8206E-11, 9.7205E-11, -26009.5755, 1.8249E-11, 9.3031E-11,
    -23643.6162, 1.8292E-11, 8.8857E-11, -21301.7038, 1.8335E-11, 8.4684E-11,
    -18986.5103, 1.8378E-11, 8.051E-11, -16700.7075, 1.8421E-11, 7.6337E-11,
    -14446.9673, 1.8464E-11, 7.2163E-11, -12227.9615, 1.8507E-11, 6.799E-11,
    -10046.362, 1.855E-11, 6.3816E-11, 2.0581E-11, -51516.5995, 4.5042E-10,
    1.9912E-11, -50008.028, 4.307E-10, 1.9244E-11, -48499.4565, 4.1098E-10,
    1.8575E-11, -46990.885, 3.9126E-10, 1.7906E-11, -65482.3136, 3.7155E-10,
    1.7237E-11, -43020.1665, 3.5183E-10, 1.6569E-11, -40560.6914, 3.3211E-10,
    1.59E-11, -38106.56, 3.124E-10, 1.5231E-11, -35660.4443, 2.9268E-10,
    1.4563E-11, -33225.0161, 2.7296E-10, 1.3894E-11, -30802.9473, 2.5325E-10,
    1.3225E-11, -28396.9098, 2.3353E-10, 1.2557E-11, -26009.5755, 2.1381E-10,
    1.1888E-11, -23643.6162, 1.9409E-10, 1.1219E-11, -21301.7038, 1.7438E-10,
    1.055E-11, -18986.5103, 1.5466E-10, 9.8817E-12, -16700.7075, 1.3494E-10,
    9.213E-12, -14446.9673, 1.1523E-10, 8.5443E-12, -12227.9615, 9.551E-11,
    7.8756E-12, -10046.362, 7.5793E-11, 1.5057E-10, 4.5065E-10, -51516.5995,
    1.4453E-10, 4.3087E-10, -50008.028, 1.3849E-10, 4.111E-10, -48499.4565,
    1.3246E-10, 3.9132E-10, -46990.885, 1.2642E-10, 3.7155E-10, -65482.3136,
    1.2038E-10, 3.5177E-10, -43020.1665, 1.1435E-10, 3.32E-10, -40560.6914,
    1.0831E-10, 3.1222E-10, -38106.56, 1.0227E-10, 2.9245E-10, -35660.4443,
    9.6234E-11, 2.7267E-10, -33225.0161, 9.0197E-11, 2.529E-10, -30802.9473,
    8.416E-11, 2.3313E-10, -28396.9098, 7.8123E-11, 2.1335E-10, -26009.5755,
    7.2086E-11, 1.9358E-10, -23643.6162, 6.6049E-11, 1.738E-10, -21301.7038,
    6.0012E-11, 1.5403E-10, -18986.5103, 5.3975E-11, 1.3425E-10, -16700.7075,
    4.7938E-11, 1.1448E-10, -14446.9673, 4.1901E-11, 9.4703E-11, -12227.9615,
    3.5864E-11, 7.4928E-11, -10046.362, -48647.756, 1.7057E-11, 1.3624E-10,
    -47240.8586, 1.7102E-11, 1.3227E-10, -45833.9613, 1.7147E-11, 1.2831E-10,
    -44427.0639, 1.7192E-11, 1.2435E-10, -43020.1665, 1.7237E-11, 1.2038E-10,
    -61613.2692, 1.7283E-11, 1.1642E-10, -39252.7963, 1.7328E-11, 1.1246E-10,
    -36894.9953, 1.7373E-11, 1.0849E-10, -34542.538, 1.7418E-11, 1.0453E-10,
    -32198.0964, 1.7463E-11, 1.0057E-10, -29864.3424, 1.7508E-11, 9.6602E-11,
    -27543.9477, 1.7554E-11, 9.2639E-11, -25239.5843, 1.7599E-11, 8.8675E-11,
    -22953.9241, 1.7644E-11, 8.4712E-11, -20689.639, 1.7689E-11, 8.0748E-11,
    -18449.4008, 1.7734E-11, 7.6785E-11, -16235.8814, 1.7779E-11, 7.2822E-11,
    -14051.7527, 1.7825E-11, 6.8858E-11, -11899.6866, 1.787E-11, 6.4895E-11,
    -9782.3549, 1.7915E-11, 6.0932E-11, 2.0615E-11, -48647.756, 4.2651E-10,
    1.9948E-11, -47240.8586, 4.0783E-10, 1.9282E-11, -45833.9613, 3.8914E-10,
    1.8615E-11, -44427.0639, 3.7046E-10, 1.7949E-11, -43020.1665, 3.5177E-10,
    1.7283E-11, -61613.2692, 3.3309E-10, 1.6616E-11, -39252.7963, 3.144E-10,
    1.595E-11, -36894.9953, 2.9572E-10, 1.5283E-11, -34542.538, 2.7703E-10,
    1.4617E-11, -32198.0964, 2.5835E-10, 1.395E-11, -29864.3424, 2.3966E-10,
    1.3284E-11, -27543.9477, 2.2098E-10, 1.2618E-11, -25239.5843, 2.0229E-10,
    1.1951E-11, -22953.9241, 1.8361E-10, 1.1285E-11, -20689.639, 1.6493E-10,
    1.0618E-11, -18449.4008, 1.4624E-10, 9.9518E-12, -16235.8814, 1.2756E-10,
    9.2854E-12, -14051.7527, 1.0887E-10, 8.6189E-12, -11899.6866, 9.0186E-11,
    7.9525E-12, -9782.3549, 7.1501E-11, 1.4555E-10, 4.268E-10, -48647.756,
    1.3973E-10, 4.0806E-10, -47240.8586, 1.339E-10, 3.8932E-10, -45833.9613,
    1.2807E-10, 3.7057E-10, -44427.0639, 1.2225E-10, 3.5183E-10, -43020.1665,
    1.1642E-10, 3.3309E-10, -61613.2692, 1.1059E-10, 3.1435E-10, -39252.7963,
    1.0477E-10, 2.956E-10, -36894.9953, 9.8938E-11, 2.7686E-10, -34542.538,
    9.3111E-11, 2.5812E-10, -32198.0964, 8.7284E-11, 2.3938E-10, -29864.3424,
    8.1458E-11, 2.2063E-10, -27543.9477, 7.5631E-11, 2.0189E-10, -25239.5843,
    6.9804E-11, 1.8315E-10, -22953.9241, 6.3977E-11, 1.6441E-10, -20689.639,
    5.815E-11, 1.4566E-10, -18449.4008, 5.2323E-11, 1.2692E-10, -16235.8814,
    4.6496E-11, 1.0818E-10, -14051.7527, 4.0669E-11, 8.9437E-11, -11899.6866,
    3.4843E-11, 7.0694E-11, -9782.3549, -45792.2718, 1.6379E-11, 1.2936E-10,
    -44484.3767, 1.6426E-11, 1.256E-10, -43176.4816, 1.6474E-11, 1.2185E-10,
    -41868.5865, 1.6521E-11, 1.181E-10, -40560.6914, 1.6569E-11, 1.1435E-10,
    -39252.7963, 1.6616E-11, 1.1059E-10, -57944.9012, 1.6664E-11, 1.0684E-10,
    -35683.4305, 1.6711E-11, 1.0309E-10, -33424.6318, 1.6758E-11, 9.9332E-11,
    -31171.1768, 1.6806E-11, 9.5579E-11, -28925.7374, 1.6853E-11, 9.1826E-11,
    -26690.9856, 1.6901E-11, 8.8073E-11, -24469.5932, 1.6948E-11, 8.4319E-11,
    -22264.2321, 1.6996E-11, 8.0566E-11, -20077.5741, 1.7043E-11, 7.6813E-11,
    -17912.2913, 1.709E-11, 7.306E-11, -15771.0553, 1.7138E-11, 6.9307E-11,
    -13656.5382, 1.7185E-11, 6.5553E-11, -11571.4117, 1.7233E-11, 6.18E-11,
    -9518.3479, 1.728E-11, 5.8047E-11, 2.0649E-11, -45792.2718, 4.0261E-10,
    1.9984E-11, -44484.3767, 3.8496E-10, 1.932E-11, -43176.4816, 3.673E-10,
    1.8656E-11, -41868.5865, 3.4965E-10, 1.7992E-11, -40560.6914, 3.32E-10,
    1.7328E-11, -39252.7963, 3.1435E-10, 1.6664E-11, -57944.9012, 2.9669E-10,
    1.5999E-11, -35683.4305, 2.7904E-10, 1.5335E-11, -33424.6318, 2.6139E-10,
    1.4671E-11, -31171.1768, 2.4374E-10, 1.4007E-11, -28925.7374, 2.2608E-10,
    1.3343E-11, -26690.9856, 2.0843E-10, 1.2679E-11, -24469.5932, 1.9078E-10,
    1.2014E-11, -22264.2321, 1.7313E-10, 1.135E-11, -20077.5741, 1.5547E-10,
    1.0686E-11, -17912.2913, 1.3782E-10, 1.0022E-11, -15771.0553, 1.2017E-10,
    9.3577E-12, -13656.5382, 1.0252E-10, 8.6936E-12, -11571.4117, 8.4863E-11,
    8.0294E-12, -9518.3479, 6.721E-11, 1.4054E-10, 4.0295E-10, -45792.2718,
    1.3492E-10, 3.8524E-10, -44484.3767, 1.2931E-10, 3.6753E-10, -43176.4816,
    1.2369E-10, 3.4982E-10, -41868.5865, 1.1807E-10, 3.3211E-10, -40560.6914,
    1.1246E-10, 3.144E-10, -39252.7963, 1.0684E-10, 2.9669E-10, -57944.9012,
    1.0122E-10, 2.7898E-10, -35683.4305, 9.5605E-11, 2.6127E-10, -33424.6318,
    8.9989E-11, 2.4356E-10, -31171.1768, 8.4372E-11, 2.2585E-10, -28925.7374,
    7.8755E-11, 2.0814E-10, -26690.9856, 7.3138E-11, 1.9043E-10, -24469.5932,
    6.7522E-11, 1.7272E-10, -22264.2321, 6.1905E-11, 1.5501E-10, -20077.5741,
    5.6288E-11, 1.373E-10, -17912.2913, 5.0672E-11, 1.1959E-10, -15771.0553,
    4.5055E-11, 1.0188E-10, -13656.5382, 3.9438E-11, 8.4171E-11, -11571.4117,
    3.3821E-11, 6.6461E-11, -9518.3479, -42952.819, 1.5701E-11, 1.2248E-10,
    -41741.2542, 1.5751E-11, 1.1894E-10, -40529.6895, 1.5801E-11, 1.1539E-10,
    -39318.1248, 1.585E-11, 1.1185E-10, -38106.56, 1.59E-11, 1.0831E-10,
    -36894.9953, 1.595E-11, 1.0477E-10, -35683.4305, 1.5999E-11, 1.0122E-10,
    -54471.8658, 1.6049E-11, 9.7679E-11, -32306.7255, 1.6099E-11, 9.4136E-11,
    -30144.2571, 1.6148E-11, 9.0593E-11, -27987.1325, 1.6198E-11, 8.705E-11,
    -25838.0235, 1.6248E-11, 8.3507E-11, -23699.6021, 1.6298E-11, 7.9964E-11,
    -21574.54, 1.6347E-11, 7.6421E-11, -19465.5093, 1.6397E-11, 7.2878E-11,
    -17375.1817, 1.6447E-11, 6.9335E-11, -15306.2292, 1.6496E-11, 6.5792E-11,
    -13261.3236, 1.6546E-11, 6.2248E-11, -11243.1368, 1.6596E-11, 5.8705E-11,
    -9254.3408, 1.6645E-11, 5.5162E-11, 2.0682E-11, -42952.819, 3.7871E-10,
    2.0021E-11, -41741.2542, 3.6208E-10, 1.9359E-11, -40529.6895, 3.4546E-10,
    1.8697E-11, -39318.1248, 3.2884E-10, 1.8035E-11, -38106.56, 3.1222E-10,
    1.7373E-11, -36894.9953, 2.956E-10, 1.6711E-11, -35683.4305, 2.7898E-10,
    1.6049E-11, -54471.8658, 2.6236E-10, 1.5387E-11, -32306.7255, 2.4574E-10,
    1.4725E-11, -30144.2571, 2.2912E-10, 1.4063E-11, -27987.1325, 2.125E-10,
    1.3402E-11, -25838.0235, 1.9588E-10, 1.274E-11, -23699.6021, 1.7926E-10,
    1.2078E-11, -21574.54, 1.6264E-10, 1.1416E-11, -19465.5093, 1.4602E-10,
    1.0754E-11, -17375.1817, 1.294E-10, 1.0092E-11, -15306.2292, 1.1278E-10,
    9.4301E-12, -13261.3236, 9.6159E-11, 8.7682E-12, -11243.1368, 7.9539E-11,
    8.1063E-12, -9254.3408, 6.2919E-11, 1.3552E-10, 3.7911E-10, -42952.819,
    1.3012E-10, 3.6243E-10, -41741.2542, 1.2471E-10, 3.4575E-10, -40529.6895,
    1.1931E-10, 3.2907E-10, -39318.1248, 1.139E-10, 3.124E-10, -38106.56,
    1.0849E-10, 2.9572E-10, -36894.9953, 1.0309E-10, 2.7904E-10, -35683.4305,
    9.7679E-11, 2.6236E-10, -54471.8658, 9.2272E-11, 2.4568E-10, -32306.7255,
    8.6866E-11, 2.2901E-10, -30144.2571, 8.1459E-11, 2.1233E-10, -27987.1325,
    7.6053E-11, 1.9565E-10, -25838.0235, 7.0646E-11, 1.7897E-10, -23699.6021,
    6.524E-11, 1.623E-10, -21574.54, 5.9833E-11, 1.4562E-10, -19465.5093,
    5.4427E-11, 1.2894E-10, -17375.1817, 4.902E-11, 1.1226E-10, -15306.2292,
    4.3613E-11, 9.5583E-11, -13261.3236, 3.8207E-11, 7.8905E-11, -11243.1368,
    3.28E-11, 6.2227E-11, -9254.3408, -40132.0693, 1.5023E-11, 1.156E-10,
    -39014.1631, 1.5075E-11, 1.1227E-10, -37896.2568, 1.5127E-11, 1.0894E-10,
    -36778.3505, 1.5179E-11, 1.056E-10, -35660.4443, 1.5231E-11, 1.0227E-10,
    -34542.538, 1.5283E-11, 9.8938E-11, -33424.6318, 1.5335E-11, 9.5605E-11,
    -32306.7255, 1.5387E-11, 9.2272E-11, -51188.8193, 1.5439E-11, 8.894E-11,
    -29117.3375, 1.5491E-11, 8.5607E-11, -27048.5276, 1.5543E-11, 8.2274E-11,
    -24985.0615, 1.5595E-11, 7.8941E-11, -22929.611, 1.5647E-11, 7.5608E-11,
    -20884.848, 1.5699E-11, 7.2275E-11, -18853.4444, 1.5751E-11, 6.8942E-11,
    -16838.0722, 1.5803E-11, 6.5609E-11, -14841.4031, 1.5855E-11, 6.2276E-11,
    -12866.1091, 1.5907E-11, 5.8944E-11, -10914.862, 1.5959E-11, 5.5611E-11,
    -8990.3337, 1.6011E-11, 5.2278E-11, 2.0716E-11, -40132.0693, 3.548E-10,
    2.0057E-11, -39014.1631, 3.3921E-10, 1.9397E-11, -37896.2568, 3.2363E-10,
    1.8737E-11, -36778.3505, 3.0804E-10, 1.8078E-11, -35660.4443, 2.9245E-10,
    1.7418E-11, -34542.538, 2.7686E-10, 1.6758E-11, -33424.6318, 2.6127E-10,
    1.6099E-11, -32306.7255, 2.4568E-10, 1.5439E-11, -51188.8193, 2.301E-10,
    1.478E-11, -29117.3375, 2.1451E-10, 1.412E-11, -27048.5276, 1.9892E-10,
    1.346E-11, -24985.0615, 1.8333E-10, 1.2801E-11, -22929.611, 1.6774E-10,
    1.2141E-11, -20884.848, 1.5216E-10, 1.1481E-11, -18853.4444, 1.3657E-10,
    1.0822E-11, -16838.0722, 1.2098E-10, 1.0162E-11, -14841.4031, 1.0539E-10,
    9.5025E-12, -12866.1091, 8.9804E-11, 8.8428E-12, -10914.862, 7.4216E-11,
    8.1832E-12, -8990.3337, 5.8627E-11, 1.3051E-10, 3.5526E-10, -40132.0693,
    1.2531E-10, 3.3962E-10, -39014.1631, 1.2012E-10, 3.2397E-10, -37896.2568,
    1.1492E-10, 3.0833E-10, -36778.3505, 1.0973E-10, 2.9268E-10, -35660.4443,
    1.0453E-10, 2.7703E-10, -34542.538, 9.9332E-11, 2.6139E-10, -33424.6318,
    9.4136E-11, 2.4574E-10, -32306.7255, 8.894E-11, 2.301E-10, -51188.8193,
    8.3743E-11, 2.1445E-10, -29117.3375, 7.8547E-11, 1.9881E-10, -27048.5276,
    7.335E-11, 1.8316E-10, -24985.0615, 6.8154E-11, 1.6751E-10, -22929.611,
    6.2958E-11, 1.5187E-10, -20884.848, 5.7761E-11, 1.3622E-10, -18853.4444,
    5.2565E-11, 1.2058E-10, -16838.0722, 4.7368E-11, 1.0493E-10, -14841.4031,
    4.2172E-11, 8.9285E-11, -12866.1091, 3.6976E-11, 7.3639E-11, -10914.862,
    3.1779E-11, 5.7993E-11, -8990.3337, -37332.6947, 1.4346E-11, 1.0872E-10,
    -36305.775, 1.44E-11, 1.056E-10, -35278.8554, 1.4454E-11, 1.0248E-10,
    -34251.9357, 1.4508E-11, 9.9357E-11, -33225.0161, 1.4563E-11, 9.6234E-11,
    -32198.0964, 1.4617E-11, 9.3111E-11, -31171.1768, 1.4671E-11, 8.9989E-11,
    -30144.2571, 1.4725E-11, 8.6866E-11, -29117.3375, 1.478E-11, 8.3743E-11,
    -48090.4178, 1.4834E-11, 8.062E-11, -26109.9227, 1.4888E-11, 7.7498E-11,
    -24132.0994, 1.4942E-11, 7.4375E-11, -22159.6198, 1.4996E-11, 7.1252E-11,
    -20195.156, 1.5051E-11, 6.813E-11, -18241.3796, 1.5105E-11, 6.5007E-11,
    -16300.9627, 1.5159E-11, 6.1884E-11, -14376.577, 1.5213E-11, 5.8761E-11,
    -12470.8945, 1.5268E-11, 5.5639E-11, -10586.5871, 1.5322E-11, 5.2516E-11,
    -8726.3266, 1.5376E-11, 4.9393E-11, 2.075E-11, -37332.6947, 3.309E-10,
    2.0093E-11, -36305.775, 3.1634E-10, 1.9435E-11, -35278.8554, 3.0179E-10,
    1.8778E-11, -34251.9357, 2.8723E-10, 1.8121E-11, -33225.0161, 2.7267E-10,
    1.7463E-11, -32198.0964, 2.5812E-10, 1.6806E-11, -31171.1768, 2.4356E-10,
    1.6148E-11, -30144.2571, 2.2901E-10, 1.5491E-11, -29117.3375, 2.1445E-10,
    1.4834E-11, -48090.4178, 1.999E-10, 1.4176E-11, -26109.9227, 1.8534E-10,
    1.3519E-11, -24132.0994, 1.7078E-10, 1.2862E-11, -22159.6198, 1.5623E-10,
    1.2204E-11, -20195.156, 1.4167E-10, 1.1547E-11, -18241.3796, 1.2712E-10,
    1.089E-11, -16300.9627, 1.1256E-10, 1.0232E-11, -14376.577, 9.8004E-11,
    9.5748E-12, -12470.8945, 8.3448E-11, 8.9175E-12, -10586.5871, 6.8892E-11,
    8.2601E-12, -8726.3266, 5.4336E-11, 1.255E-10, 3.3142E-10, -37332.6947,
    1.2051E-10, 3.168E-10, -36305.775, 1.1552E-10, 3.0219E-10, -35278.8554,
    1.1054E-10, 2.8758E-10, -34251.9357, 1.0555E-10, 2.7296E-10, -33225.0161,
    1.0057E-10, 2.5835E-10, -32198.0964, 9.5579E-11, 2.4374E-10, -31171.1768,
    9.0593E-11, 2.2912E-10, -30144.2571, 8.5607E-11, 2.1451E-10, -29117.3375,
    8.062E-11, 1.999E-10, -48090.4178, 7.5634E-11, 1.8528E-10, -26109.9227,
    7.0648E-11, 1.7067E-10, -24132.0994, 6.5662E-11, 1.5605E-10, -22159.6198,
    6.0676E-11, 1.4144E-10, -20195.156, 5.5689E-11, 1.2683E-10, -18241.3796,
    5.0703E-11, 1.1221E-10, -16300.9627, 4.5717E-11, 9.76E-11, -14376.577,
    4.0731E-11, 8.2987E-11, -12470.8945, 3.5744E-11, 6.8373E-11, -10586.5871,
    3.0758E-11, 5.376E-11, -8726.3266, -34557.367, 1.3668E-11, 1.0185E-10,
    -33618.7621, 1.3724E-11, 9.8935E-11, -32680.1571, 1.3781E-11, 9.6022E-11,
    -31741.5522, 1.3837E-11, 9.311E-11, -30802.9473, 1.3894E-11, 9.0197E-11,
    -29864.3424, 1.395E-11, 8.7284E-11, -28925.7374, 1.4007E-11, 8.4372E-11,
    -27987.1325, 1.4063E-11, 8.1459E-11, -27048.5276, 1.412E-11, 7.8547E-11,
    -26109.9227, 1.4176E-11, 7.5634E-11, -45171.3177, 1.4233E-11, 7.2722E-11,
    -23279.1373, 1.4289E-11, 6.9809E-11, -21389.6287, 1.4346E-11, 6.6897E-11,
    -19505.4639, 1.4402E-11, 6.3984E-11, -17629.3148, 1.4459E-11, 6.1071E-11,
    -15763.8531, 1.4515E-11, 5.8159E-11, -13911.7509, 1.4572E-11, 5.5246E-11,
    -12075.68, 1.4628E-11, 5.2334E-11, -10258.3122, 1.4685E-11, 4.9421E-11,
    -8462.3195, 1.4741E-11, 4.6509E-11, 2.0784E-11, -34557.367, 3.0699E-10,
    2.0129E-11, -33618.7621, 2.9347E-10, 1.9474E-11, -32680.1571, 2.7995E-10,
    1.8819E-11, -31741.5522, 2.6642E-10, 1.8163E-11, -30802.9473, 2.529E-10,
    1.7508E-11, -29864.3424, 2.3938E-10, 1.6853E-11, -28925.7374, 2.2585E-10,
    1.6198E-11, -27987.1325, 2.1233E-10, 1.5543E-11, -27048.5276, 1.9881E-10,
    1.4888E-11, -26109.9227, 1.8528E-10, 1.4233E-11, -45171.3177, 1.7176E-10,
    1.3578E-11, -23279.1373, 1.5823E-10, 1.2923E-11, -21389.6287, 1.4471E-10,
    1.2268E-11, -19505.4639, 1.3119E-10, 1.1612E-11, -17629.3148, 1.1766E-10,
    1.0957E-11, -15763.8531, 1.0414E-10, 1.0302E-11, -13911.7509, 9.0616E-11,
    9.6472E-12, -12075.68, 7.7092E-11, 8.9921E-12, -10258.3122, 6.3569E-11,
    8.337E-12, -8462.3195, 5.0045E-11, 1.2048E-10, 3.0757E-10, -34557.367,
    1.1571E-10, 2.9399E-10, -33618.7621, 1.1093E-10, 2.8041E-10, -32680.1571,
    1.0615E-10, 2.6683E-10, -31741.5522, 1.0138E-10, 2.5325E-10, -30802.9473,
    9.6602E-11, 2.3966E-10, -29864.3424, 9.1826E-11, 2.2608E-10, -28925.7374,
    8.705E-11, 2.125E-10, -27987.1325, 8.2274E-11, 1.9892E-10, -27048.5276,
    7.7498E-11, 1.8534E-10, -26109.9227, 7.2722E-11, 1.7176E-10, -45171.3177,
    6.7946E-11, 1.5818E-10, -23279.1373, 6.317E-11, 1.446E-10, -21389.6287,
    5.8393E-11, 1.3101E-10, -19505.4639, 5.3617E-11, 1.1743E-10, -17629.3148,
    4.8841E-11, 1.0385E-10, -15763.8531, 4.4065E-11, 9.027E-11, -13911.7509,
    3.9289E-11, 7.6689E-11, -12075.68, 3.4513E-11, 6.3107E-11, -10258.3122,
    2.9737E-11, 4.9526E-11, -8462.3195, -31808.7581, 1.299E-11, 9.497E-11,
    -30955.796, 1.3049E-11, 9.2267E-11, -30102.834, 1.3108E-11, 8.9565E-11,
    -29249.8719, 1.3166E-11, 8.6862E-11, -28396.9098, 1.3225E-11, 8.416E-11,
    -27543.9477, 1.3284E-11, 8.1458E-11, -26690.9856, 1.3343E-11, 7.8755E-11,
    -25838.0235, 1.3402E-11, 7.6053E-11, -24985.0615, 1.346E-11, 7.335E-11,
    -24132.0994, 1.3519E-11, 7.0648E-11, -23279.1373, 1.3578E-11, 6.7946E-11,
    -42426.1752, 1.3637E-11, 6.5243E-11, -20619.6376, 1.3695E-11, 6.2541E-11,
    -18815.7719, 1.3754E-11, 5.9838E-11, -17017.2499, 1.3813E-11, 5.7136E-11,
    -15226.7436, 1.3872E-11, 5.4434E-11, -13446.9248, 1.393E-11, 5.1731E-11,
    -11680.4654, 1.3989E-11, 4.9029E-11, -9930.0373, 1.4048E-11, 4.6326E-11,
    -8198.3124, 1.4107E-11, 4.3624E-11, 2.0818E-11, -31808.7581, 2.8309E-10,
    2.0165E-11, -30955.796, 2.706E-10, 1.9512E-11, -30102.834, 2.5811E-10,
    1.8859E-11, -29249.8719, 2.4562E-10, 1.8206E-11, -28396.9098, 2.3313E-10,
    1.7554E-11, -27543.9477, 2.2063E-10, 1.6901E-11, -26690.9856, 2.0814E-10,
    1.6248E-11, -25838.0235, 1.9565E-10, 1.5595E-11, -24985.0615, 1.8316E-10,
    1.4942E-11, -24132.0994, 1.7067E-10, 1.4289E-11, -23279.1373, 1.5818E-10,
    1.3637E-11, -42426.1752, 1.4569E-10, 1.2984E-11, -20619.6376, 1.3319E-10,
    1.2331E-11, -18815.7719, 1.207E-10, 1.1678E-11, -17017.2499, 1.0821E-10,
    1.1025E-11, -15226.7436, 9.5719E-11, 1.0372E-11, -13446.9248, 8.3228E-11,
    9.7196E-12, -11680.4654, 7.0736E-11, 9.0667E-12, -9930.0373, 5.8245E-11,
    8.4139E-12, -8198.3124, 4.5754E-11, 1.1547E-10, 2.8373E-10, -31808.7581,
    1.109E-10, 2.7118E-10, -30955.796, 1.0634E-10, 2.5863E-10, -30102.834,
    1.0177E-10, 2.4608E-10, -29249.8719, 9.7205E-11, 2.3353E-10, -28396.9098,
    9.2639E-11, 2.2098E-10, -27543.9477, 8.8073E-11, 2.0843E-10, -26690.9856,
    8.3507E-11, 1.9588E-10, -25838.0235, 7.8941E-11, 1.8333E-10, -24985.0615,
    7.4375E-11, 1.7078E-10, -24132.0994, 6.9809E-11, 1.5823E-10, -23279.1373,
    6.5243E-11, 1.4569E-10, -42426.1752, 6.0677E-11, 1.3314E-10, -20619.6376,
    5.6111E-11, 1.2059E-10, -18815.7719, 5.1545E-11, 1.0804E-10, -17017.2499,
    4.698E-11, 9.5489E-11, -15226.7436, 4.2414E-11, 8.294E-11, -13446.9248,
    3.7848E-11, 7.0391E-11, -11680.4654, 3.3282E-11, 5.7842E-11, -9930.0373,
    2.8716E-11, 4.5293E-11, -8198.3124, -29089.5399, 1.2312E-11, 8.8092E-11,
    -28319.5488, 1.2373E-11, 8.56E-11, -27549.5577, 1.2434E-11, 8.3107E-11,
    -26779.5666, 1.2495E-11, 8.0615E-11, -26009.5755, 1.2557E-11, 7.8123E-11,
    -25239.5843, 1.2618E-11, 7.5631E-11, -24469.5932, 1.2679E-11, 7.3138E-11,
    -23699.6021, 1.274E-11, 7.0646E-11, -22929.611, 1.2801E-11, 6.8154E-11,
    -22159.6198, 1.2862E-11, 6.5662E-11, -21389.6287, 1.2923E-11, 6.317E-11,
    -20619.6376, 1.2984E-11, 6.0677E-11, -39849.6465, 1.3045E-11, 5.8185E-11,
    -18126.0798, 1.3106E-11, 5.5693E-11, -16405.1851, 1.3167E-11, 5.3201E-11,
    -14689.6341, 1.3228E-11, 5.0708E-11, -12982.0987, 1.3289E-11, 4.8216E-11,
    -11285.2509, 1.335E-11, 4.5724E-11, -9601.7625, 1.3411E-11, 4.3232E-11,
    -7934.3053, 1.3472E-11, 4.0739E-11, 2.0852E-11, -29089.5399, 2.5919E-10,
    2.0201E-11, -28319.5488, 2.4773E-10, 1.955E-11, -27549.5577, 2.3627E-10,
    1.89E-11, -26779.5666, 2.2481E-10, 1.8249E-11, -26009.5755, 2.1335E-10,
    1.7599E-11, -25239.5843, 2.0189E-10, 1.6948E-11, -24469.5932, 1.9043E-10,
    1.6298E-11, -23699.6021, 1.7897E-10, 1.5647E-11, -22929.611, 1.6751E-10,
    1.4996E-11, -22159.6198, 1.5605E-10, 1.4346E-11, -21389.6287, 1.446E-10,
    1.3695E-11, -20619.6376, 1.3314E-10, 1.3045E-11, -39849.6465, 1.2168E-10,
    1.2394E-11, -18126.0798, 1.1022E-10, 1.1744E-11, -16405.1851, 9.8758E-11,
    1.1093E-11, -14689.6341, 8.7299E-11, 1.0443E-11, -12982.0987, 7.584E-11,
    9.7919E-12, -11285.2509, 6.4381E-11, 9.1414E-12, -9601.7625, 5.2922E-11,
    8.4908E-12, -7934.3053, 4.1462E-11, 1.1045E-10, 2.5988E-10, -29089.5399,
    1.061E-10, 2.4836E-10, -28319.5488, 1.0174E-10, 2.3685E-10, -27549.5577,
    9.7387E-11, 2.2533E-10, -26779.5666, 9.3031E-11, 2.1381E-10, -26009.5755,
    8.8675E-11, 2.0229E-10, -25239.5843, 8.4319E-11, 1.9078E-10, -24469.5932,
    7.9964E-11, 1.7926E-10, -23699.6021, 7.5608E-11, 1.6774E-10, -22929.611,
    7.1252E-11, 1.5623E-10, -22159.6198, 6.6897E-11, 1.4471E-10, -21389.6287,
    6.2541E-11, 1.3319E-10, -20619.6376, 5.8185E-11, 1.2168E-10, -39849.6465,
    5.3829E-11, 1.1016E-10, -18126.0798, 4.9474E-11, 9.8643E-11, -16405.1851,
    4.5118E-11, 8.7126E-11, -14689.6341, 4.0762E-11, 7.5609E-11, -12982.0987,
    3.6406E-11, 6.4093E-11, -11285.2509, 3.2051E-11, 5.2576E-11, -9601.7625,
    2.7695E-11, 4.1059E-11, -7934.3053, -26402.3844, 1.1635E-11, 8.1214E-11,
    -25712.6923, 1.1698E-11, 7.8932E-11, -25023.0003, 1.1761E-11, 7.665E-11,
    -24333.3082, 1.1825E-11, 7.4368E-11, -23643.6162, 1.1888E-11, 7.2086E-11,
    -22953.9241, 1.1951E-11, 6.9804E-11, -22264.2321, 1.2014E-11, 6.7522E-11,
    -21574.54, 1.2078E-11, 6.524E-11, -20884.848, 1.2141E-11, 6.2958E-11,
    -20195.156, 1.2204E-11, 6.0676E-11, -19505.4639, 1.2268E-11, 5.8393E-11,
    -18815.7719, 1.2331E-11, 5.6111E-11, -18126.0798, 1.2394E-11, 5.3829E-11,
    -37436.3878, 1.2457E-11, 5.1547E-11, -15793.1202, 1.2521E-11, 4.9265E-11,
    -14152.5245, 1.2584E-11, 4.6983E-11, -12517.2726, 1.2647E-11, 4.4701E-11,
    -10890.0363, 1.2711E-11, 4.2419E-11, -9273.4876, 1.2774E-11, 4.0137E-11,
    -7670.2982, 1.2837E-11, 3.7855E-11, 2.0885E-11, -26402.3844, 2.3528E-10,
    2.0237E-11, -25712.6923, 2.2486E-10, 1.9589E-11, -25023.0003, 2.1443E-10,
    1.894E-11, -24333.3082, 2.04E-10, 1.8292E-11, -23643.6162, 1.9358E-10,
    1.7644E-11, -22953.9241, 1.8315E-10, 1.6996E-11, -22264.2321, 1.7272E-10,
    1.6347E-11, -21574.54, 1.623E-10, 1.5699E-11, -20884.848, 1.5187E-10,
    1.5051E-11, -20195.156, 1.4144E-10, 1.4402E-11, -19505.4639, 1.3101E-10,
    1.3754E-11, -18815.7719, 1.2059E-10, 1.3106E-11, -18126.0798, 1.1016E-10,
    1.2457E-11, -37436.3878, 9.9733E-11, 1.1809E-11, -15793.1202, 8.9306E-11,
    1.1161E-11, -14152.5245, 7.8879E-11, 1.0513E-11, -12517.2726, 6.8452E-11,
    9.8643E-12, -10890.0363, 5.8025E-11, 9.216E-12, -9273.4876, 4.7598E-11,
    8.5677E-12, -7670.2982, 3.7171E-11, 1.0544E-10, 2.3603E-10, -26402.3844,
    1.0129E-10, 2.2555E-10, -25712.6923, 9.7149E-11, 2.1506E-10, -25023.0003,
    9.3003E-11, 2.0458E-10, -24333.3082, 8.8857E-11, 1.9409E-10, -23643.6162,
    8.4712E-11, 1.8361E-10, -22953.9241, 8.0566E-11, 1.7313E-10, -22264.2321,
    7.6421E-11, 1.6264E-10, -21574.54, 7.2275E-11, 1.5216E-10, -20884.848,
    6.813E-11, 1.4167E-10, -20195.156, 6.3984E-11, 1.3119E-10, -19505.4639,
    5.9838E-11, 1.207E-10, -18815.7719, 5.5693E-11, 1.1022E-10, -18126.0798,
    5.1547E-11, 9.9733E-11, -37436.3878, 4.7402E-11, 8.9248E-11, -15793.1202,
    4.3256E-11, 7.8764E-11, -14152.5245, 3.911E-11, 6.8279E-11, -12517.2726,
    3.4965E-11, 5.7794E-11, -10890.0363, 3.0819E-11, 4.731E-11, -9273.4876,
    2.6674E-11, 3.6825E-11, -7670.2982, -23749.9632, 1.0957E-11, 7.4336E-11,
    -23137.8984, 1.1022E-11, 7.2265E-11, -22525.8335, 1.1088E-11, 7.0193E-11,
    -21913.7687, 1.1154E-11, 6.8121E-11, -21301.7038, 1.1219E-11, 6.6049E-11,
    -20689.639, 1.1285E-11, 6.3977E-11, -20077.5741, 1.135E-11, 6.1905E-11,
    -19465.5093, 1.1416E-11, 5.9833E-11, -18853.4444, 1.1481E-11, 5.7761E-11,
    -18241.3796, 1.1547E-11, 5.5689E-11, -17629.3148, 1.1612E-11, 5.3617E-11,
    -17017.2499, 1.1678E-11, 5.1545E-11, -16405.1851, 1.1744E-11, 4.9474E-11,
    -15793.1202, 1.1809E-11, 4.7402E-11, -35181.0554, 1.1875E-11, 4.533E-11,
    -13615.415, 1.194E-11, 4.3258E-11, -12052.4465, 1.2006E-11, 4.1186E-11,
    -10494.8218, 1.2071E-11, 3.9114E-11, -8945.2127, 1.2137E-11, 3.7042E-11,
    -7406.2911, 1.2203E-11, 3.497E-11, 2.0919E-11, -23749.9632, 2.1138E-10,
    2.0273E-11, -23137.8984, 2.0199E-10, 1.9627E-11, -22525.8335, 1.9259E-10,
    1.8981E-11, -21913.7687, 1.832E-10, 1.8335E-11, -21301.7038, 1.738E-10,
    1.7689E-11, -20689.639, 1.6441E-10, 1.7043E-11, -20077.5741, 1.5501E-10,
    1.6397E-11, -19465.5093, 1.4562E-10, 1.5751E-11, -18853.4444, 1.3622E-10,
    1.5105E-11, -18241.3796, 1.2683E-10, 1.4459E-11, -17629.3148, 1.1743E-10,
    1.3813E-11, -17017.2499, 1.0804E-10, 1.3167E-11, -16405.1851, 9.8643E-11,
    1.2521E-11, -15793.1202, 8.9248E-11, 1.1875E-11, -35181.0554, 7.9854E-11,
    1.1229E-11, -13615.415, 7.0459E-11, 1.0583E-11, -12052.4465, 6.1064E-11,
    9.9367E-12, -10494.8218, 5.1669E-11, 9.2906E-12, -8945.2127, 4.2275E-11,
    8.6446E-12, -7406.2911, 3.288E-11, 1.0043E-10, 2.1219E-10, -23749.9632,
    9.649E-11, 2.0273E-10, -23137.8984, 9.2555E-11, 1.9328E-10, -22525.8335,
    8.8619E-11, 1.8383E-10, -21913.7687, 8.4684E-11, 1.7438E-10, -21301.7038,
    8.0748E-11, 1.6493E-10, -20689.639, 7.6813E-11, 1.5547E-10, -20077.5741,
    7.2878E-11, 1.4602E-10, -19465.5093, 6.8942E-11, 1.3657E-10, -18853.4444,
    6.5007E-11, 1.2712E-10, -18241.3796, 6.1071E-11, 1.1766E-10, -17629.3148,
    5.7136E-11, 1.0821E-10, -17017.2499, 5.3201E-11, 9.8758E-11, -16405.1851,
    4.9265E-11, 8.9306E-11, -15793.1202, 4.533E-11, 7.9854E-11, -35181.0554,
    4.1394E-11, 7.0401E-11, -13615.415, 3.7459E-11, 6.0949E-11, -12052.4465,
    3.3523E-11, 5.1496E-11, -10494.8218, 2.9588E-11, 4.2044E-11, -8945.2127,
    2.5653E-11, 3.2592E-11, -7406.2911, -21134.9484, 1.0279E-11, 6.7459E-11,
    -20597.8389, 1.0347E-11, 6.5597E-11, -20060.7294, 1.0415E-11, 6.3735E-11,
    -19523.6198, 1.0483E-11, 6.1874E-11, -18986.5103, 1.055E-11, 6.0012E-11,
    -18449.4008, 1.0618E-11, 5.815E-11, -17912.2913, 1.0686E-11, 5.6288E-11,
    -17375.1817, 1.0754E-11, 5.4427E-11, -16838.0722, 1.0822E-11, 5.2565E-11,
    -16300.9627, 1.089E-11, 5.0703E-11, -15763.8531, 1.0957E-11, 4.8841E-11,
    -15226.7436, 1.1025E-11, 4.698E-11, -14689.6341, 1.1093E-11, 4.5118E-11,
    -14152.5245, 1.1161E-11, 4.3256E-11, -13615.415, 1.1229E-11, 4.1394E-11,
    -33078.3055, 1.1297E-11, 3.9533E-11, -11587.6204, 1.1364E-11, 3.7671E-11,
    -10099.6072, 1.1432E-11, 3.5809E-11, -8616.9378, 1.15E-11, 3.3947E-11,
    -7142.2841, 1.1568E-11, 3.2086E-11, 2.0953E-11, -21134.9484, 1.8748E-10,
    2.0309E-11, -20597.8389, 1.7911E-10, 1.9665E-11, -20060.7294, 1.7075E-10,
    1.9022E-11, -19523.6198, 1.6239E-10, 1.8378E-11, -18986.5103, 1.5403E-10,
    1.7734E-11, -18449.4008, 1.4566E-10, 1.709E-11, -17912.2913, 1.373E-10,
    1.6447E-11, -17375.1817, 1.2894E-10, 1.5803E-11, -16838.0722, 1.2058E-10,
    1.5159E-11, -16300.9627, 1.1221E-10, 1.4515E-11, -15763.8531, 1.0385E-10,
    1.3872E-11, -15226.7436, 9.5489E-11, 1.3228E-11, -14689.6341, 8.7126E-11,
    1.2584E-11, -14152.5245, 7.8764E-11, 1.194E-11, -13615.415, 7.0401E-11,
    1.1297E-11, -33078.3055, 6.2039E-11, 1.0653E-11, -11587.6204, 5.3676E-11,
    1.0009E-11, -10099.6072, 4.5314E-11, 9.3653E-12, -8616.9378, 3.6951E-11,
    8.7215E-12, -7142.2841, 2.8588E-11, 9.5411E-11, 1.8834E-10, -21134.9484,
    9.1686E-11, 1.7992E-10, -20597.8389, 8.7961E-11, 1.715E-10, -20060.7294,
    8.4236E-11, 1.6308E-10, -19523.6198, 8.051E-11, 1.5466E-10, -18986.5103,
    7.6785E-11, 1.4624E-10, -18449.4008, 7.306E-11, 1.3782E-10, -17912.2913,
    6.9335E-11, 1.294E-10, -17375.1817, 6.5609E-11, 1.2098E-10, -16838.0722,
    6.1884E-11, 1.1256E-10, -16300.9627, 5.8159E-11, 1.0414E-10, -15763.8531,
    5.4434E-11, 9.5719E-11, -15226.7436, 5.0708E-11, 8.7299E-11, -14689.6341,
    4.6983E-11, 7.8879E-11, -14152.5245, 4.3258E-11, 7.0459E-11, -13615.415,
    3.9533E-11, 6.2039E-11, -33078.3055, 3.5807E-11, 5.3618E-11, -11587.6204,
    3.2082E-11, 4.5198E-11, -10099.6072, 2.8357E-11, 3.6778E-11, -8616.9378,
    2.4632E-11, 2.8358E-11, -7142.2841, -18560.0119, 9.6013E-12, 6.0581E-11,
    -18095.1858, 9.6714E-12, 5.893E-11, -17630.3597, 9.7415E-12, 5.7278E-11,
    -17165.5336, 9.8116E-12, 5.5626E-11, -16700.7075, 9.8817E-12, 5.3975E-11,
    -16235.8814, 9.9518E-12, 5.2323E-11, -15771.0553, 1.0022E-11, 5.0672E-11,
    -15306.2292, 1.0092E-11, 4.902E-11, -14841.4031, 1.0162E-11, 4.7368E-11,
    -14376.577, 1.0232E-11, 4.5717E-11, -13911.7509, 1.0302E-11, 4.4065E-11,
    -13446.9248, 1.0372E-11, 4.2414E-11, -12982.0987, 1.0443E-11, 4.0762E-11,
    -12517.2726, 1.0513E-11, 3.911E-11, -12052.4465, 1.0583E-11, 3.7459E-11,
    -11587.6204, 1.0653E-11, 3.5807E-11, -31122.7943, 1.0723E-11, 3.4156E-11,
    -9704.3927, 1.0793E-11, 3.2504E-11, -8288.6629, 1.0863E-11, 3.0853E-11,
    -6878.277, 1.0933E-11, 2.9201E-11, 2.0987E-11, -18560.0119, 1.6357E-10,
    2.0345E-11, -18095.1858, 1.5624E-10, 1.9704E-11, -17630.3597, 1.4891E-10,
    1.9062E-11, -17165.5336, 1.4158E-10, 1.8421E-11, -16700.7075, 1.3425E-10,
    1.7779E-11, -16235.8814, 1.2692E-10, 1.7138E-11, -15771.0553, 1.1959E-10,
    1.6496E-11, -15306.2292, 1.1226E-10, 1.5855E-11, -14841.4031, 1.0493E-10,
    1.5213E-11, -14376.577, 9.76E-11, 1.4572E-11, -13911.7509, 9.027E-11,
    1.393E-11, -13446.9248, 8.294E-11, 1.3289E-11, -12982.0987, 7.5609E-11,
    1.2647E-11, -12517.2726, 6.8279E-11, 1.2006E-11, -12052.4465, 6.0949E-11,
    1.1364E-11, -11587.6204, 5.3618E-11, 1.0723E-11, -31122.7943, 4.6288E-11,
    1.0081E-11, -9704.3927, 3.8958E-11, 9.4399E-12, -8288.6629, 3.1627E-11,
    8.7984E-12, -6878.277, 2.4297E-11, 9.0397E-11, 1.645E-10, -18560.0119,
    8.6882E-11, 1.5711E-10, -18095.1858, 8.3367E-11, 1.4972E-10, -17630.3597,
    7.9852E-11, 1.4233E-10, -17165.5336, 7.6337E-11, 1.3494E-10, -16700.7075,
    7.2822E-11, 1.2756E-10, -16235.8814, 6.9307E-11, 1.2017E-10, -15771.0553,
    6.5792E-11, 1.1278E-10, -15306.2292, 6.2276E-11, 1.0539E-10, -14841.4031,
    5.8761E-11, 9.8004E-11, -14376.577, 5.5246E-11, 9.0616E-11, -13911.7509,
    5.1731E-11, 8.3228E-11, -13446.9248, 4.8216E-11, 7.584E-11, -12982.0987,
    4.4701E-11, 6.8452E-11, -12517.2726, 4.1186E-11, 6.1064E-11, -12052.4465,
    3.7671E-11, 5.3676E-11, -11587.6204, 3.4156E-11, 4.6288E-11, -31122.7943,
    3.0641E-11, 3.89E-11, -9704.3927, 2.7126E-11, 3.1512E-11, -8288.6629,
    2.361E-11, 2.4124E-11, -6878.277, -16027.8254, 8.9235E-12, 5.3703E-11,
    -15632.6109, 8.9959E-12, 5.2262E-11, -15237.3963, 9.0683E-12, 5.0821E-11,
    -14842.1818, 9.1406E-12, 4.9379E-11, -14446.9673, 9.213E-12, 4.7938E-11,
    -14051.7527, 9.2854E-12, 4.6496E-11, -13656.5382, 9.3577E-12, 4.5055E-11,
    -13261.3236, 9.4301E-12, 4.3613E-11, -12866.1091, 9.5025E-12, 4.2172E-11,
    -12470.8945, 9.5748E-12, 4.0731E-11, -12075.68, 9.6472E-12, 3.9289E-11,
    -11680.4654, 9.7196E-12, 3.7848E-11, -11285.2509, 9.7919E-12, 3.6406E-11,
    -10890.0363, 9.8643E-12, 3.4965E-11, -10494.8218, 9.9367E-12, 3.3523E-11,
    -10099.6072, 1.0009E-11, 3.2082E-11, -9704.3927, 1.0081E-11, 3.0641E-11,
    -29309.1781, 1.0154E-11, 2.9199E-11, -7960.3881, 1.0226E-11, 2.7758E-11,
    -6614.2699, 1.0299E-11, 2.6316E-11, 2.1021E-11, -16027.8254, 1.3967E-10,
    2.0381E-11, -15632.6109, 1.3337E-10, 1.9742E-11, -15237.3963, 1.2707E-10,
    1.9103E-11, -14842.1818, 1.2078E-10, 1.8464E-11, -14446.9673, 1.1448E-10,
    1.7825E-11, -14051.7527, 1.0818E-10, 1.7185E-11, -13656.5382, 1.0188E-10,
    1.6546E-11, -13261.3236, 9.5583E-11, 1.5907E-11, -12866.1091, 8.9285E-11,
    1.5268E-11, -12470.8945, 8.2987E-11, 1.4628E-11, -12075.68, 7.6689E-11,
    1.3989E-11, -11680.4654, 7.0391E-11, 1.335E-11, -11285.2509, 6.4093E-11,
    1.2711E-11, -10890.0363, 5.7794E-11, 1.2071E-11, -10494.8218, 5.1496E-11,
    1.1432E-11, -10099.6072, 4.5198E-11, 1.0793E-11, -9704.3927, 3.89E-11,
    1.0154E-11, -29309.1781, 3.2602E-11, 9.5145E-12, -7960.3881, 2.6304E-11,
    8.8753E-12, -6614.2699, 2.0006E-11, 8.5383E-11, 1.4065E-10, -16027.8254,
    8.2078E-11, 1.3429E-10, -15632.6109, 7.8773E-11, 1.2794E-10, -15237.3963,
    7.5468E-11, 1.2158E-10, -14842.1818, 7.2163E-11, 1.1523E-10, -14446.9673,
    6.8858E-11, 1.0887E-10, -14051.7527, 6.5553E-11, 1.0252E-10, -13656.5382,
    6.2248E-11, 9.6159E-11, -13261.3236, 5.8944E-11, 8.9804E-11, -12866.1091,
    5.5639E-11, 8.3448E-11, -12470.8945, 5.2334E-11, 7.7092E-11, -12075.68,
    4.9029E-11, 7.0736E-11, -11680.4654, 4.5724E-11, 6.4381E-11, -11285.2509,
    4.2419E-11, 5.8025E-11, -10890.0363, 3.9114E-11, 5.1669E-11, -10494.8218,
    3.5809E-11, 4.5314E-11, -10099.6072, 3.2504E-11, 3.8958E-11, -9704.3927,
    2.9199E-11, 3.2602E-11, -29309.1781, 2.5894E-11, 2.6246E-11, -7960.3881,
    2.2589E-11, 1.9891E-11, -6614.2699, -13541.061, 8.2457E-12, 4.6826E-11,
    -13212.7861, 8.3204E-12, 4.5594E-11, -12884.5112, 8.395E-12, 4.4363E-11,
    -12556.2363, 8.4697E-12, 4.3132E-11, -12227.9615, 8.5443E-12, 4.1901E-11,
    -11899.6866, 8.6189E-12, 4.0669E-11, -11571.4117, 8.6936E-12, 3.9438E-11,
    -11243.1368, 8.7682E-12, 3.8207E-11, -10914.862, 8.8428E-12, 3.6976E-11,
    -10586.5871, 8.9175E-12, 3.5744E-11, -10258.3122, 8.9921E-12, 3.4513E-11,
    -9930.0373, 9.0667E-12, 3.3282E-11, -9601.7625, 9.1414E-12, 3.2051E-11,
    -9273.4876, 9.216E-12, 3.0819E-11, -8945.2127, 9.2906E-12, 2.9588E-11,
    -8616.9378, 9.3653E-12, 2.8357E-11, -8288.6629, 9.4399E-12, 2.7126E-11,
    -7960.3881, 9.5145E-12, 2.5894E-11, -27632.1132, 9.5892E-12, 2.4663E-11,
    -6350.2628, 9.6638E-12, 2.3432E-11, 2.1054E-11, -13541.061, 1.1577E-10,
    2.0418E-11, -13212.7861, 1.105E-10, 1.9781E-11, -12884.5112, 1.0523E-10,
    1.9144E-11, -12556.2363, 9.9969E-11, 1.8507E-11, -12227.9615, 9.4703E-11,
    1.787E-11, -11899.6866, 8.9437E-11, 1.7233E-11, -11571.4117, 8.4171E-11,
    1.6596E-11, -11243.1368, 7.8905E-11, 1.5959E-11, -10914.862, 7.3639E-11,
    1.5322E-11, -10586.5871, 6.8373E-11, 1.4685E-11, -10258.3122, 6.3107E-11,
    1.4048E-11, -9930.0373, 5.7842E-11, 1.3411E-11, -9601.7625, 5.2576E-11,
    1.2774E-11, -9273.4876, 4.731E-11, 1.2137E-11, -8945.2127, 4.2044E-11,
    1.15E-11, -8616.9378, 3.6778E-11, 1.0863E-11, -8288.6629, 3.1512E-11,
    1.0226E-11, -7960.3881, 2.6246E-11, 9.5892E-12, -27632.1132, 2.098E-11,
    8.9522E-12, -6350.2628, 1.5715E-11, 8.0369E-11, 1.168E-10, -13541.061,
    7.7274E-11, 1.1148E-10, -13212.7861, 7.4179E-11, 1.0616E-10, -12884.5112,
    7.1085E-11, 1.0083E-10, -12556.2363, 6.799E-11, 9.551E-11, -12227.9615,
    6.4895E-11, 9.0186E-11, -11899.6866, 6.18E-11, 8.4863E-11, -11571.4117,
    5.8705E-11, 7.9539E-11, -11243.1368, 5.5611E-11, 7.4216E-11, -10914.862,
    5.2516E-11, 6.8892E-11, -10586.5871, 4.9421E-11, 6.3569E-11, -10258.3122,
    4.6326E-11, 5.8245E-11, -9930.0373, 4.3232E-11, 5.2922E-11, -9601.7625,
    4.0137E-11, 4.7598E-11, -9273.4876, 3.7042E-11, 4.2275E-11, -8945.2127,
    3.3947E-11, 3.6951E-11, -8616.9378, 3.0853E-11, 3.1627E-11, -8288.6629,
    2.7758E-11, 2.6304E-11, -7960.3881, 2.4663E-11, 2.098E-11, -27632.1132,
    2.1568E-11, 1.5657E-11, -6350.2628, -11102.3904, 7.568E-12, 3.9948E-11,
    -10838.3833, 7.6449E-12, 3.8927E-11, -10574.3762, 7.7218E-12, 3.7906E-11,
    -10310.3691, 7.7987E-12, 3.6885E-11, -10046.362, 7.8756E-12, 3.5864E-11,
    -9782.3549, 7.9525E-12, 3.4843E-11, -9518.3479, 8.0294E-12, 3.3821E-11,
    -9254.3408, 8.1063E-12, 3.28E-11, -8990.3337, 8.1832E-12, 3.1779E-11,
    -8726.3266, 8.2601E-12, 3.0758E-11, -8462.3195, 8.337E-12, 2.9737E-11,
    -8198.3124, 8.4139E-12, 2.8716E-11, -7934.3053, 8.4908E-12, 2.7695E-11,
    -7670.2982, 8.5677E-12, 2.6674E-11, -7406.2911, 8.6446E-12, 2.5653E-11,
    -7142.2841, 8.7215E-12, 2.4632E-11, -6878.277, 8.7984E-12, 2.361E-11,
    -6614.2699, 8.8753E-12, 2.2589E-11, -6350.2628, 8.9522E-12, 2.1568E-11,
    -26086.2557, 9.0291E-12, 2.0547E-11, 2.1088E-11, -11102.3904, 9.1863E-11,
    2.0454E-11, -10838.3833, 8.7629E-11, 1.9819E-11, -10574.3762, 8.3395E-11,
    1.9184E-11, -10310.3691, 7.9162E-11, 1.855E-11, -10046.362, 7.4928E-11,
    1.7915E-11, -9782.3549, 7.0694E-11, 1.728E-11, -9518.3479, 6.6461E-11,
    1.6645E-11, -9254.3408, 6.2227E-11, 1.6011E-11, -8990.3337, 5.7993E-11,
    1.5376E-11, -8726.3266, 5.376E-11, 1.4741E-11, -8462.3195, 4.9526E-11,
    1.4107E-11, -8198.3124, 4.5293E-11, 1.3472E-11, -7934.3053, 4.1059E-11,
    1.2837E-11, -7670.2982, 3.6825E-11, 1.2203E-11, -7406.2911, 3.2592E-11,
    1.1568E-11, -7142.2841, 2.8358E-11, 1.0933E-11, -6878.277, 2.4124E-11,
    1.0299E-11, -6614.2699, 1.9891E-11, 9.6638E-12, -6350.2628, 1.5657E-11,
    9.0291E-12, -26086.2557, 1.1423E-11, 7.5355E-11, 9.2958E-11, -11102.3904,
    7.247E-11, 8.8667E-11, -10838.3833, 6.9585E-11, 8.4375E-11, -10574.3762,
    6.6701E-11, 8.0084E-11, -10310.3691, 6.3816E-11, 7.5793E-11, -10046.362,
    6.0932E-11, 7.1501E-11, -9782.3549, 5.8047E-11, 6.721E-11, -9518.3479,
    5.5162E-11, 6.2919E-11, -9254.3408, 5.2278E-11, 5.8627E-11, -8990.3337,
    4.9393E-11, 5.4336E-11, -8726.3266, 4.6509E-11, 5.0045E-11, -8462.3195,
    4.3624E-11, 4.5754E-11, -8198.3124, 4.0739E-11, 4.1462E-11, -7934.3053,
    3.7855E-11, 3.7171E-11, -7670.2982, 3.497E-11, 3.288E-11, -7406.2911,
    3.2086E-11, 2.8588E-11, -7142.2841, 2.9201E-11, 2.4297E-11, -6878.277,
    2.6316E-11, 2.0006E-11, -6614.2699, 2.3432E-11, 1.5715E-11, -6350.2628,
    2.0547E-11, 1.1423E-11, -26086.2557 };

  static const real_T dv[60] = { 0.073259, 2.0176E-16, -0.073259, 0.055425,
    1.6355E-16, -0.055425, 0.041191, 1.3311E-16, -0.041191, 0.029862, 1.0894E-16,
    -0.029862, 0.020873, 8.982E-17, -0.020873, 0.013771, 7.4759E-17, -0.013771,
    0.0081856, 6.2947E-17, -0.0081856, 0.0038209, 5.3715E-17, -0.0038209,
    0.0004359, 4.6514E-17, -0.0004359, -0.0021638, 4.0876E-17, 0.0021638,
    -0.004135, 3.6419E-17, 0.004135, -0.0056045, 3.2791E-17, 0.0056045,
    -0.0066744, 2.9685E-17, 0.0066744, -0.0074272, 2.6814E-17, 0.0074272,
    -0.0079291, 2.3899E-17, 0.0079291, -0.0082337, 2.0651E-17, 0.0082337,
    -0.0083839, 1.6767E-17, 0.0083839, -0.0084142, 1.1909E-17, 0.0084142,
    -0.0083523, 5.6883E-18, 0.0083523, -0.0082203, -2.3484E-18, 0.0082203 };

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
  //  P_MPC=[1567.5765 1.5888e-12 2.7019e-12 8536.6551 1.0896e-11 3.5199e-11 ;1.5888e-12 1567.5765 9.8891e-12 6.7324e-11 8536.6551 6.1982e-11 ;2.7019e-12 9.8891e-12 1567.5765 1.4677e-11 9.7557e-11 8536.6551 ;8536.6551 6.7324e-11 1.4677e-11 423802.2092 1.1408e-09 6.5756e-10 ;1.0896e-11 8536.6551 9.7557e-11 1.1408e-09 423802.2092 6.8455e-10 ;3.5199e-11 6.1982e-11 8536.6551 6.5756e-10 6.8455e-10 423802.2092 ;]; 
  //  R_MPC=[100 0 0 ;0 100 0 ;0 0 100];
  //  V_TRMPC=[0.25752 ;0.25752 ;0.25752 ;0.25752 ;0.25752 ;0.25752];
  //  Hp = 20;
  //   X_QP=[0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ]; 
  //  X_QP = zeros(3*Hp,1);
  //  solution matrices x_traj = A_tilde x_0 + B_tilde u_traj
  //  n=6;
  //  n = size(A,1);
  //  m=3;
  //  m= size(B,2);
  //  A_tilde = zeros((Hp)*n,n);
  //  A_tilde = [1 0 0 2 0 0 ;0 1 0 0 2 0 ;0 0 1 0 0 2 ;0 0 0 1 0 0 ;0 0 0 0 1 0 ;0 0 0 0 0 1 ;1 0 0 4 0 0 ;0 1 0 0 4 0 ;0 0 1 0 0 4 ;0 0 0 1 0 0 ;0 0 0 0 1 0 ;0 0 0 0 0 1 ;1 0 0 6 0 0 ;0 1 0 0 6 0 ;0 0 1 0 0 6 ;0 0 0 1 0 0 ;0 0 0 0 1 0 ;0 0 0 0 0 1 ;1 0 0 8 0 0 ;0 1 0 0 8 0 ;0 0 1 0 0 8 ;0 0 0 1 0 0 ;0 0 0 0 1 0 ;0 0 0 0 0 1 ;1 0 0 10 0 0 ;0 1 0 0 10 0 ;0 0 1 0 0 10 ;0 0 0 1 0 0 ;0 0 0 0 1 0 ;0 0 0 0 0 1 ;1 0 0 12 0 0 ;0 1 0 0 12 0 ;0 0 1 0 0 12 ;0 0 0 1 0 0 ;0 0 0 0 1 0 ;0 0 0 0 0 1 ;1 0 0 14 0 0 ;0 1 0 0 14 0 ;0 0 1 0 0 14 ;0 0 0 1 0 0 ;0 0 0 0 1 0 ;0 0 0 0 0 1 ;1 0 0 16 0 0 ;0 1 0 0 16 0 ;0 0 1 0 0 16 ;0 0 0 1 0 0 ;0 0 0 0 1 0 ;0 0 0 0 0 1 ;1 0 0 18 0 0 ;0 1 0 0 18 0 ;0 0 1 0 0 18 ;0 0 0 1 0 0 ;0 0 0 0 1 0 ;0 0 0 0 0 1 ;1 0 0 20 0 0 ;0 1 0 0 20 0 ;0 0 1 0 0 20 ;0 0 0 1 0 0 ;0 0 0 0 1 0 ;0 0 0 0 0 1 ;1 0 0 22 0 0 ;0 1 0 0 22 0 ;0 0 1 0 0 22 ;0 0 0 1 0 0 ;0 0 0 0 1 0 ;0 0 0 0 0 1 ;1 0 0 24 0 0 ;0 1 0 0 24 0 ;0 0 1 0 0 24 ;0 0 0 1 0 0 ;0 0 0 0 1 0 ;0 0 0 0 0 1 ;1 0 0 26 0 0 ;0 1 0 0 26 0 ;0 0 1 0 0 26 ;0 0 0 1 0 0 ;0 0 0 0 1 0 ;0 0 0 0 0 1 ;1 0 0 28 0 0 ;0 1 0 0 28 0 ;0 0 1 0 0 28 ;0 0 0 1 0 0 ;0 0 0 0 1 0 ;0 0 0 0 0 1 ;1 0 0 30 0 0 ;0 1 0 0 30 0 ;0 0 1 0 0 30 ;0 0 0 1 0 0 ;0 0 0 0 1 0 ;0 0 0 0 0 1 ;1 0 0 32 0 0 ;0 1 0 0 32 0 ;0 0 1 0 0 32 ;0 0 0 1 0 0 ;0 0 0 0 1 0 ;0 0 0 0 0 1 ;1 0 0 34 0 0 ;0 1 0 0 34 0 ;0 0 1 0 0 34 ;0 0 0 1 0 0 ;0 0 0 0 1 0 ;0 0 0 0 0 1 ;1 0 0 36 0 0 ;0 1 0 0 36 0 ;0 0 1 0 0 36 ;0 0 0 1 0 0 ;0 0 0 0 1 0 ;0 0 0 0 0 1 ;1 0 0 38 0 0 ;0 1 0 0 38 0 ;0 0 1 0 0 38 ;0 0 0 1 0 0 ;0 0 0 0 1 0 ;0 0 0 0 0 1 ;1 0 0 40 0 0 ;0 1 0 0 40 0 ;0 0 1 0 0 40 ;0 0 0 1 0 0 ;0 0 0 0 1 0 ;0 0 0 0 0 1 ;]; 
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
  //  B_tilde=[0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;2.4802 0 0 2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 2.4802 0 0 2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 2.4802 0 0 2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;2.6786 0 0 2.4802 0 0 2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 2.6786 0 0 2.4802 0 0 2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 2.6786 0 0 2.4802 0 0 2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;2.877 0 0 2.6786 0 0 2.4802 0 0 2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 2.877 0 0 2.6786 0 0 2.4802 0 0 2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 2.877 0 0 2.6786 0 0 2.4802 0 0 2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;3.0754 0 0 2.877 0 0 2.6786 0 0 2.4802 0 0 2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 3.0754 0 0 2.877 0 0 2.6786 0 0 2.4802 0 0 2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 3.0754 0 0 2.877 0 0 2.6786 0 0 2.4802 0 0 2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 ;0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 0 ;3.2738 0 0 3.0754 0 0 2.877 0 0 2.6786 0 0 2.4802 0 0 2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 ;0 3.2738 0 0 3.0754 0 0 2.877 0 0 2.6786 0 0 2.4802 0 0 2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 ;0 0 3.2738 0 0 3.0754 0 0 2.877 0 0 2.6786 0 0 2.4802 0 0 2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 0 ;0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 0 ;0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 0 ;0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 0 ;3.4722 0 0 3.2738 0 0 3.0754 0 0 2.877 0 0 2.6786 0 0 2.4802 0 0 2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 0 ;0 3.4722 0 0 3.2738 0 0 3.0754 0 0 2.877 0 0 2.6786 0 0 2.4802 0 0 2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 0 ;0 0 3.4722 0 0 3.2738 0 0 3.0754 0 0 2.877 0 0 2.6786 0 0 2.4802 0 0 2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 0 ;0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 0 ;0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 0 ;0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 0 ;3.6706 0 0 3.4722 0 0 3.2738 0 0 3.0754 0 0 2.877 0 0 2.6786 0 0 2.4802 0 0 2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 0 ;0 3.6706 0 0 3.4722 0 0 3.2738 0 0 3.0754 0 0 2.877 0 0 2.6786 0 0 2.4802 0 0 2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 0 ;0 0 3.6706 0 0 3.4722 0 0 3.2738 0 0 3.0754 0 0 2.877 0 0 2.6786 0 0 2.4802 0 0 2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 0 ;0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 0 ;0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 0 ;0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0 ;3.869 0 0 3.6706 0 0 3.4722 0 0 3.2738 0 0 3.0754 0 0 2.877 0 0 2.6786 0 0 2.4802 0 0 2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 0 ;0 3.869 0 0 3.6706 0 0 3.4722 0 0 3.2738 0 0 3.0754 0 0 2.877 0 0 2.6786 0 0 2.4802 0 0 2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 0 ;0 0 3.869 0 0 3.6706 0 0 3.4722 0 0 3.2738 0 0 3.0754 0 0 2.877 0 0 2.6786 0 0 2.4802 0 0 2.2817 0 0 2.0833 0 0 1.8849 0 0 1.6865 0 0 1.4881 0 0 1.2897 0 0 1.0913 0 0 0.89285 0 0 0.69444 0 0 0.49603 0 0 0.29762 0 0 0.099206 ;0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 ;0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 ;0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 0 0 0.099206 ;]; 
  //  cost
  //  Q_tilde = blkdiag(kron(eye(Hp-1),Q_MPC),P_MPC);
  //  Q_tilde = [70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1567.5765 1.5888e-12 2.7019e-12 8536.6551 1.0896e-11 3.5199e-11 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.5888e-12 1567.5765 9.8891e-12 6.7324e-11 8536.6551 6.1982e-11 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2.7019e-12 9.8891e-12 1567.5765 1.4677e-11 9.7557e-11 8536.6551 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8536.6551 6.7324e-11 1.4677e-11 423802.2092 1.1408e-09 6.5756e-10 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.0896e-11 8536.6551 9.7557e-11 1.1408e-09 423802.2092 6.8455e-10 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3.5199e-11 6.1982e-11 8536.6551 6.5756e-10 6.8455e-10 423802.2092 ;]; 
  //  R_tilde = kron(eye(Hp),R_MPC);
  //  R_tilde=[100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 ]; 
  //  F_tilde = kron([ones(m-1,1);zeros(m*(Hp-1)+1,1)],F_MPC);
  //  L3 = blkdiag(kron(eye(Hp-1),Q_MPC),zeros(n,n));
  //  L3=[70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100000 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;]; 
  //  L4 = zeros(n,Hp*n);
  //  L4(:,end-(n-1):end) = eye(n);
  //  L4 = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 ;]; 
  //  PHI = blkdiag(kron(eye(Hp),eye(6,6)));
  //  PHI = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 ;]; 
  //  constraints
  //  CONTROL INPUTS
  //  Define and Solve QP
  //  GAMMA = zeros(n*Hp,1);%Dimension of GAMMA is n*horizon x n
  //  for i = 1:Hp
  //      GAMMA(i*n-5:i*n,1) = xfinal;
  //  end
  //  GAMMA = [0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;0 ;]; 
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
  //  L = [-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 ;1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 ;]; 
  //  b = [kron(ones(Hp,1), V_TRMPC(1:m,1)); kron(ones(Hp,1), V_TRMPC(m+1:end,1))]; 
  //  b=[0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;0.22901 ;]; 
  //  Call QP Solver
  //  Call QP Solver
  //  QQ = (QQ+QQ')/2;
  //  QQ = [138371.1425 1.3007e-10 1.3212e-10 112666.5366 1.2498e-10 1.274e-10 106967.442 1.1989e-10 1.2267e-10 101279.3704 1.148e-10 1.1794e-10 95607.833 1.0971e-10 1.1322e-10 89958.3413 1.0462e-10 1.0849e-10 84336.4067 9.953e-11 1.0377e-10 78747.5407 9.4441e-11 9.904e-11 73197.2546 8.9351e-11 9.4314e-11 67691.06 8.4261e-11 8.9588e-11 62234.4682 7.9172e-11 8.4862e-11 56832.9906 7.4082e-11 8.0136e-11 51492.1387 6.8992e-11 7.541e-11 46217.4239 6.3903e-11 7.0684e-11 41014.3576 5.8813e-11 6.5958e-11 35888.4513 5.3723e-11 6.1232e-11 30845.2164 4.8633e-11 5.6505e-11 25890.1642 4.3544e-11 5.1779e-11 21028.8063 3.8454e-11 4.7053e-11 16266.654 3.3364e-11 4.2327e-11 ;1.3007e-10 138371.1425 4.3201e-10 1.272e-10 112666.5366 4.1299e-10 1.2433e-10 106967.442 3.9397e-10 1.2146e-10 101279.3704 3.7494e-10 1.186e-10 95607.833 3.5592e-10 1.1573e-10 89958.3413 3.369e-10 1.1286e-10 84336.4067 3.1787e-10 1.0999e-10 78747.5407 2.9885e-10 1.0712e-10 73197.2546 2.7983e-10 1.0425e-10 67691.06 2.608e-10 1.0139e-10 62234.4682 2.4178e-10 9.8518e-11 56832.9906 2.2276e-10 9.5649e-11 51492.1387 2.0373e-10 9.2781e-11 46217.4239 1.8471e-10 8.9913e-11 41014.3576 1.6569e-10 8.7044e-11 35888.4513 1.4666e-10 8.4176e-11 30845.2164 1.2764e-10 8.1308e-11 25890.1642 1.0861e-10 7.844e-11 21028.8063 8.9591e-11 7.5571e-11 16266.654 7.0567e-11 ;1.3212e-10 4.3201e-10 138371.1425 1.2659e-10 4.1439e-10 112666.5366 1.2105e-10 3.9677e-10 106967.442 1.1552e-10 3.7915e-10 101279.3704 1.0999e-10 3.6152e-10 95607.833 1.0445e-10 3.439e-10 89958.3413 9.8919e-11 3.2628e-10 84336.4067 9.3385e-11 3.0865e-10 78747.5407 8.7851e-11 2.9103e-10 73197.2546 8.2317e-11 2.7341e-10 67691.06 7.6783e-11 2.5578e-10 62234.4682 7.1249e-11 2.3816e-10 56832.9906 6.5715e-11 2.2054e-10 51492.1387 6.0181e-11 2.0292e-10 46217.4239 5.4647e-11 1.8529e-10 41014.3576 4.9113e-11 1.6767e-10 35888.4513 4.3579e-11 1.5005e-10 30845.2164 3.8045e-11 1.3242e-10 25890.1642 3.2511e-11 1.148e-10 21028.8063 2.6977e-11 9.7177e-11 16266.654 ;112666.5366 1.272e-10 1.2659e-10 129154.3028 1.2224e-10 1.2208e-10 103675.0808 1.1727e-10 1.1756e-10 98201.3702 1.1231e-10 1.1305e-10 92738.6825 1.0734e-10 1.0854e-10 87292.529 1.0238e-10 1.0402e-10 81868.4212 9.7413e-11 9.9509e-11 76471.8706 9.2448e-11 9.4995e-11 71108.3885 8.7483e-11 9.0482e-11 65783.4864 8.2519e-11 8.5969e-11 60502.6757 7.7554e-11 8.1455e-11 55271.4678 7.259e-11 7.6942e-11 50095.3742 6.7625e-11 7.2429e-11 44979.9063 6.2661e-11 6.7915e-11 39930.5754 5.7696e-11 6.3402e-11 34952.8931 5.2731e-11 5.8889e-11 30052.3707 4.7767e-11 5.4375e-11 25234.5197 4.2802e-11 4.9862e-11 20504.8515 3.7838e-11 4.5348e-11 15868.8775 3.2873e-11 4.0835e-11 ;1.2498e-10 112666.5366 4.1439e-10 1.2224e-10 129154.3028 3.9615e-10 1.1949e-10 103675.0808 3.779e-10 1.1675e-10 98201.3702 3.5966e-10 1.1401e-10 92738.6825 3.4141e-10 1.1126e-10 87292.529 3.2317e-10 1.0852e-10 81868.4212 3.0492e-10 1.0578e-10 76471.8706 2.8668e-10 1.0303e-10 71108.3885 2.6843e-10 1.0029e-10 65783.4864 2.5019e-10 9.7547e-11 60502.6757 2.3194e-10 9.4804e-11 55271.4678 2.137e-10 9.2061e-11 50095.3742 1.9545e-10 8.9318e-11 44979.9063 1.7721e-10 8.6574e-11 39930.5754 1.5896e-10 8.3831e-11 34952.8931 1.4072e-10 8.1088e-11 30052.3707 1.2247e-10 7.8345e-11 25234.5197 1.0423e-10 7.5602e-11 20504.8515 8.5983e-11 7.2858e-11 15868.8775 6.7738e-11 ;1.274e-10 4.1299e-10 112666.5366 1.2208e-10 3.9615e-10 129154.3028 1.1675e-10 3.793e-10 103675.0808 1.1143e-10 3.6246e-10 98201.3702 1.0611e-10 3.4561e-10 92738.6825 1.0079e-10 3.2877e-10 87292.529 9.5469e-11 3.1192e-10 81868.4212 9.0148e-11 2.9508e-10 76471.8706 8.4827e-11 2.7824e-10 71108.3885 7.9506e-11 2.6139e-10 65783.4864 7.4184e-11 2.4455e-10 60502.6757 6.8863e-11 2.277e-10 55271.4678 6.3542e-11 2.1086e-10 50095.3742 5.8221e-11 1.9401e-10 44979.9063 5.2899e-11 1.7717e-10 39930.5754 4.7578e-11 1.6032e-10 34952.8931 4.2257e-11 1.4348e-10 30052.3707 3.6936e-11 1.2664e-10 25234.5197 3.1614e-11 1.0979e-10 20504.8515 2.6293e-11 9.2947e-11 15868.8775 ;106967.442 1.2433e-10 1.2105e-10 103675.0808 1.1949e-10 1.1675e-10 120382.7195 1.1465e-10 1.1245e-10 95123.37 1.0981e-10 1.0815e-10 89869.532 1.0497e-10 1.0385e-10 84626.7167 1.0013e-10 9.9552e-11 79400.4358 9.5295e-11 9.5251e-11 74196.2006 9.0455e-11 9.0951e-11 69019.5224 8.5616e-11 8.665e-11 63875.9129 8.0776e-11 8.2349e-11 58770.8833 7.5937e-11 7.8049e-11 53709.9451 7.1097e-11 7.3748e-11 48698.6097 6.6258e-11 6.9447e-11 43742.3886 6.1418e-11 6.5147e-11 38846.7932 5.6579e-11 6.0846e-11 34017.3348 5.1739e-11 5.6546e-11 29259.525 4.69e-11 5.2245e-11 24578.8751 4.206e-11 4.7944e-11 19980.8967 3.7221e-11 4.3644e-11 15471.101 3.2381e-11 3.9343e-11 ;1.1989e-10 106967.442 3.9677e-10 1.1727e-10 103675.0808 3.793e-10 1.1465e-10 120382.7195 3.6184e-10 1.1203e-10 95123.37 3.4437e-10 1.0942e-10 89869.532 3.269e-10 1.068e-10 84626.7167 3.0944e-10 1.0418e-10 79400.4358 2.9197e-10 1.0156e-10 74196.2006 2.745e-10 9.8944e-11 69019.5224 2.5704e-10 9.6326e-11 63875.9129 2.3957e-10 9.3708e-11 58770.8833 2.2211e-10 9.109e-11 53709.9451 2.0464e-10 8.8472e-11 48698.6097 1.8717e-10 8.5854e-11 43742.3886 1.6971e-10 8.3236e-11 38846.7932 1.5224e-10 8.0618e-11 34017.3348 1.3477e-10 7.8e-11 29259.525 1.1731e-10 7.5382e-11 24578.8751 9.9841e-11 7.2764e-11 19980.8967 8.2375e-11 7.0146e-11 15471.101 6.4909e-11 ;1.2267e-10 3.9397e-10 106967.442 1.1756e-10 3.779e-10 103675.0808 1.1245e-10 3.6184e-10 120382.7195 1.0735e-10 3.4577e-10 95123.37 1.0224e-10 3.297e-10 89869.532 9.7128e-11 3.1364e-10 84626.7167 9.202e-11 2.9757e-10 79400.4358 8.6911e-11 2.8151e-10 74196.2006 8.1803e-11 2.6544e-10 69019.5224 7.6694e-11 2.4938e-10 63875.9129 7.1586e-11 2.3331e-10 58770.8833 6.6477e-11 2.1724e-10 53709.9451 6.1368e-11 2.0118e-10 48698.6097 5.626e-11 1.8511e-10 43742.3886 5.1151e-11 1.6905e-10 38846.7932 4.6043e-11 1.5298e-10 34017.3348 4.0934e-11 1.3691e-10 29259.525 3.5826e-11 1.2085e-10 24578.8751 3.0717e-11 1.0478e-10 19980.8967 2.5609e-11 8.8717e-11 15471.101 ;101279.3704 1.2146e-10 1.1552e-10 98201.3702 1.1675e-10 1.1143e-10 95123.37 1.1203e-10 1.0735e-10 112045.3699 1.0732e-10 1.0326e-10 87000.3814 1.0261e-10 9.9169e-11 81960.9045 9.7892e-11 9.5081e-11 76932.4503 9.3177e-11 9.0994e-11 71920.5305 8.8463e-11 8.6906e-11 66930.6563 8.3748e-11 8.2818e-11 61968.3393 7.9034e-11 7.873e-11 57039.0908 7.432e-11 7.4642e-11 52148.4223 6.9605e-11 7.0554e-11 47301.8452 6.4891e-11 6.6466e-11 42504.871 6.0176e-11 6.2378e-11 37763.0109 5.5462e-11 5.829e-11 33081.7766 5.0748e-11 5.4203e-11 28466.6793 4.6033e-11 5.0115e-11 23923.2306 4.1319e-11 4.6027e-11 19456.9418 3.6604e-11 4.1939e-11 15073.3244 3.189e-11 3.7851e-11 ;1.148e-10 101279.3704 3.7915e-10 1.1231e-10 98201.3702 3.6246e-10 1.0981e-10 95123.37 3.4577e-10 1.0732e-10 112045.3699 3.2908e-10 1.0483e-10 87000.3814 3.1239e-10 1.0233e-10 81960.9045 2.9571e-10 9.9842e-11 76932.4503 2.7902e-10 9.7349e-11 71920.5305 2.6233e-10 9.4856e-11 66930.6563 2.4564e-10 9.2363e-11 61968.3393 2.2896e-10 8.987e-11 57039.0908 2.1227e-10 8.7377e-11 52148.4223 1.9558e-10 8.4884e-11 47301.8452 1.7889e-10 8.2391e-11 42504.871 1.6221e-10 7.9898e-11 37763.0109 1.4552e-10 7.7405e-11 33081.7766 1.2883e-10 7.4912e-11 28466.6793 1.1214e-10 7.2419e-11 23923.2306 9.5455e-11 6.9926e-11 19456.9418 7.8767e-11 6.7433e-11 15073.3244 6.2079e-11 ;1.1794e-10 3.7494e-10 101279.3704 1.1305e-10 3.5966e-10 98201.3702 1.0815e-10 3.4437e-10 95123.37 1.0326e-10 3.2908e-10 112045.3699 9.8361e-11 3.138e-10 87000.3814 9.3466e-11 2.9851e-10 81960.9045 8.857e-11 2.8322e-10 76932.4503 8.3674e-11 2.6793e-10 71920.5305 7.8778e-11 2.5265e-10 66930.6563 7.3883e-11 2.3736e-10 61968.3393 6.8987e-11 2.2207e-10 57039.0908 6.4091e-11 2.0678e-10 52148.4223 5.9195e-11 1.915e-10 47301.8452 5.4299e-11 1.7621e-10 42504.871 4.9404e-11 1.6092e-10 37763.0109 4.4508e-11 1.4564e-10 33081.7766 3.9612e-11 1.3035e-10 28466.6793 3.4716e-11 1.1506e-10 23923.2306 2.982e-11 9.9774e-11 19456.9418 2.4925e-11 8.4487e-11 15073.3244 ;95607.833 1.186e-10 1.0999e-10 92738.6825 1.1401e-10 1.0611e-10 89869.532 1.0942e-10 1.0224e-10 87000.3814 1.0483e-10 9.8361e-11 104131.2309 1.0024e-10 9.4486e-11 79295.0922 9.5649e-11 9.0611e-11 74464.4649 9.106e-11 8.6736e-11 69644.8604 8.647e-11 8.2861e-11 64841.7902 8.1881e-11 7.8986e-11 60060.7658 7.7292e-11 7.5111e-11 55307.2984 7.2702e-11 7.1235e-11 50586.8996 6.8113e-11 6.736e-11 45905.0807 6.3524e-11 6.3485e-11 41267.3533 5.8934e-11 5.961e-11 36679.2287 5.4345e-11 5.5735e-11 32146.2183 4.9756e-11 5.186e-11 27673.8337 4.5166e-11 4.7984e-11 23267.5861 4.0577e-11 4.4109e-11 18932.987 3.5988e-11 4.0234e-11 14675.5479 3.1398e-11 3.6359e-11 ;1.0971e-10 95607.833 3.6152e-10 1.0734e-10 92738.6825 3.4561e-10 1.0497e-10 89869.532 3.297e-10 1.0261e-10 87000.3814 3.138e-10 1.0024e-10 104131.2309 2.9789e-10 9.787e-11 79295.0922 2.8198e-10 9.5502e-11 74464.4649 2.6607e-10 9.3135e-11 69644.8604 2.5016e-10 9.0767e-11 64841.7902 2.3425e-10 8.8399e-11 60060.7658 2.1834e-10 8.6031e-11 55307.2984 2.0243e-10 8.3663e-11 50586.8996 1.8652e-10 8.1295e-11 45905.0807 1.7061e-10 7.8927e-11 41267.3533 1.547e-10 7.6559e-11 36679.2287 1.388e-10 7.4191e-11 32146.2183 1.2289e-10 7.1823e-11 27673.8337 1.0698e-10 6.9456e-11 23267.5861 9.1068e-11 6.7088e-11 18932.987 7.5159e-11 6.472e-11 14675.5479 5.925e-11 ;1.1322e-10 3.5592e-10 95607.833 1.0854e-10 3.4141e-10 92738.6825 1.0385e-10 3.269e-10 89869.532 9.9169e-11 3.1239e-10 87000.3814 9.4486e-11 2.9789e-10 104131.2309 8.9803e-11 2.8338e-10 79295.0922 8.512e-11 2.6887e-10 74464.4649 8.0437e-11 2.5436e-10 69644.8604 7.5754e-11 2.3985e-10 64841.7902 7.1071e-11 2.2534e-10 60060.7658 6.6388e-11 2.1083e-10 55307.2984 6.1705e-11 1.9633e-10 50586.8996 5.7022e-11 1.8182e-10 45905.0807 5.2339e-11 1.6731e-10 41267.3533 4.7656e-11 1.528e-10 36679.2287 4.2973e-11 1.3829e-10 32146.2183 3.829e-11 1.2378e-10 27673.8337 3.3607e-11 1.0927e-10 23267.5861 2.8924e-11 9.4766e-11 18932.987 2.4241e-11 8.0257e-11 14675.5479 ;89958.3413 1.1573e-10 1.0445e-10 87292.529 1.1126e-10 1.0079e-10 84626.7167 1.068e-10 9.7128e-11 81960.9045 1.0233e-10 9.3466e-11 79295.0922 9.787e-11 8.9803e-11 96629.2799 9.3406e-11 8.6141e-11 71996.4794 8.8942e-11 8.2478e-11 67369.1904 8.4478e-11 7.8816e-11 62752.9241 8.0013e-11 7.5154e-11 58153.1922 7.5549e-11 7.1491e-11 53575.5059 7.1085e-11 6.7829e-11 49025.3768 6.6621e-11 6.4166e-11 44508.3163 6.2157e-11 6.0504e-11 40029.8357 5.7692e-11 5.6841e-11 35595.4465 5.3228e-11 5.3179e-11 31210.6601 4.8764e-11 4.9517e-11 26880.988 4.43e-11 4.5854e-11 22611.9415 3.9835e-11 4.2192e-11 18409.0322 3.5371e-11 3.8529e-11 14277.7714 3.0907e-11 3.4867e-11 ;1.0462e-10 89958.3413 3.439e-10 1.0238e-10 87292.529 3.2877e-10 1.0013e-10 84626.7167 3.1364e-10 9.7892e-11 81960.9045 2.9851e-10 9.5649e-11 79295.0922 2.8338e-10 9.3406e-11 96629.2799 2.6825e-10 9.1163e-11 71996.4794 2.5312e-10 8.892e-11 67369.1904 2.3799e-10 8.6678e-11 62752.9241 2.2286e-10 8.4435e-11 58153.1922 2.0773e-10 8.2192e-11 53575.5059 1.9259e-10 7.9949e-11 49025.3768 1.7746e-10 7.7706e-11 44508.3163 1.6233e-10 7.5464e-11 40029.8357 1.472e-10 7.3221e-11 35595.4465 1.3207e-10 7.0978e-11 31210.6601 1.1694e-10 6.8735e-11 26880.988 1.0181e-10 6.6492e-11 22611.9415 8.6682e-11 6.425e-11 18409.0322 7.1551e-11 6.2007e-11 14277.7714 5.6421e-11 ;1.0849e-10 3.369e-10 89958.3413 1.0402e-10 3.2317e-10 87292.529 9.9552e-11 3.0944e-10 84626.7167 9.5081e-11 2.9571e-10 81960.9045 9.0611e-11 2.8198e-10 79295.0922 8.6141e-11 2.6825e-10 96629.2799 8.1671e-11 2.5452e-10 71996.4794 7.72e-11 2.4079e-10 67369.1904 7.273e-11 2.2706e-10 62752.9241 6.826e-11 2.1333e-10 58153.1922 6.3789e-11 1.996e-10 53575.5059 5.9319e-11 1.8587e-10 49025.3768 5.4849e-11 1.7214e-10 44508.3163 5.0378e-11 1.5841e-10 40029.8357 4.5908e-11 1.4468e-10 35595.4465 4.1438e-11 1.3095e-10 31210.6601 3.6967e-11 1.1722e-10 26880.988 3.2497e-11 1.0349e-10 22611.9415 2.8027e-11 8.9758e-11 18409.0322 2.3556e-11 7.6028e-11 14277.7714 ;84336.4067 1.1286e-10 9.8919e-11 81868.4212 1.0852e-10 9.5469e-11 79400.4358 1.0418e-10 9.202e-11 76932.4503 9.9842e-11 8.857e-11 74464.4649 9.5502e-11 8.512e-11 71996.4794 9.1163e-11 8.1671e-11 89528.494 8.6824e-11 7.8221e-11 65093.5203 8.2485e-11 7.4771e-11 60664.058 7.8146e-11 7.1321e-11 56245.6186 7.3807e-11 6.7872e-11 51843.7135 6.9468e-11 6.4422e-11 47463.8541 6.5128e-11 6.0972e-11 43111.5518 6.0789e-11 5.7523e-11 38792.318 5.645e-11 5.4073e-11 34511.6642 5.2111e-11 5.0623e-11 30275.1018 4.7772e-11 4.7174e-11 26088.1423 4.3433e-11 4.3724e-11 21956.297 3.9094e-11 4.0274e-11 17885.0774 3.4755e-11 3.6825e-11 13879.9948 3.0415e-11 3.3375e-11 ;9.953e-11 84336.4067 3.2628e-10 9.7413e-11 81868.4212 3.1192e-10 9.5295e-11 79400.4358 2.9757e-10 9.3177e-11 76932.4503 2.8322e-10 9.106e-11 74464.4649 2.6887e-10 8.8942e-11 71996.4794 2.5452e-10 8.6824e-11 89528.494 2.4017e-10 8.4706e-11 65093.5203 2.2581e-10 8.2589e-11 60664.058 2.1146e-10 8.0471e-11 56245.6186 1.9711e-10 7.8353e-11 51843.7135 1.8276e-10 7.6236e-11 47463.8541 1.6841e-10 7.4118e-11 43111.5518 1.5405e-10 7.2e-11 38792.318 1.397e-10 6.9882e-11 34511.6642 1.2535e-10 6.7765e-11 30275.1018 1.11e-10 6.5647e-11 26088.1423 9.6647e-11 6.3529e-11 21956.297 8.2295e-11 6.1412e-11 17885.0774 6.7943e-11 5.9294e-11 13879.9948 5.3591e-11 ;1.0377e-10 3.1787e-10 84336.4067 9.9509e-11 3.0492e-10 81868.4212 9.5251e-11 2.9197e-10 79400.4358 9.0994e-11 2.7902e-10 76932.4503 8.6736e-11 2.6607e-10 74464.4649 8.2478e-11 2.5312e-10 71996.4794 7.8221e-11 2.4017e-10 89528.494 7.3963e-11 2.2721e-10 65093.5203 6.9706e-11 2.1426e-10 60664.058 6.5448e-11 2.0131e-10 56245.6186 6.119e-11 1.8836e-10 51843.7135 5.6933e-11 1.7541e-10 47463.8541 5.2675e-11 1.6246e-10 43111.5518 4.8418e-11 1.4951e-10 38792.318 4.416e-11 1.3655e-10 34511.6642 3.9903e-11 1.236e-10 30275.1018 3.5645e-11 1.1065e-10 26088.1423 3.1387e-11 9.77e-11 21956.297 2.713e-11 8.4749e-11 17885.0774 2.2872e-11 7.1798e-11 13879.9948 ;78747.5407 1.0999e-10 9.3385e-11 76471.8706 1.0578e-10 9.0148e-11 74196.2006 1.0156e-10 8.6911e-11 71920.5305 9.7349e-11 8.3674e-11 69644.8604 9.3135e-11 8.0437e-11 67369.1904 8.892e-11 7.72e-11 65093.5203 8.4706e-11 7.3963e-11 82817.8502 8.0492e-11 7.0726e-11 58575.192 7.6278e-11 6.7489e-11 54338.0451 7.2064e-11 6.4252e-11 50111.9211 6.785e-11 6.1015e-11 45902.3313 6.3636e-11 5.7778e-11 41714.7873 5.9422e-11 5.4542e-11 37554.8004 5.5208e-11 5.1305e-11 33427.882 5.0994e-11 4.8068e-11 29339.5436 4.678e-11 4.4831e-11 25295.2966 4.2566e-11 4.1594e-11 21300.6524 3.8352e-11 3.8357e-11 17361.1225 3.4138e-11 3.512e-11 13482.2183 2.9924e-11 3.1883e-11 ;9.4441e-11 78747.5407 3.0865e-10 9.2448e-11 76471.8706 2.9508e-10 9.0455e-11 74196.2006 2.8151e-10 8.8463e-11 71920.5305 2.6793e-10 8.647e-11 69644.8604 2.5436e-10 8.4478e-11 67369.1904 2.4079e-10 8.2485e-11 65093.5203 2.2721e-10 8.0492e-11 82817.8502 2.1364e-10 7.85e-11 58575.192 2.0007e-10 7.6507e-11 54338.0451 1.8649e-10 7.4515e-11 50111.9211 1.7292e-10 7.2522e-11 45902.3313 1.5935e-10 7.0529e-11 41714.7873 1.4577e-10 6.8537e-11 37554.8004 1.322e-10 6.6544e-11 33427.882 1.1863e-10 6.4551e-11 29339.5436 1.0505e-10 6.2559e-11 25295.2966 9.1482e-11 6.0566e-11 21300.6524 7.7908e-11 5.8574e-11 17361.1225 6.4335e-11 5.6581e-11 13482.2183 5.0762e-11 ;9.904e-11 2.9885e-10 78747.5407 9.4995e-11 2.8668e-10 76471.8706 9.0951e-11 2.745e-10 74196.2006 8.6906e-11 2.6233e-10 71920.5305 8.2861e-11 2.5016e-10 69644.8604 7.8816e-11 2.3799e-10 67369.1904 7.4771e-11 2.2581e-10 65093.5203 7.0726e-11 2.1364e-10 82817.8502 6.6681e-11 2.0147e-10 58575.192 6.2637e-11 1.893e-10 54338.0451 5.8592e-11 1.7712e-10 50111.9211 5.4547e-11 1.6495e-10 45902.3313 5.0502e-11 1.5278e-10 41714.7873 4.6457e-11 1.406e-10 37554.8004 4.2412e-11 1.2843e-10 33427.882 3.8367e-11 1.1626e-10 29339.5436 3.4323e-11 1.0409e-10 25295.2966 3.0278e-11 9.1913e-11 21300.6524 2.6233e-11 7.9741e-11 17361.1225 2.2188e-11 6.7568e-11 13482.2183 ;73197.2546 1.0712e-10 8.7851e-11 71108.3885 1.0303e-10 8.4827e-11 69019.5224 9.8944e-11 8.1803e-11 66930.6563 9.4856e-11 7.8778e-11 64841.7902 9.0767e-11 7.5754e-11 62752.9241 8.6678e-11 7.273e-11 60664.058 8.2589e-11 6.9706e-11 58575.192 7.85e-11 6.6681e-11 76486.3259 7.4411e-11 6.3657e-11 52430.4715 7.0322e-11 6.0633e-11 48380.1286 6.6233e-11 5.7609e-11 44340.8086 6.2144e-11 5.4585e-11 40318.0228 5.8055e-11 5.156e-11 36317.2827 5.3966e-11 4.8536e-11 32344.0998 4.9877e-11 4.5512e-11 28403.9854 4.5788e-11 4.2488e-11 24502.4509 4.1699e-11 3.9463e-11 20645.0079 3.761e-11 3.6439e-11 16837.1677 3.3521e-11 3.3415e-11 13084.4418 2.9432e-11 3.0391e-11 ;8.9351e-11 73197.2546 2.9103e-10 8.7483e-11 71108.3885 2.7824e-10 8.5616e-11 69019.5224 2.6544e-10 8.3748e-11 66930.6563 2.5265e-10 8.1881e-11 64841.7902 2.3985e-10 8.0013e-11 62752.9241 2.2706e-10 7.8146e-11 60664.058 2.1426e-10 7.6278e-11 58575.192 2.0147e-10 7.4411e-11 76486.3259 1.8867e-10 7.2543e-11 52430.4715 1.7588e-10 7.0676e-11 48380.1286 1.6308e-10 6.8808e-11 44340.8086 1.5029e-10 6.6941e-11 40318.0228 1.3749e-10 6.5073e-11 36317.2827 1.247e-10 6.3206e-11 32344.0998 1.1191e-10 6.1338e-11 28403.9854 9.9111e-11 5.9471e-11 24502.4509 8.6316e-11 5.7603e-11 20645.0079 7.3522e-11 5.5736e-11 16837.1677 6.0727e-11 5.3868e-11 13084.4418 4.7933e-11 ;9.4314e-11 2.7983e-10 73197.2546 9.0482e-11 2.6843e-10 71108.3885 8.665e-11 2.5704e-10 69019.5224 8.2818e-11 2.4564e-10 66930.6563 7.8986e-11 2.3425e-10 64841.7902 7.5154e-11 2.2286e-10 62752.9241 7.1321e-11 2.1146e-10 60664.058 6.7489e-11 2.0007e-10 58575.192 6.3657e-11 1.8867e-10 76486.3259 5.9825e-11 1.7728e-10 52430.4715 5.5993e-11 1.6589e-10 48380.1286 5.2161e-11 1.5449e-10 44340.8086 4.8329e-11 1.431e-10 40318.0228 4.4497e-11 1.317e-10 36317.2827 4.0664e-11 1.2031e-10 32344.0998 3.6832e-11 1.0891e-10 28403.9854 3.3e-11 9.752e-11 24502.4509 2.9168e-11 8.6126e-11 20645.0079 2.5336e-11 7.4732e-11 16837.1677 2.1504e-11 6.3338e-11 13084.4418 ;67691.06 1.0425e-10 8.2317e-11 65783.4864 1.0029e-10 7.9506e-11 63875.9129 9.6326e-11 7.6694e-11 61968.3393 9.2363e-11 7.3883e-11 60060.7658 8.8399e-11 7.1071e-11 58153.1922 8.4435e-11 6.826e-11 56245.6186 8.0471e-11 6.5448e-11 54338.0451 7.6507e-11 6.2637e-11 52430.4715 7.2543e-11 5.9825e-11 70522.898 6.8579e-11 5.7014e-11 46648.3362 6.4616e-11 5.4202e-11 42779.2858 6.0652e-11 5.1391e-11 38921.2583 5.6688e-11 4.8579e-11 35079.7651 5.2724e-11 4.5768e-11 31260.3175 4.876e-11 4.2956e-11 27468.4271 4.4796e-11 4.0145e-11 23709.6053 4.0833e-11 3.7333e-11 19989.3634 3.6869e-11 3.4522e-11 16313.2129 3.2905e-11 3.171e-11 12686.6652 2.8941e-11 2.8899e-11 ;8.4261e-11 67691.06 2.7341e-10 8.2519e-11 65783.4864 2.6139e-10 8.0776e-11 63875.9129 2.4938e-10 7.9034e-11 61968.3393 2.3736e-10 7.7292e-11 60060.7658 2.2534e-10 7.5549e-11 58153.1922 2.1333e-10 7.3807e-11 56245.6186 2.0131e-10 7.2064e-11 54338.0451 1.893e-10 7.0322e-11 52430.4715 1.7728e-10 6.8579e-11 70522.898 1.6526e-10 6.6837e-11 46648.3362 1.5325e-10 6.5095e-11 42779.2858 1.4123e-10 6.3352e-11 38921.2583 1.2922e-10 6.161e-11 35079.7651 1.172e-10 5.9867e-11 31260.3175 1.0518e-10 5.8125e-11 27468.4271 9.3167e-11 5.6382e-11 23709.6053 8.1151e-11 5.464e-11 19989.3634 6.9135e-11 5.2898e-11 16313.2129 5.7119e-11 5.1155e-11 12686.6652 4.5103e-11 ;8.9588e-11 2.608e-10 67691.06 8.5969e-11 2.5019e-10 65783.4864 8.2349e-11 2.3957e-10 63875.9129 7.873e-11 2.2896e-10 61968.3393 7.5111e-11 2.1834e-10 60060.7658 7.1491e-11 2.0773e-10 58153.1922 6.7872e-11 1.9711e-10 56245.6186 6.4252e-11 1.8649e-10 54338.0451 6.0633e-11 1.7588e-10 52430.4715 5.7014e-11 1.6526e-10 70522.898 5.3394e-11 1.5465e-10 46648.3362 4.9775e-11 1.4403e-10 42779.2858 4.6155e-11 1.3342e-10 38921.2583 4.2536e-11 1.228e-10 35079.7651 3.8917e-11 1.1219e-10 31260.3175 3.5297e-11 1.0157e-10 27468.4271 3.1678e-11 9.0955e-11 23709.6053 2.8058e-11 8.0339e-11 19989.3634 2.4439e-11 6.9724e-11 16313.2129 2.082e-11 5.9108e-11 12686.6652 ;62234.4682 1.0139e-10 7.6783e-11 60502.6757 9.7547e-11 7.4184e-11 58770.8833 9.3708e-11 7.1586e-11 57039.0908 8.987e-11 6.8987e-11 55307.2984 8.6031e-11 6.6388e-11 53575.5059 8.2192e-11 6.3789e-11 51843.7135 7.8353e-11 6.119e-11 50111.9211 7.4515e-11 5.8592e-11 48380.1286 7.0676e-11 5.5993e-11 46648.3362 6.6837e-11 5.3394e-11 64916.5437 6.2998e-11 5.0795e-11 41217.7631 5.916e-11 4.8197e-11 37524.4938 5.5321e-11 4.5598e-11 33842.2474 5.1482e-11 4.2999e-11 30176.5353 4.7643e-11 4.04e-11 26532.8689 4.3805e-11 3.7802e-11 22916.7596 3.9966e-11 3.5203e-11 19333.7188 3.6127e-11 3.2604e-11 15789.2581 3.2288e-11 3.0005e-11 12288.8887 2.8449e-11 2.7407e-11 ;7.9172e-11 62234.4682 2.5578e-10 7.7554e-11 60502.6757 2.4455e-10 7.5937e-11 58770.8833 2.3331e-10 7.432e-11 57039.0908 2.2207e-10 7.2702e-11 55307.2984 2.1083e-10 7.1085e-11 53575.5059 1.996e-10 6.9468e-11 51843.7135 1.8836e-10 6.785e-11 50111.9211 1.7712e-10 6.6233e-11 48380.1286 1.6589e-10 6.4616e-11 46648.3362 1.5465e-10 6.2998e-11 64916.5437 1.4341e-10 6.1381e-11 41217.7631 1.3217e-10 5.9764e-11 37524.4938 1.2094e-10 5.8146e-11 33842.2474 1.097e-10 5.6529e-11 30176.5353 9.8461e-11 5.4912e-11 26532.8689 8.7223e-11 5.3294e-11 22916.7596 7.5986e-11 5.1677e-11 19333.7188 6.4749e-11 5.006e-11 15789.2581 5.3511e-11 4.8442e-11 12288.8887 4.2274e-11 ;8.4862e-11 2.4178e-10 62234.4682 8.1455e-11 2.3194e-10 60502.6757 7.8049e-11 2.2211e-10 58770.8833 7.4642e-11 2.1227e-10 57039.0908 7.1235e-11 2.0243e-10 55307.2984 6.7829e-11 1.9259e-10 53575.5059 6.4422e-11 1.8276e-10 51843.7135 6.1015e-11 1.7292e-10 50111.9211 5.7609e-11 1.6308e-10 48380.1286 5.4202e-11 1.5325e-10 46648.3362 5.0795e-11 1.4341e-10 64916.5437 4.7389e-11 1.3357e-10 41217.7631 4.3982e-11 1.2374e-10 37524.4938 4.0575e-11 1.139e-10 33842.2474 3.7169e-11 1.0406e-10 30176.5353 3.3762e-11 9.4226e-11 26532.8689 3.0356e-11 8.4389e-11 22916.7596 2.6949e-11 7.4552e-11 19333.7188 2.3542e-11 6.4715e-11 15789.2581 2.0136e-11 5.4878e-11 12288.8887 ;56832.9906 9.8518e-11 7.1249e-11 55271.4678 9.4804e-11 6.8863e-11 53709.9451 9.109e-11 6.6477e-11 52148.4223 8.7377e-11 6.4091e-11 50586.8996 8.3663e-11 6.1705e-11 49025.3768 7.9949e-11 5.9319e-11 47463.8541 7.6236e-11 5.6933e-11 45902.3313 7.2522e-11 5.4547e-11 44340.8086 6.8808e-11 5.2161e-11 42779.2858 6.5095e-11 4.9775e-11 41217.7631 6.1381e-11 4.7389e-11 59656.2403 5.7667e-11 4.5003e-11 36127.7293 5.3954e-11 4.2617e-11 32604.7298 5.024e-11 4.0231e-11 29092.7531 4.6526e-11 3.7845e-11 25597.3106 4.2813e-11 3.5459e-11 22123.9139 3.9099e-11 3.3073e-11 18678.0743 3.5385e-11 3.0687e-11 15265.3032 3.1672e-11 2.8301e-11 11891.1121 2.7958e-11 2.5915e-11 ;7.4082e-11 56832.9906 2.3816e-10 7.259e-11 55271.4678 2.277e-10 7.1097e-11 53709.9451 2.1724e-10 6.9605e-11 52148.4223 2.0678e-10 6.8113e-11 50586.8996 1.9633e-10 6.6621e-11 49025.3768 1.8587e-10 6.5128e-11 47463.8541 1.7541e-10 6.3636e-11 45902.3313 1.6495e-10 6.2144e-11 44340.8086 1.5449e-10 6.0652e-11 42779.2858 1.4403e-10 5.916e-11 41217.7631 1.3357e-10 5.7667e-11 59656.2403 1.2311e-10 5.6175e-11 36127.7293 1.1266e-10 5.4683e-11 32604.7298 1.022e-10 5.3191e-11 29092.7531 9.1738e-11 5.1698e-11 25597.3106 8.128e-11 5.0206e-11 22123.9139 7.0821e-11 4.8714e-11 18678.0743 6.0362e-11 4.7222e-11 15265.3032 4.9903e-11 4.5729e-11 11891.1121 3.9444e-11 ;8.0136e-11 2.2276e-10 56832.9906 7.6942e-11 2.137e-10 55271.4678 7.3748e-11 2.0464e-10 53709.9451 7.0554e-11 1.9558e-10 52148.4223 6.736e-11 1.8652e-10 50586.8996 6.4166e-11 1.7746e-10 49025.3768 6.0972e-11 1.6841e-10 47463.8541 5.7778e-11 1.5935e-10 45902.3313 5.4585e-11 1.5029e-10 44340.8086 5.1391e-11 1.4123e-10 42779.2858 4.8197e-11 1.3217e-10 41217.7631 4.5003e-11 1.2311e-10 59656.2403 4.1809e-11 1.1406e-10 36127.7293 3.8615e-11 1.05e-10 32604.7298 3.5421e-11 9.594e-11 29092.7531 3.2227e-11 8.6882e-11 25597.3106 2.9033e-11 7.7823e-11 22123.9139 2.5839e-11 6.8765e-11 18678.0743 2.2645e-11 5.9707e-11 15265.3032 1.9451e-11 5.0648e-11 11891.1121 ;51492.1387 9.5649e-11 6.5715e-11 50095.3742 9.2061e-11 6.3542e-11 48698.6097 8.8472e-11 6.1368e-11 47301.8452 8.4884e-11 5.9195e-11 45905.0807 8.1295e-11 5.7022e-11 44508.3163 7.7706e-11 5.4849e-11 43111.5518 7.4118e-11 5.2675e-11 41714.7873 7.0529e-11 5.0502e-11 40318.0228 6.6941e-11 4.8329e-11 38921.2583 6.3352e-11 4.6155e-11 37524.4938 5.9764e-11 4.3982e-11 36127.7293 5.6175e-11 4.1809e-11 54730.9648 5.2586e-11 3.9636e-11 31367.2121 4.8998e-11 3.7462e-11 28008.9708 4.5409e-11 3.5289e-11 24661.7524 4.1821e-11 3.3116e-11 21331.0682 3.8232e-11 3.0942e-11 18022.4297 3.4644e-11 2.8769e-11 14741.3484 3.1055e-11 2.6596e-11 11493.3356 2.7467e-11 2.4423e-11 ;6.8992e-11 51492.1387 2.2054e-10 6.7625e-11 50095.3742 2.1086e-10 6.6258e-11 48698.6097 2.0118e-10 6.4891e-11 47301.8452 1.915e-10 6.3524e-11 45905.0807 1.8182e-10 6.2157e-11 44508.3163 1.7214e-10 6.0789e-11 43111.5518 1.6246e-10 5.9422e-11 41714.7873 1.5278e-10 5.8055e-11 40318.0228 1.431e-10 5.6688e-11 38921.2583 1.3342e-10 5.5321e-11 37524.4938 1.2374e-10 5.3954e-11 36127.7293 1.1406e-10 5.2586e-11 54730.9648 1.0438e-10 5.1219e-11 31367.2121 9.4696e-11 4.9852e-11 28008.9708 8.5016e-11 4.8485e-11 24661.7524 7.5336e-11 4.7118e-11 21331.0682 6.5656e-11 4.5751e-11 18022.4297 5.5975e-11 4.4384e-11 14741.3484 4.6295e-11 4.3016e-11 11493.3356 3.6615e-11 ;7.541e-11 2.0373e-10 51492.1387 7.2429e-11 1.9545e-10 50095.3742 6.9447e-11 1.8717e-10 48698.6097 6.6466e-11 1.7889e-10 47301.8452 6.3485e-11 1.7061e-10 45905.0807 6.0504e-11 1.6233e-10 44508.3163 5.7523e-11 1.5405e-10 43111.5518 5.4542e-11 1.4577e-10 41714.7873 5.156e-11 1.3749e-10 40318.0228 4.8579e-11 1.2922e-10 38921.2583 4.5598e-11 1.2094e-10 37524.4938 4.2617e-11 1.1266e-10 36127.7293 3.9636e-11 1.0438e-10 54730.9648 3.6654e-11 9.6097e-11 31367.2121 3.3673e-11 8.7817e-11 28008.9708 3.0692e-11 7.9537e-11 24661.7524 2.7711e-11 7.1257e-11 21331.0682 2.473e-11 6.2978e-11 18022.4297 2.1748e-11 5.4698e-11 14741.3484 1.8767e-11 4.6418e-11 11493.3356 ;46217.4239 9.2781e-11 6.0181e-11 44979.9063 8.9318e-11 5.8221e-11 43742.3886 8.5854e-11 5.626e-11 42504.871 8.2391e-11 5.4299e-11 41267.3533 7.8927e-11 5.2339e-11 40029.8357 7.5464e-11 5.0378e-11 38792.318 7.2e-11 4.8418e-11 37554.8004 6.8537e-11 4.6457e-11 36317.2827 6.5073e-11 4.4497e-11 35079.7651 6.161e-11 4.2536e-11 33842.2474 5.8146e-11 4.0575e-11 32604.7298 5.4683e-11 3.8615e-11 31367.2121 5.1219e-11 3.6654e-11 50129.6945 4.7756e-11 3.4694e-11 26925.1886 4.4292e-11 3.2733e-11 23726.1941 4.0829e-11 3.0773e-11 20538.2225 3.7365e-11 2.8812e-11 17366.7852 3.3902e-11 2.6852e-11 14217.3936 3.0438e-11 2.4891e-11 11095.5591 2.6975e-11 2.293e-11 ;6.3903e-11 46217.4239 2.0292e-10 6.2661e-11 44979.9063 1.9401e-10 6.1418e-11 43742.3886 1.8511e-10 6.0176e-11 42504.871 1.7621e-10 5.8934e-11 41267.3533 1.6731e-10 5.7692e-11 40029.8357 1.5841e-10 5.645e-11 38792.318 1.4951e-10 5.5208e-11 37554.8004 1.406e-10 5.3966e-11 36317.2827 1.317e-10 5.2724e-11 35079.7651 1.228e-10 5.1482e-11 33842.2474 1.139e-10 5.024e-11 32604.7298 1.05e-10 4.8998e-11 31367.2121 9.6097e-11 4.7756e-11 50129.6945 8.7195e-11 4.6514e-11 26925.1886 7.8293e-11 4.5272e-11 23726.1941 6.9392e-11 4.403e-11 20538.2225 6.049e-11 4.2788e-11 17366.7852 5.1589e-11 4.1546e-11 14217.3936 4.2687e-11 4.0304e-11 11095.5591 3.3786e-11 ;7.0684e-11 1.8471e-10 46217.4239 6.7915e-11 1.7721e-10 44979.9063 6.5147e-11 1.6971e-10 43742.3886 6.2378e-11 1.6221e-10 42504.871 5.961e-11 1.547e-10 41267.3533 5.6841e-11 1.472e-10 40029.8357 5.4073e-11 1.397e-10 38792.318 5.1305e-11 1.322e-10 37554.8004 4.8536e-11 1.247e-10 36317.2827 4.5768e-11 1.172e-10 35079.7651 4.2999e-11 1.097e-10 33842.2474 4.0231e-11 1.022e-10 32604.7298 3.7462e-11 9.4696e-11 31367.2121 3.4694e-11 8.7195e-11 50129.6945 3.1925e-11 7.9694e-11 26925.1886 2.9157e-11 7.2193e-11 23726.1941 2.6388e-11 6.4692e-11 20538.2225 2.362e-11 5.7191e-11 17366.7852 2.0852e-11 4.969e-11 14217.3936 1.8083e-11 4.2189e-11 11095.5591 ;41014.3576 8.9913e-11 5.4647e-11 39930.5754 8.6574e-11 5.2899e-11 38846.7932 8.3236e-11 5.1151e-11 37763.0109 7.9898e-11 4.9404e-11 36679.2287 7.6559e-11 4.7656e-11 35595.4465 7.3221e-11 4.5908e-11 34511.6642 6.9882e-11 4.416e-11 33427.882 6.6544e-11 4.2412e-11 32344.0998 6.3206e-11 4.0664e-11 31260.3175 5.9867e-11 3.8917e-11 30176.5353 5.6529e-11 3.7169e-11 29092.7531 5.3191e-11 3.5421e-11 28008.9708 4.9852e-11 3.3673e-11 26925.1886 4.6514e-11 3.1925e-11 45841.4064 4.3175e-11 3.0178e-11 22790.6359 3.9837e-11 2.843e-11 19745.3769 3.6499e-11 2.6682e-11 16711.1407 3.316e-11 2.4934e-11 13693.4388 2.9822e-11 2.3186e-11 10697.7825 2.6484e-11 2.1438e-11 ;5.8813e-11 41014.3576 1.8529e-10 5.7696e-11 39930.5754 1.7717e-10 5.6579e-11 38846.7932 1.6905e-10 5.5462e-11 37763.0109 1.6092e-10 5.4345e-11 36679.2287 1.528e-10 5.3228e-11 35595.4465 1.4468e-10 5.2111e-11 34511.6642 1.3655e-10 5.0994e-11 33427.882 1.2843e-10 4.9877e-11 32344.0998 1.2031e-10 4.876e-11 31260.3175 1.1219e-10 4.7643e-11 30176.5353 1.0406e-10 4.6526e-11 29092.7531 9.594e-11 4.5409e-11 28008.9708 8.7817e-11 4.4292e-11 26925.1886 7.9694e-11 4.3175e-11 45841.4064 7.1571e-11 4.2058e-11 22790.6359 6.3448e-11 4.0941e-11 19745.3769 5.5325e-11 3.9825e-11 16711.1407 4.7202e-11 3.8708e-11 13693.4388 3.9079e-11 3.7591e-11 10697.7825 3.0956e-11 ;6.5958e-11 1.6569e-10 41014.3576 6.3402e-11 1.5896e-10 39930.5754 6.0846e-11 1.5224e-10 38846.7932 5.829e-11 1.4552e-10 37763.0109 5.5735e-11 1.388e-10 36679.2287 5.3179e-11 1.3207e-10 35595.4465 5.0623e-11 1.2535e-10 34511.6642 4.8068e-11 1.1863e-10 33427.882 4.5512e-11 1.1191e-10 32344.0998 4.2956e-11 1.0518e-10 31260.3175 4.04e-11 9.8461e-11 30176.5353 3.7845e-11 9.1738e-11 29092.7531 3.5289e-11 8.5016e-11 28008.9708 3.2733e-11 7.8293e-11 26925.1886 3.0178e-11 7.1571e-11 45841.4064 2.7622e-11 6.4849e-11 22790.6359 2.5066e-11 5.8126e-11 19745.3769 2.251e-11 5.1404e-11 16711.1407 1.9955e-11 4.4681e-11 13693.4388 1.7399e-11 3.7959e-11 10697.7825 ;35888.4513 8.7044e-11 4.9113e-11 34952.8931 8.3831e-11 4.7578e-11 34017.3348 8.0618e-11 4.6043e-11 33081.7766 7.7405e-11 4.4508e-11 32146.2183 7.4191e-11 4.2973e-11 31210.6601 7.0978e-11 4.1438e-11 30275.1018 6.7765e-11 3.9903e-11 29339.5436 6.4551e-11 3.8367e-11 28403.9854 6.1338e-11 3.6832e-11 27468.4271 5.8125e-11 3.5297e-11 26532.8689 5.4912e-11 3.3762e-11 25597.3106 5.1698e-11 3.2227e-11 24661.7524 4.8485e-11 3.0692e-11 23726.1941 4.5272e-11 2.9157e-11 22790.6359 4.2058e-11 2.7622e-11 41855.0777 3.8845e-11 2.6087e-11 18952.5312 3.5632e-11 2.4552e-11 16055.4961 3.2419e-11 2.3017e-11 13169.4839 2.9205e-11 2.1481e-11 10300.006 2.5992e-11 1.9946e-11 ;5.3723e-11 35888.4513 1.6767e-10 5.2731e-11 34952.8931 1.6032e-10 5.1739e-11 34017.3348 1.5298e-10 5.0748e-11 33081.7766 1.4564e-10 4.9756e-11 32146.2183 1.3829e-10 4.8764e-11 31210.6601 1.3095e-10 4.7772e-11 30275.1018 1.236e-10 4.678e-11 29339.5436 1.1626e-10 4.5788e-11 28403.9854 1.0891e-10 4.4796e-11 27468.4271 1.0157e-10 4.3805e-11 26532.8689 9.4226e-11 4.2813e-11 25597.3106 8.6882e-11 4.1821e-11 24661.7524 7.9537e-11 4.0829e-11 23726.1941 7.2193e-11 3.9837e-11 22790.6359 6.4849e-11 3.8845e-11 41855.0777 5.7504e-11 3.7853e-11 18952.5312 5.016e-11 3.6861e-11 16055.4961 4.2816e-11 3.587e-11 13169.4839 3.5471e-11 3.4878e-11 10300.006 2.8127e-11 ;6.1232e-11 1.4666e-10 35888.4513 5.8889e-11 1.4072e-10 34952.8931 5.6546e-11 1.3477e-10 34017.3348 5.4203e-11 1.2883e-10 33081.7766 5.186e-11 1.2289e-10 32146.2183 4.9517e-11 1.1694e-10 31210.6601 4.7174e-11 1.11e-10 30275.1018 4.4831e-11 1.0505e-10 29339.5436 4.2488e-11 9.9111e-11 28403.9854 4.0145e-11 9.3167e-11 27468.4271 3.7802e-11 8.7223e-11 26532.8689 3.5459e-11 8.128e-11 25597.3106 3.3116e-11 7.5336e-11 24661.7524 3.0773e-11 6.9392e-11 23726.1941 2.843e-11 6.3448e-11 22790.6359 2.6087e-11 5.7504e-11 41855.0777 2.3744e-11 5.156e-11 18952.5312 2.1401e-11 4.5617e-11 16055.4961 1.9058e-11 3.9673e-11 13169.4839 1.6715e-11 3.3729e-11 10300.006 ;30845.2164 8.4176e-11 4.3579e-11 30052.3707 8.1088e-11 4.2257e-11 29259.525 7.8e-11 4.0934e-11 28466.6793 7.4912e-11 3.9612e-11 27673.8337 7.1823e-11 3.829e-11 26880.988 6.8735e-11 3.6967e-11 26088.1423 6.5647e-11 3.5645e-11 25295.2966 6.2559e-11 3.4323e-11 24502.4509 5.9471e-11 3.3e-11 23709.6053 5.6382e-11 3.1678e-11 22916.7596 5.3294e-11 3.0356e-11 22123.9139 5.0206e-11 2.9033e-11 21331.0682 4.7118e-11 2.7711e-11 20538.2225 4.403e-11 2.6388e-11 19745.3769 4.0941e-11 2.5066e-11 18952.5312 3.7853e-11 2.3744e-11 38159.6855 3.4765e-11 2.2421e-11 15399.8516 3.1677e-11 2.1099e-11 12645.5291 2.8589e-11 1.9777e-11 9902.2295 2.5501e-11 1.8454e-11 ;4.8633e-11 30845.2164 1.5005e-10 4.7767e-11 30052.3707 1.4348e-10 4.69e-11 29259.525 1.3691e-10 4.6033e-11 28466.6793 1.3035e-10 4.5166e-11 27673.8337 1.2378e-10 4.43e-11 26880.988 1.1722e-10 4.3433e-11 26088.1423 1.1065e-10 4.2566e-11 25295.2966 1.0409e-10 4.1699e-11 24502.4509 9.752e-11 4.0833e-11 23709.6053 9.0955e-11 3.9966e-11 22916.7596 8.4389e-11 3.9099e-11 22123.9139 7.7823e-11 3.8232e-11 21331.0682 7.1257e-11 3.7365e-11 20538.2225 6.4692e-11 3.6499e-11 19745.3769 5.8126e-11 3.5632e-11 18952.5312 5.156e-11 3.4765e-11 38159.6855 4.4995e-11 3.3898e-11 15399.8516 3.8429e-11 3.3032e-11 12645.5291 3.1863e-11 3.2165e-11 9902.2295 2.5298e-11 ;5.6505e-11 1.2764e-10 30845.2164 5.4375e-11 1.2247e-10 30052.3707 5.2245e-11 1.1731e-10 29259.525 5.0115e-11 1.1214e-10 28466.6793 4.7984e-11 1.0698e-10 27673.8337 4.5854e-11 1.0181e-10 26880.988 4.3724e-11 9.6647e-11 26088.1423 4.1594e-11 9.1482e-11 25295.2966 3.9463e-11 8.6316e-11 24502.4509 3.7333e-11 8.1151e-11 23709.6053 3.5203e-11 7.5986e-11 22916.7596 3.3073e-11 7.0821e-11 22123.9139 3.0942e-11 6.5656e-11 21331.0682 2.8812e-11 6.049e-11 20538.2225 2.6682e-11 5.5325e-11 19745.3769 2.4552e-11 5.016e-11 18952.5312 2.2421e-11 4.4995e-11 38159.6855 2.0291e-11 3.9829e-11 15399.8516 1.8161e-11 3.4664e-11 12645.5291 1.6031e-11 2.9499e-11 9902.2295 ;25890.1642 8.1308e-11 3.8045e-11 25234.5197 7.8345e-11 3.6936e-11 24578.8751 7.5382e-11 3.5826e-11 23923.2306 7.2419e-11 3.4716e-11 23267.5861 6.9456e-11 3.3607e-11 22611.9415 6.6492e-11 3.2497e-11 21956.297 6.3529e-11 3.1387e-11 21300.6524 6.0566e-11 3.0278e-11 20645.0079 5.7603e-11 2.9168e-11 19989.3634 5.464e-11 2.8058e-11 19333.7188 5.1677e-11 2.6949e-11 18678.0743 4.8714e-11 2.5839e-11 18022.4297 4.5751e-11 2.473e-11 17366.7852 4.2788e-11 2.362e-11 16711.1407 3.9825e-11 2.251e-11 16055.4961 3.6861e-11 2.1401e-11 15399.8516 3.3898e-11 2.0291e-11 34744.207 3.0935e-11 1.9181e-11 12121.5743 2.7972e-11 1.8072e-11 9504.4529 2.5009e-11 1.6962e-11 ;4.3544e-11 25890.1642 1.3242e-10 4.2802e-11 25234.5197 1.2664e-10 4.206e-11 24578.8751 1.2085e-10 4.1319e-11 23923.2306 1.1506e-10 4.0577e-11 23267.5861 1.0927e-10 3.9835e-11 22611.9415 1.0349e-10 3.9094e-11 21956.297 9.77e-11 3.8352e-11 21300.6524 9.1913e-11 3.761e-11 20645.0079 8.6126e-11 3.6869e-11 19989.3634 8.0339e-11 3.6127e-11 19333.7188 7.4552e-11 3.5385e-11 18678.0743 6.8765e-11 3.4644e-11 18022.4297 6.2978e-11 3.3902e-11 17366.7852 5.7191e-11 3.316e-11 16711.1407 5.1404e-11 3.2419e-11 16055.4961 4.5617e-11 3.1677e-11 15399.8516 3.9829e-11 3.0935e-11 34744.207 3.4042e-11 3.0194e-11 12121.5743 2.8255e-11 2.9452e-11 9504.4529 2.2468e-11 ;5.1779e-11 1.0861e-10 25890.1642 4.9862e-11 1.0423e-10 25234.5197 4.7944e-11 9.9841e-11 24578.8751 4.6027e-11 9.5455e-11 23923.2306 4.4109e-11 9.1068e-11 23267.5861 4.2192e-11 8.6682e-11 22611.9415 4.0274e-11 8.2295e-11 21956.297 3.8357e-11 7.7908e-11 21300.6524 3.6439e-11 7.3522e-11 20645.0079 3.4522e-11 6.9135e-11 19989.3634 3.2604e-11 6.4749e-11 19333.7188 3.0687e-11 6.0362e-11 18678.0743 2.8769e-11 5.5975e-11 18022.4297 2.6852e-11 5.1589e-11 17366.7852 2.4934e-11 4.7202e-11 16711.1407 2.3017e-11 4.2816e-11 16055.4961 2.1099e-11 3.8429e-11 15399.8516 1.9181e-11 3.4042e-11 34744.207 1.7264e-11 2.9656e-11 12121.5743 1.5346e-11 2.5269e-11 9504.4529 ;21028.8063 7.844e-11 3.2511e-11 20504.8515 7.5602e-11 3.1614e-11 19980.8967 7.2764e-11 3.0717e-11 19456.9418 6.9926e-11 2.982e-11 18932.987 6.7088e-11 2.8924e-11 18409.0322 6.425e-11 2.8027e-11 17885.0774 6.1412e-11 2.713e-11 17361.1225 5.8574e-11 2.6233e-11 16837.1677 5.5736e-11 2.5336e-11 16313.2129 5.2898e-11 2.4439e-11 15789.2581 5.006e-11 2.3542e-11 15265.3032 4.7222e-11 2.2645e-11 14741.3484 4.4384e-11 2.1748e-11 14217.3936 4.1546e-11 2.0852e-11 13693.4388 3.8708e-11 1.9955e-11 13169.4839 3.587e-11 1.9058e-11 12645.5291 3.3032e-11 1.8161e-11 12121.5743 3.0194e-11 1.7264e-11 31597.6194 2.7356e-11 1.6367e-11 9106.6764 2.4518e-11 1.547e-11 ;3.8454e-11 21028.8063 1.148e-10 3.7838e-11 20504.8515 1.0979e-10 3.7221e-11 19980.8967 1.0478e-10 3.6604e-11 19456.9418 9.9774e-11 3.5988e-11 18932.987 9.4766e-11 3.5371e-11 18409.0322 8.9758e-11 3.4755e-11 17885.0774 8.4749e-11 3.4138e-11 17361.1225 7.9741e-11 3.3521e-11 16837.1677 7.4732e-11 3.2905e-11 16313.2129 6.9724e-11 3.2288e-11 15789.2581 6.4715e-11 3.1672e-11 15265.3032 5.9707e-11 3.1055e-11 14741.3484 5.4698e-11 3.0438e-11 14217.3936 4.969e-11 2.9822e-11 13693.4388 4.4681e-11 2.9205e-11 13169.4839 3.9673e-11 2.8589e-11 12645.5291 3.4664e-11 2.7972e-11 12121.5743 2.9656e-11 2.7356e-11 31597.6194 2.4647e-11 2.6739e-11 9106.6764 1.9639e-11 ;4.7053e-11 8.9591e-11 21028.8063 4.5348e-11 8.5983e-11 20504.8515 4.3644e-11 8.2375e-11 19980.8967 4.1939e-11 7.8767e-11 19456.9418 4.0234e-11 7.5159e-11 18932.987 3.8529e-11 7.1551e-11 18409.0322 3.6825e-11 6.7943e-11 17885.0774 3.512e-11 6.4335e-11 17361.1225 3.3415e-11 6.0727e-11 16837.1677 3.171e-11 5.7119e-11 16313.2129 3.0005e-11 5.3511e-11 15789.2581 2.8301e-11 4.9903e-11 15265.3032 2.6596e-11 4.6295e-11 14741.3484 2.4891e-11 4.2687e-11 14217.3936 2.3186e-11 3.9079e-11 13693.4388 2.1481e-11 3.5471e-11 13169.4839 1.9777e-11 3.1863e-11 12645.5291 1.8072e-11 2.8255e-11 12121.5743 1.6367e-11 2.4647e-11 31597.6194 1.4662e-11 2.1039e-11 9106.6764 ;16266.654 7.5571e-11 2.6977e-11 15868.8775 7.2858e-11 2.6293e-11 15471.101 7.0146e-11 2.5609e-11 15073.3244 6.7433e-11 2.4925e-11 14675.5479 6.472e-11 2.4241e-11 14277.7714 6.2007e-11 2.3556e-11 13879.9948 5.9294e-11 2.2872e-11 13482.2183 5.6581e-11 2.2188e-11 13084.4418 5.3868e-11 2.1504e-11 12686.6652 5.1155e-11 2.082e-11 12288.8887 4.8442e-11 2.0136e-11 11891.1121 4.5729e-11 1.9451e-11 11493.3356 4.3016e-11 1.8767e-11 11095.5591 4.0304e-11 1.8083e-11 10697.7825 3.7591e-11 1.7399e-11 10300.006 3.4878e-11 1.6715e-11 9902.2295 3.2165e-11 1.6031e-11 9504.4529 2.9452e-11 1.5346e-11 9106.6764 2.6739e-11 1.4662e-11 28708.8999 2.4026e-11 1.3978e-11 ;3.3364e-11 16266.654 9.7177e-11 3.2873e-11 15868.8775 9.2947e-11 3.2381e-11 15471.101 8.8717e-11 3.189e-11 15073.3244 8.4487e-11 3.1398e-11 14675.5479 8.0257e-11 3.0907e-11 14277.7714 7.6028e-11 3.0415e-11 13879.9948 7.1798e-11 2.9924e-11 13482.2183 6.7568e-11 2.9432e-11 13084.4418 6.3338e-11 2.8941e-11 12686.6652 5.9108e-11 2.8449e-11 12288.8887 5.4878e-11 2.7958e-11 11891.1121 5.0648e-11 2.7467e-11 11493.3356 4.6418e-11 2.6975e-11 11095.5591 4.2189e-11 2.6484e-11 10697.7825 3.7959e-11 2.5992e-11 10300.006 3.3729e-11 2.5501e-11 9902.2295 2.9499e-11 2.5009e-11 9504.4529 2.5269e-11 2.4518e-11 9106.6764 2.1039e-11 2.4026e-11 28708.8999 1.6809e-11 ;4.2327e-11 7.0567e-11 16266.654 4.0835e-11 6.7738e-11 15868.8775 3.9343e-11 6.4909e-11 15471.101 3.7851e-11 6.2079e-11 15073.3244 3.6359e-11 5.925e-11 14675.5479 3.4867e-11 5.6421e-11 14277.7714 3.3375e-11 5.3591e-11 13879.9948 3.1883e-11 5.0762e-11 13482.2183 3.0391e-11 4.7933e-11 13084.4418 2.8899e-11 4.5103e-11 12686.6652 2.7407e-11 4.2274e-11 12288.8887 2.5915e-11 3.9444e-11 11891.1121 2.4423e-11 3.6615e-11 11493.3356 2.293e-11 3.3786e-11 11095.5591 2.1438e-11 3.0956e-11 10697.7825 1.9946e-11 2.8127e-11 10300.006 1.8454e-11 2.5298e-11 9902.2295 1.6962e-11 2.2468e-11 9504.4529 1.547e-11 1.9639e-11 9106.6764 1.3978e-11 1.6809e-11 28708.8999 ;]; 
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

      igr2[i] = nlamt[i] * ((d + 0.20531) - mesil[i]);
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

