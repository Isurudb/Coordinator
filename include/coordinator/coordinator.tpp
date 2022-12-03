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
    double r=0, p=0, y=0;  // Rotate the previous pose by 180* about Z

        q_ref.setRPY(r, p, y);
        tf2::convert(attitude,attitude_);
        q_ref_inv=q_ref.inverse();//
  q_e= q_ref_inv*attitude_;  // Calculate the new orientation
  q_e.normalize();


    position_.x = px;
    position_.y = py;
    position_.z = pz;

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

   double a[18] = { 0.0992, 0.0, 0.0, 0.0992, 0.0, 0.0, 0.0,
    0.0992, 0.0, 0.0, 0.0992, 0.0, 0.0, 0.0, 0.0992, 0.0, 0.0,
    0.0992 };

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

 double a_1[18] = { -0.003068559900000001, -0.0, -0.0, -0.0,
    -0.07183028100000001, -0.0, -0.0, -0.0, -0.072222055180000008, -0.153427995,
    -0.0, -0.0, -0.0, -3.59151405, -0.0, -0.0, -0.0, -3.611102759 };

  double dv[9];
  double b_qx[6];
  double U_body[3];
  double b_fx[3];
  double d;
  double d1;
  double d2;
  double d3;
  double d4;

  //  Distance between two Astrobee com
  //  =2*0.153427995 Ixx of compund object
  // =2*0.14271405+m*d^2 Iyy of compund object
  // =2*0.162302759+m*d^2 Izz of compund object
  // 1
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
  b_qx[0] = qx;
  b_qx[1] = qy;
  b_qx[2] = qz;
  b_qx[3] = arg_omegax;
  b_qx[4] = arg_omegay;
  b_qx[5] = arg_omegaz;
  for (int32_T i = 0; i < 3; i++) {
    U_body[i] = ((dv[i] * fx) + (dv[i + 3] * fy)) + (dv[i + 6] * fz);
    d = 0.0;
    for (int32_T i1 = 0; i1 < 6; i1++) {
      d += a_1[i + (3 * i1)] * b_qx[i1];
    }

    b_fx[i] = d;
  }

  /* *tau_x = b_fx[0];
  *tau_y = b_fx[1] - (-0.3 * U_body[2]);
  *tau_z = b_fx[2] - (0.3 * U_body[1]); */

 arg_tau_x = b_fx[0];
 arg_tau_y = b_fx[1] - (-0.3 * U_body[2]);
 arg_tau_z = b_fx[2] - (0.3 * U_body[1]);


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
    0.0, 0.0, 0.0, 0.0, 1567.5765, 1.5888E-12, 2.7019E-12, 8536.6551, 1.0896E-11,
    3.5199E-11, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.5888E-12, 1567.5765,
    9.8891E-12, 6.7324E-11, 8536.6551, 6.1982E-11, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 2.7019E-12, 9.8891E-12, 1567.5765, 1.4677E-11, 9.7557E-11,
    8536.6551, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8536.6551, 6.7324E-11,
    1.4677E-11, 423802.2092, 1.1408E-9, 6.5756E-10, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0896E-11, 8536.6551, 9.7557E-11, 1.1408E-9, 423802.2092,
    6.8455E-10, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.5199E-11,
    6.1982E-11, 8536.6551, 6.5756E-10, 6.8455E-10, 423802.2092 };

  static const real_T b_b[7200] = { 0.099206, 0.0, 0.0, 0.099206, 0.0, 0.0,
    0.29762, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.49603, 0.0, 0.0, 0.099206, 0.0, 0.0,
    0.69444, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.89285, 0.0, 0.0, 0.099206, 0.0, 0.0,
    1.0913, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.2897, 0.0, 0.0, 0.099206, 0.0, 0.0,
    1.4881, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.6865, 0.0, 0.0, 0.099206, 0.0, 0.0,
    1.8849, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.0833, 0.0, 0.0, 0.099206, 0.0, 0.0,
    2.2817, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.4802, 0.0, 0.0, 0.099206, 0.0, 0.0,
    2.6786, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.877, 0.0, 0.0, 0.099206, 0.0, 0.0,
    3.0754, 0.0, 0.0, 0.099206, 0.0, 0.0, 3.2738, 0.0, 0.0, 0.099206, 0.0, 0.0,
    3.4722, 0.0, 0.0, 0.099206, 0.0, 0.0, 3.6706, 0.0, 0.0, 0.099206, 0.0, 0.0,
    3.869, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.099206, 0.0,
    0.0, 0.29762, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.49603, 0.0, 0.0, 0.099206, 0.0,
    0.0, 0.69444, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.89285, 0.0, 0.0, 0.099206, 0.0,
    0.0, 1.0913, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.2897, 0.0, 0.0, 0.099206, 0.0,
    0.0, 1.4881, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.6865, 0.0, 0.0, 0.099206, 0.0,
    0.0, 1.8849, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.0833, 0.0, 0.0, 0.099206, 0.0,
    0.0, 2.2817, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.4802, 0.0, 0.0, 0.099206, 0.0,
    0.0, 2.6786, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.877, 0.0, 0.0, 0.099206, 0.0,
    0.0, 3.0754, 0.0, 0.0, 0.099206, 0.0, 0.0, 3.2738, 0.0, 0.0, 0.099206, 0.0,
    0.0, 3.4722, 0.0, 0.0, 0.099206, 0.0, 0.0, 3.6706, 0.0, 0.0, 0.099206, 0.0,
    0.0, 3.869, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.099206,
    0.0, 0.0, 0.29762, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.49603, 0.0, 0.0, 0.099206,
    0.0, 0.0, 0.69444, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.89285, 0.0, 0.0, 0.099206,
    0.0, 0.0, 1.0913, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.2897, 0.0, 0.0, 0.099206,
    0.0, 0.0, 1.4881, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.6865, 0.0, 0.0, 0.099206,
    0.0, 0.0, 1.8849, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.0833, 0.0, 0.0, 0.099206,
    0.0, 0.0, 2.2817, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.4802, 0.0, 0.0, 0.099206,
    0.0, 0.0, 2.6786, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.877, 0.0, 0.0, 0.099206,
    0.0, 0.0, 3.0754, 0.0, 0.0, 0.099206, 0.0, 0.0, 3.2738, 0.0, 0.0, 0.099206,
    0.0, 0.0, 3.4722, 0.0, 0.0, 0.099206, 0.0, 0.0, 3.6706, 0.0, 0.0, 0.099206,
    0.0, 0.0, 3.869, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.099206,
    0.0, 0.0, 0.099206, 0.0, 0.0, 0.29762, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.49603,
    0.0, 0.0, 0.099206, 0.0, 0.0, 0.69444, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.89285,
    0.0, 0.0, 0.099206, 0.0, 0.0, 1.0913, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.2897,
    0.0, 0.0, 0.099206, 0.0, 0.0, 1.4881, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.6865,
    0.0, 0.0, 0.099206, 0.0, 0.0, 1.8849, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.0833,
    0.0, 0.0, 0.099206, 0.0, 0.0, 2.2817, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.4802,
    0.0, 0.0, 0.099206, 0.0, 0.0, 2.6786, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.877,
    0.0, 0.0, 0.099206, 0.0, 0.0, 3.0754, 0.0, 0.0, 0.099206, 0.0, 0.0, 3.2738,
    0.0, 0.0, 0.099206, 0.0, 0.0, 3.4722, 0.0, 0.0, 0.099206, 0.0, 0.0, 3.6706,
    0.0, 0.0, 0.099206, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.099206,
    0.0, 0.0, 0.099206, 0.0, 0.0, 0.29762, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.49603,
    0.0, 0.0, 0.099206, 0.0, 0.0, 0.69444, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.89285,
    0.0, 0.0, 0.099206, 0.0, 0.0, 1.0913, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.2897,
    0.0, 0.0, 0.099206, 0.0, 0.0, 1.4881, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.6865,
    0.0, 0.0, 0.099206, 0.0, 0.0, 1.8849, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.0833,
    0.0, 0.0, 0.099206, 0.0, 0.0, 2.2817, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.4802,
    0.0, 0.0, 0.099206, 0.0, 0.0, 2.6786, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.877,
    0.0, 0.0, 0.099206, 0.0, 0.0, 3.0754, 0.0, 0.0, 0.099206, 0.0, 0.0, 3.2738,
    0.0, 0.0, 0.099206, 0.0, 0.0, 3.4722, 0.0, 0.0, 0.099206, 0.0, 0.0, 3.6706,
    0.0, 0.0, 0.099206, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.099206,
    0.0, 0.0, 0.099206, 0.0, 0.0, 0.29762, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.49603,
    0.0, 0.0, 0.099206, 0.0, 0.0, 0.69444, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.89285,
    0.0, 0.0, 0.099206, 0.0, 0.0, 1.0913, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.2897,
    0.0, 0.0, 0.099206, 0.0, 0.0, 1.4881, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.6865,
    0.0, 0.0, 0.099206, 0.0, 0.0, 1.8849, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.0833,
    0.0, 0.0, 0.099206, 0.0, 0.0, 2.2817, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.4802,
    0.0, 0.0, 0.099206, 0.0, 0.0, 2.6786, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.877,
    0.0, 0.0, 0.099206, 0.0, 0.0, 3.0754, 0.0, 0.0, 0.099206, 0.0, 0.0, 3.2738,
    0.0, 0.0, 0.099206, 0.0, 0.0, 3.4722, 0.0, 0.0, 0.099206, 0.0, 0.0, 3.6706,
    0.0, 0.0, 0.099206, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.099206, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.29762, 0.0, 0.0, 0.099206,
    0.0, 0.0, 0.49603, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.69444, 0.0, 0.0, 0.099206,
    0.0, 0.0, 0.89285, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.0913, 0.0, 0.0, 0.099206,
    0.0, 0.0, 1.2897, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.4881, 0.0, 0.0, 0.099206,
    0.0, 0.0, 1.6865, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.8849, 0.0, 0.0, 0.099206,
    0.0, 0.0, 2.0833, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.2817, 0.0, 0.0, 0.099206,
    0.0, 0.0, 2.4802, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.6786, 0.0, 0.0, 0.099206,
    0.0, 0.0, 2.877, 0.0, 0.0, 0.099206, 0.0, 0.0, 3.0754, 0.0, 0.0, 0.099206,
    0.0, 0.0, 3.2738, 0.0, 0.0, 0.099206, 0.0, 0.0, 3.4722, 0.0, 0.0, 0.099206,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.099206, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.29762, 0.0, 0.0, 0.099206, 0.0,
    0.0, 0.49603, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.69444, 0.0, 0.0, 0.099206, 0.0,
    0.0, 0.89285, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.0913, 0.0, 0.0, 0.099206, 0.0,
    0.0, 1.2897, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.4881, 0.0, 0.0, 0.099206, 0.0,
    0.0, 1.6865, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.8849, 0.0, 0.0, 0.099206, 0.0,
    0.0, 2.0833, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.2817, 0.0, 0.0, 0.099206, 0.0,
    0.0, 2.4802, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.6786, 0.0, 0.0, 0.099206, 0.0,
    0.0, 2.877, 0.0, 0.0, 0.099206, 0.0, 0.0, 3.0754, 0.0, 0.0, 0.099206, 0.0,
    0.0, 3.2738, 0.0, 0.0, 0.099206, 0.0, 0.0, 3.4722, 0.0, 0.0, 0.099206, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.099206, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.29762, 0.0, 0.0, 0.099206, 0.0,
    0.0, 0.49603, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.69444, 0.0, 0.0, 0.099206, 0.0,
    0.0, 0.89285, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.0913, 0.0, 0.0, 0.099206, 0.0,
    0.0, 1.2897, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.4881, 0.0, 0.0, 0.099206, 0.0,
    0.0, 1.6865, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.8849, 0.0, 0.0, 0.099206, 0.0,
    0.0, 2.0833, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.2817, 0.0, 0.0, 0.099206, 0.0,
    0.0, 2.4802, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.6786, 0.0, 0.0, 0.099206, 0.0,
    0.0, 2.877, 0.0, 0.0, 0.099206, 0.0, 0.0, 3.0754, 0.0, 0.0, 0.099206, 0.0,
    0.0, 3.2738, 0.0, 0.0, 0.099206, 0.0, 0.0, 3.4722, 0.0, 0.0, 0.099206, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.099206, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.29762, 0.0, 0.0,
    0.099206, 0.0, 0.0, 0.49603, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.69444, 0.0, 0.0,
    0.099206, 0.0, 0.0, 0.89285, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.0913, 0.0, 0.0,
    0.099206, 0.0, 0.0, 1.2897, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.4881, 0.0, 0.0,
    0.099206, 0.0, 0.0, 1.6865, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.8849, 0.0, 0.0,
    0.099206, 0.0, 0.0, 2.0833, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.2817, 0.0, 0.0,
    0.099206, 0.0, 0.0, 2.4802, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.6786, 0.0, 0.0,
    0.099206, 0.0, 0.0, 2.877, 0.0, 0.0, 0.099206, 0.0, 0.0, 3.0754, 0.0, 0.0,
    0.099206, 0.0, 0.0, 3.2738, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.099206, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.29762, 0.0, 0.0, 0.099206, 0.0,
    0.0, 0.49603, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.69444, 0.0, 0.0, 0.099206, 0.0,
    0.0, 0.89285, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.0913, 0.0, 0.0, 0.099206, 0.0,
    0.0, 1.2897, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.4881, 0.0, 0.0, 0.099206, 0.0,
    0.0, 1.6865, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.8849, 0.0, 0.0, 0.099206, 0.0,
    0.0, 2.0833, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.2817, 0.0, 0.0, 0.099206, 0.0,
    0.0, 2.4802, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.6786, 0.0, 0.0, 0.099206, 0.0,
    0.0, 2.877, 0.0, 0.0, 0.099206, 0.0, 0.0, 3.0754, 0.0, 0.0, 0.099206, 0.0,
    0.0, 3.2738, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.099206, 0.0,
    0.0, 0.099206, 0.0, 0.0, 0.29762, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.49603, 0.0,
    0.0, 0.099206, 0.0, 0.0, 0.69444, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.89285, 0.0,
    0.0, 0.099206, 0.0, 0.0, 1.0913, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.2897, 0.0,
    0.0, 0.099206, 0.0, 0.0, 1.4881, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.6865, 0.0,
    0.0, 0.099206, 0.0, 0.0, 1.8849, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.0833, 0.0,
    0.0, 0.099206, 0.0, 0.0, 2.2817, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.4802, 0.0,
    0.0, 0.099206, 0.0, 0.0, 2.6786, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.877, 0.0,
    0.0, 0.099206, 0.0, 0.0, 3.0754, 0.0, 0.0, 0.099206, 0.0, 0.0, 3.2738, 0.0,
    0.0, 0.099206, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.099206, 0.0,
    0.0, 0.099206, 0.0, 0.0, 0.29762, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.49603, 0.0,
    0.0, 0.099206, 0.0, 0.0, 0.69444, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.89285, 0.0,
    0.0, 0.099206, 0.0, 0.0, 1.0913, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.2897, 0.0,
    0.0, 0.099206, 0.0, 0.0, 1.4881, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.6865, 0.0,
    0.0, 0.099206, 0.0, 0.0, 1.8849, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.0833, 0.0,
    0.0, 0.099206, 0.0, 0.0, 2.2817, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.4802, 0.0,
    0.0, 0.099206, 0.0, 0.0, 2.6786, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.877, 0.0,
    0.0, 0.099206, 0.0, 0.0, 3.0754, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.099206, 0.0, 0.0,
    0.29762, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.49603, 0.0, 0.0, 0.099206, 0.0, 0.0,
    0.69444, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.89285, 0.0, 0.0, 0.099206, 0.0, 0.0,
    1.0913, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.2897, 0.0, 0.0, 0.099206, 0.0, 0.0,
    1.4881, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.6865, 0.0, 0.0, 0.099206, 0.0, 0.0,
    1.8849, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.0833, 0.0, 0.0, 0.099206, 0.0, 0.0,
    2.2817, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.4802, 0.0, 0.0, 0.099206, 0.0, 0.0,
    2.6786, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.877, 0.0, 0.0, 0.099206, 0.0, 0.0,
    3.0754, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.099206, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.29762, 0.0, 0.0,
    0.099206, 0.0, 0.0, 0.49603, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.69444, 0.0, 0.0,
    0.099206, 0.0, 0.0, 0.89285, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.0913, 0.0, 0.0,
    0.099206, 0.0, 0.0, 1.2897, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.4881, 0.0, 0.0,
    0.099206, 0.0, 0.0, 1.6865, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.8849, 0.0, 0.0,
    0.099206, 0.0, 0.0, 2.0833, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.2817, 0.0, 0.0,
    0.099206, 0.0, 0.0, 2.4802, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.6786, 0.0, 0.0,
    0.099206, 0.0, 0.0, 2.877, 0.0, 0.0, 0.099206, 0.0, 0.0, 3.0754, 0.0, 0.0,
    0.099206, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.099206, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.29762, 0.0, 0.0,
    0.099206, 0.0, 0.0, 0.49603, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.69444, 0.0, 0.0,
    0.099206, 0.0, 0.0, 0.89285, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.0913, 0.0, 0.0,
    0.099206, 0.0, 0.0, 1.2897, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.4881, 0.0, 0.0,
    0.099206, 0.0, 0.0, 1.6865, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.8849, 0.0, 0.0,
    0.099206, 0.0, 0.0, 2.0833, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.2817, 0.0, 0.0,
    0.099206, 0.0, 0.0, 2.4802, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.6786, 0.0, 0.0,
    0.099206, 0.0, 0.0, 2.877, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.099206, 0.0,
    0.0, 0.099206, 0.0, 0.0, 0.29762, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.49603, 0.0,
    0.0, 0.099206, 0.0, 0.0, 0.69444, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.89285, 0.0,
    0.0, 0.099206, 0.0, 0.0, 1.0913, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.2897, 0.0,
    0.0, 0.099206, 0.0, 0.0, 1.4881, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.6865, 0.0,
    0.0, 0.099206, 0.0, 0.0, 1.8849, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.0833, 0.0,
    0.0, 0.099206, 0.0, 0.0, 2.2817, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.4802, 0.0,
    0.0, 0.099206, 0.0, 0.0, 2.6786, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.877, 0.0,
    0.0, 0.099206, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.099206, 0.0, 0.0,
    0.29762, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.49603, 0.0, 0.0, 0.099206, 0.0, 0.0,
    0.69444, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.89285, 0.0, 0.0, 0.099206, 0.0, 0.0,
    1.0913, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.2897, 0.0, 0.0, 0.099206, 0.0, 0.0,
    1.4881, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.6865, 0.0, 0.0, 0.099206, 0.0, 0.0,
    1.8849, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.0833, 0.0, 0.0, 0.099206, 0.0, 0.0,
    2.2817, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.4802, 0.0, 0.0, 0.099206, 0.0, 0.0,
    2.6786, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.877, 0.0, 0.0, 0.099206, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.29762, 0.0,
    0.0, 0.099206, 0.0, 0.0, 0.49603, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.69444, 0.0,
    0.0, 0.099206, 0.0, 0.0, 0.89285, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.0913, 0.0,
    0.0, 0.099206, 0.0, 0.0, 1.2897, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.4881, 0.0,
    0.0, 0.099206, 0.0, 0.0, 1.6865, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.8849, 0.0,
    0.0, 0.099206, 0.0, 0.0, 2.0833, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.2817, 0.0,
    0.0, 0.099206, 0.0, 0.0, 2.4802, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.6786, 0.0,
    0.0, 0.099206, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.099206, 0.0,
    0.0, 0.099206, 0.0, 0.0, 0.29762, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.49603, 0.0,
    0.0, 0.099206, 0.0, 0.0, 0.69444, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.89285, 0.0,
    0.0, 0.099206, 0.0, 0.0, 1.0913, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.2897, 0.0,
    0.0, 0.099206, 0.0, 0.0, 1.4881, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.6865, 0.0,
    0.0, 0.099206, 0.0, 0.0, 1.8849, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.0833, 0.0,
    0.0, 0.099206, 0.0, 0.0, 2.2817, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.4802, 0.0,
    0.0, 0.099206, 0.0, 0.0, 2.6786, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.29762, 0.0,
    0.0, 0.099206, 0.0, 0.0, 0.49603, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.69444, 0.0,
    0.0, 0.099206, 0.0, 0.0, 0.89285, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.0913, 0.0,
    0.0, 0.099206, 0.0, 0.0, 1.2897, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.4881, 0.0,
    0.0, 0.099206, 0.0, 0.0, 1.6865, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.8849, 0.0,
    0.0, 0.099206, 0.0, 0.0, 2.0833, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.2817, 0.0,
    0.0, 0.099206, 0.0, 0.0, 2.4802, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.6786, 0.0,
    0.0, 0.099206, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.099206, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.29762, 0.0, 0.0, 0.099206, 0.0,
    0.0, 0.49603, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.69444, 0.0, 0.0, 0.099206, 0.0,
    0.0, 0.89285, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.0913, 0.0, 0.0, 0.099206, 0.0,
    0.0, 1.2897, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.4881, 0.0, 0.0, 0.099206, 0.0,
    0.0, 1.6865, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.8849, 0.0, 0.0, 0.099206, 0.0,
    0.0, 2.0833, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.2817, 0.0, 0.0, 0.099206, 0.0,
    0.0, 2.4802, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.099206, 0.0, 0.0,
    0.29762, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.49603, 0.0, 0.0, 0.099206, 0.0, 0.0,
    0.69444, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.89285, 0.0, 0.0, 0.099206, 0.0, 0.0,
    1.0913, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.2897, 0.0, 0.0, 0.099206, 0.0, 0.0,
    1.4881, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.6865, 0.0, 0.0, 0.099206, 0.0, 0.0,
    1.8849, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.0833, 0.0, 0.0, 0.099206, 0.0, 0.0,
    2.2817, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.4802, 0.0, 0.0, 0.099206, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.099206,
    0.0, 0.0, 0.099206, 0.0, 0.0, 0.29762, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.49603,
    0.0, 0.0, 0.099206, 0.0, 0.0, 0.69444, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.89285,
    0.0, 0.0, 0.099206, 0.0, 0.0, 1.0913, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.2897,
    0.0, 0.0, 0.099206, 0.0, 0.0, 1.4881, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.6865,
    0.0, 0.0, 0.099206, 0.0, 0.0, 1.8849, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.0833,
    0.0, 0.0, 0.099206, 0.0, 0.0, 2.2817, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.4802,
    0.0, 0.0, 0.099206, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.099206, 0.0, 0.0,
    0.29762, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.49603, 0.0, 0.0, 0.099206, 0.0, 0.0,
    0.69444, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.89285, 0.0, 0.0, 0.099206, 0.0, 0.0,
    1.0913, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.2897, 0.0, 0.0, 0.099206, 0.0, 0.0,
    1.4881, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.6865, 0.0, 0.0, 0.099206, 0.0, 0.0,
    1.8849, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.0833, 0.0, 0.0, 0.099206, 0.0, 0.0,
    2.2817, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.099206, 0.0, 0.0,
    0.099206, 0.0, 0.0, 0.29762, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.49603, 0.0, 0.0,
    0.099206, 0.0, 0.0, 0.69444, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.89285, 0.0, 0.0,
    0.099206, 0.0, 0.0, 1.0913, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.2897, 0.0, 0.0,
    0.099206, 0.0, 0.0, 1.4881, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.6865, 0.0, 0.0,
    0.099206, 0.0, 0.0, 1.8849, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.0833, 0.0, 0.0,
    0.099206, 0.0, 0.0, 2.2817, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.099206, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.29762, 0.0, 0.0, 0.099206, 0.0,
    0.0, 0.49603, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.69444, 0.0, 0.0, 0.099206, 0.0,
    0.0, 0.89285, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.0913, 0.0, 0.0, 0.099206, 0.0,
    0.0, 1.2897, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.4881, 0.0, 0.0, 0.099206, 0.0,
    0.0, 1.6865, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.8849, 0.0, 0.0, 0.099206, 0.0,
    0.0, 2.0833, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.2817, 0.0, 0.0, 0.099206, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.099206, 0.0,
    0.0, 0.29762, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.49603, 0.0, 0.0, 0.099206, 0.0,
    0.0, 0.69444, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.89285, 0.0, 0.0, 0.099206, 0.0,
    0.0, 1.0913, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.2897, 0.0, 0.0, 0.099206, 0.0,
    0.0, 1.4881, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.6865, 0.0, 0.0, 0.099206, 0.0,
    0.0, 1.8849, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.0833, 0.0, 0.0, 0.099206, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.099206, 0.0, 0.0,
    0.099206, 0.0, 0.0, 0.29762, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.49603, 0.0, 0.0,
    0.099206, 0.0, 0.0, 0.69444, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.89285, 0.0, 0.0,
    0.099206, 0.0, 0.0, 1.0913, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.2897, 0.0, 0.0,
    0.099206, 0.0, 0.0, 1.4881, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.6865, 0.0, 0.0,
    0.099206, 0.0, 0.0, 1.8849, 0.0, 0.0, 0.099206, 0.0, 0.0, 2.0833, 0.0, 0.0,
    0.099206, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.099206, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.29762, 0.0, 0.0, 0.099206, 0.0,
    0.0, 0.49603, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.69444, 0.0, 0.0, 0.099206, 0.0,
    0.0, 0.89285, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.0913, 0.0, 0.0, 0.099206, 0.0,
    0.0, 1.2897, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.4881, 0.0, 0.0, 0.099206, 0.0,
    0.0, 1.6865, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.8849, 0.0, 0.0, 0.099206, 0.0,
    0.0, 2.0833, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.099206, 0.0, 0.0,
    0.29762, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.49603, 0.0, 0.0, 0.099206, 0.0, 0.0,
    0.69444, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.89285, 0.0, 0.0, 0.099206, 0.0, 0.0,
    1.0913, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.2897, 0.0, 0.0, 0.099206, 0.0, 0.0,
    1.4881, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.6865, 0.0, 0.0, 0.099206, 0.0, 0.0,
    1.8849, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.099206, 0.0,
    0.0, 0.29762, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.49603, 0.0, 0.0, 0.099206, 0.0,
    0.0, 0.69444, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.89285, 0.0, 0.0, 0.099206, 0.0,
    0.0, 1.0913, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.2897, 0.0, 0.0, 0.099206, 0.0,
    0.0, 1.4881, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.6865, 0.0, 0.0, 0.099206, 0.0,
    0.0, 1.8849, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.099206,
    0.0, 0.0, 0.29762, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.49603, 0.0, 0.0, 0.099206,
    0.0, 0.0, 0.69444, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.89285, 0.0, 0.0, 0.099206,
    0.0, 0.0, 1.0913, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.2897, 0.0, 0.0, 0.099206,
    0.0, 0.0, 1.4881, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.6865, 0.0, 0.0, 0.099206,
    0.0, 0.0, 1.8849, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.099206,
    0.0, 0.0, 0.099206, 0.0, 0.0, 0.29762, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.49603,
    0.0, 0.0, 0.099206, 0.0, 0.0, 0.69444, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.89285,
    0.0, 0.0, 0.099206, 0.0, 0.0, 1.0913, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.2897,
    0.0, 0.0, 0.099206, 0.0, 0.0, 1.4881, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.6865,
    0.0, 0.0, 0.099206, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.099206,
    0.0, 0.0, 0.099206, 0.0, 0.0, 0.29762, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.49603,
    0.0, 0.0, 0.099206, 0.0, 0.0, 0.69444, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.89285,
    0.0, 0.0, 0.099206, 0.0, 0.0, 1.0913, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.2897,
    0.0, 0.0, 0.099206, 0.0, 0.0, 1.4881, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.6865,
    0.0, 0.0, 0.099206, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.099206,
    0.0, 0.0, 0.099206, 0.0, 0.0, 0.29762, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.49603,
    0.0, 0.0, 0.099206, 0.0, 0.0, 0.69444, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.89285,
    0.0, 0.0, 0.099206, 0.0, 0.0, 1.0913, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.2897,
    0.0, 0.0, 0.099206, 0.0, 0.0, 1.4881, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.6865,
    0.0, 0.0, 0.099206, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.099206, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.29762, 0.0, 0.0, 0.099206,
    0.0, 0.0, 0.49603, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.69444, 0.0, 0.0, 0.099206,
    0.0, 0.0, 0.89285, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.0913, 0.0, 0.0, 0.099206,
    0.0, 0.0, 1.2897, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.4881, 0.0, 0.0, 0.099206,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.099206, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.29762, 0.0, 0.0, 0.099206, 0.0,
    0.0, 0.49603, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.69444, 0.0, 0.0, 0.099206, 0.0,
    0.0, 0.89285, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.0913, 0.0, 0.0, 0.099206, 0.0,
    0.0, 1.2897, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.4881, 0.0, 0.0, 0.099206, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.099206, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.29762, 0.0, 0.0, 0.099206, 0.0,
    0.0, 0.49603, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.69444, 0.0, 0.0, 0.099206, 0.0,
    0.0, 0.89285, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.0913, 0.0, 0.0, 0.099206, 0.0,
    0.0, 1.2897, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.4881, 0.0, 0.0, 0.099206, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.099206, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.29762, 0.0, 0.0,
    0.099206, 0.0, 0.0, 0.49603, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.69444, 0.0, 0.0,
    0.099206, 0.0, 0.0, 0.89285, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.0913, 0.0, 0.0,
    0.099206, 0.0, 0.0, 1.2897, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.099206, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.29762, 0.0, 0.0, 0.099206, 0.0,
    0.0, 0.49603, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.69444, 0.0, 0.0, 0.099206, 0.0,
    0.0, 0.89285, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.0913, 0.0, 0.0, 0.099206, 0.0,
    0.0, 1.2897, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.099206, 0.0,
    0.0, 0.099206, 0.0, 0.0, 0.29762, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.49603, 0.0,
    0.0, 0.099206, 0.0, 0.0, 0.69444, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.89285, 0.0,
    0.0, 0.099206, 0.0, 0.0, 1.0913, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.2897, 0.0,
    0.0, 0.099206, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.099206, 0.0,
    0.0, 0.099206, 0.0, 0.0, 0.29762, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.49603, 0.0,
    0.0, 0.099206, 0.0, 0.0, 0.69444, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.89285, 0.0,
    0.0, 0.099206, 0.0, 0.0, 1.0913, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.099206, 0.0, 0.0,
    0.29762, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.49603, 0.0, 0.0, 0.099206, 0.0, 0.0,
    0.69444, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.89285, 0.0, 0.0, 0.099206, 0.0, 0.0,
    1.0913, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.099206, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.29762, 0.0, 0.0,
    0.099206, 0.0, 0.0, 0.49603, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.69444, 0.0, 0.0,
    0.099206, 0.0, 0.0, 0.89285, 0.0, 0.0, 0.099206, 0.0, 0.0, 1.0913, 0.0, 0.0,
    0.099206, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.099206, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.29762, 0.0, 0.0,
    0.099206, 0.0, 0.0, 0.49603, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.69444, 0.0, 0.0,
    0.099206, 0.0, 0.0, 0.89285, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.099206,
    0.0, 0.0, 0.099206, 0.0, 0.0, 0.29762, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.49603,
    0.0, 0.0, 0.099206, 0.0, 0.0, 0.69444, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.89285,
    0.0, 0.0, 0.099206, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.099206, 0.0, 0.0,
    0.29762, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.49603, 0.0, 0.0, 0.099206, 0.0, 0.0,
    0.69444, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.89285, 0.0, 0.0, 0.099206, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.29762, 0.0,
    0.0, 0.099206, 0.0, 0.0, 0.49603, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.69444, 0.0,
    0.0, 0.099206, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.099206, 0.0,
    0.0, 0.099206, 0.0, 0.0, 0.29762, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.49603, 0.0,
    0.0, 0.099206, 0.0, 0.0, 0.69444, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.29762,
    0.0, 0.0, 0.099206, 0.0, 0.0, 0.49603, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.69444,
    0.0, 0.0, 0.099206, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.099206, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.29762, 0.0, 0.0, 0.099206,
    0.0, 0.0, 0.49603, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.099206, 0.0,
    0.0, 0.29762, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.49603, 0.0, 0.0, 0.099206, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.099206, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.29762, 0.0, 0.0, 0.099206, 0.0,
    0.0, 0.49603, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.099206, 0.0, 0.0,
    0.099206, 0.0, 0.0, 0.29762, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.099206, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.29762, 0.0, 0.0, 0.099206,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.099206, 0.0, 0.0,
    0.29762, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.099206, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.099206, 0.0, 0.0, 0.099206 };

  static const real_T b_H[3600] = { 138371.1425, 1.3007E-10, 1.3212E-10,
    112666.5366, 1.2498E-10, 1.274E-10, 106967.442, 1.1989E-10, 1.2267E-10,
    101279.3704, 1.148E-10, 1.1794E-10, 95607.833, 1.0971E-10, 1.1322E-10,
    89958.3413, 1.0462E-10, 1.0849E-10, 84336.4067, 9.953E-11, 1.0377E-10,
    78747.5407, 9.4441E-11, 9.904E-11, 73197.2546, 8.9351E-11, 9.4314E-11,
    67691.06, 8.4261E-11, 8.9588E-11, 62234.4682, 7.9172E-11, 8.4862E-11,
    56832.9906, 7.4082E-11, 8.0136E-11, 51492.1387, 6.8992E-11, 7.541E-11,
    46217.4239, 6.3903E-11, 7.0684E-11, 41014.3576, 5.8813E-11, 6.5958E-11,
    35888.4513, 5.3723E-11, 6.1232E-11, 30845.2164, 4.8633E-11, 5.6505E-11,
    25890.1642, 4.3544E-11, 5.1779E-11, 21028.8063, 3.8454E-11, 4.7053E-11,
    16266.654, 3.3364E-11, 4.2327E-11, 1.3007E-10, 138371.1425, 4.3201E-10,
    1.272E-10, 112666.5366, 4.1299E-10, 1.2433E-10, 106967.442, 3.9397E-10,
    1.2146E-10, 101279.3704, 3.7494E-10, 1.186E-10, 95607.833, 3.5592E-10,
    1.1573E-10, 89958.3413, 3.369E-10, 1.1286E-10, 84336.4067, 3.1787E-10,
    1.0999E-10, 78747.5407, 2.9885E-10, 1.0712E-10, 73197.2546, 2.7983E-10,
    1.0425E-10, 67691.06, 2.608E-10, 1.0139E-10, 62234.4682, 2.4178E-10,
    9.8518E-11, 56832.9906, 2.2276E-10, 9.5649E-11, 51492.1387, 2.0373E-10,
    9.2781E-11, 46217.4239, 1.8471E-10, 8.9913E-11, 41014.3576, 1.6569E-10,
    8.7044E-11, 35888.4513, 1.4666E-10, 8.4176E-11, 30845.2164, 1.2764E-10,
    8.1308E-11, 25890.1642, 1.0861E-10, 7.844E-11, 21028.8063, 8.9591E-11,
    7.5571E-11, 16266.654, 7.0567E-11, 1.3212E-10, 4.3201E-10, 138371.1425,
    1.2659E-10, 4.1439E-10, 112666.5366, 1.2105E-10, 3.9677E-10, 106967.442,
    1.1552E-10, 3.7915E-10, 101279.3704, 1.0999E-10, 3.6152E-10, 95607.833,
    1.0445E-10, 3.439E-10, 89958.3413, 9.8919E-11, 3.2628E-10, 84336.4067,
    9.3385E-11, 3.0865E-10, 78747.5407, 8.7851E-11, 2.9103E-10, 73197.2546,
    8.2317E-11, 2.7341E-10, 67691.06, 7.6783E-11, 2.5578E-10, 62234.4682,
    7.1249E-11, 2.3816E-10, 56832.9906, 6.5715E-11, 2.2054E-10, 51492.1387,
    6.0181E-11, 2.0292E-10, 46217.4239, 5.4647E-11, 1.8529E-10, 41014.3576,
    4.9113E-11, 1.6767E-10, 35888.4513, 4.3579E-11, 1.5005E-10, 30845.2164,
    3.8045E-11, 1.3242E-10, 25890.1642, 3.2511E-11, 1.148E-10, 21028.8063,
    2.6977E-11, 9.7177E-11, 16266.654, 112666.5366, 1.272E-10, 1.2659E-10,
    129154.3028, 1.2224E-10, 1.2208E-10, 103675.0808, 1.1727E-10, 1.1756E-10,
    98201.3702, 1.1231E-10, 1.1305E-10, 92738.6825, 1.0734E-10, 1.0854E-10,
    87292.529, 1.0238E-10, 1.0402E-10, 81868.4212, 9.7413E-11, 9.9509E-11,
    76471.8706, 9.2448E-11, 9.4995E-11, 71108.3885, 8.7483E-11, 9.0482E-11,
    65783.4864, 8.2519E-11, 8.5969E-11, 60502.6757, 7.7554E-11, 8.1455E-11,
    55271.4678, 7.259E-11, 7.6942E-11, 50095.3742, 6.7625E-11, 7.2429E-11,
    44979.9063, 6.2661E-11, 6.7915E-11, 39930.5754, 5.7696E-11, 6.3402E-11,
    34952.8931, 5.2731E-11, 5.8889E-11, 30052.3707, 4.7767E-11, 5.4375E-11,
    25234.5197, 4.2802E-11, 4.9862E-11, 20504.8515, 3.7838E-11, 4.5348E-11,
    15868.8775, 3.2873E-11, 4.0835E-11, 1.2498E-10, 112666.5366, 4.1439E-10,
    1.2224E-10, 129154.3028, 3.9615E-10, 1.1949E-10, 103675.0808, 3.779E-10,
    1.1675E-10, 98201.3702, 3.5966E-10, 1.1401E-10, 92738.6825, 3.4141E-10,
    1.1126E-10, 87292.529, 3.2317E-10, 1.0852E-10, 81868.4212, 3.0492E-10,
    1.0578E-10, 76471.8706, 2.8668E-10, 1.0303E-10, 71108.3885, 2.6843E-10,
    1.0029E-10, 65783.4864, 2.5019E-10, 9.7547E-11, 60502.6757, 2.3194E-10,
    9.4804E-11, 55271.4678, 2.137E-10, 9.2061E-11, 50095.3742, 1.9545E-10,
    8.9318E-11, 44979.9063, 1.7721E-10, 8.6574E-11, 39930.5754, 1.5896E-10,
    8.3831E-11, 34952.8931, 1.4072E-10, 8.1088E-11, 30052.3707, 1.2247E-10,
    7.8345E-11, 25234.5197, 1.0423E-10, 7.5602E-11, 20504.8515, 8.5983E-11,
    7.2858E-11, 15868.8775, 6.7738E-11, 1.274E-10, 4.1299E-10, 112666.5366,
    1.2208E-10, 3.9615E-10, 129154.3028, 1.1675E-10, 3.793E-10, 103675.0808,
    1.1143E-10, 3.6246E-10, 98201.3702, 1.0611E-10, 3.4561E-10, 92738.6825,
    1.0079E-10, 3.2877E-10, 87292.529, 9.5469E-11, 3.1192E-10, 81868.4212,
    9.0148E-11, 2.9508E-10, 76471.8706, 8.4827E-11, 2.7824E-10, 71108.3885,
    7.9506E-11, 2.6139E-10, 65783.4864, 7.4184E-11, 2.4455E-10, 60502.6757,
    6.8863E-11, 2.277E-10, 55271.4678, 6.3542E-11, 2.1086E-10, 50095.3742,
    5.8221E-11, 1.9401E-10, 44979.9063, 5.2899E-11, 1.7717E-10, 39930.5754,
    4.7578E-11, 1.6032E-10, 34952.8931, 4.2257E-11, 1.4348E-10, 30052.3707,
    3.6936E-11, 1.2664E-10, 25234.5197, 3.1614E-11, 1.0979E-10, 20504.8515,
    2.6293E-11, 9.2947E-11, 15868.8775, 106967.442, 1.2433E-10, 1.2105E-10,
    103675.0808, 1.1949E-10, 1.1675E-10, 120382.7195, 1.1465E-10, 1.1245E-10,
    95123.37, 1.0981E-10, 1.0815E-10, 89869.532, 1.0497E-10, 1.0385E-10,
    84626.7167, 1.0013E-10, 9.9552E-11, 79400.4358, 9.5295E-11, 9.5251E-11,
    74196.2006, 9.0455E-11, 9.0951E-11, 69019.5224, 8.5616E-11, 8.665E-11,
    63875.9129, 8.0776E-11, 8.2349E-11, 58770.8833, 7.5937E-11, 7.8049E-11,
    53709.9451, 7.1097E-11, 7.3748E-11, 48698.6097, 6.6258E-11, 6.9447E-11,
    43742.3886, 6.1418E-11, 6.5147E-11, 38846.7932, 5.6579E-11, 6.0846E-11,
    34017.3348, 5.1739E-11, 5.6546E-11, 29259.525, 4.69E-11, 5.2245E-11,
    24578.8751, 4.206E-11, 4.7944E-11, 19980.8967, 3.7221E-11, 4.3644E-11,
    15471.101, 3.2381E-11, 3.9343E-11, 1.1989E-10, 106967.442, 3.9677E-10,
    1.1727E-10, 103675.0808, 3.793E-10, 1.1465E-10, 120382.7195, 3.6184E-10,
    1.1203E-10, 95123.37, 3.4437E-10, 1.0942E-10, 89869.532, 3.269E-10,
    1.068E-10, 84626.7167, 3.0944E-10, 1.0418E-10, 79400.4358, 2.9197E-10,
    1.0156E-10, 74196.2006, 2.745E-10, 9.8944E-11, 69019.5224, 2.5704E-10,
    9.6326E-11, 63875.9129, 2.3957E-10, 9.3708E-11, 58770.8833, 2.2211E-10,
    9.109E-11, 53709.9451, 2.0464E-10, 8.8472E-11, 48698.6097, 1.8717E-10,
    8.5854E-11, 43742.3886, 1.6971E-10, 8.3236E-11, 38846.7932, 1.5224E-10,
    8.0618E-11, 34017.3348, 1.3477E-10, 7.8E-11, 29259.525, 1.1731E-10,
    7.5382E-11, 24578.8751, 9.9841E-11, 7.2764E-11, 19980.8967, 8.2375E-11,
    7.0146E-11, 15471.101, 6.4909E-11, 1.2267E-10, 3.9397E-10, 106967.442,
    1.1756E-10, 3.779E-10, 103675.0808, 1.1245E-10, 3.6184E-10, 120382.7195,
    1.0735E-10, 3.4577E-10, 95123.37, 1.0224E-10, 3.297E-10, 89869.532,
    9.7128E-11, 3.1364E-10, 84626.7167, 9.202E-11, 2.9757E-10, 79400.4358,
    8.6911E-11, 2.8151E-10, 74196.2006, 8.1803E-11, 2.6544E-10, 69019.5224,
    7.6694E-11, 2.4938E-10, 63875.9129, 7.1586E-11, 2.3331E-10, 58770.8833,
    6.6477E-11, 2.1724E-10, 53709.9451, 6.1368E-11, 2.0118E-10, 48698.6097,
    5.626E-11, 1.8511E-10, 43742.3886, 5.1151E-11, 1.6905E-10, 38846.7932,
    4.6043E-11, 1.5298E-10, 34017.3348, 4.0934E-11, 1.3691E-10, 29259.525,
    3.5826E-11, 1.2085E-10, 24578.8751, 3.0717E-11, 1.0478E-10, 19980.8967,
    2.5609E-11, 8.8717E-11, 15471.101, 101279.3704, 1.2146E-10, 1.1552E-10,
    98201.3702, 1.1675E-10, 1.1143E-10, 95123.37, 1.1203E-10, 1.0735E-10,
    112045.3699, 1.0732E-10, 1.0326E-10, 87000.3814, 1.0261E-10, 9.9169E-11,
    81960.9045, 9.7892E-11, 9.5081E-11, 76932.4503, 9.3177E-11, 9.0994E-11,
    71920.5305, 8.8463E-11, 8.6906E-11, 66930.6563, 8.3748E-11, 8.2818E-11,
    61968.3393, 7.9034E-11, 7.873E-11, 57039.0908, 7.432E-11, 7.4642E-11,
    52148.4223, 6.9605E-11, 7.0554E-11, 47301.8452, 6.4891E-11, 6.6466E-11,
    42504.871, 6.0176E-11, 6.2378E-11, 37763.0109, 5.5462E-11, 5.829E-11,
    33081.7766, 5.0748E-11, 5.4203E-11, 28466.6793, 4.6033E-11, 5.0115E-11,
    23923.2306, 4.1319E-11, 4.6027E-11, 19456.9418, 3.6604E-11, 4.1939E-11,
    15073.3244, 3.189E-11, 3.7851E-11, 1.148E-10, 101279.3704, 3.7915E-10,
    1.1231E-10, 98201.3702, 3.6246E-10, 1.0981E-10, 95123.37, 3.4577E-10,
    1.0732E-10, 112045.3699, 3.2908E-10, 1.0483E-10, 87000.3814, 3.1239E-10,
    1.0233E-10, 81960.9045, 2.9571E-10, 9.9842E-11, 76932.4503, 2.7902E-10,
    9.7349E-11, 71920.5305, 2.6233E-10, 9.4856E-11, 66930.6563, 2.4564E-10,
    9.2363E-11, 61968.3393, 2.2896E-10, 8.987E-11, 57039.0908, 2.1227E-10,
    8.7377E-11, 52148.4223, 1.9558E-10, 8.4884E-11, 47301.8452, 1.7889E-10,
    8.2391E-11, 42504.871, 1.6221E-10, 7.9898E-11, 37763.0109, 1.4552E-10,
    7.7405E-11, 33081.7766, 1.2883E-10, 7.4912E-11, 28466.6793, 1.1214E-10,
    7.2419E-11, 23923.2306, 9.5455E-11, 6.9926E-11, 19456.9418, 7.8767E-11,
    6.7433E-11, 15073.3244, 6.2079E-11, 1.1794E-10, 3.7494E-10, 101279.3704,
    1.1305E-10, 3.5966E-10, 98201.3702, 1.0815E-10, 3.4437E-10, 95123.37,
    1.0326E-10, 3.2908E-10, 112045.3699, 9.8361E-11, 3.138E-10, 87000.3814,
    9.3466E-11, 2.9851E-10, 81960.9045, 8.857E-11, 2.8322E-10, 76932.4503,
    8.3674E-11, 2.6793E-10, 71920.5305, 7.8778E-11, 2.5265E-10, 66930.6563,
    7.3883E-11, 2.3736E-10, 61968.3393, 6.8987E-11, 2.2207E-10, 57039.0908,
    6.4091E-11, 2.0678E-10, 52148.4223, 5.9195E-11, 1.915E-10, 47301.8452,
    5.4299E-11, 1.7621E-10, 42504.871, 4.9404E-11, 1.6092E-10, 37763.0109,
    4.4508E-11, 1.4564E-10, 33081.7766, 3.9612E-11, 1.3035E-10, 28466.6793,
    3.4716E-11, 1.1506E-10, 23923.2306, 2.982E-11, 9.9774E-11, 19456.9418,
    2.4925E-11, 8.4487E-11, 15073.3244, 95607.833, 1.186E-10, 1.0999E-10,
    92738.6825, 1.1401E-10, 1.0611E-10, 89869.532, 1.0942E-10, 1.0224E-10,
    87000.3814, 1.0483E-10, 9.8361E-11, 104131.2309, 1.0024E-10, 9.4486E-11,
    79295.0922, 9.5649E-11, 9.0611E-11, 74464.4649, 9.106E-11, 8.6736E-11,
    69644.8604, 8.647E-11, 8.2861E-11, 64841.7902, 8.1881E-11, 7.8986E-11,
    60060.7658, 7.7292E-11, 7.5111E-11, 55307.2984, 7.2702E-11, 7.1235E-11,
    50586.8996, 6.8113E-11, 6.736E-11, 45905.0807, 6.3524E-11, 6.3485E-11,
    41267.3533, 5.8934E-11, 5.961E-11, 36679.2287, 5.4345E-11, 5.5735E-11,
    32146.2183, 4.9756E-11, 5.186E-11, 27673.8337, 4.5166E-11, 4.7984E-11,
    23267.5861, 4.0577E-11, 4.4109E-11, 18932.987, 3.5988E-11, 4.0234E-11,
    14675.5479, 3.1398E-11, 3.6359E-11, 1.0971E-10, 95607.833, 3.6152E-10,
    1.0734E-10, 92738.6825, 3.4561E-10, 1.0497E-10, 89869.532, 3.297E-10,
    1.0261E-10, 87000.3814, 3.138E-10, 1.0024E-10, 104131.2309, 2.9789E-10,
    9.787E-11, 79295.0922, 2.8198E-10, 9.5502E-11, 74464.4649, 2.6607E-10,
    9.3135E-11, 69644.8604, 2.5016E-10, 9.0767E-11, 64841.7902, 2.3425E-10,
    8.8399E-11, 60060.7658, 2.1834E-10, 8.6031E-11, 55307.2984, 2.0243E-10,
    8.3663E-11, 50586.8996, 1.8652E-10, 8.1295E-11, 45905.0807, 1.7061E-10,
    7.8927E-11, 41267.3533, 1.547E-10, 7.6559E-11, 36679.2287, 1.388E-10,
    7.4191E-11, 32146.2183, 1.2289E-10, 7.1823E-11, 27673.8337, 1.0698E-10,
    6.9456E-11, 23267.5861, 9.1068E-11, 6.7088E-11, 18932.987, 7.5159E-11,
    6.472E-11, 14675.5479, 5.925E-11, 1.1322E-10, 3.5592E-10, 95607.833,
    1.0854E-10, 3.4141E-10, 92738.6825, 1.0385E-10, 3.269E-10, 89869.532,
    9.9169E-11, 3.1239E-10, 87000.3814, 9.4486E-11, 2.9789E-10, 104131.2309,
    8.9803E-11, 2.8338E-10, 79295.0922, 8.512E-11, 2.6887E-10, 74464.4649,
    8.0437E-11, 2.5436E-10, 69644.8604, 7.5754E-11, 2.3985E-10, 64841.7902,
    7.1071E-11, 2.2534E-10, 60060.7658, 6.6388E-11, 2.1083E-10, 55307.2984,
    6.1705E-11, 1.9633E-10, 50586.8996, 5.7022E-11, 1.8182E-10, 45905.0807,
    5.2339E-11, 1.6731E-10, 41267.3533, 4.7656E-11, 1.528E-10, 36679.2287,
    4.2973E-11, 1.3829E-10, 32146.2183, 3.829E-11, 1.2378E-10, 27673.8337,
    3.3607E-11, 1.0927E-10, 23267.5861, 2.8924E-11, 9.4766E-11, 18932.987,
    2.4241E-11, 8.0257E-11, 14675.5479, 89958.3413, 1.1573E-10, 1.0445E-10,
    87292.529, 1.1126E-10, 1.0079E-10, 84626.7167, 1.068E-10, 9.7128E-11,
    81960.9045, 1.0233E-10, 9.3466E-11, 79295.0922, 9.787E-11, 8.9803E-11,
    96629.2799, 9.3406E-11, 8.6141E-11, 71996.4794, 8.8942E-11, 8.2478E-11,
    67369.1904, 8.4478E-11, 7.8816E-11, 62752.9241, 8.0013E-11, 7.5154E-11,
    58153.1922, 7.5549E-11, 7.1491E-11, 53575.5059, 7.1085E-11, 6.7829E-11,
    49025.3768, 6.6621E-11, 6.4166E-11, 44508.3163, 6.2157E-11, 6.0504E-11,
    40029.8357, 5.7692E-11, 5.6841E-11, 35595.4465, 5.3228E-11, 5.3179E-11,
    31210.6601, 4.8764E-11, 4.9517E-11, 26880.988, 4.43E-11, 4.5854E-11,
    22611.9415, 3.9835E-11, 4.2192E-11, 18409.0322, 3.5371E-11, 3.8529E-11,
    14277.7714, 3.0907E-11, 3.4867E-11, 1.0462E-10, 89958.3413, 3.439E-10,
    1.0238E-10, 87292.529, 3.2877E-10, 1.0013E-10, 84626.7167, 3.1364E-10,
    9.7892E-11, 81960.9045, 2.9851E-10, 9.5649E-11, 79295.0922, 2.8338E-10,
    9.3406E-11, 96629.2799, 2.6825E-10, 9.1163E-11, 71996.4794, 2.5312E-10,
    8.892E-11, 67369.1904, 2.3799E-10, 8.6678E-11, 62752.9241, 2.2286E-10,
    8.4435E-11, 58153.1922, 2.0773E-10, 8.2192E-11, 53575.5059, 1.9259E-10,
    7.9949E-11, 49025.3768, 1.7746E-10, 7.7706E-11, 44508.3163, 1.6233E-10,
    7.5464E-11, 40029.8357, 1.472E-10, 7.3221E-11, 35595.4465, 1.3207E-10,
    7.0978E-11, 31210.6601, 1.1694E-10, 6.8735E-11, 26880.988, 1.0181E-10,
    6.6492E-11, 22611.9415, 8.6682E-11, 6.425E-11, 18409.0322, 7.1551E-11,
    6.2007E-11, 14277.7714, 5.6421E-11, 1.0849E-10, 3.369E-10, 89958.3413,
    1.0402E-10, 3.2317E-10, 87292.529, 9.9552E-11, 3.0944E-10, 84626.7167,
    9.5081E-11, 2.9571E-10, 81960.9045, 9.0611E-11, 2.8198E-10, 79295.0922,
    8.6141E-11, 2.6825E-10, 96629.2799, 8.1671E-11, 2.5452E-10, 71996.4794,
    7.72E-11, 2.4079E-10, 67369.1904, 7.273E-11, 2.2706E-10, 62752.9241,
    6.826E-11, 2.1333E-10, 58153.1922, 6.3789E-11, 1.996E-10, 53575.5059,
    5.9319E-11, 1.8587E-10, 49025.3768, 5.4849E-11, 1.7214E-10, 44508.3163,
    5.0378E-11, 1.5841E-10, 40029.8357, 4.5908E-11, 1.4468E-10, 35595.4465,
    4.1438E-11, 1.3095E-10, 31210.6601, 3.6967E-11, 1.1722E-10, 26880.988,
    3.2497E-11, 1.0349E-10, 22611.9415, 2.8027E-11, 8.9758E-11, 18409.0322,
    2.3556E-11, 7.6028E-11, 14277.7714, 84336.4067, 1.1286E-10, 9.8919E-11,
    81868.4212, 1.0852E-10, 9.5469E-11, 79400.4358, 1.0418E-10, 9.202E-11,
    76932.4503, 9.9842E-11, 8.857E-11, 74464.4649, 9.5502E-11, 8.512E-11,
    71996.4794, 9.1163E-11, 8.1671E-11, 89528.494, 8.6824E-11, 7.8221E-11,
    65093.5203, 8.2485E-11, 7.4771E-11, 60664.058, 7.8146E-11, 7.1321E-11,
    56245.6186, 7.3807E-11, 6.7872E-11, 51843.7135, 6.9468E-11, 6.4422E-11,
    47463.8541, 6.5128E-11, 6.0972E-11, 43111.5518, 6.0789E-11, 5.7523E-11,
    38792.318, 5.645E-11, 5.4073E-11, 34511.6642, 5.2111E-11, 5.0623E-11,
    30275.1018, 4.7772E-11, 4.7174E-11, 26088.1423, 4.3433E-11, 4.3724E-11,
    21956.297, 3.9094E-11, 4.0274E-11, 17885.0774, 3.4755E-11, 3.6825E-11,
    13879.9948, 3.0415E-11, 3.3375E-11, 9.953E-11, 84336.4067, 3.2628E-10,
    9.7413E-11, 81868.4212, 3.1192E-10, 9.5295E-11, 79400.4358, 2.9757E-10,
    9.3177E-11, 76932.4503, 2.8322E-10, 9.106E-11, 74464.4649, 2.6887E-10,
    8.8942E-11, 71996.4794, 2.5452E-10, 8.6824E-11, 89528.494, 2.4017E-10,
    8.4706E-11, 65093.5203, 2.2581E-10, 8.2589E-11, 60664.058, 2.1146E-10,
    8.0471E-11, 56245.6186, 1.9711E-10, 7.8353E-11, 51843.7135, 1.8276E-10,
    7.6236E-11, 47463.8541, 1.6841E-10, 7.4118E-11, 43111.5518, 1.5405E-10,
    7.2E-11, 38792.318, 1.397E-10, 6.9882E-11, 34511.6642, 1.2535E-10,
    6.7765E-11, 30275.1018, 1.11E-10, 6.5647E-11, 26088.1423, 9.6647E-11,
    6.3529E-11, 21956.297, 8.2295E-11, 6.1412E-11, 17885.0774, 6.7943E-11,
    5.9294E-11, 13879.9948, 5.3591E-11, 1.0377E-10, 3.1787E-10, 84336.4067,
    9.9509E-11, 3.0492E-10, 81868.4212, 9.5251E-11, 2.9197E-10, 79400.4358,
    9.0994E-11, 2.7902E-10, 76932.4503, 8.6736E-11, 2.6607E-10, 74464.4649,
    8.2478E-11, 2.5312E-10, 71996.4794, 7.8221E-11, 2.4017E-10, 89528.494,
    7.3963E-11, 2.2721E-10, 65093.5203, 6.9706E-11, 2.1426E-10, 60664.058,
    6.5448E-11, 2.0131E-10, 56245.6186, 6.119E-11, 1.8836E-10, 51843.7135,
    5.6933E-11, 1.7541E-10, 47463.8541, 5.2675E-11, 1.6246E-10, 43111.5518,
    4.8418E-11, 1.4951E-10, 38792.318, 4.416E-11, 1.3655E-10, 34511.6642,
    3.9903E-11, 1.236E-10, 30275.1018, 3.5645E-11, 1.1065E-10, 26088.1423,
    3.1387E-11, 9.77E-11, 21956.297, 2.713E-11, 8.4749E-11, 17885.0774,
    2.2872E-11, 7.1798E-11, 13879.9948, 78747.5407, 1.0999E-10, 9.3385E-11,
    76471.8706, 1.0578E-10, 9.0148E-11, 74196.2006, 1.0156E-10, 8.6911E-11,
    71920.5305, 9.7349E-11, 8.3674E-11, 69644.8604, 9.3135E-11, 8.0437E-11,
    67369.1904, 8.892E-11, 7.72E-11, 65093.5203, 8.4706E-11, 7.3963E-11,
    82817.8502, 8.0492E-11, 7.0726E-11, 58575.192, 7.6278E-11, 6.7489E-11,
    54338.0451, 7.2064E-11, 6.4252E-11, 50111.9211, 6.785E-11, 6.1015E-11,
    45902.3313, 6.3636E-11, 5.7778E-11, 41714.7873, 5.9422E-11, 5.4542E-11,
    37554.8004, 5.5208E-11, 5.1305E-11, 33427.882, 5.0994E-11, 4.8068E-11,
    29339.5436, 4.678E-11, 4.4831E-11, 25295.2966, 4.2566E-11, 4.1594E-11,
    21300.6524, 3.8352E-11, 3.8357E-11, 17361.1225, 3.4138E-11, 3.512E-11,
    13482.2183, 2.9924E-11, 3.1883E-11, 9.4441E-11, 78747.5407, 3.0865E-10,
    9.2448E-11, 76471.8706, 2.9508E-10, 9.0455E-11, 74196.2006, 2.8151E-10,
    8.8463E-11, 71920.5305, 2.6793E-10, 8.647E-11, 69644.8604, 2.5436E-10,
    8.4478E-11, 67369.1904, 2.4079E-10, 8.2485E-11, 65093.5203, 2.2721E-10,
    8.0492E-11, 82817.8502, 2.1364E-10, 7.85E-11, 58575.192, 2.0007E-10,
    7.6507E-11, 54338.0451, 1.8649E-10, 7.4515E-11, 50111.9211, 1.7292E-10,
    7.2522E-11, 45902.3313, 1.5935E-10, 7.0529E-11, 41714.7873, 1.4577E-10,
    6.8537E-11, 37554.8004, 1.322E-10, 6.6544E-11, 33427.882, 1.1863E-10,
    6.4551E-11, 29339.5436, 1.0505E-10, 6.2559E-11, 25295.2966, 9.1482E-11,
    6.0566E-11, 21300.6524, 7.7908E-11, 5.8574E-11, 17361.1225, 6.4335E-11,
    5.6581E-11, 13482.2183, 5.0762E-11, 9.904E-11, 2.9885E-10, 78747.5407,
    9.4995E-11, 2.8668E-10, 76471.8706, 9.0951E-11, 2.745E-10, 74196.2006,
    8.6906E-11, 2.6233E-10, 71920.5305, 8.2861E-11, 2.5016E-10, 69644.8604,
    7.8816E-11, 2.3799E-10, 67369.1904, 7.4771E-11, 2.2581E-10, 65093.5203,
    7.0726E-11, 2.1364E-10, 82817.8502, 6.6681E-11, 2.0147E-10, 58575.192,
    6.2637E-11, 1.893E-10, 54338.0451, 5.8592E-11, 1.7712E-10, 50111.9211,
    5.4547E-11, 1.6495E-10, 45902.3313, 5.0502E-11, 1.5278E-10, 41714.7873,
    4.6457E-11, 1.406E-10, 37554.8004, 4.2412E-11, 1.2843E-10, 33427.882,
    3.8367E-11, 1.1626E-10, 29339.5436, 3.4323E-11, 1.0409E-10, 25295.2966,
    3.0278E-11, 9.1913E-11, 21300.6524, 2.6233E-11, 7.9741E-11, 17361.1225,
    2.2188E-11, 6.7568E-11, 13482.2183, 73197.2546, 1.0712E-10, 8.7851E-11,
    71108.3885, 1.0303E-10, 8.4827E-11, 69019.5224, 9.8944E-11, 8.1803E-11,
    66930.6563, 9.4856E-11, 7.8778E-11, 64841.7902, 9.0767E-11, 7.5754E-11,
    62752.9241, 8.6678E-11, 7.273E-11, 60664.058, 8.2589E-11, 6.9706E-11,
    58575.192, 7.85E-11, 6.6681E-11, 76486.3259, 7.4411E-11, 6.3657E-11,
    52430.4715, 7.0322E-11, 6.0633E-11, 48380.1286, 6.6233E-11, 5.7609E-11,
    44340.8086, 6.2144E-11, 5.4585E-11, 40318.0228, 5.8055E-11, 5.156E-11,
    36317.2827, 5.3966E-11, 4.8536E-11, 32344.0998, 4.9877E-11, 4.5512E-11,
    28403.9854, 4.5788E-11, 4.2488E-11, 24502.4509, 4.1699E-11, 3.9463E-11,
    20645.0079, 3.761E-11, 3.6439E-11, 16837.1677, 3.3521E-11, 3.3415E-11,
    13084.4418, 2.9432E-11, 3.0391E-11, 8.9351E-11, 73197.2546, 2.9103E-10,
    8.7483E-11, 71108.3885, 2.7824E-10, 8.5616E-11, 69019.5224, 2.6544E-10,
    8.3748E-11, 66930.6563, 2.5265E-10, 8.1881E-11, 64841.7902, 2.3985E-10,
    8.0013E-11, 62752.9241, 2.2706E-10, 7.8146E-11, 60664.058, 2.1426E-10,
    7.6278E-11, 58575.192, 2.0147E-10, 7.4411E-11, 76486.3259, 1.8867E-10,
    7.2543E-11, 52430.4715, 1.7588E-10, 7.0676E-11, 48380.1286, 1.6308E-10,
    6.8808E-11, 44340.8086, 1.5029E-10, 6.6941E-11, 40318.0228, 1.3749E-10,
    6.5073E-11, 36317.2827, 1.247E-10, 6.3206E-11, 32344.0998, 1.1191E-10,
    6.1338E-11, 28403.9854, 9.9111E-11, 5.9471E-11, 24502.4509, 8.6316E-11,
    5.7603E-11, 20645.0079, 7.3522E-11, 5.5736E-11, 16837.1677, 6.0727E-11,
    5.3868E-11, 13084.4418, 4.7933E-11, 9.4314E-11, 2.7983E-10, 73197.2546,
    9.0482E-11, 2.6843E-10, 71108.3885, 8.665E-11, 2.5704E-10, 69019.5224,
    8.2818E-11, 2.4564E-10, 66930.6563, 7.8986E-11, 2.3425E-10, 64841.7902,
    7.5154E-11, 2.2286E-10, 62752.9241, 7.1321E-11, 2.1146E-10, 60664.058,
    6.7489E-11, 2.0007E-10, 58575.192, 6.3657E-11, 1.8867E-10, 76486.3259,
    5.9825E-11, 1.7728E-10, 52430.4715, 5.5993E-11, 1.6589E-10, 48380.1286,
    5.2161E-11, 1.5449E-10, 44340.8086, 4.8329E-11, 1.431E-10, 40318.0228,
    4.4497E-11, 1.317E-10, 36317.2827, 4.0664E-11, 1.2031E-10, 32344.0998,
    3.6832E-11, 1.0891E-10, 28403.9854, 3.3E-11, 9.752E-11, 24502.4509,
    2.9168E-11, 8.6126E-11, 20645.0079, 2.5336E-11, 7.4732E-11, 16837.1677,
    2.1504E-11, 6.3338E-11, 13084.4418, 67691.06, 1.0425E-10, 8.2317E-11,
    65783.4864, 1.0029E-10, 7.9506E-11, 63875.9129, 9.6326E-11, 7.6694E-11,
    61968.3393, 9.2363E-11, 7.3883E-11, 60060.7658, 8.8399E-11, 7.1071E-11,
    58153.1922, 8.4435E-11, 6.826E-11, 56245.6186, 8.0471E-11, 6.5448E-11,
    54338.0451, 7.6507E-11, 6.2637E-11, 52430.4715, 7.2543E-11, 5.9825E-11,
    70522.898, 6.8579E-11, 5.7014E-11, 46648.3362, 6.4616E-11, 5.4202E-11,
    42779.2858, 6.0652E-11, 5.1391E-11, 38921.2583, 5.6688E-11, 4.8579E-11,
    35079.7651, 5.2724E-11, 4.5768E-11, 31260.3175, 4.876E-11, 4.2956E-11,
    27468.4271, 4.4796E-11, 4.0145E-11, 23709.6053, 4.0833E-11, 3.7333E-11,
    19989.3634, 3.6869E-11, 3.4522E-11, 16313.2129, 3.2905E-11, 3.171E-11,
    12686.6652, 2.8941E-11, 2.8899E-11, 8.4261E-11, 67691.06, 2.7341E-10,
    8.2519E-11, 65783.4864, 2.6139E-10, 8.0776E-11, 63875.9129, 2.4938E-10,
    7.9034E-11, 61968.3393, 2.3736E-10, 7.7292E-11, 60060.7658, 2.2534E-10,
    7.5549E-11, 58153.1922, 2.1333E-10, 7.3807E-11, 56245.6186, 2.0131E-10,
    7.2064E-11, 54338.0451, 1.893E-10, 7.0322E-11, 52430.4715, 1.7728E-10,
    6.8579E-11, 70522.898, 1.6526E-10, 6.6837E-11, 46648.3362, 1.5325E-10,
    6.5095E-11, 42779.2858, 1.4123E-10, 6.3352E-11, 38921.2583, 1.2922E-10,
    6.161E-11, 35079.7651, 1.172E-10, 5.9867E-11, 31260.3175, 1.0518E-10,
    5.8125E-11, 27468.4271, 9.3167E-11, 5.6382E-11, 23709.6053, 8.1151E-11,
    5.464E-11, 19989.3634, 6.9135E-11, 5.2898E-11, 16313.2129, 5.7119E-11,
    5.1155E-11, 12686.6652, 4.5103E-11, 8.9588E-11, 2.608E-10, 67691.06,
    8.5969E-11, 2.5019E-10, 65783.4864, 8.2349E-11, 2.3957E-10, 63875.9129,
    7.873E-11, 2.2896E-10, 61968.3393, 7.5111E-11, 2.1834E-10, 60060.7658,
    7.1491E-11, 2.0773E-10, 58153.1922, 6.7872E-11, 1.9711E-10, 56245.6186,
    6.4252E-11, 1.8649E-10, 54338.0451, 6.0633E-11, 1.7588E-10, 52430.4715,
    5.7014E-11, 1.6526E-10, 70522.898, 5.3394E-11, 1.5465E-10, 46648.3362,
    4.9775E-11, 1.4403E-10, 42779.2858, 4.6155E-11, 1.3342E-10, 38921.2583,
    4.2536E-11, 1.228E-10, 35079.7651, 3.8917E-11, 1.1219E-10, 31260.3175,
    3.5297E-11, 1.0157E-10, 27468.4271, 3.1678E-11, 9.0955E-11, 23709.6053,
    2.8058E-11, 8.0339E-11, 19989.3634, 2.4439E-11, 6.9724E-11, 16313.2129,
    2.082E-11, 5.9108E-11, 12686.6652, 62234.4682, 1.0139E-10, 7.6783E-11,
    60502.6757, 9.7547E-11, 7.4184E-11, 58770.8833, 9.3708E-11, 7.1586E-11,
    57039.0908, 8.987E-11, 6.8987E-11, 55307.2984, 8.6031E-11, 6.6388E-11,
    53575.5059, 8.2192E-11, 6.3789E-11, 51843.7135, 7.8353E-11, 6.119E-11,
    50111.9211, 7.4515E-11, 5.8592E-11, 48380.1286, 7.0676E-11, 5.5993E-11,
    46648.3362, 6.6837E-11, 5.3394E-11, 64916.5437, 6.2998E-11, 5.0795E-11,
    41217.7631, 5.916E-11, 4.8197E-11, 37524.4938, 5.5321E-11, 4.5598E-11,
    33842.2474, 5.1482E-11, 4.2999E-11, 30176.5353, 4.7643E-11, 4.04E-11,
    26532.8689, 4.3805E-11, 3.7802E-11, 22916.7596, 3.9966E-11, 3.5203E-11,
    19333.7188, 3.6127E-11, 3.2604E-11, 15789.2581, 3.2288E-11, 3.0005E-11,
    12288.8887, 2.8449E-11, 2.7407E-11, 7.9172E-11, 62234.4682, 2.5578E-10,
    7.7554E-11, 60502.6757, 2.4455E-10, 7.5937E-11, 58770.8833, 2.3331E-10,
    7.432E-11, 57039.0908, 2.2207E-10, 7.2702E-11, 55307.2984, 2.1083E-10,
    7.1085E-11, 53575.5059, 1.996E-10, 6.9468E-11, 51843.7135, 1.8836E-10,
    6.785E-11, 50111.9211, 1.7712E-10, 6.6233E-11, 48380.1286, 1.6589E-10,
    6.4616E-11, 46648.3362, 1.5465E-10, 6.2998E-11, 64916.5437, 1.4341E-10,
    6.1381E-11, 41217.7631, 1.3217E-10, 5.9764E-11, 37524.4938, 1.2094E-10,
    5.8146E-11, 33842.2474, 1.097E-10, 5.6529E-11, 30176.5353, 9.8461E-11,
    5.4912E-11, 26532.8689, 8.7223E-11, 5.3294E-11, 22916.7596, 7.5986E-11,
    5.1677E-11, 19333.7188, 6.4749E-11, 5.006E-11, 15789.2581, 5.3511E-11,
    4.8442E-11, 12288.8887, 4.2274E-11, 8.4862E-11, 2.4178E-10, 62234.4682,
    8.1455E-11, 2.3194E-10, 60502.6757, 7.8049E-11, 2.2211E-10, 58770.8833,
    7.4642E-11, 2.1227E-10, 57039.0908, 7.1235E-11, 2.0243E-10, 55307.2984,
    6.7829E-11, 1.9259E-10, 53575.5059, 6.4422E-11, 1.8276E-10, 51843.7135,
    6.1015E-11, 1.7292E-10, 50111.9211, 5.7609E-11, 1.6308E-10, 48380.1286,
    5.4202E-11, 1.5325E-10, 46648.3362, 5.0795E-11, 1.4341E-10, 64916.5437,
    4.7389E-11, 1.3357E-10, 41217.7631, 4.3982E-11, 1.2374E-10, 37524.4938,
    4.0575E-11, 1.139E-10, 33842.2474, 3.7169E-11, 1.0406E-10, 30176.5353,
    3.3762E-11, 9.4226E-11, 26532.8689, 3.0356E-11, 8.4389E-11, 22916.7596,
    2.6949E-11, 7.4552E-11, 19333.7188, 2.3542E-11, 6.4715E-11, 15789.2581,
    2.0136E-11, 5.4878E-11, 12288.8887, 56832.9906, 9.8518E-11, 7.1249E-11,
    55271.4678, 9.4804E-11, 6.8863E-11, 53709.9451, 9.109E-11, 6.6477E-11,
    52148.4223, 8.7377E-11, 6.4091E-11, 50586.8996, 8.3663E-11, 6.1705E-11,
    49025.3768, 7.9949E-11, 5.9319E-11, 47463.8541, 7.6236E-11, 5.6933E-11,
    45902.3313, 7.2522E-11, 5.4547E-11, 44340.8086, 6.8808E-11, 5.2161E-11,
    42779.2858, 6.5095E-11, 4.9775E-11, 41217.7631, 6.1381E-11, 4.7389E-11,
    59656.2403, 5.7667E-11, 4.5003E-11, 36127.7293, 5.3954E-11, 4.2617E-11,
    32604.7298, 5.024E-11, 4.0231E-11, 29092.7531, 4.6526E-11, 3.7845E-11,
    25597.3106, 4.2813E-11, 3.5459E-11, 22123.9139, 3.9099E-11, 3.3073E-11,
    18678.0743, 3.5385E-11, 3.0687E-11, 15265.3032, 3.1672E-11, 2.8301E-11,
    11891.1121, 2.7958E-11, 2.5915E-11, 7.4082E-11, 56832.9906, 2.3816E-10,
    7.259E-11, 55271.4678, 2.277E-10, 7.1097E-11, 53709.9451, 2.1724E-10,
    6.9605E-11, 52148.4223, 2.0678E-10, 6.8113E-11, 50586.8996, 1.9633E-10,
    6.6621E-11, 49025.3768, 1.8587E-10, 6.5128E-11, 47463.8541, 1.7541E-10,
    6.3636E-11, 45902.3313, 1.6495E-10, 6.2144E-11, 44340.8086, 1.5449E-10,
    6.0652E-11, 42779.2858, 1.4403E-10, 5.916E-11, 41217.7631, 1.3357E-10,
    5.7667E-11, 59656.2403, 1.2311E-10, 5.6175E-11, 36127.7293, 1.1266E-10,
    5.4683E-11, 32604.7298, 1.022E-10, 5.3191E-11, 29092.7531, 9.1738E-11,
    5.1698E-11, 25597.3106, 8.128E-11, 5.0206E-11, 22123.9139, 7.0821E-11,
    4.8714E-11, 18678.0743, 6.0362E-11, 4.7222E-11, 15265.3032, 4.9903E-11,
    4.5729E-11, 11891.1121, 3.9444E-11, 8.0136E-11, 2.2276E-10, 56832.9906,
    7.6942E-11, 2.137E-10, 55271.4678, 7.3748E-11, 2.0464E-10, 53709.9451,
    7.0554E-11, 1.9558E-10, 52148.4223, 6.736E-11, 1.8652E-10, 50586.8996,
    6.4166E-11, 1.7746E-10, 49025.3768, 6.0972E-11, 1.6841E-10, 47463.8541,
    5.7778E-11, 1.5935E-10, 45902.3313, 5.4585E-11, 1.5029E-10, 44340.8086,
    5.1391E-11, 1.4123E-10, 42779.2858, 4.8197E-11, 1.3217E-10, 41217.7631,
    4.5003E-11, 1.2311E-10, 59656.2403, 4.1809E-11, 1.1406E-10, 36127.7293,
    3.8615E-11, 1.05E-10, 32604.7298, 3.5421E-11, 9.594E-11, 29092.7531,
    3.2227E-11, 8.6882E-11, 25597.3106, 2.9033E-11, 7.7823E-11, 22123.9139,
    2.5839E-11, 6.8765E-11, 18678.0743, 2.2645E-11, 5.9707E-11, 15265.3032,
    1.9451E-11, 5.0648E-11, 11891.1121, 51492.1387, 9.5649E-11, 6.5715E-11,
    50095.3742, 9.2061E-11, 6.3542E-11, 48698.6097, 8.8472E-11, 6.1368E-11,
    47301.8452, 8.4884E-11, 5.9195E-11, 45905.0807, 8.1295E-11, 5.7022E-11,
    44508.3163, 7.7706E-11, 5.4849E-11, 43111.5518, 7.4118E-11, 5.2675E-11,
    41714.7873, 7.0529E-11, 5.0502E-11, 40318.0228, 6.6941E-11, 4.8329E-11,
    38921.2583, 6.3352E-11, 4.6155E-11, 37524.4938, 5.9764E-11, 4.3982E-11,
    36127.7293, 5.6175E-11, 4.1809E-11, 54730.9648, 5.2586E-11, 3.9636E-11,
    31367.2121, 4.8998E-11, 3.7462E-11, 28008.9708, 4.5409E-11, 3.5289E-11,
    24661.7524, 4.1821E-11, 3.3116E-11, 21331.0682, 3.8232E-11, 3.0942E-11,
    18022.4297, 3.4644E-11, 2.8769E-11, 14741.3484, 3.1055E-11, 2.6596E-11,
    11493.3356, 2.7467E-11, 2.4423E-11, 6.8992E-11, 51492.1387, 2.2054E-10,
    6.7625E-11, 50095.3742, 2.1086E-10, 6.6258E-11, 48698.6097, 2.0118E-10,
    6.4891E-11, 47301.8452, 1.915E-10, 6.3524E-11, 45905.0807, 1.8182E-10,
    6.2157E-11, 44508.3163, 1.7214E-10, 6.0789E-11, 43111.5518, 1.6246E-10,
    5.9422E-11, 41714.7873, 1.5278E-10, 5.8055E-11, 40318.0228, 1.431E-10,
    5.6688E-11, 38921.2583, 1.3342E-10, 5.5321E-11, 37524.4938, 1.2374E-10,
    5.3954E-11, 36127.7293, 1.1406E-10, 5.2586E-11, 54730.9648, 1.0438E-10,
    5.1219E-11, 31367.2121, 9.4696E-11, 4.9852E-11, 28008.9708, 8.5016E-11,
    4.8485E-11, 24661.7524, 7.5336E-11, 4.7118E-11, 21331.0682, 6.5656E-11,
    4.5751E-11, 18022.4297, 5.5975E-11, 4.4384E-11, 14741.3484, 4.6295E-11,
    4.3016E-11, 11493.3356, 3.6615E-11, 7.541E-11, 2.0373E-10, 51492.1387,
    7.2429E-11, 1.9545E-10, 50095.3742, 6.9447E-11, 1.8717E-10, 48698.6097,
    6.6466E-11, 1.7889E-10, 47301.8452, 6.3485E-11, 1.7061E-10, 45905.0807,
    6.0504E-11, 1.6233E-10, 44508.3163, 5.7523E-11, 1.5405E-10, 43111.5518,
    5.4542E-11, 1.4577E-10, 41714.7873, 5.156E-11, 1.3749E-10, 40318.0228,
    4.8579E-11, 1.2922E-10, 38921.2583, 4.5598E-11, 1.2094E-10, 37524.4938,
    4.2617E-11, 1.1266E-10, 36127.7293, 3.9636E-11, 1.0438E-10, 54730.9648,
    3.6654E-11, 9.6097E-11, 31367.2121, 3.3673E-11, 8.7817E-11, 28008.9708,
    3.0692E-11, 7.9537E-11, 24661.7524, 2.7711E-11, 7.1257E-11, 21331.0682,
    2.473E-11, 6.2978E-11, 18022.4297, 2.1748E-11, 5.4698E-11, 14741.3484,
    1.8767E-11, 4.6418E-11, 11493.3356, 46217.4239, 9.2781E-11, 6.0181E-11,
    44979.9063, 8.9318E-11, 5.8221E-11, 43742.3886, 8.5854E-11, 5.626E-11,
    42504.871, 8.2391E-11, 5.4299E-11, 41267.3533, 7.8927E-11, 5.2339E-11,
    40029.8357, 7.5464E-11, 5.0378E-11, 38792.318, 7.2E-11, 4.8418E-11,
    37554.8004, 6.8537E-11, 4.6457E-11, 36317.2827, 6.5073E-11, 4.4497E-11,
    35079.7651, 6.161E-11, 4.2536E-11, 33842.2474, 5.8146E-11, 4.0575E-11,
    32604.7298, 5.4683E-11, 3.8615E-11, 31367.2121, 5.1219E-11, 3.6654E-11,
    50129.6945, 4.7756E-11, 3.4694E-11, 26925.1886, 4.4292E-11, 3.2733E-11,
    23726.1941, 4.0829E-11, 3.0773E-11, 20538.2225, 3.7365E-11, 2.8812E-11,
    17366.7852, 3.3902E-11, 2.6852E-11, 14217.3936, 3.0438E-11, 2.4891E-11,
    11095.5591, 2.6975E-11, 2.293E-11, 6.3903E-11, 46217.4239, 2.0292E-10,
    6.2661E-11, 44979.9063, 1.9401E-10, 6.1418E-11, 43742.3886, 1.8511E-10,
    6.0176E-11, 42504.871, 1.7621E-10, 5.8934E-11, 41267.3533, 1.6731E-10,
    5.7692E-11, 40029.8357, 1.5841E-10, 5.645E-11, 38792.318, 1.4951E-10,
    5.5208E-11, 37554.8004, 1.406E-10, 5.3966E-11, 36317.2827, 1.317E-10,
    5.2724E-11, 35079.7651, 1.228E-10, 5.1482E-11, 33842.2474, 1.139E-10,
    5.024E-11, 32604.7298, 1.05E-10, 4.8998E-11, 31367.2121, 9.6097E-11,
    4.7756E-11, 50129.6945, 8.7195E-11, 4.6514E-11, 26925.1886, 7.8293E-11,
    4.5272E-11, 23726.1941, 6.9392E-11, 4.403E-11, 20538.2225, 6.049E-11,
    4.2788E-11, 17366.7852, 5.1589E-11, 4.1546E-11, 14217.3936, 4.2687E-11,
    4.0304E-11, 11095.5591, 3.3786E-11, 7.0684E-11, 1.8471E-10, 46217.4239,
    6.7915E-11, 1.7721E-10, 44979.9063, 6.5147E-11, 1.6971E-10, 43742.3886,
    6.2378E-11, 1.6221E-10, 42504.871, 5.961E-11, 1.547E-10, 41267.3533,
    5.6841E-11, 1.472E-10, 40029.8357, 5.4073E-11, 1.397E-10, 38792.318,
    5.1305E-11, 1.322E-10, 37554.8004, 4.8536E-11, 1.247E-10, 36317.2827,
    4.5768E-11, 1.172E-10, 35079.7651, 4.2999E-11, 1.097E-10, 33842.2474,
    4.0231E-11, 1.022E-10, 32604.7298, 3.7462E-11, 9.4696E-11, 31367.2121,
    3.4694E-11, 8.7195E-11, 50129.6945, 3.1925E-11, 7.9694E-11, 26925.1886,
    2.9157E-11, 7.2193E-11, 23726.1941, 2.6388E-11, 6.4692E-11, 20538.2225,
    2.362E-11, 5.7191E-11, 17366.7852, 2.0852E-11, 4.969E-11, 14217.3936,
    1.8083E-11, 4.2189E-11, 11095.5591, 41014.3576, 8.9913E-11, 5.4647E-11,
    39930.5754, 8.6574E-11, 5.2899E-11, 38846.7932, 8.3236E-11, 5.1151E-11,
    37763.0109, 7.9898E-11, 4.9404E-11, 36679.2287, 7.6559E-11, 4.7656E-11,
    35595.4465, 7.3221E-11, 4.5908E-11, 34511.6642, 6.9882E-11, 4.416E-11,
    33427.882, 6.6544E-11, 4.2412E-11, 32344.0998, 6.3206E-11, 4.0664E-11,
    31260.3175, 5.9867E-11, 3.8917E-11, 30176.5353, 5.6529E-11, 3.7169E-11,
    29092.7531, 5.3191E-11, 3.5421E-11, 28008.9708, 4.9852E-11, 3.3673E-11,
    26925.1886, 4.6514E-11, 3.1925E-11, 45841.4064, 4.3175E-11, 3.0178E-11,
    22790.6359, 3.9837E-11, 2.843E-11, 19745.3769, 3.6499E-11, 2.6682E-11,
    16711.1407, 3.316E-11, 2.4934E-11, 13693.4388, 2.9822E-11, 2.3186E-11,
    10697.7825, 2.6484E-11, 2.1438E-11, 5.8813E-11, 41014.3576, 1.8529E-10,
    5.7696E-11, 39930.5754, 1.7717E-10, 5.6579E-11, 38846.7932, 1.6905E-10,
    5.5462E-11, 37763.0109, 1.6092E-10, 5.4345E-11, 36679.2287, 1.528E-10,
    5.3228E-11, 35595.4465, 1.4468E-10, 5.2111E-11, 34511.6642, 1.3655E-10,
    5.0994E-11, 33427.882, 1.2843E-10, 4.9877E-11, 32344.0998, 1.2031E-10,
    4.876E-11, 31260.3175, 1.1219E-10, 4.7643E-11, 30176.5353, 1.0406E-10,
    4.6526E-11, 29092.7531, 9.594E-11, 4.5409E-11, 28008.9708, 8.7817E-11,
    4.4292E-11, 26925.1886, 7.9694E-11, 4.3175E-11, 45841.4064, 7.1571E-11,
    4.2058E-11, 22790.6359, 6.3448E-11, 4.0941E-11, 19745.3769, 5.5325E-11,
    3.9825E-11, 16711.1407, 4.7202E-11, 3.8708E-11, 13693.4388, 3.9079E-11,
    3.7591E-11, 10697.7825, 3.0956E-11, 6.5958E-11, 1.6569E-10, 41014.3576,
    6.3402E-11, 1.5896E-10, 39930.5754, 6.0846E-11, 1.5224E-10, 38846.7932,
    5.829E-11, 1.4552E-10, 37763.0109, 5.5735E-11, 1.388E-10, 36679.2287,
    5.3179E-11, 1.3207E-10, 35595.4465, 5.0623E-11, 1.2535E-10, 34511.6642,
    4.8068E-11, 1.1863E-10, 33427.882, 4.5512E-11, 1.1191E-10, 32344.0998,
    4.2956E-11, 1.0518E-10, 31260.3175, 4.04E-11, 9.8461E-11, 30176.5353,
    3.7845E-11, 9.1738E-11, 29092.7531, 3.5289E-11, 8.5016E-11, 28008.9708,
    3.2733E-11, 7.8293E-11, 26925.1886, 3.0178E-11, 7.1571E-11, 45841.4064,
    2.7622E-11, 6.4849E-11, 22790.6359, 2.5066E-11, 5.8126E-11, 19745.3769,
    2.251E-11, 5.1404E-11, 16711.1407, 1.9955E-11, 4.4681E-11, 13693.4388,
    1.7399E-11, 3.7959E-11, 10697.7825, 35888.4513, 8.7044E-11, 4.9113E-11,
    34952.8931, 8.3831E-11, 4.7578E-11, 34017.3348, 8.0618E-11, 4.6043E-11,
    33081.7766, 7.7405E-11, 4.4508E-11, 32146.2183, 7.4191E-11, 4.2973E-11,
    31210.6601, 7.0978E-11, 4.1438E-11, 30275.1018, 6.7765E-11, 3.9903E-11,
    29339.5436, 6.4551E-11, 3.8367E-11, 28403.9854, 6.1338E-11, 3.6832E-11,
    27468.4271, 5.8125E-11, 3.5297E-11, 26532.8689, 5.4912E-11, 3.3762E-11,
    25597.3106, 5.1698E-11, 3.2227E-11, 24661.7524, 4.8485E-11, 3.0692E-11,
    23726.1941, 4.5272E-11, 2.9157E-11, 22790.6359, 4.2058E-11, 2.7622E-11,
    41855.0777, 3.8845E-11, 2.6087E-11, 18952.5312, 3.5632E-11, 2.4552E-11,
    16055.4961, 3.2419E-11, 2.3017E-11, 13169.4839, 2.9205E-11, 2.1481E-11,
    10300.006, 2.5992E-11, 1.9946E-11, 5.3723E-11, 35888.4513, 1.6767E-10,
    5.2731E-11, 34952.8931, 1.6032E-10, 5.1739E-11, 34017.3348, 1.5298E-10,
    5.0748E-11, 33081.7766, 1.4564E-10, 4.9756E-11, 32146.2183, 1.3829E-10,
    4.8764E-11, 31210.6601, 1.3095E-10, 4.7772E-11, 30275.1018, 1.236E-10,
    4.678E-11, 29339.5436, 1.1626E-10, 4.5788E-11, 28403.9854, 1.0891E-10,
    4.4796E-11, 27468.4271, 1.0157E-10, 4.3805E-11, 26532.8689, 9.4226E-11,
    4.2813E-11, 25597.3106, 8.6882E-11, 4.1821E-11, 24661.7524, 7.9537E-11,
    4.0829E-11, 23726.1941, 7.2193E-11, 3.9837E-11, 22790.6359, 6.4849E-11,
    3.8845E-11, 41855.0777, 5.7504E-11, 3.7853E-11, 18952.5312, 5.016E-11,
    3.6861E-11, 16055.4961, 4.2816E-11, 3.587E-11, 13169.4839, 3.5471E-11,
    3.4878E-11, 10300.006, 2.8127E-11, 6.1232E-11, 1.4666E-10, 35888.4513,
    5.8889E-11, 1.4072E-10, 34952.8931, 5.6546E-11, 1.3477E-10, 34017.3348,
    5.4203E-11, 1.2883E-10, 33081.7766, 5.186E-11, 1.2289E-10, 32146.2183,
    4.9517E-11, 1.1694E-10, 31210.6601, 4.7174E-11, 1.11E-10, 30275.1018,
    4.4831E-11, 1.0505E-10, 29339.5436, 4.2488E-11, 9.9111E-11, 28403.9854,
    4.0145E-11, 9.3167E-11, 27468.4271, 3.7802E-11, 8.7223E-11, 26532.8689,
    3.5459E-11, 8.128E-11, 25597.3106, 3.3116E-11, 7.5336E-11, 24661.7524,
    3.0773E-11, 6.9392E-11, 23726.1941, 2.843E-11, 6.3448E-11, 22790.6359,
    2.6087E-11, 5.7504E-11, 41855.0777, 2.3744E-11, 5.156E-11, 18952.5312,
    2.1401E-11, 4.5617E-11, 16055.4961, 1.9058E-11, 3.9673E-11, 13169.4839,
    1.6715E-11, 3.3729E-11, 10300.006, 30845.2164, 8.4176E-11, 4.3579E-11,
    30052.3707, 8.1088E-11, 4.2257E-11, 29259.525, 7.8E-11, 4.0934E-11,
    28466.6793, 7.4912E-11, 3.9612E-11, 27673.8337, 7.1823E-11, 3.829E-11,
    26880.988, 6.8735E-11, 3.6967E-11, 26088.1423, 6.5647E-11, 3.5645E-11,
    25295.2966, 6.2559E-11, 3.4323E-11, 24502.4509, 5.9471E-11, 3.3E-11,
    23709.6053, 5.6382E-11, 3.1678E-11, 22916.7596, 5.3294E-11, 3.0356E-11,
    22123.9139, 5.0206E-11, 2.9033E-11, 21331.0682, 4.7118E-11, 2.7711E-11,
    20538.2225, 4.403E-11, 2.6388E-11, 19745.3769, 4.0941E-11, 2.5066E-11,
    18952.5312, 3.7853E-11, 2.3744E-11, 38159.6855, 3.4765E-11, 2.2421E-11,
    15399.8516, 3.1677E-11, 2.1099E-11, 12645.5291, 2.8589E-11, 1.9777E-11,
    9902.2295, 2.5501E-11, 1.8454E-11, 4.8633E-11, 30845.2164, 1.5005E-10,
    4.7767E-11, 30052.3707, 1.4348E-10, 4.69E-11, 29259.525, 1.3691E-10,
    4.6033E-11, 28466.6793, 1.3035E-10, 4.5166E-11, 27673.8337, 1.2378E-10,
    4.43E-11, 26880.988, 1.1722E-10, 4.3433E-11, 26088.1423, 1.1065E-10,
    4.2566E-11, 25295.2966, 1.0409E-10, 4.1699E-11, 24502.4509, 9.752E-11,
    4.0833E-11, 23709.6053, 9.0955E-11, 3.9966E-11, 22916.7596, 8.4389E-11,
    3.9099E-11, 22123.9139, 7.7823E-11, 3.8232E-11, 21331.0682, 7.1257E-11,
    3.7365E-11, 20538.2225, 6.4692E-11, 3.6499E-11, 19745.3769, 5.8126E-11,
    3.5632E-11, 18952.5312, 5.156E-11, 3.4765E-11, 38159.6855, 4.4995E-11,
    3.3898E-11, 15399.8516, 3.8429E-11, 3.3032E-11, 12645.5291, 3.1863E-11,
    3.2165E-11, 9902.2295, 2.5298E-11, 5.6505E-11, 1.2764E-10, 30845.2164,
    5.4375E-11, 1.2247E-10, 30052.3707, 5.2245E-11, 1.1731E-10, 29259.525,
    5.0115E-11, 1.1214E-10, 28466.6793, 4.7984E-11, 1.0698E-10, 27673.8337,
    4.5854E-11, 1.0181E-10, 26880.988, 4.3724E-11, 9.6647E-11, 26088.1423,
    4.1594E-11, 9.1482E-11, 25295.2966, 3.9463E-11, 8.6316E-11, 24502.4509,
    3.7333E-11, 8.1151E-11, 23709.6053, 3.5203E-11, 7.5986E-11, 22916.7596,
    3.3073E-11, 7.0821E-11, 22123.9139, 3.0942E-11, 6.5656E-11, 21331.0682,
    2.8812E-11, 6.049E-11, 20538.2225, 2.6682E-11, 5.5325E-11, 19745.3769,
    2.4552E-11, 5.016E-11, 18952.5312, 2.2421E-11, 4.4995E-11, 38159.6855,
    2.0291E-11, 3.9829E-11, 15399.8516, 1.8161E-11, 3.4664E-11, 12645.5291,
    1.6031E-11, 2.9499E-11, 9902.2295, 25890.1642, 8.1308E-11, 3.8045E-11,
    25234.5197, 7.8345E-11, 3.6936E-11, 24578.8751, 7.5382E-11, 3.5826E-11,
    23923.2306, 7.2419E-11, 3.4716E-11, 23267.5861, 6.9456E-11, 3.3607E-11,
    22611.9415, 6.6492E-11, 3.2497E-11, 21956.297, 6.3529E-11, 3.1387E-11,
    21300.6524, 6.0566E-11, 3.0278E-11, 20645.0079, 5.7603E-11, 2.9168E-11,
    19989.3634, 5.464E-11, 2.8058E-11, 19333.7188, 5.1677E-11, 2.6949E-11,
    18678.0743, 4.8714E-11, 2.5839E-11, 18022.4297, 4.5751E-11, 2.473E-11,
    17366.7852, 4.2788E-11, 2.362E-11, 16711.1407, 3.9825E-11, 2.251E-11,
    16055.4961, 3.6861E-11, 2.1401E-11, 15399.8516, 3.3898E-11, 2.0291E-11,
    34744.207, 3.0935E-11, 1.9181E-11, 12121.5743, 2.7972E-11, 1.8072E-11,
    9504.4529, 2.5009E-11, 1.6962E-11, 4.3544E-11, 25890.1642, 1.3242E-10,
    4.2802E-11, 25234.5197, 1.2664E-10, 4.206E-11, 24578.8751, 1.2085E-10,
    4.1319E-11, 23923.2306, 1.1506E-10, 4.0577E-11, 23267.5861, 1.0927E-10,
    3.9835E-11, 22611.9415, 1.0349E-10, 3.9094E-11, 21956.297, 9.77E-11,
    3.8352E-11, 21300.6524, 9.1913E-11, 3.761E-11, 20645.0079, 8.6126E-11,
    3.6869E-11, 19989.3634, 8.0339E-11, 3.6127E-11, 19333.7188, 7.4552E-11,
    3.5385E-11, 18678.0743, 6.8765E-11, 3.4644E-11, 18022.4297, 6.2978E-11,
    3.3902E-11, 17366.7852, 5.7191E-11, 3.316E-11, 16711.1407, 5.1404E-11,
    3.2419E-11, 16055.4961, 4.5617E-11, 3.1677E-11, 15399.8516, 3.9829E-11,
    3.0935E-11, 34744.207, 3.4042E-11, 3.0194E-11, 12121.5743, 2.8255E-11,
    2.9452E-11, 9504.4529, 2.2468E-11, 5.1779E-11, 1.0861E-10, 25890.1642,
    4.9862E-11, 1.0423E-10, 25234.5197, 4.7944E-11, 9.9841E-11, 24578.8751,
    4.6027E-11, 9.5455E-11, 23923.2306, 4.4109E-11, 9.1068E-11, 23267.5861,
    4.2192E-11, 8.6682E-11, 22611.9415, 4.0274E-11, 8.2295E-11, 21956.297,
    3.8357E-11, 7.7908E-11, 21300.6524, 3.6439E-11, 7.3522E-11, 20645.0079,
    3.4522E-11, 6.9135E-11, 19989.3634, 3.2604E-11, 6.4749E-11, 19333.7188,
    3.0687E-11, 6.0362E-11, 18678.0743, 2.8769E-11, 5.5975E-11, 18022.4297,
    2.6852E-11, 5.1589E-11, 17366.7852, 2.4934E-11, 4.7202E-11, 16711.1407,
    2.3017E-11, 4.2816E-11, 16055.4961, 2.1099E-11, 3.8429E-11, 15399.8516,
    1.9181E-11, 3.4042E-11, 34744.207, 1.7264E-11, 2.9656E-11, 12121.5743,
    1.5346E-11, 2.5269E-11, 9504.4529, 21028.8063, 7.844E-11, 3.2511E-11,
    20504.8515, 7.5602E-11, 3.1614E-11, 19980.8967, 7.2764E-11, 3.0717E-11,
    19456.9418, 6.9926E-11, 2.982E-11, 18932.987, 6.7088E-11, 2.8924E-11,
    18409.0322, 6.425E-11, 2.8027E-11, 17885.0774, 6.1412E-11, 2.713E-11,
    17361.1225, 5.8574E-11, 2.6233E-11, 16837.1677, 5.5736E-11, 2.5336E-11,
    16313.2129, 5.2898E-11, 2.4439E-11, 15789.2581, 5.006E-11, 2.3542E-11,
    15265.3032, 4.7222E-11, 2.2645E-11, 14741.3484, 4.4384E-11, 2.1748E-11,
    14217.3936, 4.1546E-11, 2.0852E-11, 13693.4388, 3.8708E-11, 1.9955E-11,
    13169.4839, 3.587E-11, 1.9058E-11, 12645.5291, 3.3032E-11, 1.8161E-11,
    12121.5743, 3.0194E-11, 1.7264E-11, 31597.6194, 2.7356E-11, 1.6367E-11,
    9106.6764, 2.4518E-11, 1.547E-11, 3.8454E-11, 21028.8063, 1.148E-10,
    3.7838E-11, 20504.8515, 1.0979E-10, 3.7221E-11, 19980.8967, 1.0478E-10,
    3.6604E-11, 19456.9418, 9.9774E-11, 3.5988E-11, 18932.987, 9.4766E-11,
    3.5371E-11, 18409.0322, 8.9758E-11, 3.4755E-11, 17885.0774, 8.4749E-11,
    3.4138E-11, 17361.1225, 7.9741E-11, 3.3521E-11, 16837.1677, 7.4732E-11,
    3.2905E-11, 16313.2129, 6.9724E-11, 3.2288E-11, 15789.2581, 6.4715E-11,
    3.1672E-11, 15265.3032, 5.9707E-11, 3.1055E-11, 14741.3484, 5.4698E-11,
    3.0438E-11, 14217.3936, 4.969E-11, 2.9822E-11, 13693.4388, 4.4681E-11,
    2.9205E-11, 13169.4839, 3.9673E-11, 2.8589E-11, 12645.5291, 3.4664E-11,
    2.7972E-11, 12121.5743, 2.9656E-11, 2.7356E-11, 31597.6194, 2.4647E-11,
    2.6739E-11, 9106.6764, 1.9639E-11, 4.7053E-11, 8.9591E-11, 21028.8063,
    4.5348E-11, 8.5983E-11, 20504.8515, 4.3644E-11, 8.2375E-11, 19980.8967,
    4.1939E-11, 7.8767E-11, 19456.9418, 4.0234E-11, 7.5159E-11, 18932.987,
    3.8529E-11, 7.1551E-11, 18409.0322, 3.6825E-11, 6.7943E-11, 17885.0774,
    3.512E-11, 6.4335E-11, 17361.1225, 3.3415E-11, 6.0727E-11, 16837.1677,
    3.171E-11, 5.7119E-11, 16313.2129, 3.0005E-11, 5.3511E-11, 15789.2581,
    2.8301E-11, 4.9903E-11, 15265.3032, 2.6596E-11, 4.6295E-11, 14741.3484,
    2.4891E-11, 4.2687E-11, 14217.3936, 2.3186E-11, 3.9079E-11, 13693.4388,
    2.1481E-11, 3.5471E-11, 13169.4839, 1.9777E-11, 3.1863E-11, 12645.5291,
    1.8072E-11, 2.8255E-11, 12121.5743, 1.6367E-11, 2.4647E-11, 31597.6194,
    1.4662E-11, 2.1039E-11, 9106.6764, 16266.654, 7.5571E-11, 2.6977E-11,
    15868.8775, 7.2858E-11, 2.6293E-11, 15471.101, 7.0146E-11, 2.5609E-11,
    15073.3244, 6.7433E-11, 2.4925E-11, 14675.5479, 6.472E-11, 2.4241E-11,
    14277.7714, 6.2007E-11, 2.3556E-11, 13879.9948, 5.9294E-11, 2.2872E-11,
    13482.2183, 5.6581E-11, 2.2188E-11, 13084.4418, 5.3868E-11, 2.1504E-11,
    12686.6652, 5.1155E-11, 2.082E-11, 12288.8887, 4.8442E-11, 2.0136E-11,
    11891.1121, 4.5729E-11, 1.9451E-11, 11493.3356, 4.3016E-11, 1.8767E-11,
    11095.5591, 4.0304E-11, 1.8083E-11, 10697.7825, 3.7591E-11, 1.7399E-11,
    10300.006, 3.4878E-11, 1.6715E-11, 9902.2295, 3.2165E-11, 1.6031E-11,
    9504.4529, 2.9452E-11, 1.5346E-11, 9106.6764, 2.6739E-11, 1.4662E-11,
    28708.8999, 2.4026E-11, 1.3978E-11, 3.3364E-11, 16266.654, 9.7177E-11,
    3.2873E-11, 15868.8775, 9.2947E-11, 3.2381E-11, 15471.101, 8.8717E-11,
    3.189E-11, 15073.3244, 8.4487E-11, 3.1398E-11, 14675.5479, 8.0257E-11,
    3.0907E-11, 14277.7714, 7.6028E-11, 3.0415E-11, 13879.9948, 7.1798E-11,
    2.9924E-11, 13482.2183, 6.7568E-11, 2.9432E-11, 13084.4418, 6.3338E-11,
    2.8941E-11, 12686.6652, 5.9108E-11, 2.8449E-11, 12288.8887, 5.4878E-11,
    2.7958E-11, 11891.1121, 5.0648E-11, 2.7467E-11, 11493.3356, 4.6418E-11,
    2.6975E-11, 11095.5591, 4.2189E-11, 2.6484E-11, 10697.7825, 3.7959E-11,
    2.5992E-11, 10300.006, 3.3729E-11, 2.5501E-11, 9902.2295, 2.9499E-11,
    2.5009E-11, 9504.4529, 2.5269E-11, 2.4518E-11, 9106.6764, 2.1039E-11,
    2.4026E-11, 28708.8999, 1.6809E-11, 4.2327E-11, 7.0567E-11, 16266.654,
    4.0835E-11, 6.7738E-11, 15868.8775, 3.9343E-11, 6.4909E-11, 15471.101,
    3.7851E-11, 6.2079E-11, 15073.3244, 3.6359E-11, 5.925E-11, 14675.5479,
    3.4867E-11, 5.6421E-11, 14277.7714, 3.3375E-11, 5.3591E-11, 13879.9948,
    3.1883E-11, 5.0762E-11, 13482.2183, 3.0391E-11, 4.7933E-11, 13084.4418,
    2.8899E-11, 4.5103E-11, 12686.6652, 2.7407E-11, 4.2274E-11, 12288.8887,
    2.5915E-11, 3.9444E-11, 11891.1121, 2.4423E-11, 3.6615E-11, 11493.3356,
    2.293E-11, 3.3786E-11, 11095.5591, 2.1438E-11, 3.0956E-11, 10697.7825,
    1.9946E-11, 2.8127E-11, 10300.006, 1.8454E-11, 2.5298E-11, 9902.2295,
    1.6962E-11, 2.2468E-11, 9504.4529, 1.547E-11, 1.9639E-11, 9106.6764,
    1.3978E-11, 1.6809E-11, 28708.8999 };

  static const real_T b_a[3600] = { -138371.1425, -1.3007E-10, -1.3212E-10,
    -112666.5366, -1.2498E-10, -1.274E-10, -106967.442, -1.1989E-10, -1.2267E-10,
    -101279.3704, -1.148E-10, -1.1794E-10, -95607.833, -1.0971E-10, -1.1322E-10,
    -89958.3413, -1.0462E-10, -1.0849E-10, -84336.4067, -9.953E-11, -1.0377E-10,
    -78747.5407, -9.4441E-11, -9.904E-11, -73197.2546, -8.9351E-11, -9.4314E-11,
    -67691.06, -8.4261E-11, -8.9588E-11, -62234.4682, -7.9172E-11, -8.4862E-11,
    -56832.9906, -7.4082E-11, -8.0136E-11, -51492.1387, -6.8992E-11, -7.541E-11,
    -46217.4239, -6.3903E-11, -7.0684E-11, -41014.3576, -5.8813E-11, -6.5958E-11,
    -35888.4513, -5.3723E-11, -6.1232E-11, -30845.2164, -4.8633E-11, -5.6505E-11,
    -25890.1642, -4.3544E-11, -5.1779E-11, -21028.8063, -3.8454E-11, -4.7053E-11,
    -16266.654, -3.3364E-11, -4.2327E-11, -1.3007E-10, -138371.1425, -4.3201E-10,
    -1.272E-10, -112666.5366, -4.1299E-10, -1.2433E-10, -106967.442, -3.9397E-10,
    -1.2146E-10, -101279.3704, -3.7494E-10, -1.186E-10, -95607.833, -3.5592E-10,
    -1.1573E-10, -89958.3413, -3.369E-10, -1.1286E-10, -84336.4067, -3.1787E-10,
    -1.0999E-10, -78747.5407, -2.9885E-10, -1.0712E-10, -73197.2546, -2.7983E-10,
    -1.0425E-10, -67691.06, -2.608E-10, -1.0139E-10, -62234.4682, -2.4178E-10,
    -9.8518E-11, -56832.9906, -2.2276E-10, -9.5649E-11, -51492.1387, -2.0373E-10,
    -9.2781E-11, -46217.4239, -1.8471E-10, -8.9913E-11, -41014.3576, -1.6569E-10,
    -8.7044E-11, -35888.4513, -1.4666E-10, -8.4176E-11, -30845.2164, -1.2764E-10,
    -8.1308E-11, -25890.1642, -1.0861E-10, -7.844E-11, -21028.8063, -8.9591E-11,
    -7.5571E-11, -16266.654, -7.0567E-11, -1.3212E-10, -4.3201E-10, -138371.1425,
    -1.2659E-10, -4.1439E-10, -112666.5366, -1.2105E-10, -3.9677E-10,
    -106967.442, -1.1552E-10, -3.7915E-10, -101279.3704, -1.0999E-10,
    -3.6152E-10, -95607.833, -1.0445E-10, -3.439E-10, -89958.3413, -9.8919E-11,
    -3.2628E-10, -84336.4067, -9.3385E-11, -3.0865E-10, -78747.5407, -8.7851E-11,
    -2.9103E-10, -73197.2546, -8.2317E-11, -2.7341E-10, -67691.06, -7.6783E-11,
    -2.5578E-10, -62234.4682, -7.1249E-11, -2.3816E-10, -56832.9906, -6.5715E-11,
    -2.2054E-10, -51492.1387, -6.0181E-11, -2.0292E-10, -46217.4239, -5.4647E-11,
    -1.8529E-10, -41014.3576, -4.9113E-11, -1.6767E-10, -35888.4513, -4.3579E-11,
    -1.5005E-10, -30845.2164, -3.8045E-11, -1.3242E-10, -25890.1642, -3.2511E-11,
    -1.148E-10, -21028.8063, -2.6977E-11, -9.7177E-11, -16266.654, -112666.5366,
    -1.272E-10, -1.2659E-10, -129154.3028, -1.2224E-10, -1.2208E-10,
    -103675.0808, -1.1727E-10, -1.1756E-10, -98201.3702, -1.1231E-10,
    -1.1305E-10, -92738.6825, -1.0734E-10, -1.0854E-10, -87292.529, -1.0238E-10,
    -1.0402E-10, -81868.4212, -9.7413E-11, -9.9509E-11, -76471.8706, -9.2448E-11,
    -9.4995E-11, -71108.3885, -8.7483E-11, -9.0482E-11, -65783.4864, -8.2519E-11,
    -8.5969E-11, -60502.6757, -7.7554E-11, -8.1455E-11, -55271.4678, -7.259E-11,
    -7.6942E-11, -50095.3742, -6.7625E-11, -7.2429E-11, -44979.9063, -6.2661E-11,
    -6.7915E-11, -39930.5754, -5.7696E-11, -6.3402E-11, -34952.8931, -5.2731E-11,
    -5.8889E-11, -30052.3707, -4.7767E-11, -5.4375E-11, -25234.5197, -4.2802E-11,
    -4.9862E-11, -20504.8515, -3.7838E-11, -4.5348E-11, -15868.8775, -3.2873E-11,
    -4.0835E-11, -1.2498E-10, -112666.5366, -4.1439E-10, -1.2224E-10,
    -129154.3028, -3.9615E-10, -1.1949E-10, -103675.0808, -3.779E-10,
    -1.1675E-10, -98201.3702, -3.5966E-10, -1.1401E-10, -92738.6825, -3.4141E-10,
    -1.1126E-10, -87292.529, -3.2317E-10, -1.0852E-10, -81868.4212, -3.0492E-10,
    -1.0578E-10, -76471.8706, -2.8668E-10, -1.0303E-10, -71108.3885, -2.6843E-10,
    -1.0029E-10, -65783.4864, -2.5019E-10, -9.7547E-11, -60502.6757, -2.3194E-10,
    -9.4804E-11, -55271.4678, -2.137E-10, -9.2061E-11, -50095.3742, -1.9545E-10,
    -8.9318E-11, -44979.9063, -1.7721E-10, -8.6574E-11, -39930.5754, -1.5896E-10,
    -8.3831E-11, -34952.8931, -1.4072E-10, -8.1088E-11, -30052.3707, -1.2247E-10,
    -7.8345E-11, -25234.5197, -1.0423E-10, -7.5602E-11, -20504.8515, -8.5983E-11,
    -7.2858E-11, -15868.8775, -6.7738E-11, -1.274E-10, -4.1299E-10, -112666.5366,
    -1.2208E-10, -3.9615E-10, -129154.3028, -1.1675E-10, -3.793E-10,
    -103675.0808, -1.1143E-10, -3.6246E-10, -98201.3702, -1.0611E-10,
    -3.4561E-10, -92738.6825, -1.0079E-10, -3.2877E-10, -87292.529, -9.5469E-11,
    -3.1192E-10, -81868.4212, -9.0148E-11, -2.9508E-10, -76471.8706, -8.4827E-11,
    -2.7824E-10, -71108.3885, -7.9506E-11, -2.6139E-10, -65783.4864, -7.4184E-11,
    -2.4455E-10, -60502.6757, -6.8863E-11, -2.277E-10, -55271.4678, -6.3542E-11,
    -2.1086E-10, -50095.3742, -5.8221E-11, -1.9401E-10, -44979.9063, -5.2899E-11,
    -1.7717E-10, -39930.5754, -4.7578E-11, -1.6032E-10, -34952.8931, -4.2257E-11,
    -1.4348E-10, -30052.3707, -3.6936E-11, -1.2664E-10, -25234.5197, -3.1614E-11,
    -1.0979E-10, -20504.8515, -2.6293E-11, -9.2947E-11, -15868.8775, -106967.442,
    -1.2433E-10, -1.2105E-10, -103675.0808, -1.1949E-10, -1.1675E-10,
    -120382.7195, -1.1465E-10, -1.1245E-10, -95123.37, -1.0981E-10, -1.0815E-10,
    -89869.532, -1.0497E-10, -1.0385E-10, -84626.7167, -1.0013E-10, -9.9552E-11,
    -79400.4358, -9.5295E-11, -9.5251E-11, -74196.2006, -9.0455E-11, -9.0951E-11,
    -69019.5224, -8.5616E-11, -8.665E-11, -63875.9129, -8.0776E-11, -8.2349E-11,
    -58770.8833, -7.5937E-11, -7.8049E-11, -53709.9451, -7.1097E-11, -7.3748E-11,
    -48698.6097, -6.6258E-11, -6.9447E-11, -43742.3886, -6.1418E-11, -6.5147E-11,
    -38846.7932, -5.6579E-11, -6.0846E-11, -34017.3348, -5.1739E-11, -5.6546E-11,
    -29259.525, -4.69E-11, -5.2245E-11, -24578.8751, -4.206E-11, -4.7944E-11,
    -19980.8967, -3.7221E-11, -4.3644E-11, -15471.101, -3.2381E-11, -3.9343E-11,
    -1.1989E-10, -106967.442, -3.9677E-10, -1.1727E-10, -103675.0808, -3.793E-10,
    -1.1465E-10, -120382.7195, -3.6184E-10, -1.1203E-10, -95123.37, -3.4437E-10,
    -1.0942E-10, -89869.532, -3.269E-10, -1.068E-10, -84626.7167, -3.0944E-10,
    -1.0418E-10, -79400.4358, -2.9197E-10, -1.0156E-10, -74196.2006, -2.745E-10,
    -9.8944E-11, -69019.5224, -2.5704E-10, -9.6326E-11, -63875.9129, -2.3957E-10,
    -9.3708E-11, -58770.8833, -2.2211E-10, -9.109E-11, -53709.9451, -2.0464E-10,
    -8.8472E-11, -48698.6097, -1.8717E-10, -8.5854E-11, -43742.3886, -1.6971E-10,
    -8.3236E-11, -38846.7932, -1.5224E-10, -8.0618E-11, -34017.3348, -1.3477E-10,
    -7.8E-11, -29259.525, -1.1731E-10, -7.5382E-11, -24578.8751, -9.9841E-11,
    -7.2764E-11, -19980.8967, -8.2375E-11, -7.0146E-11, -15471.101, -6.4909E-11,
    -1.2267E-10, -3.9397E-10, -106967.442, -1.1756E-10, -3.779E-10, -103675.0808,
    -1.1245E-10, -3.6184E-10, -120382.7195, -1.0735E-10, -3.4577E-10, -95123.37,
    -1.0224E-10, -3.297E-10, -89869.532, -9.7128E-11, -3.1364E-10, -84626.7167,
    -9.202E-11, -2.9757E-10, -79400.4358, -8.6911E-11, -2.8151E-10, -74196.2006,
    -8.1803E-11, -2.6544E-10, -69019.5224, -7.6694E-11, -2.4938E-10, -63875.9129,
    -7.1586E-11, -2.3331E-10, -58770.8833, -6.6477E-11, -2.1724E-10, -53709.9451,
    -6.1368E-11, -2.0118E-10, -48698.6097, -5.626E-11, -1.8511E-10, -43742.3886,
    -5.1151E-11, -1.6905E-10, -38846.7932, -4.6043E-11, -1.5298E-10, -34017.3348,
    -4.0934E-11, -1.3691E-10, -29259.525, -3.5826E-11, -1.2085E-10, -24578.8751,
    -3.0717E-11, -1.0478E-10, -19980.8967, -2.5609E-11, -8.8717E-11, -15471.101,
    -101279.3704, -1.2146E-10, -1.1552E-10, -98201.3702, -1.1675E-10,
    -1.1143E-10, -95123.37, -1.1203E-10, -1.0735E-10, -112045.3699, -1.0732E-10,
    -1.0326E-10, -87000.3814, -1.0261E-10, -9.9169E-11, -81960.9045, -9.7892E-11,
    -9.5081E-11, -76932.4503, -9.3177E-11, -9.0994E-11, -71920.5305, -8.8463E-11,
    -8.6906E-11, -66930.6563, -8.3748E-11, -8.2818E-11, -61968.3393, -7.9034E-11,
    -7.873E-11, -57039.0908, -7.432E-11, -7.4642E-11, -52148.4223, -6.9605E-11,
    -7.0554E-11, -47301.8452, -6.4891E-11, -6.6466E-11, -42504.871, -6.0176E-11,
    -6.2378E-11, -37763.0109, -5.5462E-11, -5.829E-11, -33081.7766, -5.0748E-11,
    -5.4203E-11, -28466.6793, -4.6033E-11, -5.0115E-11, -23923.2306, -4.1319E-11,
    -4.6027E-11, -19456.9418, -3.6604E-11, -4.1939E-11, -15073.3244, -3.189E-11,
    -3.7851E-11, -1.148E-10, -101279.3704, -3.7915E-10, -1.1231E-10, -98201.3702,
    -3.6246E-10, -1.0981E-10, -95123.37, -3.4577E-10, -1.0732E-10, -112045.3699,
    -3.2908E-10, -1.0483E-10, -87000.3814, -3.1239E-10, -1.0233E-10, -81960.9045,
    -2.9571E-10, -9.9842E-11, -76932.4503, -2.7902E-10, -9.7349E-11, -71920.5305,
    -2.6233E-10, -9.4856E-11, -66930.6563, -2.4564E-10, -9.2363E-11, -61968.3393,
    -2.2896E-10, -8.987E-11, -57039.0908, -2.1227E-10, -8.7377E-11, -52148.4223,
    -1.9558E-10, -8.4884E-11, -47301.8452, -1.7889E-10, -8.2391E-11, -42504.871,
    -1.6221E-10, -7.9898E-11, -37763.0109, -1.4552E-10, -7.7405E-11, -33081.7766,
    -1.2883E-10, -7.4912E-11, -28466.6793, -1.1214E-10, -7.2419E-11, -23923.2306,
    -9.5455E-11, -6.9926E-11, -19456.9418, -7.8767E-11, -6.7433E-11, -15073.3244,
    -6.2079E-11, -1.1794E-10, -3.7494E-10, -101279.3704, -1.1305E-10,
    -3.5966E-10, -98201.3702, -1.0815E-10, -3.4437E-10, -95123.37, -1.0326E-10,
    -3.2908E-10, -112045.3699, -9.8361E-11, -3.138E-10, -87000.3814, -9.3466E-11,
    -2.9851E-10, -81960.9045, -8.857E-11, -2.8322E-10, -76932.4503, -8.3674E-11,
    -2.6793E-10, -71920.5305, -7.8778E-11, -2.5265E-10, -66930.6563, -7.3883E-11,
    -2.3736E-10, -61968.3393, -6.8987E-11, -2.2207E-10, -57039.0908, -6.4091E-11,
    -2.0678E-10, -52148.4223, -5.9195E-11, -1.915E-10, -47301.8452, -5.4299E-11,
    -1.7621E-10, -42504.871, -4.9404E-11, -1.6092E-10, -37763.0109, -4.4508E-11,
    -1.4564E-10, -33081.7766, -3.9612E-11, -1.3035E-10, -28466.6793, -3.4716E-11,
    -1.1506E-10, -23923.2306, -2.982E-11, -9.9774E-11, -19456.9418, -2.4925E-11,
    -8.4487E-11, -15073.3244, -95607.833, -1.186E-10, -1.0999E-10, -92738.6825,
    -1.1401E-10, -1.0611E-10, -89869.532, -1.0942E-10, -1.0224E-10, -87000.3814,
    -1.0483E-10, -9.8361E-11, -104131.2309, -1.0024E-10, -9.4486E-11,
    -79295.0922, -9.5649E-11, -9.0611E-11, -74464.4649, -9.106E-11, -8.6736E-11,
    -69644.8604, -8.647E-11, -8.2861E-11, -64841.7902, -8.1881E-11, -7.8986E-11,
    -60060.7658, -7.7292E-11, -7.5111E-11, -55307.2984, -7.2702E-11, -7.1235E-11,
    -50586.8996, -6.8113E-11, -6.736E-11, -45905.0807, -6.3524E-11, -6.3485E-11,
    -41267.3533, -5.8934E-11, -5.961E-11, -36679.2287, -5.4345E-11, -5.5735E-11,
    -32146.2183, -4.9756E-11, -5.186E-11, -27673.8337, -4.5166E-11, -4.7984E-11,
    -23267.5861, -4.0577E-11, -4.4109E-11, -18932.987, -3.5988E-11, -4.0234E-11,
    -14675.5479, -3.1398E-11, -3.6359E-11, -1.0971E-10, -95607.833, -3.6152E-10,
    -1.0734E-10, -92738.6825, -3.4561E-10, -1.0497E-10, -89869.532, -3.297E-10,
    -1.0261E-10, -87000.3814, -3.138E-10, -1.0024E-10, -104131.2309, -2.9789E-10,
    -9.787E-11, -79295.0922, -2.8198E-10, -9.5502E-11, -74464.4649, -2.6607E-10,
    -9.3135E-11, -69644.8604, -2.5016E-10, -9.0767E-11, -64841.7902, -2.3425E-10,
    -8.8399E-11, -60060.7658, -2.1834E-10, -8.6031E-11, -55307.2984, -2.0243E-10,
    -8.3663E-11, -50586.8996, -1.8652E-10, -8.1295E-11, -45905.0807, -1.7061E-10,
    -7.8927E-11, -41267.3533, -1.547E-10, -7.6559E-11, -36679.2287, -1.388E-10,
    -7.4191E-11, -32146.2183, -1.2289E-10, -7.1823E-11, -27673.8337, -1.0698E-10,
    -6.9456E-11, -23267.5861, -9.1068E-11, -6.7088E-11, -18932.987, -7.5159E-11,
    -6.472E-11, -14675.5479, -5.925E-11, -1.1322E-10, -3.5592E-10, -95607.833,
    -1.0854E-10, -3.4141E-10, -92738.6825, -1.0385E-10, -3.269E-10, -89869.532,
    -9.9169E-11, -3.1239E-10, -87000.3814, -9.4486E-11, -2.9789E-10,
    -104131.2309, -8.9803E-11, -2.8338E-10, -79295.0922, -8.512E-11, -2.6887E-10,
    -74464.4649, -8.0437E-11, -2.5436E-10, -69644.8604, -7.5754E-11, -2.3985E-10,
    -64841.7902, -7.1071E-11, -2.2534E-10, -60060.7658, -6.6388E-11, -2.1083E-10,
    -55307.2984, -6.1705E-11, -1.9633E-10, -50586.8996, -5.7022E-11, -1.8182E-10,
    -45905.0807, -5.2339E-11, -1.6731E-10, -41267.3533, -4.7656E-11, -1.528E-10,
    -36679.2287, -4.2973E-11, -1.3829E-10, -32146.2183, -3.829E-11, -1.2378E-10,
    -27673.8337, -3.3607E-11, -1.0927E-10, -23267.5861, -2.8924E-11, -9.4766E-11,
    -18932.987, -2.4241E-11, -8.0257E-11, -14675.5479, -89958.3413, -1.1573E-10,
    -1.0445E-10, -87292.529, -1.1126E-10, -1.0079E-10, -84626.7167, -1.068E-10,
    -9.7128E-11, -81960.9045, -1.0233E-10, -9.3466E-11, -79295.0922, -9.787E-11,
    -8.9803E-11, -96629.2799, -9.3406E-11, -8.6141E-11, -71996.4794, -8.8942E-11,
    -8.2478E-11, -67369.1904, -8.4478E-11, -7.8816E-11, -62752.9241, -8.0013E-11,
    -7.5154E-11, -58153.1922, -7.5549E-11, -7.1491E-11, -53575.5059, -7.1085E-11,
    -6.7829E-11, -49025.3768, -6.6621E-11, -6.4166E-11, -44508.3163, -6.2157E-11,
    -6.0504E-11, -40029.8357, -5.7692E-11, -5.6841E-11, -35595.4465, -5.3228E-11,
    -5.3179E-11, -31210.6601, -4.8764E-11, -4.9517E-11, -26880.988, -4.43E-11,
    -4.5854E-11, -22611.9415, -3.9835E-11, -4.2192E-11, -18409.0322, -3.5371E-11,
    -3.8529E-11, -14277.7714, -3.0907E-11, -3.4867E-11, -1.0462E-10, -89958.3413,
    -3.439E-10, -1.0238E-10, -87292.529, -3.2877E-10, -1.0013E-10, -84626.7167,
    -3.1364E-10, -9.7892E-11, -81960.9045, -2.9851E-10, -9.5649E-11, -79295.0922,
    -2.8338E-10, -9.3406E-11, -96629.2799, -2.6825E-10, -9.1163E-11, -71996.4794,
    -2.5312E-10, -8.892E-11, -67369.1904, -2.3799E-10, -8.6678E-11, -62752.9241,
    -2.2286E-10, -8.4435E-11, -58153.1922, -2.0773E-10, -8.2192E-11, -53575.5059,
    -1.9259E-10, -7.9949E-11, -49025.3768, -1.7746E-10, -7.7706E-11, -44508.3163,
    -1.6233E-10, -7.5464E-11, -40029.8357, -1.472E-10, -7.3221E-11, -35595.4465,
    -1.3207E-10, -7.0978E-11, -31210.6601, -1.1694E-10, -6.8735E-11, -26880.988,
    -1.0181E-10, -6.6492E-11, -22611.9415, -8.6682E-11, -6.425E-11, -18409.0322,
    -7.1551E-11, -6.2007E-11, -14277.7714, -5.6421E-11, -1.0849E-10, -3.369E-10,
    -89958.3413, -1.0402E-10, -3.2317E-10, -87292.529, -9.9552E-11, -3.0944E-10,
    -84626.7167, -9.5081E-11, -2.9571E-10, -81960.9045, -9.0611E-11, -2.8198E-10,
    -79295.0922, -8.6141E-11, -2.6825E-10, -96629.2799, -8.1671E-11, -2.5452E-10,
    -71996.4794, -7.72E-11, -2.4079E-10, -67369.1904, -7.273E-11, -2.2706E-10,
    -62752.9241, -6.826E-11, -2.1333E-10, -58153.1922, -6.3789E-11, -1.996E-10,
    -53575.5059, -5.9319E-11, -1.8587E-10, -49025.3768, -5.4849E-11, -1.7214E-10,
    -44508.3163, -5.0378E-11, -1.5841E-10, -40029.8357, -4.5908E-11, -1.4468E-10,
    -35595.4465, -4.1438E-11, -1.3095E-10, -31210.6601, -3.6967E-11, -1.1722E-10,
    -26880.988, -3.2497E-11, -1.0349E-10, -22611.9415, -2.8027E-11, -8.9758E-11,
    -18409.0322, -2.3556E-11, -7.6028E-11, -14277.7714, -84336.4067, -1.1286E-10,
    -9.8919E-11, -81868.4212, -1.0852E-10, -9.5469E-11, -79400.4358, -1.0418E-10,
    -9.202E-11, -76932.4503, -9.9842E-11, -8.857E-11, -74464.4649, -9.5502E-11,
    -8.512E-11, -71996.4794, -9.1163E-11, -8.1671E-11, -89528.494, -8.6824E-11,
    -7.8221E-11, -65093.5203, -8.2485E-11, -7.4771E-11, -60664.058, -7.8146E-11,
    -7.1321E-11, -56245.6186, -7.3807E-11, -6.7872E-11, -51843.7135, -6.9468E-11,
    -6.4422E-11, -47463.8541, -6.5128E-11, -6.0972E-11, -43111.5518, -6.0789E-11,
    -5.7523E-11, -38792.318, -5.645E-11, -5.4073E-11, -34511.6642, -5.2111E-11,
    -5.0623E-11, -30275.1018, -4.7772E-11, -4.7174E-11, -26088.1423, -4.3433E-11,
    -4.3724E-11, -21956.297, -3.9094E-11, -4.0274E-11, -17885.0774, -3.4755E-11,
    -3.6825E-11, -13879.9948, -3.0415E-11, -3.3375E-11, -9.953E-11, -84336.4067,
    -3.2628E-10, -9.7413E-11, -81868.4212, -3.1192E-10, -9.5295E-11, -79400.4358,
    -2.9757E-10, -9.3177E-11, -76932.4503, -2.8322E-10, -9.106E-11, -74464.4649,
    -2.6887E-10, -8.8942E-11, -71996.4794, -2.5452E-10, -8.6824E-11, -89528.494,
    -2.4017E-10, -8.4706E-11, -65093.5203, -2.2581E-10, -8.2589E-11, -60664.058,
    -2.1146E-10, -8.0471E-11, -56245.6186, -1.9711E-10, -7.8353E-11, -51843.7135,
    -1.8276E-10, -7.6236E-11, -47463.8541, -1.6841E-10, -7.4118E-11, -43111.5518,
    -1.5405E-10, -7.2E-11, -38792.318, -1.397E-10, -6.9882E-11, -34511.6642,
    -1.2535E-10, -6.7765E-11, -30275.1018, -1.11E-10, -6.5647E-11, -26088.1423,
    -9.6647E-11, -6.3529E-11, -21956.297, -8.2295E-11, -6.1412E-11, -17885.0774,
    -6.7943E-11, -5.9294E-11, -13879.9948, -5.3591E-11, -1.0377E-10, -3.1787E-10,
    -84336.4067, -9.9509E-11, -3.0492E-10, -81868.4212, -9.5251E-11, -2.9197E-10,
    -79400.4358, -9.0994E-11, -2.7902E-10, -76932.4503, -8.6736E-11, -2.6607E-10,
    -74464.4649, -8.2478E-11, -2.5312E-10, -71996.4794, -7.8221E-11, -2.4017E-10,
    -89528.494, -7.3963E-11, -2.2721E-10, -65093.5203, -6.9706E-11, -2.1426E-10,
    -60664.058, -6.5448E-11, -2.0131E-10, -56245.6186, -6.119E-11, -1.8836E-10,
    -51843.7135, -5.6933E-11, -1.7541E-10, -47463.8541, -5.2675E-11, -1.6246E-10,
    -43111.5518, -4.8418E-11, -1.4951E-10, -38792.318, -4.416E-11, -1.3655E-10,
    -34511.6642, -3.9903E-11, -1.236E-10, -30275.1018, -3.5645E-11, -1.1065E-10,
    -26088.1423, -3.1387E-11, -9.77E-11, -21956.297, -2.713E-11, -8.4749E-11,
    -17885.0774, -2.2872E-11, -7.1798E-11, -13879.9948, -78747.5407, -1.0999E-10,
    -9.3385E-11, -76471.8706, -1.0578E-10, -9.0148E-11, -74196.2006, -1.0156E-10,
    -8.6911E-11, -71920.5305, -9.7349E-11, -8.3674E-11, -69644.8604, -9.3135E-11,
    -8.0437E-11, -67369.1904, -8.892E-11, -7.72E-11, -65093.5203, -8.4706E-11,
    -7.3963E-11, -82817.8502, -8.0492E-11, -7.0726E-11, -58575.192, -7.6278E-11,
    -6.7489E-11, -54338.0451, -7.2064E-11, -6.4252E-11, -50111.9211, -6.785E-11,
    -6.1015E-11, -45902.3313, -6.3636E-11, -5.7778E-11, -41714.7873, -5.9422E-11,
    -5.4542E-11, -37554.8004, -5.5208E-11, -5.1305E-11, -33427.882, -5.0994E-11,
    -4.8068E-11, -29339.5436, -4.678E-11, -4.4831E-11, -25295.2966, -4.2566E-11,
    -4.1594E-11, -21300.6524, -3.8352E-11, -3.8357E-11, -17361.1225, -3.4138E-11,
    -3.512E-11, -13482.2183, -2.9924E-11, -3.1883E-11, -9.4441E-11, -78747.5407,
    -3.0865E-10, -9.2448E-11, -76471.8706, -2.9508E-10, -9.0455E-11, -74196.2006,
    -2.8151E-10, -8.8463E-11, -71920.5305, -2.6793E-10, -8.647E-11, -69644.8604,
    -2.5436E-10, -8.4478E-11, -67369.1904, -2.4079E-10, -8.2485E-11, -65093.5203,
    -2.2721E-10, -8.0492E-11, -82817.8502, -2.1364E-10, -7.85E-11, -58575.192,
    -2.0007E-10, -7.6507E-11, -54338.0451, -1.8649E-10, -7.4515E-11, -50111.9211,
    -1.7292E-10, -7.2522E-11, -45902.3313, -1.5935E-10, -7.0529E-11, -41714.7873,
    -1.4577E-10, -6.8537E-11, -37554.8004, -1.322E-10, -6.6544E-11, -33427.882,
    -1.1863E-10, -6.4551E-11, -29339.5436, -1.0505E-10, -6.2559E-11, -25295.2966,
    -9.1482E-11, -6.0566E-11, -21300.6524, -7.7908E-11, -5.8574E-11, -17361.1225,
    -6.4335E-11, -5.6581E-11, -13482.2183, -5.0762E-11, -9.904E-11, -2.9885E-10,
    -78747.5407, -9.4995E-11, -2.8668E-10, -76471.8706, -9.0951E-11, -2.745E-10,
    -74196.2006, -8.6906E-11, -2.6233E-10, -71920.5305, -8.2861E-11, -2.5016E-10,
    -69644.8604, -7.8816E-11, -2.3799E-10, -67369.1904, -7.4771E-11, -2.2581E-10,
    -65093.5203, -7.0726E-11, -2.1364E-10, -82817.8502, -6.6681E-11, -2.0147E-10,
    -58575.192, -6.2637E-11, -1.893E-10, -54338.0451, -5.8592E-11, -1.7712E-10,
    -50111.9211, -5.4547E-11, -1.6495E-10, -45902.3313, -5.0502E-11, -1.5278E-10,
    -41714.7873, -4.6457E-11, -1.406E-10, -37554.8004, -4.2412E-11, -1.2843E-10,
    -33427.882, -3.8367E-11, -1.1626E-10, -29339.5436, -3.4323E-11, -1.0409E-10,
    -25295.2966, -3.0278E-11, -9.1913E-11, -21300.6524, -2.6233E-11, -7.9741E-11,
    -17361.1225, -2.2188E-11, -6.7568E-11, -13482.2183, -73197.2546, -1.0712E-10,
    -8.7851E-11, -71108.3885, -1.0303E-10, -8.4827E-11, -69019.5224, -9.8944E-11,
    -8.1803E-11, -66930.6563, -9.4856E-11, -7.8778E-11, -64841.7902, -9.0767E-11,
    -7.5754E-11, -62752.9241, -8.6678E-11, -7.273E-11, -60664.058, -8.2589E-11,
    -6.9706E-11, -58575.192, -7.85E-11, -6.6681E-11, -76486.3259, -7.4411E-11,
    -6.3657E-11, -52430.4715, -7.0322E-11, -6.0633E-11, -48380.1286, -6.6233E-11,
    -5.7609E-11, -44340.8086, -6.2144E-11, -5.4585E-11, -40318.0228, -5.8055E-11,
    -5.156E-11, -36317.2827, -5.3966E-11, -4.8536E-11, -32344.0998, -4.9877E-11,
    -4.5512E-11, -28403.9854, -4.5788E-11, -4.2488E-11, -24502.4509, -4.1699E-11,
    -3.9463E-11, -20645.0079, -3.761E-11, -3.6439E-11, -16837.1677, -3.3521E-11,
    -3.3415E-11, -13084.4418, -2.9432E-11, -3.0391E-11, -8.9351E-11, -73197.2546,
    -2.9103E-10, -8.7483E-11, -71108.3885, -2.7824E-10, -8.5616E-11, -69019.5224,
    -2.6544E-10, -8.3748E-11, -66930.6563, -2.5265E-10, -8.1881E-11, -64841.7902,
    -2.3985E-10, -8.0013E-11, -62752.9241, -2.2706E-10, -7.8146E-11, -60664.058,
    -2.1426E-10, -7.6278E-11, -58575.192, -2.0147E-10, -7.4411E-11, -76486.3259,
    -1.8867E-10, -7.2543E-11, -52430.4715, -1.7588E-10, -7.0676E-11, -48380.1286,
    -1.6308E-10, -6.8808E-11, -44340.8086, -1.5029E-10, -6.6941E-11, -40318.0228,
    -1.3749E-10, -6.5073E-11, -36317.2827, -1.247E-10, -6.3206E-11, -32344.0998,
    -1.1191E-10, -6.1338E-11, -28403.9854, -9.9111E-11, -5.9471E-11, -24502.4509,
    -8.6316E-11, -5.7603E-11, -20645.0079, -7.3522E-11, -5.5736E-11, -16837.1677,
    -6.0727E-11, -5.3868E-11, -13084.4418, -4.7933E-11, -9.4314E-11, -2.7983E-10,
    -73197.2546, -9.0482E-11, -2.6843E-10, -71108.3885, -8.665E-11, -2.5704E-10,
    -69019.5224, -8.2818E-11, -2.4564E-10, -66930.6563, -7.8986E-11, -2.3425E-10,
    -64841.7902, -7.5154E-11, -2.2286E-10, -62752.9241, -7.1321E-11, -2.1146E-10,
    -60664.058, -6.7489E-11, -2.0007E-10, -58575.192, -6.3657E-11, -1.8867E-10,
    -76486.3259, -5.9825E-11, -1.7728E-10, -52430.4715, -5.5993E-11, -1.6589E-10,
    -48380.1286, -5.2161E-11, -1.5449E-10, -44340.8086, -4.8329E-11, -1.431E-10,
    -40318.0228, -4.4497E-11, -1.317E-10, -36317.2827, -4.0664E-11, -1.2031E-10,
    -32344.0998, -3.6832E-11, -1.0891E-10, -28403.9854, -3.3E-11, -9.752E-11,
    -24502.4509, -2.9168E-11, -8.6126E-11, -20645.0079, -2.5336E-11, -7.4732E-11,
    -16837.1677, -2.1504E-11, -6.3338E-11, -13084.4418, -67691.06, -1.0425E-10,
    -8.2317E-11, -65783.4864, -1.0029E-10, -7.9506E-11, -63875.9129, -9.6326E-11,
    -7.6694E-11, -61968.3393, -9.2363E-11, -7.3883E-11, -60060.7658, -8.8399E-11,
    -7.1071E-11, -58153.1922, -8.4435E-11, -6.826E-11, -56245.6186, -8.0471E-11,
    -6.5448E-11, -54338.0451, -7.6507E-11, -6.2637E-11, -52430.4715, -7.2543E-11,
    -5.9825E-11, -70522.898, -6.8579E-11, -5.7014E-11, -46648.3362, -6.4616E-11,
    -5.4202E-11, -42779.2858, -6.0652E-11, -5.1391E-11, -38921.2583, -5.6688E-11,
    -4.8579E-11, -35079.7651, -5.2724E-11, -4.5768E-11, -31260.3175, -4.876E-11,
    -4.2956E-11, -27468.4271, -4.4796E-11, -4.0145E-11, -23709.6053, -4.0833E-11,
    -3.7333E-11, -19989.3634, -3.6869E-11, -3.4522E-11, -16313.2129, -3.2905E-11,
    -3.171E-11, -12686.6652, -2.8941E-11, -2.8899E-11, -8.4261E-11, -67691.06,
    -2.7341E-10, -8.2519E-11, -65783.4864, -2.6139E-10, -8.0776E-11, -63875.9129,
    -2.4938E-10, -7.9034E-11, -61968.3393, -2.3736E-10, -7.7292E-11, -60060.7658,
    -2.2534E-10, -7.5549E-11, -58153.1922, -2.1333E-10, -7.3807E-11, -56245.6186,
    -2.0131E-10, -7.2064E-11, -54338.0451, -1.893E-10, -7.0322E-11, -52430.4715,
    -1.7728E-10, -6.8579E-11, -70522.898, -1.6526E-10, -6.6837E-11, -46648.3362,
    -1.5325E-10, -6.5095E-11, -42779.2858, -1.4123E-10, -6.3352E-11, -38921.2583,
    -1.2922E-10, -6.161E-11, -35079.7651, -1.172E-10, -5.9867E-11, -31260.3175,
    -1.0518E-10, -5.8125E-11, -27468.4271, -9.3167E-11, -5.6382E-11, -23709.6053,
    -8.1151E-11, -5.464E-11, -19989.3634, -6.9135E-11, -5.2898E-11, -16313.2129,
    -5.7119E-11, -5.1155E-11, -12686.6652, -4.5103E-11, -8.9588E-11, -2.608E-10,
    -67691.06, -8.5969E-11, -2.5019E-10, -65783.4864, -8.2349E-11, -2.3957E-10,
    -63875.9129, -7.873E-11, -2.2896E-10, -61968.3393, -7.5111E-11, -2.1834E-10,
    -60060.7658, -7.1491E-11, -2.0773E-10, -58153.1922, -6.7872E-11, -1.9711E-10,
    -56245.6186, -6.4252E-11, -1.8649E-10, -54338.0451, -6.0633E-11, -1.7588E-10,
    -52430.4715, -5.7014E-11, -1.6526E-10, -70522.898, -5.3394E-11, -1.5465E-10,
    -46648.3362, -4.9775E-11, -1.4403E-10, -42779.2858, -4.6155E-11, -1.3342E-10,
    -38921.2583, -4.2536E-11, -1.228E-10, -35079.7651, -3.8917E-11, -1.1219E-10,
    -31260.3175, -3.5297E-11, -1.0157E-10, -27468.4271, -3.1678E-11, -9.0955E-11,
    -23709.6053, -2.8058E-11, -8.0339E-11, -19989.3634, -2.4439E-11, -6.9724E-11,
    -16313.2129, -2.082E-11, -5.9108E-11, -12686.6652, -62234.4682, -1.0139E-10,
    -7.6783E-11, -60502.6757, -9.7547E-11, -7.4184E-11, -58770.8833, -9.3708E-11,
    -7.1586E-11, -57039.0908, -8.987E-11, -6.8987E-11, -55307.2984, -8.6031E-11,
    -6.6388E-11, -53575.5059, -8.2192E-11, -6.3789E-11, -51843.7135, -7.8353E-11,
    -6.119E-11, -50111.9211, -7.4515E-11, -5.8592E-11, -48380.1286, -7.0676E-11,
    -5.5993E-11, -46648.3362, -6.6837E-11, -5.3394E-11, -64916.5437, -6.2998E-11,
    -5.0795E-11, -41217.7631, -5.916E-11, -4.8197E-11, -37524.4938, -5.5321E-11,
    -4.5598E-11, -33842.2474, -5.1482E-11, -4.2999E-11, -30176.5353, -4.7643E-11,
    -4.04E-11, -26532.8689, -4.3805E-11, -3.7802E-11, -22916.7596, -3.9966E-11,
    -3.5203E-11, -19333.7188, -3.6127E-11, -3.2604E-11, -15789.2581, -3.2288E-11,
    -3.0005E-11, -12288.8887, -2.8449E-11, -2.7407E-11, -7.9172E-11, -62234.4682,
    -2.5578E-10, -7.7554E-11, -60502.6757, -2.4455E-10, -7.5937E-11, -58770.8833,
    -2.3331E-10, -7.432E-11, -57039.0908, -2.2207E-10, -7.2702E-11, -55307.2984,
    -2.1083E-10, -7.1085E-11, -53575.5059, -1.996E-10, -6.9468E-11, -51843.7135,
    -1.8836E-10, -6.785E-11, -50111.9211, -1.7712E-10, -6.6233E-11, -48380.1286,
    -1.6589E-10, -6.4616E-11, -46648.3362, -1.5465E-10, -6.2998E-11, -64916.5437,
    -1.4341E-10, -6.1381E-11, -41217.7631, -1.3217E-10, -5.9764E-11, -37524.4938,
    -1.2094E-10, -5.8146E-11, -33842.2474, -1.097E-10, -5.6529E-11, -30176.5353,
    -9.8461E-11, -5.4912E-11, -26532.8689, -8.7223E-11, -5.3294E-11, -22916.7596,
    -7.5986E-11, -5.1677E-11, -19333.7188, -6.4749E-11, -5.006E-11, -15789.2581,
    -5.3511E-11, -4.8442E-11, -12288.8887, -4.2274E-11, -8.4862E-11, -2.4178E-10,
    -62234.4682, -8.1455E-11, -2.3194E-10, -60502.6757, -7.8049E-11, -2.2211E-10,
    -58770.8833, -7.4642E-11, -2.1227E-10, -57039.0908, -7.1235E-11, -2.0243E-10,
    -55307.2984, -6.7829E-11, -1.9259E-10, -53575.5059, -6.4422E-11, -1.8276E-10,
    -51843.7135, -6.1015E-11, -1.7292E-10, -50111.9211, -5.7609E-11, -1.6308E-10,
    -48380.1286, -5.4202E-11, -1.5325E-10, -46648.3362, -5.0795E-11, -1.4341E-10,
    -64916.5437, -4.7389E-11, -1.3357E-10, -41217.7631, -4.3982E-11, -1.2374E-10,
    -37524.4938, -4.0575E-11, -1.139E-10, -33842.2474, -3.7169E-11, -1.0406E-10,
    -30176.5353, -3.3762E-11, -9.4226E-11, -26532.8689, -3.0356E-11, -8.4389E-11,
    -22916.7596, -2.6949E-11, -7.4552E-11, -19333.7188, -2.3542E-11, -6.4715E-11,
    -15789.2581, -2.0136E-11, -5.4878E-11, -12288.8887, -56832.9906, -9.8518E-11,
    -7.1249E-11, -55271.4678, -9.4804E-11, -6.8863E-11, -53709.9451, -9.109E-11,
    -6.6477E-11, -52148.4223, -8.7377E-11, -6.4091E-11, -50586.8996, -8.3663E-11,
    -6.1705E-11, -49025.3768, -7.9949E-11, -5.9319E-11, -47463.8541, -7.6236E-11,
    -5.6933E-11, -45902.3313, -7.2522E-11, -5.4547E-11, -44340.8086, -6.8808E-11,
    -5.2161E-11, -42779.2858, -6.5095E-11, -4.9775E-11, -41217.7631, -6.1381E-11,
    -4.7389E-11, -59656.2403, -5.7667E-11, -4.5003E-11, -36127.7293, -5.3954E-11,
    -4.2617E-11, -32604.7298, -5.024E-11, -4.0231E-11, -29092.7531, -4.6526E-11,
    -3.7845E-11, -25597.3106, -4.2813E-11, -3.5459E-11, -22123.9139, -3.9099E-11,
    -3.3073E-11, -18678.0743, -3.5385E-11, -3.0687E-11, -15265.3032, -3.1672E-11,
    -2.8301E-11, -11891.1121, -2.7958E-11, -2.5915E-11, -7.4082E-11, -56832.9906,
    -2.3816E-10, -7.259E-11, -55271.4678, -2.277E-10, -7.1097E-11, -53709.9451,
    -2.1724E-10, -6.9605E-11, -52148.4223, -2.0678E-10, -6.8113E-11, -50586.8996,
    -1.9633E-10, -6.6621E-11, -49025.3768, -1.8587E-10, -6.5128E-11, -47463.8541,
    -1.7541E-10, -6.3636E-11, -45902.3313, -1.6495E-10, -6.2144E-11, -44340.8086,
    -1.5449E-10, -6.0652E-11, -42779.2858, -1.4403E-10, -5.916E-11, -41217.7631,
    -1.3357E-10, -5.7667E-11, -59656.2403, -1.2311E-10, -5.6175E-11, -36127.7293,
    -1.1266E-10, -5.4683E-11, -32604.7298, -1.022E-10, -5.3191E-11, -29092.7531,
    -9.1738E-11, -5.1698E-11, -25597.3106, -8.128E-11, -5.0206E-11, -22123.9139,
    -7.0821E-11, -4.8714E-11, -18678.0743, -6.0362E-11, -4.7222E-11, -15265.3032,
    -4.9903E-11, -4.5729E-11, -11891.1121, -3.9444E-11, -8.0136E-11, -2.2276E-10,
    -56832.9906, -7.6942E-11, -2.137E-10, -55271.4678, -7.3748E-11, -2.0464E-10,
    -53709.9451, -7.0554E-11, -1.9558E-10, -52148.4223, -6.736E-11, -1.8652E-10,
    -50586.8996, -6.4166E-11, -1.7746E-10, -49025.3768, -6.0972E-11, -1.6841E-10,
    -47463.8541, -5.7778E-11, -1.5935E-10, -45902.3313, -5.4585E-11, -1.5029E-10,
    -44340.8086, -5.1391E-11, -1.4123E-10, -42779.2858, -4.8197E-11, -1.3217E-10,
    -41217.7631, -4.5003E-11, -1.2311E-10, -59656.2403, -4.1809E-11, -1.1406E-10,
    -36127.7293, -3.8615E-11, -1.05E-10, -32604.7298, -3.5421E-11, -9.594E-11,
    -29092.7531, -3.2227E-11, -8.6882E-11, -25597.3106, -2.9033E-11, -7.7823E-11,
    -22123.9139, -2.5839E-11, -6.8765E-11, -18678.0743, -2.2645E-11, -5.9707E-11,
    -15265.3032, -1.9451E-11, -5.0648E-11, -11891.1121, -51492.1387, -9.5649E-11,
    -6.5715E-11, -50095.3742, -9.2061E-11, -6.3542E-11, -48698.6097, -8.8472E-11,
    -6.1368E-11, -47301.8452, -8.4884E-11, -5.9195E-11, -45905.0807, -8.1295E-11,
    -5.7022E-11, -44508.3163, -7.7706E-11, -5.4849E-11, -43111.5518, -7.4118E-11,
    -5.2675E-11, -41714.7873, -7.0529E-11, -5.0502E-11, -40318.0228, -6.6941E-11,
    -4.8329E-11, -38921.2583, -6.3352E-11, -4.6155E-11, -37524.4938, -5.9764E-11,
    -4.3982E-11, -36127.7293, -5.6175E-11, -4.1809E-11, -54730.9648, -5.2586E-11,
    -3.9636E-11, -31367.2121, -4.8998E-11, -3.7462E-11, -28008.9708, -4.5409E-11,
    -3.5289E-11, -24661.7524, -4.1821E-11, -3.3116E-11, -21331.0682, -3.8232E-11,
    -3.0942E-11, -18022.4297, -3.4644E-11, -2.8769E-11, -14741.3484, -3.1055E-11,
    -2.6596E-11, -11493.3356, -2.7467E-11, -2.4423E-11, -6.8992E-11, -51492.1387,
    -2.2054E-10, -6.7625E-11, -50095.3742, -2.1086E-10, -6.6258E-11, -48698.6097,
    -2.0118E-10, -6.4891E-11, -47301.8452, -1.915E-10, -6.3524E-11, -45905.0807,
    -1.8182E-10, -6.2157E-11, -44508.3163, -1.7214E-10, -6.0789E-11, -43111.5518,
    -1.6246E-10, -5.9422E-11, -41714.7873, -1.5278E-10, -5.8055E-11, -40318.0228,
    -1.431E-10, -5.6688E-11, -38921.2583, -1.3342E-10, -5.5321E-11, -37524.4938,
    -1.2374E-10, -5.3954E-11, -36127.7293, -1.1406E-10, -5.2586E-11, -54730.9648,
    -1.0438E-10, -5.1219E-11, -31367.2121, -9.4696E-11, -4.9852E-11, -28008.9708,
    -8.5016E-11, -4.8485E-11, -24661.7524, -7.5336E-11, -4.7118E-11, -21331.0682,
    -6.5656E-11, -4.5751E-11, -18022.4297, -5.5975E-11, -4.4384E-11, -14741.3484,
    -4.6295E-11, -4.3016E-11, -11493.3356, -3.6615E-11, -7.541E-11, -2.0373E-10,
    -51492.1387, -7.2429E-11, -1.9545E-10, -50095.3742, -6.9447E-11, -1.8717E-10,
    -48698.6097, -6.6466E-11, -1.7889E-10, -47301.8452, -6.3485E-11, -1.7061E-10,
    -45905.0807, -6.0504E-11, -1.6233E-10, -44508.3163, -5.7523E-11, -1.5405E-10,
    -43111.5518, -5.4542E-11, -1.4577E-10, -41714.7873, -5.156E-11, -1.3749E-10,
    -40318.0228, -4.8579E-11, -1.2922E-10, -38921.2583, -4.5598E-11, -1.2094E-10,
    -37524.4938, -4.2617E-11, -1.1266E-10, -36127.7293, -3.9636E-11, -1.0438E-10,
    -54730.9648, -3.6654E-11, -9.6097E-11, -31367.2121, -3.3673E-11, -8.7817E-11,
    -28008.9708, -3.0692E-11, -7.9537E-11, -24661.7524, -2.7711E-11, -7.1257E-11,
    -21331.0682, -2.473E-11, -6.2978E-11, -18022.4297, -2.1748E-11, -5.4698E-11,
    -14741.3484, -1.8767E-11, -4.6418E-11, -11493.3356, -46217.4239, -9.2781E-11,
    -6.0181E-11, -44979.9063, -8.9318E-11, -5.8221E-11, -43742.3886, -8.5854E-11,
    -5.626E-11, -42504.871, -8.2391E-11, -5.4299E-11, -41267.3533, -7.8927E-11,
    -5.2339E-11, -40029.8357, -7.5464E-11, -5.0378E-11, -38792.318, -7.2E-11,
    -4.8418E-11, -37554.8004, -6.8537E-11, -4.6457E-11, -36317.2827, -6.5073E-11,
    -4.4497E-11, -35079.7651, -6.161E-11, -4.2536E-11, -33842.2474, -5.8146E-11,
    -4.0575E-11, -32604.7298, -5.4683E-11, -3.8615E-11, -31367.2121, -5.1219E-11,
    -3.6654E-11, -50129.6945, -4.7756E-11, -3.4694E-11, -26925.1886, -4.4292E-11,
    -3.2733E-11, -23726.1941, -4.0829E-11, -3.0773E-11, -20538.2225, -3.7365E-11,
    -2.8812E-11, -17366.7852, -3.3902E-11, -2.6852E-11, -14217.3936, -3.0438E-11,
    -2.4891E-11, -11095.5591, -2.6975E-11, -2.293E-11, -6.3903E-11, -46217.4239,
    -2.0292E-10, -6.2661E-11, -44979.9063, -1.9401E-10, -6.1418E-11, -43742.3886,
    -1.8511E-10, -6.0176E-11, -42504.871, -1.7621E-10, -5.8934E-11, -41267.3533,
    -1.6731E-10, -5.7692E-11, -40029.8357, -1.5841E-10, -5.645E-11, -38792.318,
    -1.4951E-10, -5.5208E-11, -37554.8004, -1.406E-10, -5.3966E-11, -36317.2827,
    -1.317E-10, -5.2724E-11, -35079.7651, -1.228E-10, -5.1482E-11, -33842.2474,
    -1.139E-10, -5.024E-11, -32604.7298, -1.05E-10, -4.8998E-11, -31367.2121,
    -9.6097E-11, -4.7756E-11, -50129.6945, -8.7195E-11, -4.6514E-11, -26925.1886,
    -7.8293E-11, -4.5272E-11, -23726.1941, -6.9392E-11, -4.403E-11, -20538.2225,
    -6.049E-11, -4.2788E-11, -17366.7852, -5.1589E-11, -4.1546E-11, -14217.3936,
    -4.2687E-11, -4.0304E-11, -11095.5591, -3.3786E-11, -7.0684E-11, -1.8471E-10,
    -46217.4239, -6.7915E-11, -1.7721E-10, -44979.9063, -6.5147E-11, -1.6971E-10,
    -43742.3886, -6.2378E-11, -1.6221E-10, -42504.871, -5.961E-11, -1.547E-10,
    -41267.3533, -5.6841E-11, -1.472E-10, -40029.8357, -5.4073E-11, -1.397E-10,
    -38792.318, -5.1305E-11, -1.322E-10, -37554.8004, -4.8536E-11, -1.247E-10,
    -36317.2827, -4.5768E-11, -1.172E-10, -35079.7651, -4.2999E-11, -1.097E-10,
    -33842.2474, -4.0231E-11, -1.022E-10, -32604.7298, -3.7462E-11, -9.4696E-11,
    -31367.2121, -3.4694E-11, -8.7195E-11, -50129.6945, -3.1925E-11, -7.9694E-11,
    -26925.1886, -2.9157E-11, -7.2193E-11, -23726.1941, -2.6388E-11, -6.4692E-11,
    -20538.2225, -2.362E-11, -5.7191E-11, -17366.7852, -2.0852E-11, -4.969E-11,
    -14217.3936, -1.8083E-11, -4.2189E-11, -11095.5591, -41014.3576, -8.9913E-11,
    -5.4647E-11, -39930.5754, -8.6574E-11, -5.2899E-11, -38846.7932, -8.3236E-11,
    -5.1151E-11, -37763.0109, -7.9898E-11, -4.9404E-11, -36679.2287, -7.6559E-11,
    -4.7656E-11, -35595.4465, -7.3221E-11, -4.5908E-11, -34511.6642, -6.9882E-11,
    -4.416E-11, -33427.882, -6.6544E-11, -4.2412E-11, -32344.0998, -6.3206E-11,
    -4.0664E-11, -31260.3175, -5.9867E-11, -3.8917E-11, -30176.5353, -5.6529E-11,
    -3.7169E-11, -29092.7531, -5.3191E-11, -3.5421E-11, -28008.9708, -4.9852E-11,
    -3.3673E-11, -26925.1886, -4.6514E-11, -3.1925E-11, -45841.4064, -4.3175E-11,
    -3.0178E-11, -22790.6359, -3.9837E-11, -2.843E-11, -19745.3769, -3.6499E-11,
    -2.6682E-11, -16711.1407, -3.316E-11, -2.4934E-11, -13693.4388, -2.9822E-11,
    -2.3186E-11, -10697.7825, -2.6484E-11, -2.1438E-11, -5.8813E-11, -41014.3576,
    -1.8529E-10, -5.7696E-11, -39930.5754, -1.7717E-10, -5.6579E-11, -38846.7932,
    -1.6905E-10, -5.5462E-11, -37763.0109, -1.6092E-10, -5.4345E-11, -36679.2287,
    -1.528E-10, -5.3228E-11, -35595.4465, -1.4468E-10, -5.2111E-11, -34511.6642,
    -1.3655E-10, -5.0994E-11, -33427.882, -1.2843E-10, -4.9877E-11, -32344.0998,
    -1.2031E-10, -4.876E-11, -31260.3175, -1.1219E-10, -4.7643E-11, -30176.5353,
    -1.0406E-10, -4.6526E-11, -29092.7531, -9.594E-11, -4.5409E-11, -28008.9708,
    -8.7817E-11, -4.4292E-11, -26925.1886, -7.9694E-11, -4.3175E-11, -45841.4064,
    -7.1571E-11, -4.2058E-11, -22790.6359, -6.3448E-11, -4.0941E-11, -19745.3769,
    -5.5325E-11, -3.9825E-11, -16711.1407, -4.7202E-11, -3.8708E-11, -13693.4388,
    -3.9079E-11, -3.7591E-11, -10697.7825, -3.0956E-11, -6.5958E-11, -1.6569E-10,
    -41014.3576, -6.3402E-11, -1.5896E-10, -39930.5754, -6.0846E-11, -1.5224E-10,
    -38846.7932, -5.829E-11, -1.4552E-10, -37763.0109, -5.5735E-11, -1.388E-10,
    -36679.2287, -5.3179E-11, -1.3207E-10, -35595.4465, -5.0623E-11, -1.2535E-10,
    -34511.6642, -4.8068E-11, -1.1863E-10, -33427.882, -4.5512E-11, -1.1191E-10,
    -32344.0998, -4.2956E-11, -1.0518E-10, -31260.3175, -4.04E-11, -9.8461E-11,
    -30176.5353, -3.7845E-11, -9.1738E-11, -29092.7531, -3.5289E-11, -8.5016E-11,
    -28008.9708, -3.2733E-11, -7.8293E-11, -26925.1886, -3.0178E-11, -7.1571E-11,
    -45841.4064, -2.7622E-11, -6.4849E-11, -22790.6359, -2.5066E-11, -5.8126E-11,
    -19745.3769, -2.251E-11, -5.1404E-11, -16711.1407, -1.9955E-11, -4.4681E-11,
    -13693.4388, -1.7399E-11, -3.7959E-11, -10697.7825, -35888.4513, -8.7044E-11,
    -4.9113E-11, -34952.8931, -8.3831E-11, -4.7578E-11, -34017.3348, -8.0618E-11,
    -4.6043E-11, -33081.7766, -7.7405E-11, -4.4508E-11, -32146.2183, -7.4191E-11,
    -4.2973E-11, -31210.6601, -7.0978E-11, -4.1438E-11, -30275.1018, -6.7765E-11,
    -3.9903E-11, -29339.5436, -6.4551E-11, -3.8367E-11, -28403.9854, -6.1338E-11,
    -3.6832E-11, -27468.4271, -5.8125E-11, -3.5297E-11, -26532.8689, -5.4912E-11,
    -3.3762E-11, -25597.3106, -5.1698E-11, -3.2227E-11, -24661.7524, -4.8485E-11,
    -3.0692E-11, -23726.1941, -4.5272E-11, -2.9157E-11, -22790.6359, -4.2058E-11,
    -2.7622E-11, -41855.0777, -3.8845E-11, -2.6087E-11, -18952.5312, -3.5632E-11,
    -2.4552E-11, -16055.4961, -3.2419E-11, -2.3017E-11, -13169.4839, -2.9205E-11,
    -2.1481E-11, -10300.006, -2.5992E-11, -1.9946E-11, -5.3723E-11, -35888.4513,
    -1.6767E-10, -5.2731E-11, -34952.8931, -1.6032E-10, -5.1739E-11, -34017.3348,
    -1.5298E-10, -5.0748E-11, -33081.7766, -1.4564E-10, -4.9756E-11, -32146.2183,
    -1.3829E-10, -4.8764E-11, -31210.6601, -1.3095E-10, -4.7772E-11, -30275.1018,
    -1.236E-10, -4.678E-11, -29339.5436, -1.1626E-10, -4.5788E-11, -28403.9854,
    -1.0891E-10, -4.4796E-11, -27468.4271, -1.0157E-10, -4.3805E-11, -26532.8689,
    -9.4226E-11, -4.2813E-11, -25597.3106, -8.6882E-11, -4.1821E-11, -24661.7524,
    -7.9537E-11, -4.0829E-11, -23726.1941, -7.2193E-11, -3.9837E-11, -22790.6359,
    -6.4849E-11, -3.8845E-11, -41855.0777, -5.7504E-11, -3.7853E-11, -18952.5312,
    -5.016E-11, -3.6861E-11, -16055.4961, -4.2816E-11, -3.587E-11, -13169.4839,
    -3.5471E-11, -3.4878E-11, -10300.006, -2.8127E-11, -6.1232E-11, -1.4666E-10,
    -35888.4513, -5.8889E-11, -1.4072E-10, -34952.8931, -5.6546E-11, -1.3477E-10,
    -34017.3348, -5.4203E-11, -1.2883E-10, -33081.7766, -5.186E-11, -1.2289E-10,
    -32146.2183, -4.9517E-11, -1.1694E-10, -31210.6601, -4.7174E-11, -1.11E-10,
    -30275.1018, -4.4831E-11, -1.0505E-10, -29339.5436, -4.2488E-11, -9.9111E-11,
    -28403.9854, -4.0145E-11, -9.3167E-11, -27468.4271, -3.7802E-11, -8.7223E-11,
    -26532.8689, -3.5459E-11, -8.128E-11, -25597.3106, -3.3116E-11, -7.5336E-11,
    -24661.7524, -3.0773E-11, -6.9392E-11, -23726.1941, -2.843E-11, -6.3448E-11,
    -22790.6359, -2.6087E-11, -5.7504E-11, -41855.0777, -2.3744E-11, -5.156E-11,
    -18952.5312, -2.1401E-11, -4.5617E-11, -16055.4961, -1.9058E-11, -3.9673E-11,
    -13169.4839, -1.6715E-11, -3.3729E-11, -10300.006, -30845.2164, -8.4176E-11,
    -4.3579E-11, -30052.3707, -8.1088E-11, -4.2257E-11, -29259.525, -7.8E-11,
    -4.0934E-11, -28466.6793, -7.4912E-11, -3.9612E-11, -27673.8337, -7.1823E-11,
    -3.829E-11, -26880.988, -6.8735E-11, -3.6967E-11, -26088.1423, -6.5647E-11,
    -3.5645E-11, -25295.2966, -6.2559E-11, -3.4323E-11, -24502.4509, -5.9471E-11,
    -3.3E-11, -23709.6053, -5.6382E-11, -3.1678E-11, -22916.7596, -5.3294E-11,
    -3.0356E-11, -22123.9139, -5.0206E-11, -2.9033E-11, -21331.0682, -4.7118E-11,
    -2.7711E-11, -20538.2225, -4.403E-11, -2.6388E-11, -19745.3769, -4.0941E-11,
    -2.5066E-11, -18952.5312, -3.7853E-11, -2.3744E-11, -38159.6855, -3.4765E-11,
    -2.2421E-11, -15399.8516, -3.1677E-11, -2.1099E-11, -12645.5291, -2.8589E-11,
    -1.9777E-11, -9902.2295, -2.5501E-11, -1.8454E-11, -4.8633E-11, -30845.2164,
    -1.5005E-10, -4.7767E-11, -30052.3707, -1.4348E-10, -4.69E-11, -29259.525,
    -1.3691E-10, -4.6033E-11, -28466.6793, -1.3035E-10, -4.5166E-11, -27673.8337,
    -1.2378E-10, -4.43E-11, -26880.988, -1.1722E-10, -4.3433E-11, -26088.1423,
    -1.1065E-10, -4.2566E-11, -25295.2966, -1.0409E-10, -4.1699E-11, -24502.4509,
    -9.752E-11, -4.0833E-11, -23709.6053, -9.0955E-11, -3.9966E-11, -22916.7596,
    -8.4389E-11, -3.9099E-11, -22123.9139, -7.7823E-11, -3.8232E-11, -21331.0682,
    -7.1257E-11, -3.7365E-11, -20538.2225, -6.4692E-11, -3.6499E-11, -19745.3769,
    -5.8126E-11, -3.5632E-11, -18952.5312, -5.156E-11, -3.4765E-11, -38159.6855,
    -4.4995E-11, -3.3898E-11, -15399.8516, -3.8429E-11, -3.3032E-11, -12645.5291,
    -3.1863E-11, -3.2165E-11, -9902.2295, -2.5298E-11, -5.6505E-11, -1.2764E-10,
    -30845.2164, -5.4375E-11, -1.2247E-10, -30052.3707, -5.2245E-11, -1.1731E-10,
    -29259.525, -5.0115E-11, -1.1214E-10, -28466.6793, -4.7984E-11, -1.0698E-10,
    -27673.8337, -4.5854E-11, -1.0181E-10, -26880.988, -4.3724E-11, -9.6647E-11,
    -26088.1423, -4.1594E-11, -9.1482E-11, -25295.2966, -3.9463E-11, -8.6316E-11,
    -24502.4509, -3.7333E-11, -8.1151E-11, -23709.6053, -3.5203E-11, -7.5986E-11,
    -22916.7596, -3.3073E-11, -7.0821E-11, -22123.9139, -3.0942E-11, -6.5656E-11,
    -21331.0682, -2.8812E-11, -6.049E-11, -20538.2225, -2.6682E-11, -5.5325E-11,
    -19745.3769, -2.4552E-11, -5.016E-11, -18952.5312, -2.2421E-11, -4.4995E-11,
    -38159.6855, -2.0291E-11, -3.9829E-11, -15399.8516, -1.8161E-11, -3.4664E-11,
    -12645.5291, -1.6031E-11, -2.9499E-11, -9902.2295, -25890.1642, -8.1308E-11,
    -3.8045E-11, -25234.5197, -7.8345E-11, -3.6936E-11, -24578.8751, -7.5382E-11,
    -3.5826E-11, -23923.2306, -7.2419E-11, -3.4716E-11, -23267.5861, -6.9456E-11,
    -3.3607E-11, -22611.9415, -6.6492E-11, -3.2497E-11, -21956.297, -6.3529E-11,
    -3.1387E-11, -21300.6524, -6.0566E-11, -3.0278E-11, -20645.0079, -5.7603E-11,
    -2.9168E-11, -19989.3634, -5.464E-11, -2.8058E-11, -19333.7188, -5.1677E-11,
    -2.6949E-11, -18678.0743, -4.8714E-11, -2.5839E-11, -18022.4297, -4.5751E-11,
    -2.473E-11, -17366.7852, -4.2788E-11, -2.362E-11, -16711.1407, -3.9825E-11,
    -2.251E-11, -16055.4961, -3.6861E-11, -2.1401E-11, -15399.8516, -3.3898E-11,
    -2.0291E-11, -34744.207, -3.0935E-11, -1.9181E-11, -12121.5743, -2.7972E-11,
    -1.8072E-11, -9504.4529, -2.5009E-11, -1.6962E-11, -4.3544E-11, -25890.1642,
    -1.3242E-10, -4.2802E-11, -25234.5197, -1.2664E-10, -4.206E-11, -24578.8751,
    -1.2085E-10, -4.1319E-11, -23923.2306, -1.1506E-10, -4.0577E-11, -23267.5861,
    -1.0927E-10, -3.9835E-11, -22611.9415, -1.0349E-10, -3.9094E-11, -21956.297,
    -9.77E-11, -3.8352E-11, -21300.6524, -9.1913E-11, -3.761E-11, -20645.0079,
    -8.6126E-11, -3.6869E-11, -19989.3634, -8.0339E-11, -3.6127E-11, -19333.7188,
    -7.4552E-11, -3.5385E-11, -18678.0743, -6.8765E-11, -3.4644E-11, -18022.4297,
    -6.2978E-11, -3.3902E-11, -17366.7852, -5.7191E-11, -3.316E-11, -16711.1407,
    -5.1404E-11, -3.2419E-11, -16055.4961, -4.5617E-11, -3.1677E-11, -15399.8516,
    -3.9829E-11, -3.0935E-11, -34744.207, -3.4042E-11, -3.0194E-11, -12121.5743,
    -2.8255E-11, -2.9452E-11, -9504.4529, -2.2468E-11, -5.1779E-11, -1.0861E-10,
    -25890.1642, -4.9862E-11, -1.0423E-10, -25234.5197, -4.7944E-11, -9.9841E-11,
    -24578.8751, -4.6027E-11, -9.5455E-11, -23923.2306, -4.4109E-11, -9.1068E-11,
    -23267.5861, -4.2192E-11, -8.6682E-11, -22611.9415, -4.0274E-11, -8.2295E-11,
    -21956.297, -3.8357E-11, -7.7908E-11, -21300.6524, -3.6439E-11, -7.3522E-11,
    -20645.0079, -3.4522E-11, -6.9135E-11, -19989.3634, -3.2604E-11, -6.4749E-11,
    -19333.7188, -3.0687E-11, -6.0362E-11, -18678.0743, -2.8769E-11, -5.5975E-11,
    -18022.4297, -2.6852E-11, -5.1589E-11, -17366.7852, -2.4934E-11, -4.7202E-11,
    -16711.1407, -2.3017E-11, -4.2816E-11, -16055.4961, -2.1099E-11, -3.8429E-11,
    -15399.8516, -1.9181E-11, -3.4042E-11, -34744.207, -1.7264E-11, -2.9656E-11,
    -12121.5743, -1.5346E-11, -2.5269E-11, -9504.4529, -21028.8063, -7.844E-11,
    -3.2511E-11, -20504.8515, -7.5602E-11, -3.1614E-11, -19980.8967, -7.2764E-11,
    -3.0717E-11, -19456.9418, -6.9926E-11, -2.982E-11, -18932.987, -6.7088E-11,
    -2.8924E-11, -18409.0322, -6.425E-11, -2.8027E-11, -17885.0774, -6.1412E-11,
    -2.713E-11, -17361.1225, -5.8574E-11, -2.6233E-11, -16837.1677, -5.5736E-11,
    -2.5336E-11, -16313.2129, -5.2898E-11, -2.4439E-11, -15789.2581, -5.006E-11,
    -2.3542E-11, -15265.3032, -4.7222E-11, -2.2645E-11, -14741.3484, -4.4384E-11,
    -2.1748E-11, -14217.3936, -4.1546E-11, -2.0852E-11, -13693.4388, -3.8708E-11,
    -1.9955E-11, -13169.4839, -3.587E-11, -1.9058E-11, -12645.5291, -3.3032E-11,
    -1.8161E-11, -12121.5743, -3.0194E-11, -1.7264E-11, -31597.6194, -2.7356E-11,
    -1.6367E-11, -9106.6764, -2.4518E-11, -1.547E-11, -3.8454E-11, -21028.8063,
    -1.148E-10, -3.7838E-11, -20504.8515, -1.0979E-10, -3.7221E-11, -19980.8967,
    -1.0478E-10, -3.6604E-11, -19456.9418, -9.9774E-11, -3.5988E-11, -18932.987,
    -9.4766E-11, -3.5371E-11, -18409.0322, -8.9758E-11, -3.4755E-11, -17885.0774,
    -8.4749E-11, -3.4138E-11, -17361.1225, -7.9741E-11, -3.3521E-11, -16837.1677,
    -7.4732E-11, -3.2905E-11, -16313.2129, -6.9724E-11, -3.2288E-11, -15789.2581,
    -6.4715E-11, -3.1672E-11, -15265.3032, -5.9707E-11, -3.1055E-11, -14741.3484,
    -5.4698E-11, -3.0438E-11, -14217.3936, -4.969E-11, -2.9822E-11, -13693.4388,
    -4.4681E-11, -2.9205E-11, -13169.4839, -3.9673E-11, -2.8589E-11, -12645.5291,
    -3.4664E-11, -2.7972E-11, -12121.5743, -2.9656E-11, -2.7356E-11, -31597.6194,
    -2.4647E-11, -2.6739E-11, -9106.6764, -1.9639E-11, -4.7053E-11, -8.9591E-11,
    -21028.8063, -4.5348E-11, -8.5983E-11, -20504.8515, -4.3644E-11, -8.2375E-11,
    -19980.8967, -4.1939E-11, -7.8767E-11, -19456.9418, -4.0234E-11, -7.5159E-11,
    -18932.987, -3.8529E-11, -7.1551E-11, -18409.0322, -3.6825E-11, -6.7943E-11,
    -17885.0774, -3.512E-11, -6.4335E-11, -17361.1225, -3.3415E-11, -6.0727E-11,
    -16837.1677, -3.171E-11, -5.7119E-11, -16313.2129, -3.0005E-11, -5.3511E-11,
    -15789.2581, -2.8301E-11, -4.9903E-11, -15265.3032, -2.6596E-11, -4.6295E-11,
    -14741.3484, -2.4891E-11, -4.2687E-11, -14217.3936, -2.3186E-11, -3.9079E-11,
    -13693.4388, -2.1481E-11, -3.5471E-11, -13169.4839, -1.9777E-11, -3.1863E-11,
    -12645.5291, -1.8072E-11, -2.8255E-11, -12121.5743, -1.6367E-11, -2.4647E-11,
    -31597.6194, -1.4662E-11, -2.1039E-11, -9106.6764, -16266.654, -7.5571E-11,
    -2.6977E-11, -15868.8775, -7.2858E-11, -2.6293E-11, -15471.101, -7.0146E-11,
    -2.5609E-11, -15073.3244, -6.7433E-11, -2.4925E-11, -14675.5479, -6.472E-11,
    -2.4241E-11, -14277.7714, -6.2007E-11, -2.3556E-11, -13879.9948, -5.9294E-11,
    -2.2872E-11, -13482.2183, -5.6581E-11, -2.2188E-11, -13084.4418, -5.3868E-11,
    -2.1504E-11, -12686.6652, -5.1155E-11, -2.082E-11, -12288.8887, -4.8442E-11,
    -2.0136E-11, -11891.1121, -4.5729E-11, -1.9451E-11, -11493.3356, -4.3016E-11,
    -1.8767E-11, -11095.5591, -4.0304E-11, -1.8083E-11, -10697.7825, -3.7591E-11,
    -1.7399E-11, -10300.006, -3.4878E-11, -1.6715E-11, -9902.2295, -3.2165E-11,
    -1.6031E-11, -9504.4529, -2.9452E-11, -1.5346E-11, -9106.6764, -2.6739E-11,
    -1.4662E-11, -28708.8999, -2.4026E-11, -1.3978E-11, -3.3364E-11, -16266.654,
    -9.7177E-11, -3.2873E-11, -15868.8775, -9.2947E-11, -3.2381E-11, -15471.101,
    -8.8717E-11, -3.189E-11, -15073.3244, -8.4487E-11, -3.1398E-11, -14675.5479,
    -8.0257E-11, -3.0907E-11, -14277.7714, -7.6028E-11, -3.0415E-11, -13879.9948,
    -7.1798E-11, -2.9924E-11, -13482.2183, -6.7568E-11, -2.9432E-11, -13084.4418,
    -6.3338E-11, -2.8941E-11, -12686.6652, -5.9108E-11, -2.8449E-11, -12288.8887,
    -5.4878E-11, -2.7958E-11, -11891.1121, -5.0648E-11, -2.7467E-11, -11493.3356,
    -4.6418E-11, -2.6975E-11, -11095.5591, -4.2189E-11, -2.6484E-11, -10697.7825,
    -3.7959E-11, -2.5992E-11, -10300.006, -3.3729E-11, -2.5501E-11, -9902.2295,
    -2.9499E-11, -2.5009E-11, -9504.4529, -2.5269E-11, -2.4518E-11, -9106.6764,
    -2.1039E-11, -2.4026E-11, -28708.8999, -1.6809E-11, -4.2327E-11, -7.0567E-11,
    -16266.654, -4.0835E-11, -6.7738E-11, -15868.8775, -3.9343E-11, -6.4909E-11,
    -15471.101, -3.7851E-11, -6.2079E-11, -15073.3244, -3.6359E-11, -5.925E-11,
    -14675.5479, -3.4867E-11, -5.6421E-11, -14277.7714, -3.3375E-11, -5.3591E-11,
    -13879.9948, -3.1883E-11, -5.0762E-11, -13482.2183, -3.0391E-11, -4.7933E-11,
    -13084.4418, -2.8899E-11, -4.5103E-11, -12686.6652, -2.7407E-11, -4.2274E-11,
    -12288.8887, -2.5915E-11, -3.9444E-11, -11891.1121, -2.4423E-11, -3.6615E-11,
    -11493.3356, -2.293E-11, -3.3786E-11, -11095.5591, -2.1438E-11, -3.0956E-11,
    -10697.7825, -1.9946E-11, -2.8127E-11, -10300.006, -1.8454E-11, -2.5298E-11,
    -9902.2295, -1.6962E-11, -2.2468E-11, -9504.4529, -1.547E-11, -1.9639E-11,
    -9106.6764, -1.3978E-11, -1.6809E-11, -28708.8999 };

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

  //  R_MPC=[100 0 0 ;0 100 0 ;0 0 100];
  //  V_TRMPC=[0.25752 ;0.25752 ;0.25752 ;0.25752 ;0.25752 ;0.25752];
  //  Hp = 20;
  //  X_QP = zeros(3*Hp,1);
  //  solution matrices x_traj = A_tilde x_0 + B_tilde u_traj
  //  n=6;
  //  n = size(A,1);
  //  m=3;
  //  m= size(B,2);
  //  A_tilde = zeros((Hp)*n,n);
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
  //  cost
  //  Q_tilde = blkdiag(kron(eye(Hp-1),Q_MPC),P_MPC);
  //  R_tilde = kron(eye(Hp),R_MPC);
  //  R_tilde=[100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 100 ]; 
  //  F_tilde = kron([ones(m-1,1);zeros(m*(Hp-1)+1,1)],F_MPC);
  //  L3 = blkdiag(kron(eye(Hp-1),Q_MPC),zeros(n,n));
  //  L4 = zeros(n,Hp*n);
  //  L4(:,end-(n-1):end) = eye(n);
  //  PHI = blkdiag(kron(eye(Hp),eye(6,6)));
  //  constraints
  //  CONTROL INPUTS
  //  Define and Solve QP
  //  GAMMA = zeros(n*Hp,1);%Dimension of GAMMA is n*horizon x n
  //  for i = 1:Hp
  //      GAMMA(i*n-5:i*n,1) = xfinal;
  //  end
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
  //  b = [kron(ones(Hp,1), V_TRMPC(1:m,1)); kron(ones(Hp,1), V_TRMPC(m+1:end,1))]; 
  //  Call QP Solver
  //  Call QP Solver
  //  QQ = (QQ+QQ')/2;
  for (b_i = 0; b_i < 60; b_i++) {
    d = 0.0;
    for (i = 0; i < 120; i++) {
      d += ftol[i] * b_b[i + (120 * b_i)];
    }

    H[b_i] = d;
    X_QP[b_i] = 0.0;
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

      igr2[i] = nlamt[i] * ((d + 0.22901) - mesil[i]);
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

