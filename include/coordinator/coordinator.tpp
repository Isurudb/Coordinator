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
    double r=0, p=0, y=45;  // Rotate the previous pose by 45* about Z

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


    position_.x = px + (cm_x*R_11 + cm_y*R_21 + cm_z*R_31);
    position_.y = py + (cm_x*R_12 + cm_y*R_22 + cm_z*R_32);
    position_.z = pz + (cm_x*R_13 + cm_y*R_23 + cm_z*R_33);

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

   double a[18] = { 0.025353, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.0,
    0.025353, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.0, 0.025353, 0.0, 0.0,
    0.025353 };

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

 double a_1[18] = { -0.0032940000000000009, -0.0, -0.0, -0.0,
    -0.16079815, -0.0, -0.0, -0.0, -0.16096515, -0.1647, -0.0, -0.0, -0.0,
    -8.0399074999999982, -0.0, -0.0, -0.0, -8.0482574999999983  };

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
    0.0, 0.0, 0.0, 0.0, 2050.715, -2.1501E-11, -8.3626E-12, 33026.8856,
    -1.5911E-9, -1.5075E-10, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -2.1501E-11, 2050.715, 1.9868E-11, -6.5574E-11, 33026.8856, 5.3886E-10, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -8.3626E-12, 1.9868E-11, 2050.715,
    -4.5613E-10, 1.1167E-9, 33026.8856, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 33026.8856, -6.5574E-11, -4.5613E-10, 1.9520796866E+6, -7.9661E-9,
    -1.1146E-8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.5911E-9,
    33026.8856, 1.1167E-9, -7.9661E-9, 1.9520796866E+6, 2.6398E-8, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.5075E-10, 5.3886E-10, 33026.8856,
    -1.1146E-8, 2.6398E-8, 1.9520796866E+6 };

  static const real_T b_b[7200] = { 0.025353, 0.0, 0.0, 0.025353, 0.0, 0.0,
    0.076059, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.12677, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.17747, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.22818, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.27888, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.32959, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.38029, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.431, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.48171, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.53241, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.58312, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.63382, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.68453, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.73524, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.78594, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.83665, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.88736, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.93806, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.98877, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.0, 0.025353, 0.0, 0.0,
    0.025353, 0.0, 0.0, 0.076059, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.12677, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.17747, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.22818, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.27888, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.32959, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.38029, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.431, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.48171, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.53241, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.58312, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.63382, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.68453, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.73524, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.78594, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.83665, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.88736, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.93806, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.98877, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.0,
    0.025353, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.076059, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.12677, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.17747, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.22818, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.27888, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.32959, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.38029, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.431, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.48171, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.53241, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.58312, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.63382, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.68453, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.73524, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.78594, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.83665, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.88736, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.93806, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.98877, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.076059,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.12677, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.17747,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.22818, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.27888,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.32959, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.38029,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.431, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.48171,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.53241, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.58312,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.63382, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.68453,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.73524, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.78594,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.83665, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.88736,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.93806, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.025353, 0.0, 0.0,
    0.076059, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.12677, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.17747, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.22818, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.27888, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.32959, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.38029, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.431, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.48171, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.53241, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.58312, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.63382, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.68453, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.73524, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.78594, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.83665, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.88736, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.93806, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.076059, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.12677, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.17747, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.22818, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.27888, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.32959, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.38029, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.431, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.48171, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.53241, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.58312, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.63382, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.68453, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.73524, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.78594, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.83665, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.88736, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.93806, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.076059, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.12677,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.17747, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.22818,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.27888, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.32959,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.38029, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.431,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.48171, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.53241,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.58312, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.63382,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.68453, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.73524,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.78594, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.83665,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.88736, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.076059, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.12677,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.17747, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.22818,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.27888, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.32959,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.38029, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.431,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.48171, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.53241,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.58312, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.63382,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.68453, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.73524,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.78594, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.83665,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.88736, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.076059, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.12677,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.17747, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.22818,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.27888, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.32959,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.38029, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.431,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.48171, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.53241,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.58312, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.63382,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.68453, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.73524,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.78594, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.83665,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.88736, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.025353, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.076059, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.12677, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.17747, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.22818, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.27888, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.32959, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.38029, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.431, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.48171, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.53241, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.58312, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.63382, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.68453, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.73524, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.78594, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.83665, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.076059, 0.0, 0.0, 0.025353, 0.0, 0.0,
    0.12677, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.17747, 0.0, 0.0, 0.025353, 0.0, 0.0,
    0.22818, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.27888, 0.0, 0.0, 0.025353, 0.0, 0.0,
    0.32959, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.38029, 0.0, 0.0, 0.025353, 0.0, 0.0,
    0.431, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.48171, 0.0, 0.0, 0.025353, 0.0, 0.0,
    0.53241, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.58312, 0.0, 0.0, 0.025353, 0.0, 0.0,
    0.63382, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.68453, 0.0, 0.0, 0.025353, 0.0, 0.0,
    0.73524, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.78594, 0.0, 0.0, 0.025353, 0.0, 0.0,
    0.83665, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.076059, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.12677,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.17747, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.22818,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.27888, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.32959,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.38029, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.431,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.48171, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.53241,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.58312, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.63382,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.68453, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.73524,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.78594, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.83665,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.076059, 0.0, 0.0, 0.025353, 0.0, 0.0,
    0.12677, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.17747, 0.0, 0.0, 0.025353, 0.0, 0.0,
    0.22818, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.27888, 0.0, 0.0, 0.025353, 0.0, 0.0,
    0.32959, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.38029, 0.0, 0.0, 0.025353, 0.0, 0.0,
    0.431, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.48171, 0.0, 0.0, 0.025353, 0.0, 0.0,
    0.53241, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.58312, 0.0, 0.0, 0.025353, 0.0, 0.0,
    0.63382, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.68453, 0.0, 0.0, 0.025353, 0.0, 0.0,
    0.73524, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.78594, 0.0, 0.0, 0.025353, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.025353, 0.0, 0.0,
    0.025353, 0.0, 0.0, 0.076059, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.12677, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.17747, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.22818, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.27888, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.32959, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.38029, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.431, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.48171, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.53241, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.58312, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.63382, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.68453, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.73524, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.78594, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.076059, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.12677, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.17747, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.22818, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.27888, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.32959, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.38029, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.431, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.48171, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.53241, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.58312, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.63382, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.68453, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.73524, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.78594, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.076059, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.12677, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.17747, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.22818, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.27888, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.32959, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.38029, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.431, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.48171, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.53241, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.58312, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.63382, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.68453, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.73524, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.076059, 0.0, 0.0,
    0.025353, 0.0, 0.0, 0.12677, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.17747, 0.0, 0.0,
    0.025353, 0.0, 0.0, 0.22818, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.27888, 0.0, 0.0,
    0.025353, 0.0, 0.0, 0.32959, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.38029, 0.0, 0.0,
    0.025353, 0.0, 0.0, 0.431, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.48171, 0.0, 0.0,
    0.025353, 0.0, 0.0, 0.53241, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.58312, 0.0, 0.0,
    0.025353, 0.0, 0.0, 0.63382, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.68453, 0.0, 0.0,
    0.025353, 0.0, 0.0, 0.73524, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.076059, 0.0, 0.0, 0.025353, 0.0, 0.0,
    0.12677, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.17747, 0.0, 0.0, 0.025353, 0.0, 0.0,
    0.22818, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.27888, 0.0, 0.0, 0.025353, 0.0, 0.0,
    0.32959, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.38029, 0.0, 0.0, 0.025353, 0.0, 0.0,
    0.431, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.48171, 0.0, 0.0, 0.025353, 0.0, 0.0,
    0.53241, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.58312, 0.0, 0.0, 0.025353, 0.0, 0.0,
    0.63382, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.68453, 0.0, 0.0, 0.025353, 0.0, 0.0,
    0.73524, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.076059, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.12677,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.17747, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.22818,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.27888, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.32959,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.38029, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.431,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.48171, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.53241,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.58312, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.63382,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.68453, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.025353, 0.0, 0.0,
    0.076059, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.12677, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.17747, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.22818, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.27888, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.32959, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.38029, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.431, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.48171, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.53241, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.58312, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.63382, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.68453, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.076059, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.12677, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.17747, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.22818, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.27888, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.32959, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.38029, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.431, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.48171, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.53241, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.58312, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.63382, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.68453, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.076059, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.12677,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.17747, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.22818,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.27888, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.32959,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.38029, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.431,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.48171, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.53241,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.58312, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.63382,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.076059, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.12677, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.17747, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.22818, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.27888, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.32959, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.38029, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.431, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.48171, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.53241, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.58312, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.63382, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.025353, 0.0, 0.0,
    0.025353, 0.0, 0.0, 0.076059, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.12677, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.17747, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.22818, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.27888, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.32959, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.38029, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.431, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.48171, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.53241, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.58312, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.63382, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.025353, 0.0, 0.0,
    0.076059, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.12677, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.17747, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.22818, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.27888, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.32959, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.38029, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.431, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.48171, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.53241, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.58312, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.076059, 0.0, 0.0, 0.025353, 0.0, 0.0,
    0.12677, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.17747, 0.0, 0.0, 0.025353, 0.0, 0.0,
    0.22818, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.27888, 0.0, 0.0, 0.025353, 0.0, 0.0,
    0.32959, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.38029, 0.0, 0.0, 0.025353, 0.0, 0.0,
    0.431, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.48171, 0.0, 0.0, 0.025353, 0.0, 0.0,
    0.53241, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.58312, 0.0, 0.0, 0.025353, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.076059, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.12677, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.17747, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.22818, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.27888, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.32959, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.38029, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.431, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.48171, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.53241, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.58312, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.076059, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.12677,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.17747, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.22818,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.27888, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.32959,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.38029, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.431,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.48171, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.53241,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.076059, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.12677, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.17747, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.22818, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.27888, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.32959, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.38029, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.431, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.48171, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.53241, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.076059,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.12677, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.17747,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.22818, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.27888,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.32959, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.38029,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.431, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.48171,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.53241, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.076059, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.12677,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.17747, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.22818,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.27888, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.32959,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.38029, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.431,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.48171, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.025353, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.076059, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.12677, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.17747, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.22818, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.27888, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.32959, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.38029, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.431, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.48171, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.076059, 0.0, 0.0,
    0.025353, 0.0, 0.0, 0.12677, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.17747, 0.0, 0.0,
    0.025353, 0.0, 0.0, 0.22818, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.27888, 0.0, 0.0,
    0.025353, 0.0, 0.0, 0.32959, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.38029, 0.0, 0.0,
    0.025353, 0.0, 0.0, 0.431, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.48171, 0.0, 0.0,
    0.025353, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.076059, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.12677, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.17747, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.22818, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.27888, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.32959, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.38029, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.431, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.076059, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.12677, 0.0, 0.0,
    0.025353, 0.0, 0.0, 0.17747, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.22818, 0.0, 0.0,
    0.025353, 0.0, 0.0, 0.27888, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.32959, 0.0, 0.0,
    0.025353, 0.0, 0.0, 0.38029, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.431, 0.0, 0.0,
    0.025353, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.025353, 0.0, 0.0,
    0.025353, 0.0, 0.0, 0.076059, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.12677, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.17747, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.22818, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.27888, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.32959, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.38029, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.431, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.025353, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.076059, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.12677, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.17747, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.22818, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.27888, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.32959, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.38029, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.025353, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.076059, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.12677, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.17747, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.22818, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.27888, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.32959, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.38029, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.025353, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.076059, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.12677, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.17747, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.22818, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.27888, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.32959, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.38029, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.076059, 0.0, 0.0,
    0.025353, 0.0, 0.0, 0.12677, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.17747, 0.0, 0.0,
    0.025353, 0.0, 0.0, 0.22818, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.27888, 0.0, 0.0,
    0.025353, 0.0, 0.0, 0.32959, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.076059, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.12677, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.17747, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.22818, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.27888, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.32959, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.025353, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.076059, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.12677, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.17747, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.22818, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.27888, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.32959, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.076059, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.12677, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.17747, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.22818, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.27888, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.076059, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.12677,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.17747, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.22818,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.27888, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.076059, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.12677, 0.0, 0.0,
    0.025353, 0.0, 0.0, 0.17747, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.22818, 0.0, 0.0,
    0.025353, 0.0, 0.0, 0.27888, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.025353, 0.0, 0.0,
    0.025353, 0.0, 0.0, 0.076059, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.12677, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.17747, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.22818, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.025353, 0.0, 0.0,
    0.076059, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.12677, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.17747, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.22818, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.025353, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.076059, 0.0, 0.0,
    0.025353, 0.0, 0.0, 0.12677, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.17747, 0.0, 0.0,
    0.025353, 0.0, 0.0, 0.22818, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.076059, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.12677, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.17747, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.076059, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.12677, 0.0, 0.0,
    0.025353, 0.0, 0.0, 0.17747, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.076059, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.12677, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.17747, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.025353, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.076059, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.12677, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.025353, 0.0, 0.0,
    0.076059, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.12677, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.025353, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.076059, 0.0, 0.0, 0.025353, 0.0,
    0.0, 0.12677, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.025353, 0.0, 0.0,
    0.025353, 0.0, 0.0, 0.076059, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.025353, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.076059, 0.0, 0.0, 0.025353,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.025353, 0.0, 0.0,
    0.076059, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.025353, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.025353, 0.0, 0.0, 0.025353 };

  static const real_T b_H[3600] = { 33095.9468, -1.3534E-10, -6.1107E-11,
    12611.9638, -1.3302E-10, -5.9096E-11, 12128.3408, -1.3069E-10, -5.7084E-11,
    11645.4377, -1.2837E-10, -5.5073E-11, 11163.6144, -1.2604E-10, -5.3062E-11,
    10683.231, -1.2372E-10, -5.1051E-11, 10204.6473, -1.2139E-10, -4.9039E-11,
    9728.2234, -1.1907E-10, -4.7028E-11, 9254.3191, -1.1675E-10, -4.5017E-11,
    8783.2944, -1.1442E-10, -4.3005E-11, 8315.5094, -1.121E-10, -4.0994E-11,
    7851.3239, -1.0977E-10, -3.8983E-11, 7391.0978, -1.0745E-10, -3.6972E-11,
    6935.1912, -1.0512E-10, -3.496E-11, 6483.9641, -1.028E-10, -3.2949E-11,
    6037.7762, -1.0047E-10, -3.0938E-11, 5596.9877, -9.8149E-11, -2.8926E-11,
    5161.9584, -9.5824E-11, -2.6915E-11, 4733.0484, -9.35E-11, -2.4904E-11,
    4310.6175, -9.1175E-11, -2.2892E-11, -1.3534E-10, 33095.9468, 1.5579E-10,
    -1.291E-10, 12611.9638, 1.5092E-10, -1.2285E-10, 12128.3408, 1.4606E-10,
    -1.166E-10, 11645.4377, 1.412E-10, -1.1035E-10, 11163.6144, 1.3633E-10,
    -1.0411E-10, 10683.231, 1.3147E-10, -9.7861E-11, 10204.6473, 1.2661E-10,
    -9.1614E-11, 9728.2234, 1.2174E-10, -8.5367E-11, 9254.3191, 1.1688E-10,
    -7.912E-11, 8783.2944, 1.1202E-10, -7.2873E-11, 8315.5094, 1.0715E-10,
    -6.6627E-11, 7851.3239, 1.0229E-10, -6.038E-11, 7391.0978, 9.7428E-11,
    -5.4133E-11, 6935.1912, 9.2564E-11, -4.7886E-11, 6483.9641, 8.7701E-11,
    -4.1639E-11, 6037.7762, 8.2838E-11, -3.5392E-11, 5596.9877, 7.7974E-11,
    -2.9146E-11, 5161.9584, 7.3111E-11, -2.2899E-11, 4733.0484, 6.8247E-11,
    -1.6652E-11, 4310.6175, 6.3384E-11, -6.1107E-11, 1.5579E-10, 33095.9468,
    -5.9881E-11, 1.5241E-10, 12611.9638, -5.8655E-11, 1.4903E-10, 12128.3408,
    -5.7429E-11, 1.4566E-10, 11645.4377, -5.6203E-11, 1.4228E-10, 11163.6144,
    -5.4976E-11, 1.389E-10, 10683.231, -5.375E-11, 1.3552E-10, 10204.6473,
    -5.2524E-11, 1.3214E-10, 9728.2234, -5.1298E-11, 1.2877E-10, 9254.3191,
    -5.0072E-11, 1.2539E-10, 8783.2944, -4.8846E-11, 1.2201E-10, 8315.5094,
    -4.762E-11, 1.1863E-10, 7851.3239, -4.6393E-11, 1.1526E-10, 7391.0978,
    -4.5167E-11, 1.1188E-10, 6935.1912, -4.3941E-11, 1.085E-10, 6483.9641,
    -4.2715E-11, 1.0512E-10, 6037.7762, -4.1489E-11, 1.0175E-10, 5596.9877,
    -4.0263E-11, 9.8367E-11, 5161.9584, -3.9037E-11, 9.499E-11, 4733.0484,
    -3.7811E-11, 9.1612E-11, 4310.6175, 12611.9638, -1.291E-10, -5.9881E-11,
    32273.6501, -1.2688E-10, -5.7913E-11, 11806.8715, -1.2467E-10, -5.5944E-11,
    11340.4527, -1.2245E-10, -5.3976E-11, 10874.754, -1.2024E-10, -5.2008E-11,
    10410.135, -1.1803E-10, -5.0039E-11, 9946.9559, -1.1581E-10, -4.8071E-11,
    9485.5766, -1.136E-10, -4.6103E-11, 9026.3569, -1.1138E-10, -4.4135E-11,
    8569.657, -1.0917E-10, -4.2166E-11, 8115.8367, -1.0696E-10, -4.0198E-11,
    7665.2559, -1.0474E-10, -3.823E-11, 7218.2747, -1.0253E-10, -3.6261E-11,
    6775.253, -1.0031E-10, -3.4293E-11, 6336.5508, -9.8099E-11, -3.2325E-11,
    5902.5279, -9.5885E-11, -3.0357E-11, 5473.5444, -9.3671E-11, -2.8388E-11,
    5049.9602, -9.1457E-11, -2.642E-11, 4632.1353, -8.9243E-11, -2.4452E-11,
    4220.4295, -8.7029E-11, -2.2483E-11, -1.3302E-10, 12611.9638, 1.5241E-10,
    -1.2688E-10, 32273.6501, 1.4765E-10, -1.2074E-10, 11806.8715, 1.4289E-10,
    -1.1461E-10, 11340.4527, 1.3813E-10, -1.0847E-10, 10874.754, 1.3337E-10,
    -1.0234E-10, 10410.135, 1.286E-10, -9.62E-11, 9946.9559, 1.2384E-10,
    -9.0063E-11, 9485.5766, 1.1908E-10, -8.3927E-11, 9026.3569, 1.1432E-10,
    -7.7791E-11, 8569.657, 1.0956E-10, -7.1655E-11, 8115.8367, 1.048E-10,
    -6.5518E-11, 7665.2559, 1.0004E-10, -5.9382E-11, 7218.2747, 9.5276E-11,
    -5.3246E-11, 6775.253, 9.0515E-11, -4.7109E-11, 6336.5508, 8.5754E-11,
    -4.0973E-11, 5902.5279, 8.0992E-11, -3.4837E-11, 5473.5444, 7.6231E-11,
    -2.8701E-11, 5049.9602, 7.147E-11, -2.2564E-11, 4632.1353, 6.6709E-11,
    -1.6428E-11, 4220.4295, 6.1948E-11, -5.9096E-11, 1.5092E-10, 12611.9638,
    -5.7913E-11, 1.4765E-10, 32273.6501, -5.673E-11, 1.4437E-10, 11806.8715,
    -5.5546E-11, 1.411E-10, 11340.4527, -5.4363E-11, 1.3782E-10, 10874.754,
    -5.318E-11, 1.3455E-10, 10410.135, -5.1997E-11, 1.3127E-10, 9946.9559,
    -5.0814E-11, 1.28E-10, 9485.5766, -4.9631E-11, 1.2472E-10, 9026.3569,
    -4.8448E-11, 1.2145E-10, 8569.657, -4.7264E-11, 1.1817E-10, 8115.8367,
    -4.6081E-11, 1.1489E-10, 7665.2559, -4.4898E-11, 1.1162E-10, 7218.2747,
    -4.3715E-11, 1.0834E-10, 6775.253, -4.2532E-11, 1.0507E-10, 6336.5508,
    -4.1349E-11, 1.0179E-10, 5902.5279, -4.0166E-11, 9.8516E-11, 5473.5444,
    -3.8983E-11, 9.5241E-11, 5049.9602, -3.7799E-11, 9.1965E-11, 4632.1353,
    -3.6616E-11, 8.869E-11, 4220.4295, 12128.3408, -1.2285E-10, -5.8655E-11,
    11806.8715, -1.2074E-10, -5.673E-11, 31485.4021, -1.1864E-10, -5.4804E-11,
    11035.4678, -1.1654E-10, -5.2879E-11, 10585.8935, -1.1443E-10, -5.0954E-11,
    10137.0391, -1.1233E-10, -4.9028E-11, 9689.2645, -1.1023E-10, -4.7103E-11,
    9242.9298, -1.0812E-10, -4.5178E-11, 8798.3948, -1.0602E-10, -4.3252E-11,
    8356.0195, -1.0392E-10, -4.1327E-11, 7916.164, -1.0181E-10, -3.9402E-11,
    7479.188, -9.971E-11, -3.7477E-11, 7045.4516, -9.7607E-11, -3.5551E-11,
    6615.3148, -9.5503E-11, -3.3626E-11, 6189.1375, -9.34E-11, -3.1701E-11,
    5767.2796, -9.1297E-11, -2.9775E-11, 5350.1011, -8.9193E-11, -2.785E-11,
    4937.962, -8.709E-11, -2.5925E-11, 4531.2221, -8.4986E-11, -2.4E-11,
    4130.2416, -8.2883E-11, -2.2074E-11, -1.3069E-10, 12128.3408, 1.4903E-10,
    -1.2467E-10, 11806.8715, 1.4437E-10, -1.1864E-10, 31485.4021, 1.3971E-10,
    -1.1262E-10, 11035.4678, 1.3506E-10, -1.0659E-10, 10585.8935, 1.304E-10,
    -1.0056E-10, 10137.0391, 1.2574E-10, -9.4538E-11, 9689.2645, 1.2108E-10,
    -8.8513E-11, 9242.9298, 1.1642E-10, -8.2487E-11, 8798.3948, 1.1176E-10,
    -7.6461E-11, 8356.0195, 1.071E-10, -7.0436E-11, 7916.164, 1.0244E-10,
    -6.441E-11, 7479.188, 9.7783E-11, -5.8384E-11, 7045.4516, 9.3124E-11,
    -5.2358E-11, 6615.3148, 8.8465E-11, -4.6333E-11, 6189.1375, 8.3806E-11,
    -4.0307E-11, 5767.2796, 7.9147E-11, -3.4281E-11, 5350.1011, 7.4488E-11,
    -2.8256E-11, 4937.962, 6.9829E-11, -2.223E-11, 4531.2221, 6.517E-11,
    -1.6204E-11, 4130.2416, 6.0511E-11, -5.7084E-11, 1.4606E-10, 12128.3408,
    -5.5944E-11, 1.4289E-10, 11806.8715, -5.4804E-11, 1.3971E-10, 31485.4021,
    -5.3664E-11, 1.3654E-10, 11035.4678, -5.2524E-11, 1.3337E-10, 10585.8935,
    -5.1384E-11, 1.3019E-10, 10137.0391, -5.0244E-11, 1.2702E-10, 9689.2645,
    -4.9104E-11, 1.2385E-10, 9242.9298, -4.7963E-11, 1.2067E-10, 8798.3948,
    -4.6823E-11, 1.175E-10, 8356.0195, -4.5683E-11, 1.1433E-10, 7916.164,
    -4.4543E-11, 1.1115E-10, 7479.188, -4.3403E-11, 1.0798E-10, 7045.4516,
    -4.2263E-11, 1.0481E-10, 6615.3148, -4.1123E-11, 1.0163E-10, 6189.1375,
    -3.9983E-11, 9.8461E-11, 5767.2796, -3.8842E-11, 9.5288E-11, 5350.1011,
    -3.7702E-11, 9.2114E-11, 4937.962, -3.6562E-11, 8.8941E-11, 4531.2221,
    -3.5422E-11, 8.5768E-11, 4130.2416, 11645.4377, -1.166E-10, -5.7429E-11,
    11340.4527, -1.1461E-10, -5.5546E-11, 11035.4678, -1.1262E-10, -5.3664E-11,
    30730.4829, -1.1062E-10, -5.1782E-11, 10297.033, -1.0863E-10, -4.99E-11,
    9863.9431, -1.0664E-10, -4.8017E-11, 9431.5731, -1.0464E-10, -4.6135E-11,
    9000.283, -1.0265E-10, -4.4253E-11, 8570.4326, -1.0066E-10, -4.237E-11,
    8142.3821, -9.8665E-11, -4.0488E-11, 7716.4912, -9.6673E-11, -3.8606E-11,
    7293.1201, -9.468E-11, -3.6723E-11, 6872.6285, -9.2687E-11, -3.4841E-11,
    6455.3766, -9.0694E-11, -3.2959E-11, 6041.7242, -8.8701E-11, -3.1077E-11,
    5632.0313, -8.6708E-11, -2.9194E-11, 5226.6578, -8.4715E-11, -2.7312E-11,
    4825.9637, -8.2722E-11, -2.543E-11, 4430.309, -8.073E-11, -2.3547E-11,
    4040.0536, -7.8737E-11, -2.1665E-11, -1.2837E-10, 11645.4377, 1.4566E-10,
    -1.2245E-10, 11340.4527, 1.411E-10, -1.1654E-10, 11035.4678, 1.3654E-10,
    -1.1062E-10, 30730.4829, 1.3198E-10, -1.0471E-10, 10297.033, 1.2743E-10,
    -9.8792E-11, 9863.9431, 1.2287E-10, -9.2877E-11, 9431.5731, 1.1831E-10,
    -8.6962E-11, 9000.283, 1.1376E-10, -8.1047E-11, 8570.4326, 1.092E-10,
    -7.5132E-11, 8142.3821, 1.0464E-10, -6.9217E-11, 7716.4912, 1.0009E-10,
    -6.3301E-11, 7293.1201, 9.5529E-11, -5.7386E-11, 6872.6285, 9.0973E-11,
    -5.1471E-11, 6455.3766, 8.6416E-11, -4.5556E-11, 6041.7242, 8.1859E-11,
    -3.9641E-11, 5632.0313, 7.7302E-11, -3.3726E-11, 5226.6578, 7.2745E-11,
    -2.7811E-11, 4825.9637, 6.8188E-11, -2.1895E-11, 4430.309, 6.3631E-11,
    -1.598E-11, 4040.0536, 5.9074E-11, -5.5073E-11, 1.412E-10, 11645.4377,
    -5.3976E-11, 1.3813E-10, 11340.4527, -5.2879E-11, 1.3506E-10, 11035.4678,
    -5.1782E-11, 1.3198E-10, 30730.4829, -5.0685E-11, 1.2891E-10, 10297.033,
    -4.9588E-11, 1.2584E-10, 9863.9431, -4.849E-11, 1.2277E-10, 9431.5731,
    -4.7393E-11, 1.197E-10, 9000.283, -4.6296E-11, 1.1663E-10, 8570.4326,
    -4.5199E-11, 1.1356E-10, 8142.3821, -4.4102E-11, 1.1049E-10, 7716.4912,
    -4.3005E-11, 1.0741E-10, 7293.1201, -4.1908E-11, 1.0434E-10, 6872.6285,
    -4.0811E-11, 1.0127E-10, 6455.3766, -3.9713E-11, 9.8201E-11, 6041.7242,
    -3.8616E-11, 9.513E-11, 5632.0313, -3.7519E-11, 9.2059E-11, 5226.6578,
    -3.6422E-11, 8.8988E-11, 4825.9637, -3.5325E-11, 8.5917E-11, 4430.309,
    -3.4228E-11, 8.2845E-11, 4040.0536, 11163.6144, -1.1035E-10, -5.6203E-11,
    10874.754, -1.0847E-10, -5.4363E-11, 10585.8935, -1.0659E-10, -5.2524E-11,
    10297.033, -1.0471E-10, -5.0685E-11, 30008.1726, -1.0283E-10, -4.8845E-11,
    9590.8471, -1.0094E-10, -4.7006E-11, 9173.8817, -9.9061E-11, -4.5167E-11,
    8757.6362, -9.7178E-11, -4.3328E-11, 8342.4705, -9.5296E-11, -4.1488E-11,
    7928.7446, -9.3414E-11, -3.9649E-11, 7516.8185, -9.1531E-11, -3.781E-11,
    7107.0522, -8.9649E-11, -3.597E-11, 6699.8055, -8.7767E-11, -3.4131E-11,
    6295.4384, -8.5884E-11, -3.2292E-11, 5894.3109, -8.4002E-11, -3.0453E-11,
    5496.783, -8.212E-11, -2.8613E-11, 5103.2145, -8.0237E-11, -2.6774E-11,
    4713.9655, -7.8355E-11, -2.4935E-11, 4329.3959, -7.6473E-11, -2.3095E-11,
    3949.8656, -7.4591E-11, -2.1256E-11, -1.2604E-10, 11163.6144, 1.4228E-10,
    -1.2024E-10, 10874.754, 1.3782E-10, -1.1443E-10, 10585.8935, 1.3337E-10,
    -1.0863E-10, 10297.033, 1.2891E-10, -1.0283E-10, 30008.1726, 1.2446E-10,
    -9.7021E-11, 9590.8471, 1.2E-10, -9.1216E-11, 9173.8817, 1.1555E-10,
    -8.5411E-11, 8757.6362, 1.1109E-10, -7.9607E-11, 8342.4705, 1.0664E-10,
    -7.3802E-11, 7928.7446, 1.0219E-10, -6.7998E-11, 7516.8185, 9.773E-11,
    -6.2193E-11, 7107.0522, 9.3276E-11, -5.6389E-11, 6699.8055, 8.8821E-11,
    -5.0584E-11, 6295.4384, 8.4366E-11, -4.4779E-11, 5894.3109, 7.9911E-11,
    -3.8975E-11, 5496.783, 7.5457E-11, -3.317E-11, 5103.2145, 7.1002E-11,
    -2.7366E-11, 4713.9655, 6.6547E-11, -2.1561E-11, 4329.3959, 6.2093E-11,
    -1.5756E-11, 3949.8656, 5.7638E-11, -5.3062E-11, 1.3633E-10, 11163.6144,
    -5.2008E-11, 1.3337E-10, 10874.754, -5.0954E-11, 1.304E-10, 10585.8935,
    -4.99E-11, 1.2743E-10, 10297.033, -4.8845E-11, 1.2446E-10, 30008.1726,
    -4.7791E-11, 1.2149E-10, 9590.8471, -4.6737E-11, 1.1852E-10, 9173.8817,
    -4.5683E-11, 1.1555E-10, 8757.6362, -4.4629E-11, 1.1258E-10, 8342.4705,
    -4.3575E-11, 1.0961E-10, 7928.7446, -4.2521E-11, 1.0664E-10, 7516.8185,
    -4.1467E-11, 1.0368E-10, 7107.0522, -4.0412E-11, 1.0071E-10, 6699.8055,
    -3.9358E-11, 9.7737E-11, 6295.4384, -3.8304E-11, 9.4768E-11, 5894.3109,
    -3.725E-11, 9.1799E-11, 5496.783, -3.6196E-11, 8.883E-11, 5103.2145,
    -3.5142E-11, 8.5861E-11, 4713.9655, -3.4088E-11, 8.2892E-11, 4329.3959,
    -3.3033E-11, 7.9923E-11, 3949.8656, 10683.231, -1.0411E-10, -5.4976E-11,
    10410.135, -1.0234E-10, -5.318E-11, 10137.0391, -1.0056E-10, -5.1384E-11,
    9863.9431, -9.8792E-11, -4.9588E-11, 9590.8471, -9.7021E-11, -4.7791E-11,
    29317.7512, -9.5249E-11, -4.5995E-11, 8916.1903, -9.3477E-11, -4.4199E-11,
    8514.9894, -9.1705E-11, -4.2402E-11, 8114.5083, -8.9934E-11, -4.0606E-11,
    7715.1072, -8.8162E-11, -3.881E-11, 7317.1458, -8.639E-11, -3.7014E-11,
    6920.9842, -8.4618E-11, -3.5217E-11, 6526.9824, -8.2847E-11, -3.3421E-11,
    6135.5002, -8.1075E-11, -3.1625E-11, 5746.8976, -7.9303E-11, -2.9828E-11,
    5361.5346, -7.7531E-11, -2.8032E-11, 4979.7712, -7.576E-11, -2.6236E-11,
    4601.9673, -7.3988E-11, -2.444E-11, 4228.4828, -7.2216E-11, -2.2643E-11,
    3859.6777, -7.0444E-11, -2.0847E-11, -1.2372E-10, 10683.231, 1.389E-10,
    -1.1803E-10, 10410.135, 1.3455E-10, -1.1233E-10, 10137.0391, 1.3019E-10,
    -1.0664E-10, 9863.9431, 1.2584E-10, -1.0094E-10, 9590.8471, 1.2149E-10,
    -9.5249E-11, 29317.7512, 1.1714E-10, -8.9555E-11, 8916.1903, 1.1278E-10,
    -8.3861E-11, 8514.9894, 1.0843E-10, -7.8167E-11, 8114.5083, 1.0408E-10,
    -7.2473E-11, 7715.1072, 9.9727E-11, -6.6779E-11, 7317.1458, 9.5374E-11,
    -6.1085E-11, 6920.9842, 9.1022E-11, -5.5391E-11, 6526.9824, 8.6669E-11,
    -4.9697E-11, 6135.5002, 8.2317E-11, -4.4003E-11, 5746.8976, 7.7964E-11,
    -3.8309E-11, 5361.5346, 7.3612E-11, -3.2615E-11, 4979.7712, 6.9259E-11,
    -2.6921E-11, 4601.9673, 6.4906E-11, -2.1227E-11, 4228.4828, 6.0554E-11,
    -1.5532E-11, 3859.6777, 5.6201E-11, -5.1051E-11, 1.3147E-10, 10683.231,
    -5.0039E-11, 1.286E-10, 10410.135, -4.9028E-11, 1.2574E-10, 10137.0391,
    -4.8017E-11, 1.2287E-10, 9863.9431, -4.7006E-11, 1.2E-10, 9590.8471,
    -4.5995E-11, 1.1714E-10, 29317.7512, -4.4984E-11, 1.1427E-10, 8916.1903,
    -4.3973E-11, 1.114E-10, 8514.9894, -4.2962E-11, 1.0854E-10, 8114.5083,
    -4.195E-11, 1.0567E-10, 7715.1072, -4.0939E-11, 1.028E-10, 7317.1458,
    -3.9928E-11, 9.9936E-11, 6920.9842, -3.8917E-11, 9.7069E-11, 6526.9824,
    -3.7906E-11, 9.4202E-11, 6135.5002, -3.6895E-11, 9.1335E-11, 5746.8976,
    -3.5884E-11, 8.8468E-11, 5361.5346, -3.4873E-11, 8.5602E-11, 4979.7712,
    -3.3861E-11, 8.2735E-11, 4601.9673, -3.285E-11, 7.9868E-11, 4228.4828,
    -3.1839E-11, 7.7001E-11, 3859.6777, 10204.6473, -9.7861E-11, -5.375E-11,
    9946.9559, -9.62E-11, -5.1997E-11, 9689.2645, -9.4538E-11, -5.0244E-11,
    9431.5731, -9.2877E-11, -4.849E-11, 9173.8817, -9.1216E-11, -4.6737E-11,
    8916.1903, -8.9555E-11, -4.4984E-11, 28658.4989, -8.7894E-11, -4.3231E-11,
    8272.3426, -8.6232E-11, -4.1477E-11, 7886.5462, -8.4571E-11, -3.9724E-11,
    7501.4697, -8.291E-11, -3.7971E-11, 7117.4731, -8.1249E-11, -3.6217E-11,
    6734.9163, -7.9588E-11, -3.4464E-11, 6354.1593, -7.7927E-11, -3.2711E-11,
    5975.562, -7.6265E-11, -3.0958E-11, 5599.4843, -7.4604E-11, -2.9204E-11,
    5226.2863, -7.2943E-11, -2.7451E-11, 4856.3279, -7.1282E-11, -2.5698E-11,
    4489.969, -6.9621E-11, -2.3944E-11, 4127.5697, -6.7959E-11, -2.2191E-11,
    3769.4897, -6.6298E-11, -2.0438E-11, -1.2139E-10, 10204.6473, 1.3552E-10,
    -1.1581E-10, 9946.9559, 1.3127E-10, -1.1023E-10, 9689.2645, 1.2702E-10,
    -1.0464E-10, 9431.5731, 1.2277E-10, -9.9061E-11, 9173.8817, 1.1852E-10,
    -9.3477E-11, 8916.1903, 1.1427E-10, -8.7894E-11, 28658.4989, 1.1002E-10,
    -8.231E-11, 8272.3426, 1.0577E-10, -7.6727E-11, 7886.5462, 1.0152E-10,
    -7.1143E-11, 7501.4697, 9.7269E-11, -6.556E-11, 7117.4731, 9.3018E-11,
    -5.9976E-11, 6734.9163, 8.8768E-11, -5.4393E-11, 6354.1593, 8.4517E-11,
    -4.8809E-11, 5975.562, 8.0267E-11, -4.3226E-11, 5599.4843, 7.6017E-11,
    -3.7642E-11, 5226.2863, 7.1766E-11, -3.2059E-11, 4856.3279, 6.7516E-11,
    -2.6476E-11, 4489.969, 6.3266E-11, -2.0892E-11, 4127.5697, 5.9015E-11,
    -1.5309E-11, 3769.4897, 5.4765E-11, -4.9039E-11, 1.2661E-10, 10204.6473,
    -4.8071E-11, 1.2384E-10, 9946.9559, -4.7103E-11, 1.2108E-10, 9689.2645,
    -4.6135E-11, 1.1831E-10, 9431.5731, -4.5167E-11, 1.1555E-10, 9173.8817,
    -4.4199E-11, 1.1278E-10, 8916.1903, -4.3231E-11, 1.1002E-10, 28658.4989,
    -4.2262E-11, 1.0726E-10, 8272.3426, -4.1294E-11, 1.0449E-10, 7886.5462,
    -4.0326E-11, 1.0173E-10, 7501.4697, -3.9358E-11, 9.8961E-11, 7117.4731,
    -3.839E-11, 9.6196E-11, 6734.9163, -3.7422E-11, 9.3432E-11, 6354.1593,
    -3.6454E-11, 9.0667E-11, 5975.562, -3.5486E-11, 8.7902E-11, 5599.4843,
    -3.4517E-11, 8.5138E-11, 5226.2863, -3.3549E-11, 8.2373E-11, 4856.3279,
    -3.2581E-11, 7.9608E-11, 4489.969, -3.1613E-11, 7.6843E-11, 4127.5697,
    -3.0645E-11, 7.4079E-11, 3769.4897, 9728.2234, -9.1614E-11, -5.2524E-11,
    9485.5766, -9.0063E-11, -5.0814E-11, 9242.9298, -8.8513E-11, -4.9104E-11,
    9000.283, -8.6962E-11, -4.7393E-11, 8757.6362, -8.5411E-11, -4.5683E-11,
    8514.9894, -8.3861E-11, -4.3973E-11, 8272.3426, -8.231E-11, -4.2262E-11,
    28029.6958, -8.076E-11, -4.0552E-11, 7658.584, -7.9209E-11, -3.8842E-11,
    7287.8323, -7.7658E-11, -3.7132E-11, 6917.8004, -7.6108E-11, -3.5421E-11,
    6548.8484, -7.4557E-11, -3.3711E-11, 6181.3362, -7.3006E-11, -3.2001E-11,
    5815.6237, -7.1456E-11, -3.029E-11, 5452.071, -6.9905E-11, -2.858E-11,
    5091.038, -6.8355E-11, -2.687E-11, 4732.8846, -6.6804E-11, -2.516E-11,
    4377.9708, -6.5253E-11, -2.3449E-11, 4026.6565, -6.3703E-11, -2.1739E-11,
    3679.3018, -6.2152E-11, -2.0029E-11, -1.1907E-10, 9728.2234, 1.3214E-10,
    -1.136E-10, 9485.5766, 1.28E-10, -1.0812E-10, 9242.9298, 1.2385E-10,
    -1.0265E-10, 9000.283, 1.197E-10, -9.7178E-11, 8757.6362, 1.1555E-10,
    -9.1705E-11, 8514.9894, 1.114E-10, -8.6232E-11, 8272.3426, 1.0726E-10,
    -8.076E-11, 28029.6958, 1.0311E-10, -7.5287E-11, 7658.584, 9.8959E-11,
    -6.9814E-11, 7287.8323, 9.481E-11, -6.4341E-11, 6917.8004, 9.0662E-11,
    -5.8868E-11, 6548.8484, 8.6514E-11, -5.3395E-11, 6181.3362, 8.2366E-11,
    -4.7922E-11, 5815.6237, 7.8218E-11, -4.2449E-11, 5452.071, 7.4069E-11,
    -3.6976E-11, 5091.038, 6.9921E-11, -3.1503E-11, 4732.8846, 6.5773E-11,
    -2.6031E-11, 4377.9708, 6.1625E-11, -2.0558E-11, 4026.6565, 5.7477E-11,
    -1.5085E-11, 3679.3018, 5.3328E-11, -4.7028E-11, 1.2174E-10, 9728.2234,
    -4.6103E-11, 1.1908E-10, 9485.5766, -4.5178E-11, 1.1642E-10, 9242.9298,
    -4.4253E-11, 1.1376E-10, 9000.283, -4.3328E-11, 1.1109E-10, 8757.6362,
    -4.2402E-11, 1.0843E-10, 8514.9894, -4.1477E-11, 1.0577E-10, 8272.3426,
    -4.0552E-11, 1.0311E-10, 28029.6958, -3.9627E-11, 1.0044E-10, 7658.584,
    -3.8702E-11, 9.7782E-11, 7287.8323, -3.7777E-11, 9.5119E-11, 6917.8004,
    -3.6852E-11, 9.2457E-11, 6548.8484, -3.5927E-11, 8.9794E-11, 6181.3362,
    -3.5001E-11, 8.7132E-11, 5815.6237, -3.4076E-11, 8.4469E-11, 5452.071,
    -3.3151E-11, 8.1807E-11, 5091.038, -3.2226E-11, 7.9144E-11, 4732.8846,
    -3.1301E-11, 7.6482E-11, 4377.9708, -3.0376E-11, 7.3819E-11, 4026.6565,
    -2.9451E-11, 7.1157E-11, 3679.3018, 9254.3191, -8.5367E-11, -5.1298E-11,
    9026.3569, -8.3927E-11, -4.9631E-11, 8798.3948, -8.2487E-11, -4.7963E-11,
    8570.4326, -8.1047E-11, -4.6296E-11, 8342.4705, -7.9607E-11, -4.4629E-11,
    8114.5083, -7.8167E-11, -4.2962E-11, 7886.5462, -7.6727E-11, -4.1294E-11,
    7658.584, -7.5287E-11, -3.9627E-11, 27430.6219, -7.3847E-11, -3.796E-11,
    7074.1948, -7.2407E-11, -3.6292E-11, 6718.1277, -7.0966E-11, -3.4625E-11,
    6362.7804, -6.9526E-11, -3.2958E-11, 6008.5131, -6.8086E-11, -3.1291E-11,
    5655.6855, -6.6646E-11, -2.9623E-11, 5304.6578, -6.5206E-11, -2.7956E-11,
    4955.7897, -6.3766E-11, -2.6289E-11, 4609.4413, -6.2326E-11, -2.4622E-11,
    4265.9726, -6.0886E-11, -2.2954E-11, 3925.7434, -5.9446E-11, -2.1287E-11,
    3589.1138, -5.8006E-11, -1.962E-11, -1.1675E-10, 9254.3191, 1.2877E-10,
    -1.1138E-10, 9026.3569, 1.2472E-10, -1.0602E-10, 8798.3948, 1.2067E-10,
    -1.0066E-10, 8570.4326, 1.1663E-10, -9.5296E-11, 8342.4705, 1.1258E-10,
    -8.9934E-11, 8114.5083, 1.0854E-10, -8.4571E-11, 7886.5462, 1.0449E-10,
    -7.9209E-11, 7658.584, 1.0044E-10, -7.3847E-11, 27430.6219, 9.6398E-11,
    -6.8484E-11, 7074.1948, 9.2352E-11, -6.3122E-11, 6718.1277, 8.8306E-11,
    -5.776E-11, 6362.7804, 8.426E-11, -5.2397E-11, 6008.5131, 8.0214E-11,
    -4.7035E-11, 5655.6855, 7.6168E-11, -4.1673E-11, 5304.6578, 7.2122E-11,
    -3.631E-11, 4955.7897, 6.8076E-11, -3.0948E-11, 4609.4413, 6.403E-11,
    -2.5586E-11, 4265.9726, 5.9984E-11, -2.0223E-11, 3925.7434, 5.5938E-11,
    -1.4861E-11, 3589.1138, 5.1892E-11, -4.5017E-11, 1.1688E-10, 9254.3191,
    -4.4135E-11, 1.1432E-10, 9026.3569, -4.3252E-11, 1.1176E-10, 8798.3948,
    -4.237E-11, 1.092E-10, 8570.4326, -4.1488E-11, 1.0664E-10, 8342.4705,
    -4.0606E-11, 1.0408E-10, 8114.5083, -3.9724E-11, 1.0152E-10, 7886.5462,
    -3.8842E-11, 9.8959E-11, 7658.584, -3.796E-11, 9.6398E-11, 27430.6219,
    -3.7078E-11, 9.3838E-11, 7074.1948, -3.6196E-11, 9.1278E-11, 6718.1277,
    -3.5313E-11, 8.8717E-11, 6362.7804, -3.4431E-11, 8.6157E-11, 6008.5131,
    -3.3549E-11, 8.3596E-11, 5655.6855, -3.2667E-11, 8.1036E-11, 5304.6578,
    -3.1785E-11, 7.8476E-11, 4955.7897, -3.0903E-11, 7.5915E-11, 4609.4413,
    -3.0021E-11, 7.3355E-11, 4265.9726, -2.9139E-11, 7.0795E-11, 3925.7434,
    -2.8256E-11, 6.8234E-11, 3589.1138, 8783.2944, -7.912E-11, -5.0072E-11,
    8569.657, -7.7791E-11, -4.8448E-11, 8356.0195, -7.6461E-11, -4.6823E-11,
    8142.3821, -7.5132E-11, -4.5199E-11, 7928.7446, -7.3802E-11, -4.3575E-11,
    7715.1072, -7.2473E-11, -4.195E-11, 7501.4697, -7.1143E-11, -4.0326E-11,
    7287.8323, -6.9814E-11, -3.8702E-11, 7074.1948, -6.8484E-11, -3.7078E-11,
    26860.5573, -6.7155E-11, -3.5453E-11, 6518.455, -6.5825E-11, -3.3829E-11,
    6176.7125, -6.4496E-11, -3.2205E-11, 5835.69, -6.3166E-11, -3.0581E-11,
    5495.7473, -6.1837E-11, -2.8956E-11, 5157.2445, -6.0507E-11, -2.7332E-11,
    4820.5414, -5.9178E-11, -2.5708E-11, 4485.998, -5.7848E-11, -2.4083E-11,
    4153.9743, -5.6519E-11, -2.2459E-11, 3824.8303, -5.5189E-11, -2.0835E-11,
    3498.9258, -5.386E-11, -1.9211E-11, -1.1442E-10, 8783.2944, 1.2539E-10,
    -1.0917E-10, 8569.657, 1.2145E-10, -1.0392E-10, 8356.0195, 1.175E-10,
    -9.8665E-11, 8142.3821, 1.1356E-10, -9.3414E-11, 7928.7446, 1.0961E-10,
    -8.8162E-11, 7715.1072, 1.0567E-10, -8.291E-11, 7501.4697, 1.0173E-10,
    -7.7658E-11, 7287.8323, 9.7782E-11, -7.2407E-11, 7074.1948, 9.3838E-11,
    -6.7155E-11, 26860.5573, 8.9894E-11, -6.1903E-11, 6518.455, 8.595E-11,
    -5.6651E-11, 6176.7125, 8.2006E-11, -5.1399E-11, 5835.69, 7.8062E-11,
    -4.6148E-11, 5495.7473, 7.4118E-11, -4.0896E-11, 5157.2445, 7.0175E-11,
    -3.5644E-11, 4820.5414, 6.6231E-11, -3.0392E-11, 4485.998, 6.2287E-11,
    -2.5141E-11, 4153.9743, 5.8343E-11, -1.9889E-11, 3824.8303, 5.4399E-11,
    -1.4637E-11, 3498.9258, 5.0455E-11, -4.3005E-11, 1.1202E-10, 8783.2944,
    -4.2166E-11, 1.0956E-10, 8569.657, -4.1327E-11, 1.071E-10, 8356.0195,
    -4.0488E-11, 1.0464E-10, 8142.3821, -3.9649E-11, 1.0219E-10, 7928.7446,
    -3.881E-11, 9.9727E-11, 7715.1072, -3.7971E-11, 9.7269E-11, 7501.4697,
    -3.7132E-11, 9.481E-11, 7287.8323, -3.6292E-11, 9.2352E-11, 7074.1948,
    -3.5453E-11, 8.9894E-11, 26860.5573, -3.4614E-11, 8.7436E-11, 6518.455,
    -3.3775E-11, 8.4978E-11, 6176.7125, -3.2936E-11, 8.2519E-11, 5835.69,
    -3.2097E-11, 8.0061E-11, 5495.7473, -3.1258E-11, 7.7603E-11, 5157.2445,
    -3.0419E-11, 7.5145E-11, 4820.5414, -2.958E-11, 7.2687E-11, 4485.998,
    -2.874E-11, 7.0228E-11, 4153.9743, -2.7901E-11, 6.777E-11, 3824.8303,
    -2.7062E-11, 6.5312E-11, 3498.9258, 8315.5094, -7.2873E-11, -4.8846E-11,
    8115.8367, -7.1655E-11, -4.7264E-11, 7916.164, -7.0436E-11, -4.5683E-11,
    7716.4912, -6.9217E-11, -4.4102E-11, 7516.8185, -6.7998E-11, -4.2521E-11,
    7317.1458, -6.6779E-11, -4.0939E-11, 7117.4731, -6.556E-11, -3.9358E-11,
    6917.8004, -6.4341E-11, -3.7777E-11, 6718.1277, -6.3122E-11, -3.6196E-11,
    6518.455, -6.1903E-11, -3.4614E-11, 26318.7822, -6.0684E-11, -3.3033E-11,
    5990.6446, -5.9465E-11, -3.1452E-11, 5662.8669, -5.8246E-11, -2.987E-11,
    5335.8091, -5.7027E-11, -2.8289E-11, 5009.8312, -5.5808E-11, -2.6708E-11,
    4685.2931, -5.4589E-11, -2.5127E-11, 4362.5547, -5.337E-11, -2.3545E-11,
    4041.9761, -5.2151E-11, -2.1964E-11, 3723.9172, -5.0933E-11, -2.0383E-11,
    3408.7379, -4.9714E-11, -1.8802E-11, -1.121E-10, 8315.5094, 1.2201E-10,
    -1.0696E-10, 8115.8367, 1.1817E-10, -1.0181E-10, 7916.164, 1.1433E-10,
    -9.6673E-11, 7716.4912, 1.1049E-10, -9.1531E-11, 7516.8185, 1.0664E-10,
    -8.639E-11, 7317.1458, 1.028E-10, -8.1249E-11, 7117.4731, 9.8961E-11,
    -7.6108E-11, 6917.8004, 9.5119E-11, -7.0966E-11, 6718.1277, 9.1278E-11,
    -6.5825E-11, 6518.455, 8.7436E-11, -6.0684E-11, 26318.7822, 8.3594E-11,
    -5.5543E-11, 5990.6446, 7.9752E-11, -5.0402E-11, 5662.8669, 7.5911E-11,
    -4.526E-11, 5335.8091, 7.2069E-11, -4.0119E-11, 5009.8312, 6.8227E-11,
    -3.4978E-11, 4685.2931, 6.4386E-11, -2.9837E-11, 4362.5547, 6.0544E-11,
    -2.4696E-11, 4041.9761, 5.6702E-11, -1.9554E-11, 3723.9172, 5.286E-11,
    -1.4413E-11, 3408.7379, 4.9019E-11, -4.0994E-11, 1.0715E-10, 8315.5094,
    -4.0198E-11, 1.048E-10, 8115.8367, -3.9402E-11, 1.0244E-10, 7916.164,
    -3.8606E-11, 1.0009E-10, 7716.4912, -3.781E-11, 9.773E-11, 7516.8185,
    -3.7014E-11, 9.5374E-11, 7317.1458, -3.6217E-11, 9.3018E-11, 7117.4731,
    -3.5421E-11, 9.0662E-11, 6917.8004, -3.4625E-11, 8.8306E-11, 6718.1277,
    -3.3829E-11, 8.595E-11, 6518.455, -3.3033E-11, 8.3594E-11, 26318.7822,
    -3.2237E-11, 8.1238E-11, 5990.6446, -3.1441E-11, 7.8882E-11, 5662.8669,
    -3.0645E-11, 7.6526E-11, 5335.8091, -2.9849E-11, 7.417E-11, 5009.8312,
    -2.9052E-11, 7.1814E-11, 4685.2931, -2.8256E-11, 6.9458E-11, 4362.5547,
    -2.746E-11, 6.7102E-11, 4041.9761, -2.6664E-11, 6.4746E-11, 3723.9172,
    -2.5868E-11, 6.239E-11, 3408.7379, 7851.3239, -6.6627E-11, -4.762E-11,
    7665.2559, -6.5518E-11, -4.6081E-11, 7479.188, -6.441E-11, -4.4543E-11,
    7293.1201, -6.3301E-11, -4.3005E-11, 7107.0522, -6.2193E-11, -4.1467E-11,
    6920.9842, -6.1085E-11, -3.9928E-11, 6734.9163, -5.9976E-11, -3.839E-11,
    6548.8484, -5.8868E-11, -3.6852E-11, 6362.7804, -5.776E-11, -3.5313E-11,
    6176.7125, -5.6651E-11, -3.3775E-11, 5990.6446, -5.5543E-11, -3.2237E-11,
    25804.5767, -5.4434E-11, -3.0699E-11, 5490.0438, -5.3326E-11, -2.916E-11,
    5175.8709, -5.2218E-11, -2.7622E-11, 4862.4179, -5.1109E-11, -2.6084E-11,
    4550.0448, -5.0001E-11, -2.4546E-11, 4239.1114, -4.8893E-11, -2.3007E-11,
    3929.9779, -4.7784E-11, -2.1469E-11, 3623.0041, -4.6676E-11, -1.9931E-11,
    3318.5499, -4.5567E-11, -1.8392E-11, -1.0977E-10, 7851.3239, 1.1863E-10,
    -1.0474E-10, 7665.2559, 1.1489E-10, -9.971E-11, 7479.188, 1.1115E-10,
    -9.468E-11, 7293.1201, 1.0741E-10, -8.9649E-11, 7107.0522, 1.0368E-10,
    -8.4618E-11, 6920.9842, 9.9936E-11, -7.9588E-11, 6734.9163, 9.6196E-11,
    -7.4557E-11, 6548.8484, 9.2457E-11, -6.9526E-11, 6362.7804, 8.8717E-11,
    -6.4496E-11, 6176.7125, 8.4978E-11, -5.9465E-11, 5990.6446, 8.1238E-11,
    -5.4434E-11, 25804.5767, 7.7499E-11, -4.9404E-11, 5490.0438, 7.3759E-11,
    -4.4373E-11, 5175.8709, 7.0019E-11, -3.9342E-11, 4862.4179, 6.628E-11,
    -3.4312E-11, 4550.0448, 6.254E-11, -2.9281E-11, 4239.1114, 5.8801E-11,
    -2.4251E-11, 3929.9779, 5.5061E-11, -1.922E-11, 3623.0041, 5.1322E-11,
    -1.4189E-11, 3318.5499, 4.7582E-11, -3.8983E-11, 1.0229E-10, 7851.3239,
    -3.823E-11, 1.0004E-10, 7665.2559, -3.7477E-11, 9.7783E-11, 7479.188,
    -3.6723E-11, 9.5529E-11, 7293.1201, -3.597E-11, 9.3276E-11, 7107.0522,
    -3.5217E-11, 9.1022E-11, 6920.9842, -3.4464E-11, 8.8768E-11, 6734.9163,
    -3.3711E-11, 8.6514E-11, 6548.8484, -3.2958E-11, 8.426E-11, 6362.7804,
    -3.2205E-11, 8.2006E-11, 6176.7125, -3.1452E-11, 7.9752E-11, 5990.6446,
    -3.0699E-11, 7.7499E-11, 25804.5767, -2.9945E-11, 7.5245E-11, 5490.0438,
    -2.9192E-11, 7.2991E-11, 5175.8709, -2.8439E-11, 7.0737E-11, 4862.4179,
    -2.7686E-11, 6.8483E-11, 4550.0448, -2.6933E-11, 6.6229E-11, 4239.1114,
    -2.618E-11, 6.3975E-11, 3929.9779, -2.5427E-11, 6.1721E-11, 3623.0041,
    -2.4674E-11, 5.9468E-11, 3318.5499, 7391.0978, -6.038E-11, -4.6393E-11,
    7218.2747, -5.9382E-11, -4.4898E-11, 7045.4516, -5.8384E-11, -4.3403E-11,
    6872.6285, -5.7386E-11, -4.1908E-11, 6699.8055, -5.6389E-11, -4.0412E-11,
    6526.9824, -5.5391E-11, -3.8917E-11, 6354.1593, -5.4393E-11, -3.7422E-11,
    6181.3362, -5.3395E-11, -3.5927E-11, 6008.5131, -5.2397E-11, -3.4431E-11,
    5835.69, -5.1399E-11, -3.2936E-11, 5662.8669, -5.0402E-11, -3.1441E-11,
    5490.0438, -4.9404E-11, -2.9945E-11, 25317.2207, -4.8406E-11, -2.845E-11,
    5015.9327, -4.7408E-11, -2.6955E-11, 4715.0046, -4.641E-11, -2.546E-11,
    4414.7964, -4.5413E-11, -2.3964E-11, 4115.6681, -4.4415E-11, -2.2469E-11,
    3817.9796, -4.3417E-11, -2.0974E-11, 3522.0909, -4.2419E-11, -1.9479E-11,
    3228.3619, -4.1421E-11, -1.7983E-11, -1.0745E-10, 7391.0978, 1.1526E-10,
    -1.0253E-10, 7218.2747, 1.1162E-10, -9.7607E-11, 7045.4516, 1.0798E-10,
    -9.2687E-11, 6872.6285, 1.0434E-10, -8.7767E-11, 6699.8055, 1.0071E-10,
    -8.2847E-11, 6526.9824, 9.7069E-11, -7.7927E-11, 6354.1593, 9.3432E-11,
    -7.3006E-11, 6181.3362, 8.9794E-11, -6.8086E-11, 6008.5131, 8.6157E-11,
    -6.3166E-11, 5835.69, 8.2519E-11, -5.8246E-11, 5662.8669, 7.8882E-11,
    -5.3326E-11, 5490.0438, 7.5245E-11, -4.8406E-11, 25317.2207, 7.1607E-11,
    -4.3486E-11, 5015.9327, 6.797E-11, -3.8566E-11, 4715.0046, 6.4333E-11,
    -3.3646E-11, 4414.7964, 6.0695E-11, -2.8726E-11, 4115.6681, 5.7058E-11,
    -2.3806E-11, 3817.9796, 5.342E-11, -1.8885E-11, 3522.0909, 4.9783E-11,
    -1.3965E-11, 3228.3619, 4.6146E-11, -3.6972E-11, 9.7428E-11, 7391.0978,
    -3.6261E-11, 9.5276E-11, 7218.2747, -3.5551E-11, 9.3124E-11, 7045.4516,
    -3.4841E-11, 9.0973E-11, 6872.6285, -3.4131E-11, 8.8821E-11, 6699.8055,
    -3.3421E-11, 8.6669E-11, 6526.9824, -3.2711E-11, 8.4517E-11, 6354.1593,
    -3.2001E-11, 8.2366E-11, 6181.3362, -3.1291E-11, 8.0214E-11, 6008.5131,
    -3.0581E-11, 7.8062E-11, 5835.69, -2.987E-11, 7.5911E-11, 5662.8669,
    -2.916E-11, 7.3759E-11, 5490.0438, -2.845E-11, 7.1607E-11, 25317.2207,
    -2.774E-11, 6.9456E-11, 5015.9327, -2.703E-11, 6.7304E-11, 4715.0046,
    -2.632E-11, 6.5152E-11, 4414.7964, -2.561E-11, 6.3E-11, 4115.6681, -2.49E-11,
    6.0849E-11, 3817.9796, -2.419E-11, 5.8697E-11, 3522.0909, -2.3479E-11,
    5.6545E-11, 3228.3619, 6935.1912, -5.4133E-11, -4.5167E-11, 6775.253,
    -5.3246E-11, -4.3715E-11, 6615.3148, -5.2358E-11, -4.2263E-11, 6455.3766,
    -5.1471E-11, -4.0811E-11, 6295.4384, -5.0584E-11, -3.9358E-11, 6135.5002,
    -4.9697E-11, -3.7906E-11, 5975.562, -4.8809E-11, -3.6454E-11, 5815.6237,
    -4.7922E-11, -3.5001E-11, 5655.6855, -4.7035E-11, -3.3549E-11, 5495.7473,
    -4.6148E-11, -3.2097E-11, 5335.8091, -4.526E-11, -3.0645E-11, 5175.8709,
    -4.4373E-11, -2.9192E-11, 5015.9327, -4.3486E-11, -2.774E-11, 24855.9945,
    -4.2599E-11, -2.6288E-11, 4567.5913, -4.1711E-11, -2.4836E-11, 4279.5481,
    -4.0824E-11, -2.3383E-11, 3992.2248, -3.9937E-11, -2.1931E-11, 3705.9814,
    -3.905E-11, -2.0479E-11, 3421.1778, -3.8162E-11, -1.9026E-11, 3138.174,
    -3.7275E-11, -1.7574E-11, -1.0512E-10, 6935.1912, 1.1188E-10, -1.0031E-10,
    6775.253, 1.0834E-10, -9.5503E-11, 6615.3148, 1.0481E-10, -9.0694E-11,
    6455.3766, 1.0127E-10, -8.5884E-11, 6295.4384, 9.7737E-11, -8.1075E-11,
    6135.5002, 9.4202E-11, -7.6265E-11, 5975.562, 9.0667E-11, -7.1456E-11,
    5815.6237, 8.7132E-11, -6.6646E-11, 5655.6855, 8.3596E-11, -6.1837E-11,
    5495.7473, 8.0061E-11, -5.7027E-11, 5335.8091, 7.6526E-11, -5.2218E-11,
    5175.8709, 7.2991E-11, -4.7408E-11, 5015.9327, 6.9456E-11, -4.2599E-11,
    24855.9945, 6.592E-11, -3.7789E-11, 4567.5913, 6.2385E-11, -3.298E-11,
    4279.5481, 5.885E-11, -2.817E-11, 3992.2248, 5.5315E-11, -2.3361E-11,
    3705.9814, 5.1779E-11, -1.8551E-11, 3421.1778, 4.8244E-11, -1.3741E-11,
    3138.174, 4.4709E-11, -3.496E-11, 9.2564E-11, 6935.1912, -3.4293E-11,
    9.0515E-11, 6775.253, -3.3626E-11, 8.8465E-11, 6615.3148, -3.2959E-11,
    8.6416E-11, 6455.3766, -3.2292E-11, 8.4366E-11, 6295.4384, -3.1625E-11,
    8.2317E-11, 6135.5002, -3.0958E-11, 8.0267E-11, 5975.562, -3.029E-11,
    7.8218E-11, 5815.6237, -2.9623E-11, 7.6168E-11, 5655.6855, -2.8956E-11,
    7.4118E-11, 5495.7473, -2.8289E-11, 7.2069E-11, 5335.8091, -2.7622E-11,
    7.0019E-11, 5175.8709, -2.6955E-11, 6.797E-11, 5015.9327, -2.6288E-11,
    6.592E-11, 24855.9945, -2.5621E-11, 6.3871E-11, 4567.5913, -2.4954E-11,
    6.1821E-11, 4279.5481, -2.4287E-11, 5.9772E-11, 3992.2248, -2.3619E-11,
    5.7722E-11, 3705.9814, -2.2952E-11, 5.5673E-11, 3421.1778, -2.2285E-11,
    5.3623E-11, 3138.174, 6483.9641, -4.7886E-11, -4.3941E-11, 6336.5508,
    -4.7109E-11, -4.2532E-11, 6189.1375, -4.6333E-11, -4.1123E-11, 6041.7242,
    -4.5556E-11, -3.9713E-11, 5894.3109, -4.4779E-11, -3.8304E-11, 5746.8976,
    -4.4003E-11, -3.6895E-11, 5599.4843, -4.3226E-11, -3.5486E-11, 5452.071,
    -4.2449E-11, -3.4076E-11, 5304.6578, -4.1673E-11, -3.2667E-11, 5157.2445,
    -4.0896E-11, -3.1258E-11, 5009.8312, -4.0119E-11, -2.9849E-11, 4862.4179,
    -3.9342E-11, -2.8439E-11, 4715.0046, -3.8566E-11, -2.703E-11, 4567.5913,
    -3.7789E-11, -2.5621E-11, 24420.178, -3.7012E-11, -2.4211E-11, 4144.2998,
    -3.6236E-11, -2.2802E-11, 3868.7815, -3.5459E-11, -2.1393E-11, 3593.9832,
    -3.4682E-11, -1.9984E-11, 3320.2647, -3.3906E-11, -1.8574E-11, 3047.986,
    -3.3129E-11, -1.7165E-11, -1.028E-10, 6483.9641, 1.085E-10, -9.8099E-11,
    6336.5508, 1.0507E-10, -9.34E-11, 6189.1375, 1.0163E-10, -8.8701E-11,
    6041.7242, 9.8201E-11, -8.4002E-11, 5894.3109, 9.4768E-11, -7.9303E-11,
    5746.8976, 9.1335E-11, -7.4604E-11, 5599.4843, 8.7902E-11, -6.9905E-11,
    5452.071, 8.4469E-11, -6.5206E-11, 5304.6578, 8.1036E-11, -6.0507E-11,
    5157.2445, 7.7603E-11, -5.5808E-11, 5009.8312, 7.417E-11, -5.1109E-11,
    4862.4179, 7.0737E-11, -4.641E-11, 4715.0046, 6.7304E-11, -4.1711E-11,
    4567.5913, 6.3871E-11, -3.7012E-11, 24420.178, 6.0438E-11, -3.2313E-11,
    4144.2998, 5.7005E-11, -2.7614E-11, 3868.7815, 5.3572E-11, -2.2916E-11,
    3593.9832, 5.0139E-11, -1.8217E-11, 3320.2647, 4.6706E-11, -1.3518E-11,
    3047.986, 4.3272E-11, -3.2949E-11, 8.7701E-11, 6483.9641, -3.2325E-11,
    8.5754E-11, 6336.5508, -3.1701E-11, 8.3806E-11, 6189.1375, -3.1077E-11,
    8.1859E-11, 6041.7242, -3.0453E-11, 7.9911E-11, 5894.3109, -2.9828E-11,
    7.7964E-11, 5746.8976, -2.9204E-11, 7.6017E-11, 5599.4843, -2.858E-11,
    7.4069E-11, 5452.071, -2.7956E-11, 7.2122E-11, 5304.6578, -2.7332E-11,
    7.0175E-11, 5157.2445, -2.6708E-11, 6.8227E-11, 5009.8312, -2.6084E-11,
    6.628E-11, 4862.4179, -2.546E-11, 6.4333E-11, 4715.0046, -2.4836E-11,
    6.2385E-11, 4567.5913, -2.4211E-11, 6.0438E-11, 24420.178, -2.3587E-11,
    5.849E-11, 4144.2998, -2.2963E-11, 5.6543E-11, 3868.7815, -2.2339E-11,
    5.4596E-11, 3593.9832, -2.1715E-11, 5.2648E-11, 3320.2647, -2.1091E-11,
    5.0701E-11, 3047.986, 6037.7762, -4.1639E-11, -4.2715E-11, 5902.5279,
    -4.0973E-11, -4.1349E-11, 5767.2796, -4.0307E-11, -3.9983E-11, 5632.0313,
    -3.9641E-11, -3.8616E-11, 5496.783, -3.8975E-11, -3.725E-11, 5361.5346,
    -3.8309E-11, -3.5884E-11, 5226.2863, -3.7642E-11, -3.4517E-11, 5091.038,
    -3.6976E-11, -3.3151E-11, 4955.7897, -3.631E-11, -3.1785E-11, 4820.5414,
    -3.5644E-11, -3.0419E-11, 4685.2931, -3.4978E-11, -2.9052E-11, 4550.0448,
    -3.4312E-11, -2.7686E-11, 4414.7964, -3.3646E-11, -2.632E-11, 4279.5481,
    -3.298E-11, -2.4954E-11, 4144.2998, -3.2313E-11, -2.3587E-11, 24009.0515,
    -3.1647E-11, -2.2221E-11, 3745.3382, -3.0981E-11, -2.0855E-11, 3481.985,
    -3.0315E-11, -1.9489E-11, 3219.3516, -2.9649E-11, -1.8122E-11, 2957.798,
    -2.8983E-11, -1.6756E-11, -1.0047E-10, 6037.7762, 1.0512E-10, -9.5885E-11,
    5902.5279, 1.0179E-10, -9.1297E-11, 5767.2796, 9.8461E-11, -8.6708E-11,
    5632.0313, 9.513E-11, -8.212E-11, 5496.783, 9.1799E-11, -7.7531E-11,
    5361.5346, 8.8468E-11, -7.2943E-11, 5226.2863, 8.5138E-11, -6.8355E-11,
    5091.038, 8.1807E-11, -6.3766E-11, 4955.7897, 7.8476E-11, -5.9178E-11,
    4820.5414, 7.5145E-11, -5.4589E-11, 4685.2931, 7.1814E-11, -5.0001E-11,
    4550.0448, 6.8483E-11, -4.5413E-11, 4414.7964, 6.5152E-11, -4.0824E-11,
    4279.5481, 6.1821E-11, -3.6236E-11, 4144.2998, 5.849E-11, -3.1647E-11,
    24009.0515, 5.516E-11, -2.7059E-11, 3745.3382, 5.1829E-11, -2.2471E-11,
    3481.985, 4.8498E-11, -1.7882E-11, 3219.3516, 4.5167E-11, -1.3294E-11,
    2957.798, 4.1836E-11, -3.0938E-11, 8.2838E-11, 6037.7762, -3.0357E-11,
    8.0992E-11, 5902.5279, -2.9775E-11, 7.9147E-11, 5767.2796, -2.9194E-11,
    7.7302E-11, 5632.0313, -2.8613E-11, 7.5457E-11, 5496.783, -2.8032E-11,
    7.3612E-11, 5361.5346, -2.7451E-11, 7.1766E-11, 5226.2863, -2.687E-11,
    6.9921E-11, 5091.038, -2.6289E-11, 6.8076E-11, 4955.7897, -2.5708E-11,
    6.6231E-11, 4820.5414, -2.5127E-11, 6.4386E-11, 4685.2931, -2.4546E-11,
    6.254E-11, 4550.0448, -2.3964E-11, 6.0695E-11, 4414.7964, -2.3383E-11,
    5.885E-11, 4279.5481, -2.2802E-11, 5.7005E-11, 4144.2998, -2.2221E-11,
    5.516E-11, 24009.0515, -2.164E-11, 5.3314E-11, 3745.3382, -2.1059E-11,
    5.1469E-11, 3481.985, -2.0478E-11, 4.9624E-11, 3219.3516, -1.9897E-11,
    4.7779E-11, 2957.798, 5596.9877, -3.5392E-11, -4.1489E-11, 5473.5444,
    -3.4837E-11, -4.0166E-11, 5350.1011, -3.4281E-11, -3.8842E-11, 5226.6578,
    -3.3726E-11, -3.7519E-11, 5103.2145, -3.317E-11, -3.6196E-11, 4979.7712,
    -3.2615E-11, -3.4873E-11, 4856.3279, -3.2059E-11, -3.3549E-11, 4732.8846,
    -3.1503E-11, -3.2226E-11, 4609.4413, -3.0948E-11, -3.0903E-11, 4485.998,
    -3.0392E-11, -2.958E-11, 4362.5547, -2.9837E-11, -2.8256E-11, 4239.1114,
    -2.9281E-11, -2.6933E-11, 4115.6681, -2.8726E-11, -2.561E-11, 3992.2248,
    -2.817E-11, -2.4287E-11, 3868.7815, -2.7614E-11, -2.2963E-11, 3745.3382,
    -2.7059E-11, -2.164E-11, 23621.895, -2.6503E-11, -2.0317E-11, 3369.9867,
    -2.5948E-11, -1.8993E-11, 3118.4384, -2.5392E-11, -1.767E-11, 2867.6101,
    -2.4837E-11, -1.6347E-11, -9.8149E-11, 5596.9877, 1.0175E-10, -9.3671E-11,
    5473.5444, 9.8516E-11, -8.9193E-11, 5350.1011, 9.5288E-11, -8.4715E-11,
    5226.6578, 9.2059E-11, -8.0237E-11, 5103.2145, 8.883E-11, -7.576E-11,
    4979.7712, 8.5602E-11, -7.1282E-11, 4856.3279, 8.2373E-11, -6.6804E-11,
    4732.8846, 7.9144E-11, -6.2326E-11, 4609.4413, 7.5915E-11, -5.7848E-11,
    4485.998, 7.2687E-11, -5.337E-11, 4362.5547, 6.9458E-11, -4.8893E-11,
    4239.1114, 6.6229E-11, -4.4415E-11, 4115.6681, 6.3E-11, -3.9937E-11,
    3992.2248, 5.9772E-11, -3.5459E-11, 3868.7815, 5.6543E-11, -3.0981E-11,
    3745.3382, 5.3314E-11, -2.6503E-11, 23621.895, 5.0086E-11, -2.2026E-11,
    3369.9867, 4.6857E-11, -1.7548E-11, 3118.4384, 4.3628E-11, -1.307E-11,
    2867.6101, 4.0399E-11, -2.8926E-11, 7.7974E-11, 5596.9877, -2.8388E-11,
    7.6231E-11, 5473.5444, -2.785E-11, 7.4488E-11, 5350.1011, -2.7312E-11,
    7.2745E-11, 5226.6578, -2.6774E-11, 7.1002E-11, 5103.2145, -2.6236E-11,
    6.9259E-11, 4979.7712, -2.5698E-11, 6.7516E-11, 4856.3279, -2.516E-11,
    6.5773E-11, 4732.8846, -2.4622E-11, 6.403E-11, 4609.4413, -2.4083E-11,
    6.2287E-11, 4485.998, -2.3545E-11, 6.0544E-11, 4362.5547, -2.3007E-11,
    5.8801E-11, 4239.1114, -2.2469E-11, 5.7058E-11, 4115.6681, -2.1931E-11,
    5.5315E-11, 3992.2248, -2.1393E-11, 5.3572E-11, 3868.7815, -2.0855E-11,
    5.1829E-11, 3745.3382, -2.0317E-11, 5.0086E-11, 23621.895, -1.9779E-11,
    4.8343E-11, 3369.9867, -1.9241E-11, 4.6599E-11, 3118.4384, -1.8702E-11,
    4.4856E-11, 2867.6101, 5161.9584, -2.9146E-11, -4.0263E-11, 5049.9602,
    -2.8701E-11, -3.8983E-11, 4937.962, -2.8256E-11, -3.7702E-11, 4825.9637,
    -2.7811E-11, -3.6422E-11, 4713.9655, -2.7366E-11, -3.5142E-11, 4601.9673,
    -2.6921E-11, -3.3861E-11, 4489.969, -2.6476E-11, -3.2581E-11, 4377.9708,
    -2.6031E-11, -3.1301E-11, 4265.9726, -2.5586E-11, -3.0021E-11, 4153.9743,
    -2.5141E-11, -2.874E-11, 4041.9761, -2.4696E-11, -2.746E-11, 3929.9779,
    -2.4251E-11, -2.618E-11, 3817.9796, -2.3806E-11, -2.49E-11, 3705.9814,
    -2.3361E-11, -2.3619E-11, 3593.9832, -2.2916E-11, -2.2339E-11, 3481.985,
    -2.2471E-11, -2.1059E-11, 3369.9867, -2.2026E-11, -1.9779E-11, 23257.9885,
    -2.1581E-11, -1.8498E-11, 3017.5253, -2.1136E-11, -1.7218E-11, 2777.4221,
    -2.0691E-11, -1.5938E-11, -9.5824E-11, 5161.9584, 9.8367E-11, -9.1457E-11,
    5049.9602, 9.5241E-11, -8.709E-11, 4937.962, 9.2114E-11, -8.2722E-11,
    4825.9637, 8.8988E-11, -7.8355E-11, 4713.9655, 8.5861E-11, -7.3988E-11,
    4601.9673, 8.2735E-11, -6.9621E-11, 4489.969, 7.9608E-11, -6.5253E-11,
    4377.9708, 7.6482E-11, -6.0886E-11, 4265.9726, 7.3355E-11, -5.6519E-11,
    4153.9743, 7.0228E-11, -5.2151E-11, 4041.9761, 6.7102E-11, -4.7784E-11,
    3929.9779, 6.3975E-11, -4.3417E-11, 3817.9796, 6.0849E-11, -3.905E-11,
    3705.9814, 5.7722E-11, -3.4682E-11, 3593.9832, 5.4596E-11, -3.0315E-11,
    3481.985, 5.1469E-11, -2.5948E-11, 3369.9867, 4.8343E-11, -2.1581E-11,
    23257.9885, 4.5216E-11, -1.7213E-11, 3017.5253, 4.2089E-11, -1.2846E-11,
    2777.4221, 3.8963E-11, -2.6915E-11, 7.3111E-11, 5161.9584, -2.642E-11,
    7.147E-11, 5049.9602, -2.5925E-11, 6.9829E-11, 4937.962, -2.543E-11,
    6.8188E-11, 4825.9637, -2.4935E-11, 6.6547E-11, 4713.9655, -2.444E-11,
    6.4906E-11, 4601.9673, -2.3944E-11, 6.3266E-11, 4489.969, -2.3449E-11,
    6.1625E-11, 4377.9708, -2.2954E-11, 5.9984E-11, 4265.9726, -2.2459E-11,
    5.8343E-11, 4153.9743, -2.1964E-11, 5.6702E-11, 4041.9761, -2.1469E-11,
    5.5061E-11, 3929.9779, -2.0974E-11, 5.342E-11, 3817.9796, -2.0479E-11,
    5.1779E-11, 3705.9814, -1.9984E-11, 5.0139E-11, 3593.9832, -1.9489E-11,
    4.8498E-11, 3481.985, -1.8993E-11, 4.6857E-11, 3369.9867, -1.8498E-11,
    4.5216E-11, 23257.9885, -1.8003E-11, 4.3575E-11, 3017.5253, -1.7508E-11,
    4.1934E-11, 2777.4221, 4733.0484, -2.2899E-11, -3.9037E-11, 4632.1353,
    -2.2564E-11, -3.7799E-11, 4531.2221, -2.223E-11, -3.6562E-11, 4430.309,
    -2.1895E-11, -3.5325E-11, 4329.3959, -2.1561E-11, -3.4088E-11, 4228.4828,
    -2.1227E-11, -3.285E-11, 4127.5697, -2.0892E-11, -3.1613E-11, 4026.6565,
    -2.0558E-11, -3.0376E-11, 3925.7434, -2.0223E-11, -2.9139E-11, 3824.8303,
    -1.9889E-11, -2.7901E-11, 3723.9172, -1.9554E-11, -2.6664E-11, 3623.0041,
    -1.922E-11, -2.5427E-11, 3522.0909, -1.8885E-11, -2.419E-11, 3421.1778,
    -1.8551E-11, -2.2952E-11, 3320.2647, -1.8217E-11, -2.1715E-11, 3219.3516,
    -1.7882E-11, -2.0478E-11, 3118.4384, -1.7548E-11, -1.9241E-11, 3017.5253,
    -1.7213E-11, -1.8003E-11, 22916.6122, -1.6879E-11, -1.6766E-11, 2687.2341,
    -1.6544E-11, -1.5529E-11, -9.35E-11, 4733.0484, 9.499E-11, -8.9243E-11,
    4632.1353, 9.1965E-11, -8.4986E-11, 4531.2221, 8.8941E-11, -8.073E-11,
    4430.309, 8.5917E-11, -7.6473E-11, 4329.3959, 8.2892E-11, -7.2216E-11,
    4228.4828, 7.9868E-11, -6.7959E-11, 4127.5697, 7.6843E-11, -6.3703E-11,
    4026.6565, 7.3819E-11, -5.9446E-11, 3925.7434, 7.0795E-11, -5.5189E-11,
    3824.8303, 6.777E-11, -5.0933E-11, 3723.9172, 6.4746E-11, -4.6676E-11,
    3623.0041, 6.1721E-11, -4.2419E-11, 3522.0909, 5.8697E-11, -3.8162E-11,
    3421.1778, 5.5673E-11, -3.3906E-11, 3320.2647, 5.2648E-11, -2.9649E-11,
    3219.3516, 4.9624E-11, -2.5392E-11, 3118.4384, 4.6599E-11, -2.1136E-11,
    3017.5253, 4.3575E-11, -1.6879E-11, 22916.6122, 4.0551E-11, -1.2622E-11,
    2687.2341, 3.7526E-11, -2.4904E-11, 6.8247E-11, 4733.0484, -2.4452E-11,
    6.6709E-11, 4632.1353, -2.4E-11, 6.517E-11, 4531.2221, -2.3547E-11,
    6.3631E-11, 4430.309, -2.3095E-11, 6.2093E-11, 4329.3959, -2.2643E-11,
    6.0554E-11, 4228.4828, -2.2191E-11, 5.9015E-11, 4127.5697, -2.1739E-11,
    5.7477E-11, 4026.6565, -2.1287E-11, 5.5938E-11, 3925.7434, -2.0835E-11,
    5.4399E-11, 3824.8303, -2.0383E-11, 5.286E-11, 3723.9172, -1.9931E-11,
    5.1322E-11, 3623.0041, -1.9479E-11, 4.9783E-11, 3522.0909, -1.9026E-11,
    4.8244E-11, 3421.1778, -1.8574E-11, 4.6706E-11, 3320.2647, -1.8122E-11,
    4.5167E-11, 3219.3516, -1.767E-11, 4.3628E-11, 3118.4384, -1.7218E-11,
    4.2089E-11, 3017.5253, -1.6766E-11, 4.0551E-11, 22916.6122, -1.6314E-11,
    3.9012E-11, 2687.2341, 4310.6175, -1.6652E-11, -3.7811E-11, 4220.4295,
    -1.6428E-11, -3.6616E-11, 4130.2416, -1.6204E-11, -3.5422E-11, 4040.0536,
    -1.598E-11, -3.4228E-11, 3949.8656, -1.5756E-11, -3.3033E-11, 3859.6777,
    -1.5532E-11, -3.1839E-11, 3769.4897, -1.5309E-11, -3.0645E-11, 3679.3018,
    -1.5085E-11, -2.9451E-11, 3589.1138, -1.4861E-11, -2.8256E-11, 3498.9258,
    -1.4637E-11, -2.7062E-11, 3408.7379, -1.4413E-11, -2.5868E-11, 3318.5499,
    -1.4189E-11, -2.4674E-11, 3228.3619, -1.3965E-11, -2.3479E-11, 3138.174,
    -1.3741E-11, -2.2285E-11, 3047.986, -1.3518E-11, -2.1091E-11, 2957.798,
    -1.3294E-11, -1.9897E-11, 2867.6101, -1.307E-11, -1.8702E-11, 2777.4221,
    -1.2846E-11, -1.7508E-11, 2687.2341, -1.2622E-11, -1.6314E-11, 22597.0462,
    -1.2398E-11, -1.512E-11, -9.1175E-11, 4310.6175, 9.1612E-11, -8.7029E-11,
    4220.4295, 8.869E-11, -8.2883E-11, 4130.2416, 8.5768E-11, -7.8737E-11,
    4040.0536, 8.2845E-11, -7.4591E-11, 3949.8656, 7.9923E-11, -7.0444E-11,
    3859.6777, 7.7001E-11, -6.6298E-11, 3769.4897, 7.4079E-11, -6.2152E-11,
    3679.3018, 7.1157E-11, -5.8006E-11, 3589.1138, 6.8234E-11, -5.386E-11,
    3498.9258, 6.5312E-11, -4.9714E-11, 3408.7379, 6.239E-11, -4.5567E-11,
    3318.5499, 5.9468E-11, -4.1421E-11, 3228.3619, 5.6545E-11, -3.7275E-11,
    3138.174, 5.3623E-11, -3.3129E-11, 3047.986, 5.0701E-11, -2.8983E-11,
    2957.798, 4.7779E-11, -2.4837E-11, 2867.6101, 4.4856E-11, -2.0691E-11,
    2777.4221, 4.1934E-11, -1.6544E-11, 2687.2341, 3.9012E-11, -1.2398E-11,
    22597.0462, 3.609E-11, -2.2892E-11, 6.3384E-11, 4310.6175, -2.2483E-11,
    6.1948E-11, 4220.4295, -2.2074E-11, 6.0511E-11, 4130.2416, -2.1665E-11,
    5.9074E-11, 4040.0536, -2.1256E-11, 5.7638E-11, 3949.8656, -2.0847E-11,
    5.6201E-11, 3859.6777, -2.0438E-11, 5.4765E-11, 3769.4897, -2.0029E-11,
    5.3328E-11, 3679.3018, -1.962E-11, 5.1892E-11, 3589.1138, -1.9211E-11,
    5.0455E-11, 3498.9258, -1.8802E-11, 4.9019E-11, 3408.7379, -1.8392E-11,
    4.7582E-11, 3318.5499, -1.7983E-11, 4.6146E-11, 3228.3619, -1.7574E-11,
    4.4709E-11, 3138.174, -1.7165E-11, 4.3272E-11, 3047.986, -1.6756E-11,
    4.1836E-11, 2957.798, -1.6347E-11, 4.0399E-11, 2867.6101, -1.5938E-11,
    3.8963E-11, 2777.4221, -1.5529E-11, 3.7526E-11, 2687.2341, -1.512E-11,
    3.609E-11, 22597.0462 };

  static const real_T b_a[3600] = { -33095.9468, 1.3534E-10, 6.1107E-11,
    -12611.9638, 1.3302E-10, 5.9096E-11, -12128.3408, 1.3069E-10, 5.7084E-11,
    -11645.4377, 1.2837E-10, 5.5073E-11, -11163.6144, 1.2604E-10, 5.3062E-11,
    -10683.231, 1.2372E-10, 5.1051E-11, -10204.6473, 1.2139E-10, 4.9039E-11,
    -9728.2234, 1.1907E-10, 4.7028E-11, -9254.3191, 1.1675E-10, 4.5017E-11,
    -8783.2944, 1.1442E-10, 4.3005E-11, -8315.5094, 1.121E-10, 4.0994E-11,
    -7851.3239, 1.0977E-10, 3.8983E-11, -7391.0978, 1.0745E-10, 3.6972E-11,
    -6935.1912, 1.0512E-10, 3.496E-11, -6483.9641, 1.028E-10, 3.2949E-11,
    -6037.7762, 1.0047E-10, 3.0938E-11, -5596.9877, 9.8149E-11, 2.8926E-11,
    -5161.9584, 9.5824E-11, 2.6915E-11, -4733.0484, 9.35E-11, 2.4904E-11,
    -4310.6175, 9.1175E-11, 2.2892E-11, 1.3534E-10, -33095.9468, -1.5579E-10,
    1.291E-10, -12611.9638, -1.5092E-10, 1.2285E-10, -12128.3408, -1.4606E-10,
    1.166E-10, -11645.4377, -1.412E-10, 1.1035E-10, -11163.6144, -1.3633E-10,
    1.0411E-10, -10683.231, -1.3147E-10, 9.7861E-11, -10204.6473, -1.2661E-10,
    9.1614E-11, -9728.2234, -1.2174E-10, 8.5367E-11, -9254.3191, -1.1688E-10,
    7.912E-11, -8783.2944, -1.1202E-10, 7.2873E-11, -8315.5094, -1.0715E-10,
    6.6627E-11, -7851.3239, -1.0229E-10, 6.038E-11, -7391.0978, -9.7428E-11,
    5.4133E-11, -6935.1912, -9.2564E-11, 4.7886E-11, -6483.9641, -8.7701E-11,
    4.1639E-11, -6037.7762, -8.2838E-11, 3.5392E-11, -5596.9877, -7.7974E-11,
    2.9146E-11, -5161.9584, -7.3111E-11, 2.2899E-11, -4733.0484, -6.8247E-11,
    1.6652E-11, -4310.6175, -6.3384E-11, 6.1107E-11, -1.5579E-10, -33095.9468,
    5.9881E-11, -1.5241E-10, -12611.9638, 5.8655E-11, -1.4903E-10, -12128.3408,
    5.7429E-11, -1.4566E-10, -11645.4377, 5.6203E-11, -1.4228E-10, -11163.6144,
    5.4976E-11, -1.389E-10, -10683.231, 5.375E-11, -1.3552E-10, -10204.6473,
    5.2524E-11, -1.3214E-10, -9728.2234, 5.1298E-11, -1.2877E-10, -9254.3191,
    5.0072E-11, -1.2539E-10, -8783.2944, 4.8846E-11, -1.2201E-10, -8315.5094,
    4.762E-11, -1.1863E-10, -7851.3239, 4.6393E-11, -1.1526E-10, -7391.0978,
    4.5167E-11, -1.1188E-10, -6935.1912, 4.3941E-11, -1.085E-10, -6483.9641,
    4.2715E-11, -1.0512E-10, -6037.7762, 4.1489E-11, -1.0175E-10, -5596.9877,
    4.0263E-11, -9.8367E-11, -5161.9584, 3.9037E-11, -9.499E-11, -4733.0484,
    3.7811E-11, -9.1612E-11, -4310.6175, -12611.9638, 1.291E-10, 5.9881E-11,
    -32273.6501, 1.2688E-10, 5.7913E-11, -11806.8715, 1.2467E-10, 5.5944E-11,
    -11340.4527, 1.2245E-10, 5.3976E-11, -10874.754, 1.2024E-10, 5.2008E-11,
    -10410.135, 1.1803E-10, 5.0039E-11, -9946.9559, 1.1581E-10, 4.8071E-11,
    -9485.5766, 1.136E-10, 4.6103E-11, -9026.3569, 1.1138E-10, 4.4135E-11,
    -8569.657, 1.0917E-10, 4.2166E-11, -8115.8367, 1.0696E-10, 4.0198E-11,
    -7665.2559, 1.0474E-10, 3.823E-11, -7218.2747, 1.0253E-10, 3.6261E-11,
    -6775.253, 1.0031E-10, 3.4293E-11, -6336.5508, 9.8099E-11, 3.2325E-11,
    -5902.5279, 9.5885E-11, 3.0357E-11, -5473.5444, 9.3671E-11, 2.8388E-11,
    -5049.9602, 9.1457E-11, 2.642E-11, -4632.1353, 8.9243E-11, 2.4452E-11,
    -4220.4295, 8.7029E-11, 2.2483E-11, 1.3302E-10, -12611.9638, -1.5241E-10,
    1.2688E-10, -32273.6501, -1.4765E-10, 1.2074E-10, -11806.8715, -1.4289E-10,
    1.1461E-10, -11340.4527, -1.3813E-10, 1.0847E-10, -10874.754, -1.3337E-10,
    1.0234E-10, -10410.135, -1.286E-10, 9.62E-11, -9946.9559, -1.2384E-10,
    9.0063E-11, -9485.5766, -1.1908E-10, 8.3927E-11, -9026.3569, -1.1432E-10,
    7.7791E-11, -8569.657, -1.0956E-10, 7.1655E-11, -8115.8367, -1.048E-10,
    6.5518E-11, -7665.2559, -1.0004E-10, 5.9382E-11, -7218.2747, -9.5276E-11,
    5.3246E-11, -6775.253, -9.0515E-11, 4.7109E-11, -6336.5508, -8.5754E-11,
    4.0973E-11, -5902.5279, -8.0992E-11, 3.4837E-11, -5473.5444, -7.6231E-11,
    2.8701E-11, -5049.9602, -7.147E-11, 2.2564E-11, -4632.1353, -6.6709E-11,
    1.6428E-11, -4220.4295, -6.1948E-11, 5.9096E-11, -1.5092E-10, -12611.9638,
    5.7913E-11, -1.4765E-10, -32273.6501, 5.673E-11, -1.4437E-10, -11806.8715,
    5.5546E-11, -1.411E-10, -11340.4527, 5.4363E-11, -1.3782E-10, -10874.754,
    5.318E-11, -1.3455E-10, -10410.135, 5.1997E-11, -1.3127E-10, -9946.9559,
    5.0814E-11, -1.28E-10, -9485.5766, 4.9631E-11, -1.2472E-10, -9026.3569,
    4.8448E-11, -1.2145E-10, -8569.657, 4.7264E-11, -1.1817E-10, -8115.8367,
    4.6081E-11, -1.1489E-10, -7665.2559, 4.4898E-11, -1.1162E-10, -7218.2747,
    4.3715E-11, -1.0834E-10, -6775.253, 4.2532E-11, -1.0507E-10, -6336.5508,
    4.1349E-11, -1.0179E-10, -5902.5279, 4.0166E-11, -9.8516E-11, -5473.5444,
    3.8983E-11, -9.5241E-11, -5049.9602, 3.7799E-11, -9.1965E-11, -4632.1353,
    3.6616E-11, -8.869E-11, -4220.4295, -12128.3408, 1.2285E-10, 5.8655E-11,
    -11806.8715, 1.2074E-10, 5.673E-11, -31485.4021, 1.1864E-10, 5.4804E-11,
    -11035.4678, 1.1654E-10, 5.2879E-11, -10585.8935, 1.1443E-10, 5.0954E-11,
    -10137.0391, 1.1233E-10, 4.9028E-11, -9689.2645, 1.1023E-10, 4.7103E-11,
    -9242.9298, 1.0812E-10, 4.5178E-11, -8798.3948, 1.0602E-10, 4.3252E-11,
    -8356.0195, 1.0392E-10, 4.1327E-11, -7916.164, 1.0181E-10, 3.9402E-11,
    -7479.188, 9.971E-11, 3.7477E-11, -7045.4516, 9.7607E-11, 3.5551E-11,
    -6615.3148, 9.5503E-11, 3.3626E-11, -6189.1375, 9.34E-11, 3.1701E-11,
    -5767.2796, 9.1297E-11, 2.9775E-11, -5350.1011, 8.9193E-11, 2.785E-11,
    -4937.962, 8.709E-11, 2.5925E-11, -4531.2221, 8.4986E-11, 2.4E-11,
    -4130.2416, 8.2883E-11, 2.2074E-11, 1.3069E-10, -12128.3408, -1.4903E-10,
    1.2467E-10, -11806.8715, -1.4437E-10, 1.1864E-10, -31485.4021, -1.3971E-10,
    1.1262E-10, -11035.4678, -1.3506E-10, 1.0659E-10, -10585.8935, -1.304E-10,
    1.0056E-10, -10137.0391, -1.2574E-10, 9.4538E-11, -9689.2645, -1.2108E-10,
    8.8513E-11, -9242.9298, -1.1642E-10, 8.2487E-11, -8798.3948, -1.1176E-10,
    7.6461E-11, -8356.0195, -1.071E-10, 7.0436E-11, -7916.164, -1.0244E-10,
    6.441E-11, -7479.188, -9.7783E-11, 5.8384E-11, -7045.4516, -9.3124E-11,
    5.2358E-11, -6615.3148, -8.8465E-11, 4.6333E-11, -6189.1375, -8.3806E-11,
    4.0307E-11, -5767.2796, -7.9147E-11, 3.4281E-11, -5350.1011, -7.4488E-11,
    2.8256E-11, -4937.962, -6.9829E-11, 2.223E-11, -4531.2221, -6.517E-11,
    1.6204E-11, -4130.2416, -6.0511E-11, 5.7084E-11, -1.4606E-10, -12128.3408,
    5.5944E-11, -1.4289E-10, -11806.8715, 5.4804E-11, -1.3971E-10, -31485.4021,
    5.3664E-11, -1.3654E-10, -11035.4678, 5.2524E-11, -1.3337E-10, -10585.8935,
    5.1384E-11, -1.3019E-10, -10137.0391, 5.0244E-11, -1.2702E-10, -9689.2645,
    4.9104E-11, -1.2385E-10, -9242.9298, 4.7963E-11, -1.2067E-10, -8798.3948,
    4.6823E-11, -1.175E-10, -8356.0195, 4.5683E-11, -1.1433E-10, -7916.164,
    4.4543E-11, -1.1115E-10, -7479.188, 4.3403E-11, -1.0798E-10, -7045.4516,
    4.2263E-11, -1.0481E-10, -6615.3148, 4.1123E-11, -1.0163E-10, -6189.1375,
    3.9983E-11, -9.8461E-11, -5767.2796, 3.8842E-11, -9.5288E-11, -5350.1011,
    3.7702E-11, -9.2114E-11, -4937.962, 3.6562E-11, -8.8941E-11, -4531.2221,
    3.5422E-11, -8.5768E-11, -4130.2416, -11645.4377, 1.166E-10, 5.7429E-11,
    -11340.4527, 1.1461E-10, 5.5546E-11, -11035.4678, 1.1262E-10, 5.3664E-11,
    -30730.4829, 1.1062E-10, 5.1782E-11, -10297.033, 1.0863E-10, 4.99E-11,
    -9863.9431, 1.0664E-10, 4.8017E-11, -9431.5731, 1.0464E-10, 4.6135E-11,
    -9000.283, 1.0265E-10, 4.4253E-11, -8570.4326, 1.0066E-10, 4.237E-11,
    -8142.3821, 9.8665E-11, 4.0488E-11, -7716.4912, 9.6673E-11, 3.8606E-11,
    -7293.1201, 9.468E-11, 3.6723E-11, -6872.6285, 9.2687E-11, 3.4841E-11,
    -6455.3766, 9.0694E-11, 3.2959E-11, -6041.7242, 8.8701E-11, 3.1077E-11,
    -5632.0313, 8.6708E-11, 2.9194E-11, -5226.6578, 8.4715E-11, 2.7312E-11,
    -4825.9637, 8.2722E-11, 2.543E-11, -4430.309, 8.073E-11, 2.3547E-11,
    -4040.0536, 7.8737E-11, 2.1665E-11, 1.2837E-10, -11645.4377, -1.4566E-10,
    1.2245E-10, -11340.4527, -1.411E-10, 1.1654E-10, -11035.4678, -1.3654E-10,
    1.1062E-10, -30730.4829, -1.3198E-10, 1.0471E-10, -10297.033, -1.2743E-10,
    9.8792E-11, -9863.9431, -1.2287E-10, 9.2877E-11, -9431.5731, -1.1831E-10,
    8.6962E-11, -9000.283, -1.1376E-10, 8.1047E-11, -8570.4326, -1.092E-10,
    7.5132E-11, -8142.3821, -1.0464E-10, 6.9217E-11, -7716.4912, -1.0009E-10,
    6.3301E-11, -7293.1201, -9.5529E-11, 5.7386E-11, -6872.6285, -9.0973E-11,
    5.1471E-11, -6455.3766, -8.6416E-11, 4.5556E-11, -6041.7242, -8.1859E-11,
    3.9641E-11, -5632.0313, -7.7302E-11, 3.3726E-11, -5226.6578, -7.2745E-11,
    2.7811E-11, -4825.9637, -6.8188E-11, 2.1895E-11, -4430.309, -6.3631E-11,
    1.598E-11, -4040.0536, -5.9074E-11, 5.5073E-11, -1.412E-10, -11645.4377,
    5.3976E-11, -1.3813E-10, -11340.4527, 5.2879E-11, -1.3506E-10, -11035.4678,
    5.1782E-11, -1.3198E-10, -30730.4829, 5.0685E-11, -1.2891E-10, -10297.033,
    4.9588E-11, -1.2584E-10, -9863.9431, 4.849E-11, -1.2277E-10, -9431.5731,
    4.7393E-11, -1.197E-10, -9000.283, 4.6296E-11, -1.1663E-10, -8570.4326,
    4.5199E-11, -1.1356E-10, -8142.3821, 4.4102E-11, -1.1049E-10, -7716.4912,
    4.3005E-11, -1.0741E-10, -7293.1201, 4.1908E-11, -1.0434E-10, -6872.6285,
    4.0811E-11, -1.0127E-10, -6455.3766, 3.9713E-11, -9.8201E-11, -6041.7242,
    3.8616E-11, -9.513E-11, -5632.0313, 3.7519E-11, -9.2059E-11, -5226.6578,
    3.6422E-11, -8.8988E-11, -4825.9637, 3.5325E-11, -8.5917E-11, -4430.309,
    3.4228E-11, -8.2845E-11, -4040.0536, -11163.6144, 1.1035E-10, 5.6203E-11,
    -10874.754, 1.0847E-10, 5.4363E-11, -10585.8935, 1.0659E-10, 5.2524E-11,
    -10297.033, 1.0471E-10, 5.0685E-11, -30008.1726, 1.0283E-10, 4.8845E-11,
    -9590.8471, 1.0094E-10, 4.7006E-11, -9173.8817, 9.9061E-11, 4.5167E-11,
    -8757.6362, 9.7178E-11, 4.3328E-11, -8342.4705, 9.5296E-11, 4.1488E-11,
    -7928.7446, 9.3414E-11, 3.9649E-11, -7516.8185, 9.1531E-11, 3.781E-11,
    -7107.0522, 8.9649E-11, 3.597E-11, -6699.8055, 8.7767E-11, 3.4131E-11,
    -6295.4384, 8.5884E-11, 3.2292E-11, -5894.3109, 8.4002E-11, 3.0453E-11,
    -5496.783, 8.212E-11, 2.8613E-11, -5103.2145, 8.0237E-11, 2.6774E-11,
    -4713.9655, 7.8355E-11, 2.4935E-11, -4329.3959, 7.6473E-11, 2.3095E-11,
    -3949.8656, 7.4591E-11, 2.1256E-11, 1.2604E-10, -11163.6144, -1.4228E-10,
    1.2024E-10, -10874.754, -1.3782E-10, 1.1443E-10, -10585.8935, -1.3337E-10,
    1.0863E-10, -10297.033, -1.2891E-10, 1.0283E-10, -30008.1726, -1.2446E-10,
    9.7021E-11, -9590.8471, -1.2E-10, 9.1216E-11, -9173.8817, -1.1555E-10,
    8.5411E-11, -8757.6362, -1.1109E-10, 7.9607E-11, -8342.4705, -1.0664E-10,
    7.3802E-11, -7928.7446, -1.0219E-10, 6.7998E-11, -7516.8185, -9.773E-11,
    6.2193E-11, -7107.0522, -9.3276E-11, 5.6389E-11, -6699.8055, -8.8821E-11,
    5.0584E-11, -6295.4384, -8.4366E-11, 4.4779E-11, -5894.3109, -7.9911E-11,
    3.8975E-11, -5496.783, -7.5457E-11, 3.317E-11, -5103.2145, -7.1002E-11,
    2.7366E-11, -4713.9655, -6.6547E-11, 2.1561E-11, -4329.3959, -6.2093E-11,
    1.5756E-11, -3949.8656, -5.7638E-11, 5.3062E-11, -1.3633E-10, -11163.6144,
    5.2008E-11, -1.3337E-10, -10874.754, 5.0954E-11, -1.304E-10, -10585.8935,
    4.99E-11, -1.2743E-10, -10297.033, 4.8845E-11, -1.2446E-10, -30008.1726,
    4.7791E-11, -1.2149E-10, -9590.8471, 4.6737E-11, -1.1852E-10, -9173.8817,
    4.5683E-11, -1.1555E-10, -8757.6362, 4.4629E-11, -1.1258E-10, -8342.4705,
    4.3575E-11, -1.0961E-10, -7928.7446, 4.2521E-11, -1.0664E-10, -7516.8185,
    4.1467E-11, -1.0368E-10, -7107.0522, 4.0412E-11, -1.0071E-10, -6699.8055,
    3.9358E-11, -9.7737E-11, -6295.4384, 3.8304E-11, -9.4768E-11, -5894.3109,
    3.725E-11, -9.1799E-11, -5496.783, 3.6196E-11, -8.883E-11, -5103.2145,
    3.5142E-11, -8.5861E-11, -4713.9655, 3.4088E-11, -8.2892E-11, -4329.3959,
    3.3033E-11, -7.9923E-11, -3949.8656, -10683.231, 1.0411E-10, 5.4976E-11,
    -10410.135, 1.0234E-10, 5.318E-11, -10137.0391, 1.0056E-10, 5.1384E-11,
    -9863.9431, 9.8792E-11, 4.9588E-11, -9590.8471, 9.7021E-11, 4.7791E-11,
    -29317.7512, 9.5249E-11, 4.5995E-11, -8916.1903, 9.3477E-11, 4.4199E-11,
    -8514.9894, 9.1705E-11, 4.2402E-11, -8114.5083, 8.9934E-11, 4.0606E-11,
    -7715.1072, 8.8162E-11, 3.881E-11, -7317.1458, 8.639E-11, 3.7014E-11,
    -6920.9842, 8.4618E-11, 3.5217E-11, -6526.9824, 8.2847E-11, 3.3421E-11,
    -6135.5002, 8.1075E-11, 3.1625E-11, -5746.8976, 7.9303E-11, 2.9828E-11,
    -5361.5346, 7.7531E-11, 2.8032E-11, -4979.7712, 7.576E-11, 2.6236E-11,
    -4601.9673, 7.3988E-11, 2.444E-11, -4228.4828, 7.2216E-11, 2.2643E-11,
    -3859.6777, 7.0444E-11, 2.0847E-11, 1.2372E-10, -10683.231, -1.389E-10,
    1.1803E-10, -10410.135, -1.3455E-10, 1.1233E-10, -10137.0391, -1.3019E-10,
    1.0664E-10, -9863.9431, -1.2584E-10, 1.0094E-10, -9590.8471, -1.2149E-10,
    9.5249E-11, -29317.7512, -1.1714E-10, 8.9555E-11, -8916.1903, -1.1278E-10,
    8.3861E-11, -8514.9894, -1.0843E-10, 7.8167E-11, -8114.5083, -1.0408E-10,
    7.2473E-11, -7715.1072, -9.9727E-11, 6.6779E-11, -7317.1458, -9.5374E-11,
    6.1085E-11, -6920.9842, -9.1022E-11, 5.5391E-11, -6526.9824, -8.6669E-11,
    4.9697E-11, -6135.5002, -8.2317E-11, 4.4003E-11, -5746.8976, -7.7964E-11,
    3.8309E-11, -5361.5346, -7.3612E-11, 3.2615E-11, -4979.7712, -6.9259E-11,
    2.6921E-11, -4601.9673, -6.4906E-11, 2.1227E-11, -4228.4828, -6.0554E-11,
    1.5532E-11, -3859.6777, -5.6201E-11, 5.1051E-11, -1.3147E-10, -10683.231,
    5.0039E-11, -1.286E-10, -10410.135, 4.9028E-11, -1.2574E-10, -10137.0391,
    4.8017E-11, -1.2287E-10, -9863.9431, 4.7006E-11, -1.2E-10, -9590.8471,
    4.5995E-11, -1.1714E-10, -29317.7512, 4.4984E-11, -1.1427E-10, -8916.1903,
    4.3973E-11, -1.114E-10, -8514.9894, 4.2962E-11, -1.0854E-10, -8114.5083,
    4.195E-11, -1.0567E-10, -7715.1072, 4.0939E-11, -1.028E-10, -7317.1458,
    3.9928E-11, -9.9936E-11, -6920.9842, 3.8917E-11, -9.7069E-11, -6526.9824,
    3.7906E-11, -9.4202E-11, -6135.5002, 3.6895E-11, -9.1335E-11, -5746.8976,
    3.5884E-11, -8.8468E-11, -5361.5346, 3.4873E-11, -8.5602E-11, -4979.7712,
    3.3861E-11, -8.2735E-11, -4601.9673, 3.285E-11, -7.9868E-11, -4228.4828,
    3.1839E-11, -7.7001E-11, -3859.6777, -10204.6473, 9.7861E-11, 5.375E-11,
    -9946.9559, 9.62E-11, 5.1997E-11, -9689.2645, 9.4538E-11, 5.0244E-11,
    -9431.5731, 9.2877E-11, 4.849E-11, -9173.8817, 9.1216E-11, 4.6737E-11,
    -8916.1903, 8.9555E-11, 4.4984E-11, -28658.4989, 8.7894E-11, 4.3231E-11,
    -8272.3426, 8.6232E-11, 4.1477E-11, -7886.5462, 8.4571E-11, 3.9724E-11,
    -7501.4697, 8.291E-11, 3.7971E-11, -7117.4731, 8.1249E-11, 3.6217E-11,
    -6734.9163, 7.9588E-11, 3.4464E-11, -6354.1593, 7.7927E-11, 3.2711E-11,
    -5975.562, 7.6265E-11, 3.0958E-11, -5599.4843, 7.4604E-11, 2.9204E-11,
    -5226.2863, 7.2943E-11, 2.7451E-11, -4856.3279, 7.1282E-11, 2.5698E-11,
    -4489.969, 6.9621E-11, 2.3944E-11, -4127.5697, 6.7959E-11, 2.2191E-11,
    -3769.4897, 6.6298E-11, 2.0438E-11, 1.2139E-10, -10204.6473, -1.3552E-10,
    1.1581E-10, -9946.9559, -1.3127E-10, 1.1023E-10, -9689.2645, -1.2702E-10,
    1.0464E-10, -9431.5731, -1.2277E-10, 9.9061E-11, -9173.8817, -1.1852E-10,
    9.3477E-11, -8916.1903, -1.1427E-10, 8.7894E-11, -28658.4989, -1.1002E-10,
    8.231E-11, -8272.3426, -1.0577E-10, 7.6727E-11, -7886.5462, -1.0152E-10,
    7.1143E-11, -7501.4697, -9.7269E-11, 6.556E-11, -7117.4731, -9.3018E-11,
    5.9976E-11, -6734.9163, -8.8768E-11, 5.4393E-11, -6354.1593, -8.4517E-11,
    4.8809E-11, -5975.562, -8.0267E-11, 4.3226E-11, -5599.4843, -7.6017E-11,
    3.7642E-11, -5226.2863, -7.1766E-11, 3.2059E-11, -4856.3279, -6.7516E-11,
    2.6476E-11, -4489.969, -6.3266E-11, 2.0892E-11, -4127.5697, -5.9015E-11,
    1.5309E-11, -3769.4897, -5.4765E-11, 4.9039E-11, -1.2661E-10, -10204.6473,
    4.8071E-11, -1.2384E-10, -9946.9559, 4.7103E-11, -1.2108E-10, -9689.2645,
    4.6135E-11, -1.1831E-10, -9431.5731, 4.5167E-11, -1.1555E-10, -9173.8817,
    4.4199E-11, -1.1278E-10, -8916.1903, 4.3231E-11, -1.1002E-10, -28658.4989,
    4.2262E-11, -1.0726E-10, -8272.3426, 4.1294E-11, -1.0449E-10, -7886.5462,
    4.0326E-11, -1.0173E-10, -7501.4697, 3.9358E-11, -9.8961E-11, -7117.4731,
    3.839E-11, -9.6196E-11, -6734.9163, 3.7422E-11, -9.3432E-11, -6354.1593,
    3.6454E-11, -9.0667E-11, -5975.562, 3.5486E-11, -8.7902E-11, -5599.4843,
    3.4517E-11, -8.5138E-11, -5226.2863, 3.3549E-11, -8.2373E-11, -4856.3279,
    3.2581E-11, -7.9608E-11, -4489.969, 3.1613E-11, -7.6843E-11, -4127.5697,
    3.0645E-11, -7.4079E-11, -3769.4897, -9728.2234, 9.1614E-11, 5.2524E-11,
    -9485.5766, 9.0063E-11, 5.0814E-11, -9242.9298, 8.8513E-11, 4.9104E-11,
    -9000.283, 8.6962E-11, 4.7393E-11, -8757.6362, 8.5411E-11, 4.5683E-11,
    -8514.9894, 8.3861E-11, 4.3973E-11, -8272.3426, 8.231E-11, 4.2262E-11,
    -28029.6958, 8.076E-11, 4.0552E-11, -7658.584, 7.9209E-11, 3.8842E-11,
    -7287.8323, 7.7658E-11, 3.7132E-11, -6917.8004, 7.6108E-11, 3.5421E-11,
    -6548.8484, 7.4557E-11, 3.3711E-11, -6181.3362, 7.3006E-11, 3.2001E-11,
    -5815.6237, 7.1456E-11, 3.029E-11, -5452.071, 6.9905E-11, 2.858E-11,
    -5091.038, 6.8355E-11, 2.687E-11, -4732.8846, 6.6804E-11, 2.516E-11,
    -4377.9708, 6.5253E-11, 2.3449E-11, -4026.6565, 6.3703E-11, 2.1739E-11,
    -3679.3018, 6.2152E-11, 2.0029E-11, 1.1907E-10, -9728.2234, -1.3214E-10,
    1.136E-10, -9485.5766, -1.28E-10, 1.0812E-10, -9242.9298, -1.2385E-10,
    1.0265E-10, -9000.283, -1.197E-10, 9.7178E-11, -8757.6362, -1.1555E-10,
    9.1705E-11, -8514.9894, -1.114E-10, 8.6232E-11, -8272.3426, -1.0726E-10,
    8.076E-11, -28029.6958, -1.0311E-10, 7.5287E-11, -7658.584, -9.8959E-11,
    6.9814E-11, -7287.8323, -9.481E-11, 6.4341E-11, -6917.8004, -9.0662E-11,
    5.8868E-11, -6548.8484, -8.6514E-11, 5.3395E-11, -6181.3362, -8.2366E-11,
    4.7922E-11, -5815.6237, -7.8218E-11, 4.2449E-11, -5452.071, -7.4069E-11,
    3.6976E-11, -5091.038, -6.9921E-11, 3.1503E-11, -4732.8846, -6.5773E-11,
    2.6031E-11, -4377.9708, -6.1625E-11, 2.0558E-11, -4026.6565, -5.7477E-11,
    1.5085E-11, -3679.3018, -5.3328E-11, 4.7028E-11, -1.2174E-10, -9728.2234,
    4.6103E-11, -1.1908E-10, -9485.5766, 4.5178E-11, -1.1642E-10, -9242.9298,
    4.4253E-11, -1.1376E-10, -9000.283, 4.3328E-11, -1.1109E-10, -8757.6362,
    4.2402E-11, -1.0843E-10, -8514.9894, 4.1477E-11, -1.0577E-10, -8272.3426,
    4.0552E-11, -1.0311E-10, -28029.6958, 3.9627E-11, -1.0044E-10, -7658.584,
    3.8702E-11, -9.7782E-11, -7287.8323, 3.7777E-11, -9.5119E-11, -6917.8004,
    3.6852E-11, -9.2457E-11, -6548.8484, 3.5927E-11, -8.9794E-11, -6181.3362,
    3.5001E-11, -8.7132E-11, -5815.6237, 3.4076E-11, -8.4469E-11, -5452.071,
    3.3151E-11, -8.1807E-11, -5091.038, 3.2226E-11, -7.9144E-11, -4732.8846,
    3.1301E-11, -7.6482E-11, -4377.9708, 3.0376E-11, -7.3819E-11, -4026.6565,
    2.9451E-11, -7.1157E-11, -3679.3018, -9254.3191, 8.5367E-11, 5.1298E-11,
    -9026.3569, 8.3927E-11, 4.9631E-11, -8798.3948, 8.2487E-11, 4.7963E-11,
    -8570.4326, 8.1047E-11, 4.6296E-11, -8342.4705, 7.9607E-11, 4.4629E-11,
    -8114.5083, 7.8167E-11, 4.2962E-11, -7886.5462, 7.6727E-11, 4.1294E-11,
    -7658.584, 7.5287E-11, 3.9627E-11, -27430.6219, 7.3847E-11, 3.796E-11,
    -7074.1948, 7.2407E-11, 3.6292E-11, -6718.1277, 7.0966E-11, 3.4625E-11,
    -6362.7804, 6.9526E-11, 3.2958E-11, -6008.5131, 6.8086E-11, 3.1291E-11,
    -5655.6855, 6.6646E-11, 2.9623E-11, -5304.6578, 6.5206E-11, 2.7956E-11,
    -4955.7897, 6.3766E-11, 2.6289E-11, -4609.4413, 6.2326E-11, 2.4622E-11,
    -4265.9726, 6.0886E-11, 2.2954E-11, -3925.7434, 5.9446E-11, 2.1287E-11,
    -3589.1138, 5.8006E-11, 1.962E-11, 1.1675E-10, -9254.3191, -1.2877E-10,
    1.1138E-10, -9026.3569, -1.2472E-10, 1.0602E-10, -8798.3948, -1.2067E-10,
    1.0066E-10, -8570.4326, -1.1663E-10, 9.5296E-11, -8342.4705, -1.1258E-10,
    8.9934E-11, -8114.5083, -1.0854E-10, 8.4571E-11, -7886.5462, -1.0449E-10,
    7.9209E-11, -7658.584, -1.0044E-10, 7.3847E-11, -27430.6219, -9.6398E-11,
    6.8484E-11, -7074.1948, -9.2352E-11, 6.3122E-11, -6718.1277, -8.8306E-11,
    5.776E-11, -6362.7804, -8.426E-11, 5.2397E-11, -6008.5131, -8.0214E-11,
    4.7035E-11, -5655.6855, -7.6168E-11, 4.1673E-11, -5304.6578, -7.2122E-11,
    3.631E-11, -4955.7897, -6.8076E-11, 3.0948E-11, -4609.4413, -6.403E-11,
    2.5586E-11, -4265.9726, -5.9984E-11, 2.0223E-11, -3925.7434, -5.5938E-11,
    1.4861E-11, -3589.1138, -5.1892E-11, 4.5017E-11, -1.1688E-10, -9254.3191,
    4.4135E-11, -1.1432E-10, -9026.3569, 4.3252E-11, -1.1176E-10, -8798.3948,
    4.237E-11, -1.092E-10, -8570.4326, 4.1488E-11, -1.0664E-10, -8342.4705,
    4.0606E-11, -1.0408E-10, -8114.5083, 3.9724E-11, -1.0152E-10, -7886.5462,
    3.8842E-11, -9.8959E-11, -7658.584, 3.796E-11, -9.6398E-11, -27430.6219,
    3.7078E-11, -9.3838E-11, -7074.1948, 3.6196E-11, -9.1278E-11, -6718.1277,
    3.5313E-11, -8.8717E-11, -6362.7804, 3.4431E-11, -8.6157E-11, -6008.5131,
    3.3549E-11, -8.3596E-11, -5655.6855, 3.2667E-11, -8.1036E-11, -5304.6578,
    3.1785E-11, -7.8476E-11, -4955.7897, 3.0903E-11, -7.5915E-11, -4609.4413,
    3.0021E-11, -7.3355E-11, -4265.9726, 2.9139E-11, -7.0795E-11, -3925.7434,
    2.8256E-11, -6.8234E-11, -3589.1138, -8783.2944, 7.912E-11, 5.0072E-11,
    -8569.657, 7.7791E-11, 4.8448E-11, -8356.0195, 7.6461E-11, 4.6823E-11,
    -8142.3821, 7.5132E-11, 4.5199E-11, -7928.7446, 7.3802E-11, 4.3575E-11,
    -7715.1072, 7.2473E-11, 4.195E-11, -7501.4697, 7.1143E-11, 4.0326E-11,
    -7287.8323, 6.9814E-11, 3.8702E-11, -7074.1948, 6.8484E-11, 3.7078E-11,
    -26860.5573, 6.7155E-11, 3.5453E-11, -6518.455, 6.5825E-11, 3.3829E-11,
    -6176.7125, 6.4496E-11, 3.2205E-11, -5835.69, 6.3166E-11, 3.0581E-11,
    -5495.7473, 6.1837E-11, 2.8956E-11, -5157.2445, 6.0507E-11, 2.7332E-11,
    -4820.5414, 5.9178E-11, 2.5708E-11, -4485.998, 5.7848E-11, 2.4083E-11,
    -4153.9743, 5.6519E-11, 2.2459E-11, -3824.8303, 5.5189E-11, 2.0835E-11,
    -3498.9258, 5.386E-11, 1.9211E-11, 1.1442E-10, -8783.2944, -1.2539E-10,
    1.0917E-10, -8569.657, -1.2145E-10, 1.0392E-10, -8356.0195, -1.175E-10,
    9.8665E-11, -8142.3821, -1.1356E-10, 9.3414E-11, -7928.7446, -1.0961E-10,
    8.8162E-11, -7715.1072, -1.0567E-10, 8.291E-11, -7501.4697, -1.0173E-10,
    7.7658E-11, -7287.8323, -9.7782E-11, 7.2407E-11, -7074.1948, -9.3838E-11,
    6.7155E-11, -26860.5573, -8.9894E-11, 6.1903E-11, -6518.455, -8.595E-11,
    5.6651E-11, -6176.7125, -8.2006E-11, 5.1399E-11, -5835.69, -7.8062E-11,
    4.6148E-11, -5495.7473, -7.4118E-11, 4.0896E-11, -5157.2445, -7.0175E-11,
    3.5644E-11, -4820.5414, -6.6231E-11, 3.0392E-11, -4485.998, -6.2287E-11,
    2.5141E-11, -4153.9743, -5.8343E-11, 1.9889E-11, -3824.8303, -5.4399E-11,
    1.4637E-11, -3498.9258, -5.0455E-11, 4.3005E-11, -1.1202E-10, -8783.2944,
    4.2166E-11, -1.0956E-10, -8569.657, 4.1327E-11, -1.071E-10, -8356.0195,
    4.0488E-11, -1.0464E-10, -8142.3821, 3.9649E-11, -1.0219E-10, -7928.7446,
    3.881E-11, -9.9727E-11, -7715.1072, 3.7971E-11, -9.7269E-11, -7501.4697,
    3.7132E-11, -9.481E-11, -7287.8323, 3.6292E-11, -9.2352E-11, -7074.1948,
    3.5453E-11, -8.9894E-11, -26860.5573, 3.4614E-11, -8.7436E-11, -6518.455,
    3.3775E-11, -8.4978E-11, -6176.7125, 3.2936E-11, -8.2519E-11, -5835.69,
    3.2097E-11, -8.0061E-11, -5495.7473, 3.1258E-11, -7.7603E-11, -5157.2445,
    3.0419E-11, -7.5145E-11, -4820.5414, 2.958E-11, -7.2687E-11, -4485.998,
    2.874E-11, -7.0228E-11, -4153.9743, 2.7901E-11, -6.777E-11, -3824.8303,
    2.7062E-11, -6.5312E-11, -3498.9258, -8315.5094, 7.2873E-11, 4.8846E-11,
    -8115.8367, 7.1655E-11, 4.7264E-11, -7916.164, 7.0436E-11, 4.5683E-11,
    -7716.4912, 6.9217E-11, 4.4102E-11, -7516.8185, 6.7998E-11, 4.2521E-11,
    -7317.1458, 6.6779E-11, 4.0939E-11, -7117.4731, 6.556E-11, 3.9358E-11,
    -6917.8004, 6.4341E-11, 3.7777E-11, -6718.1277, 6.3122E-11, 3.6196E-11,
    -6518.455, 6.1903E-11, 3.4614E-11, -26318.7822, 6.0684E-11, 3.3033E-11,
    -5990.6446, 5.9465E-11, 3.1452E-11, -5662.8669, 5.8246E-11, 2.987E-11,
    -5335.8091, 5.7027E-11, 2.8289E-11, -5009.8312, 5.5808E-11, 2.6708E-11,
    -4685.2931, 5.4589E-11, 2.5127E-11, -4362.5547, 5.337E-11, 2.3545E-11,
    -4041.9761, 5.2151E-11, 2.1964E-11, -3723.9172, 5.0933E-11, 2.0383E-11,
    -3408.7379, 4.9714E-11, 1.8802E-11, 1.121E-10, -8315.5094, -1.2201E-10,
    1.0696E-10, -8115.8367, -1.1817E-10, 1.0181E-10, -7916.164, -1.1433E-10,
    9.6673E-11, -7716.4912, -1.1049E-10, 9.1531E-11, -7516.8185, -1.0664E-10,
    8.639E-11, -7317.1458, -1.028E-10, 8.1249E-11, -7117.4731, -9.8961E-11,
    7.6108E-11, -6917.8004, -9.5119E-11, 7.0966E-11, -6718.1277, -9.1278E-11,
    6.5825E-11, -6518.455, -8.7436E-11, 6.0684E-11, -26318.7822, -8.3594E-11,
    5.5543E-11, -5990.6446, -7.9752E-11, 5.0402E-11, -5662.8669, -7.5911E-11,
    4.526E-11, -5335.8091, -7.2069E-11, 4.0119E-11, -5009.8312, -6.8227E-11,
    3.4978E-11, -4685.2931, -6.4386E-11, 2.9837E-11, -4362.5547, -6.0544E-11,
    2.4696E-11, -4041.9761, -5.6702E-11, 1.9554E-11, -3723.9172, -5.286E-11,
    1.4413E-11, -3408.7379, -4.9019E-11, 4.0994E-11, -1.0715E-10, -8315.5094,
    4.0198E-11, -1.048E-10, -8115.8367, 3.9402E-11, -1.0244E-10, -7916.164,
    3.8606E-11, -1.0009E-10, -7716.4912, 3.781E-11, -9.773E-11, -7516.8185,
    3.7014E-11, -9.5374E-11, -7317.1458, 3.6217E-11, -9.3018E-11, -7117.4731,
    3.5421E-11, -9.0662E-11, -6917.8004, 3.4625E-11, -8.8306E-11, -6718.1277,
    3.3829E-11, -8.595E-11, -6518.455, 3.3033E-11, -8.3594E-11, -26318.7822,
    3.2237E-11, -8.1238E-11, -5990.6446, 3.1441E-11, -7.8882E-11, -5662.8669,
    3.0645E-11, -7.6526E-11, -5335.8091, 2.9849E-11, -7.417E-11, -5009.8312,
    2.9052E-11, -7.1814E-11, -4685.2931, 2.8256E-11, -6.9458E-11, -4362.5547,
    2.746E-11, -6.7102E-11, -4041.9761, 2.6664E-11, -6.4746E-11, -3723.9172,
    2.5868E-11, -6.239E-11, -3408.7379, -7851.3239, 6.6627E-11, 4.762E-11,
    -7665.2559, 6.5518E-11, 4.6081E-11, -7479.188, 6.441E-11, 4.4543E-11,
    -7293.1201, 6.3301E-11, 4.3005E-11, -7107.0522, 6.2193E-11, 4.1467E-11,
    -6920.9842, 6.1085E-11, 3.9928E-11, -6734.9163, 5.9976E-11, 3.839E-11,
    -6548.8484, 5.8868E-11, 3.6852E-11, -6362.7804, 5.776E-11, 3.5313E-11,
    -6176.7125, 5.6651E-11, 3.3775E-11, -5990.6446, 5.5543E-11, 3.2237E-11,
    -25804.5767, 5.4434E-11, 3.0699E-11, -5490.0438, 5.3326E-11, 2.916E-11,
    -5175.8709, 5.2218E-11, 2.7622E-11, -4862.4179, 5.1109E-11, 2.6084E-11,
    -4550.0448, 5.0001E-11, 2.4546E-11, -4239.1114, 4.8893E-11, 2.3007E-11,
    -3929.9779, 4.7784E-11, 2.1469E-11, -3623.0041, 4.6676E-11, 1.9931E-11,
    -3318.5499, 4.5567E-11, 1.8392E-11, 1.0977E-10, -7851.3239, -1.1863E-10,
    1.0474E-10, -7665.2559, -1.1489E-10, 9.971E-11, -7479.188, -1.1115E-10,
    9.468E-11, -7293.1201, -1.0741E-10, 8.9649E-11, -7107.0522, -1.0368E-10,
    8.4618E-11, -6920.9842, -9.9936E-11, 7.9588E-11, -6734.9163, -9.6196E-11,
    7.4557E-11, -6548.8484, -9.2457E-11, 6.9526E-11, -6362.7804, -8.8717E-11,
    6.4496E-11, -6176.7125, -8.4978E-11, 5.9465E-11, -5990.6446, -8.1238E-11,
    5.4434E-11, -25804.5767, -7.7499E-11, 4.9404E-11, -5490.0438, -7.3759E-11,
    4.4373E-11, -5175.8709, -7.0019E-11, 3.9342E-11, -4862.4179, -6.628E-11,
    3.4312E-11, -4550.0448, -6.254E-11, 2.9281E-11, -4239.1114, -5.8801E-11,
    2.4251E-11, -3929.9779, -5.5061E-11, 1.922E-11, -3623.0041, -5.1322E-11,
    1.4189E-11, -3318.5499, -4.7582E-11, 3.8983E-11, -1.0229E-10, -7851.3239,
    3.823E-11, -1.0004E-10, -7665.2559, 3.7477E-11, -9.7783E-11, -7479.188,
    3.6723E-11, -9.5529E-11, -7293.1201, 3.597E-11, -9.3276E-11, -7107.0522,
    3.5217E-11, -9.1022E-11, -6920.9842, 3.4464E-11, -8.8768E-11, -6734.9163,
    3.3711E-11, -8.6514E-11, -6548.8484, 3.2958E-11, -8.426E-11, -6362.7804,
    3.2205E-11, -8.2006E-11, -6176.7125, 3.1452E-11, -7.9752E-11, -5990.6446,
    3.0699E-11, -7.7499E-11, -25804.5767, 2.9945E-11, -7.5245E-11, -5490.0438,
    2.9192E-11, -7.2991E-11, -5175.8709, 2.8439E-11, -7.0737E-11, -4862.4179,
    2.7686E-11, -6.8483E-11, -4550.0448, 2.6933E-11, -6.6229E-11, -4239.1114,
    2.618E-11, -6.3975E-11, -3929.9779, 2.5427E-11, -6.1721E-11, -3623.0041,
    2.4674E-11, -5.9468E-11, -3318.5499, -7391.0978, 6.038E-11, 4.6393E-11,
    -7218.2747, 5.9382E-11, 4.4898E-11, -7045.4516, 5.8384E-11, 4.3403E-11,
    -6872.6285, 5.7386E-11, 4.1908E-11, -6699.8055, 5.6389E-11, 4.0412E-11,
    -6526.9824, 5.5391E-11, 3.8917E-11, -6354.1593, 5.4393E-11, 3.7422E-11,
    -6181.3362, 5.3395E-11, 3.5927E-11, -6008.5131, 5.2397E-11, 3.4431E-11,
    -5835.69, 5.1399E-11, 3.2936E-11, -5662.8669, 5.0402E-11, 3.1441E-11,
    -5490.0438, 4.9404E-11, 2.9945E-11, -25317.2207, 4.8406E-11, 2.845E-11,
    -5015.9327, 4.7408E-11, 2.6955E-11, -4715.0046, 4.641E-11, 2.546E-11,
    -4414.7964, 4.5413E-11, 2.3964E-11, -4115.6681, 4.4415E-11, 2.2469E-11,
    -3817.9796, 4.3417E-11, 2.0974E-11, -3522.0909, 4.2419E-11, 1.9479E-11,
    -3228.3619, 4.1421E-11, 1.7983E-11, 1.0745E-10, -7391.0978, -1.1526E-10,
    1.0253E-10, -7218.2747, -1.1162E-10, 9.7607E-11, -7045.4516, -1.0798E-10,
    9.2687E-11, -6872.6285, -1.0434E-10, 8.7767E-11, -6699.8055, -1.0071E-10,
    8.2847E-11, -6526.9824, -9.7069E-11, 7.7927E-11, -6354.1593, -9.3432E-11,
    7.3006E-11, -6181.3362, -8.9794E-11, 6.8086E-11, -6008.5131, -8.6157E-11,
    6.3166E-11, -5835.69, -8.2519E-11, 5.8246E-11, -5662.8669, -7.8882E-11,
    5.3326E-11, -5490.0438, -7.5245E-11, 4.8406E-11, -25317.2207, -7.1607E-11,
    4.3486E-11, -5015.9327, -6.797E-11, 3.8566E-11, -4715.0046, -6.4333E-11,
    3.3646E-11, -4414.7964, -6.0695E-11, 2.8726E-11, -4115.6681, -5.7058E-11,
    2.3806E-11, -3817.9796, -5.342E-11, 1.8885E-11, -3522.0909, -4.9783E-11,
    1.3965E-11, -3228.3619, -4.6146E-11, 3.6972E-11, -9.7428E-11, -7391.0978,
    3.6261E-11, -9.5276E-11, -7218.2747, 3.5551E-11, -9.3124E-11, -7045.4516,
    3.4841E-11, -9.0973E-11, -6872.6285, 3.4131E-11, -8.8821E-11, -6699.8055,
    3.3421E-11, -8.6669E-11, -6526.9824, 3.2711E-11, -8.4517E-11, -6354.1593,
    3.2001E-11, -8.2366E-11, -6181.3362, 3.1291E-11, -8.0214E-11, -6008.5131,
    3.0581E-11, -7.8062E-11, -5835.69, 2.987E-11, -7.5911E-11, -5662.8669,
    2.916E-11, -7.3759E-11, -5490.0438, 2.845E-11, -7.1607E-11, -25317.2207,
    2.774E-11, -6.9456E-11, -5015.9327, 2.703E-11, -6.7304E-11, -4715.0046,
    2.632E-11, -6.5152E-11, -4414.7964, 2.561E-11, -6.3E-11, -4115.6681,
    2.49E-11, -6.0849E-11, -3817.9796, 2.419E-11, -5.8697E-11, -3522.0909,
    2.3479E-11, -5.6545E-11, -3228.3619, -6935.1912, 5.4133E-11, 4.5167E-11,
    -6775.253, 5.3246E-11, 4.3715E-11, -6615.3148, 5.2358E-11, 4.2263E-11,
    -6455.3766, 5.1471E-11, 4.0811E-11, -6295.4384, 5.0584E-11, 3.9358E-11,
    -6135.5002, 4.9697E-11, 3.7906E-11, -5975.562, 4.8809E-11, 3.6454E-11,
    -5815.6237, 4.7922E-11, 3.5001E-11, -5655.6855, 4.7035E-11, 3.3549E-11,
    -5495.7473, 4.6148E-11, 3.2097E-11, -5335.8091, 4.526E-11, 3.0645E-11,
    -5175.8709, 4.4373E-11, 2.9192E-11, -5015.9327, 4.3486E-11, 2.774E-11,
    -24855.9945, 4.2599E-11, 2.6288E-11, -4567.5913, 4.1711E-11, 2.4836E-11,
    -4279.5481, 4.0824E-11, 2.3383E-11, -3992.2248, 3.9937E-11, 2.1931E-11,
    -3705.9814, 3.905E-11, 2.0479E-11, -3421.1778, 3.8162E-11, 1.9026E-11,
    -3138.174, 3.7275E-11, 1.7574E-11, 1.0512E-10, -6935.1912, -1.1188E-10,
    1.0031E-10, -6775.253, -1.0834E-10, 9.5503E-11, -6615.3148, -1.0481E-10,
    9.0694E-11, -6455.3766, -1.0127E-10, 8.5884E-11, -6295.4384, -9.7737E-11,
    8.1075E-11, -6135.5002, -9.4202E-11, 7.6265E-11, -5975.562, -9.0667E-11,
    7.1456E-11, -5815.6237, -8.7132E-11, 6.6646E-11, -5655.6855, -8.3596E-11,
    6.1837E-11, -5495.7473, -8.0061E-11, 5.7027E-11, -5335.8091, -7.6526E-11,
    5.2218E-11, -5175.8709, -7.2991E-11, 4.7408E-11, -5015.9327, -6.9456E-11,
    4.2599E-11, -24855.9945, -6.592E-11, 3.7789E-11, -4567.5913, -6.2385E-11,
    3.298E-11, -4279.5481, -5.885E-11, 2.817E-11, -3992.2248, -5.5315E-11,
    2.3361E-11, -3705.9814, -5.1779E-11, 1.8551E-11, -3421.1778, -4.8244E-11,
    1.3741E-11, -3138.174, -4.4709E-11, 3.496E-11, -9.2564E-11, -6935.1912,
    3.4293E-11, -9.0515E-11, -6775.253, 3.3626E-11, -8.8465E-11, -6615.3148,
    3.2959E-11, -8.6416E-11, -6455.3766, 3.2292E-11, -8.4366E-11, -6295.4384,
    3.1625E-11, -8.2317E-11, -6135.5002, 3.0958E-11, -8.0267E-11, -5975.562,
    3.029E-11, -7.8218E-11, -5815.6237, 2.9623E-11, -7.6168E-11, -5655.6855,
    2.8956E-11, -7.4118E-11, -5495.7473, 2.8289E-11, -7.2069E-11, -5335.8091,
    2.7622E-11, -7.0019E-11, -5175.8709, 2.6955E-11, -6.797E-11, -5015.9327,
    2.6288E-11, -6.592E-11, -24855.9945, 2.5621E-11, -6.3871E-11, -4567.5913,
    2.4954E-11, -6.1821E-11, -4279.5481, 2.4287E-11, -5.9772E-11, -3992.2248,
    2.3619E-11, -5.7722E-11, -3705.9814, 2.2952E-11, -5.5673E-11, -3421.1778,
    2.2285E-11, -5.3623E-11, -3138.174, -6483.9641, 4.7886E-11, 4.3941E-11,
    -6336.5508, 4.7109E-11, 4.2532E-11, -6189.1375, 4.6333E-11, 4.1123E-11,
    -6041.7242, 4.5556E-11, 3.9713E-11, -5894.3109, 4.4779E-11, 3.8304E-11,
    -5746.8976, 4.4003E-11, 3.6895E-11, -5599.4843, 4.3226E-11, 3.5486E-11,
    -5452.071, 4.2449E-11, 3.4076E-11, -5304.6578, 4.1673E-11, 3.2667E-11,
    -5157.2445, 4.0896E-11, 3.1258E-11, -5009.8312, 4.0119E-11, 2.9849E-11,
    -4862.4179, 3.9342E-11, 2.8439E-11, -4715.0046, 3.8566E-11, 2.703E-11,
    -4567.5913, 3.7789E-11, 2.5621E-11, -24420.178, 3.7012E-11, 2.4211E-11,
    -4144.2998, 3.6236E-11, 2.2802E-11, -3868.7815, 3.5459E-11, 2.1393E-11,
    -3593.9832, 3.4682E-11, 1.9984E-11, -3320.2647, 3.3906E-11, 1.8574E-11,
    -3047.986, 3.3129E-11, 1.7165E-11, 1.028E-10, -6483.9641, -1.085E-10,
    9.8099E-11, -6336.5508, -1.0507E-10, 9.34E-11, -6189.1375, -1.0163E-10,
    8.8701E-11, -6041.7242, -9.8201E-11, 8.4002E-11, -5894.3109, -9.4768E-11,
    7.9303E-11, -5746.8976, -9.1335E-11, 7.4604E-11, -5599.4843, -8.7902E-11,
    6.9905E-11, -5452.071, -8.4469E-11, 6.5206E-11, -5304.6578, -8.1036E-11,
    6.0507E-11, -5157.2445, -7.7603E-11, 5.5808E-11, -5009.8312, -7.417E-11,
    5.1109E-11, -4862.4179, -7.0737E-11, 4.641E-11, -4715.0046, -6.7304E-11,
    4.1711E-11, -4567.5913, -6.3871E-11, 3.7012E-11, -24420.178, -6.0438E-11,
    3.2313E-11, -4144.2998, -5.7005E-11, 2.7614E-11, -3868.7815, -5.3572E-11,
    2.2916E-11, -3593.9832, -5.0139E-11, 1.8217E-11, -3320.2647, -4.6706E-11,
    1.3518E-11, -3047.986, -4.3272E-11, 3.2949E-11, -8.7701E-11, -6483.9641,
    3.2325E-11, -8.5754E-11, -6336.5508, 3.1701E-11, -8.3806E-11, -6189.1375,
    3.1077E-11, -8.1859E-11, -6041.7242, 3.0453E-11, -7.9911E-11, -5894.3109,
    2.9828E-11, -7.7964E-11, -5746.8976, 2.9204E-11, -7.6017E-11, -5599.4843,
    2.858E-11, -7.4069E-11, -5452.071, 2.7956E-11, -7.2122E-11, -5304.6578,
    2.7332E-11, -7.0175E-11, -5157.2445, 2.6708E-11, -6.8227E-11, -5009.8312,
    2.6084E-11, -6.628E-11, -4862.4179, 2.546E-11, -6.4333E-11, -4715.0046,
    2.4836E-11, -6.2385E-11, -4567.5913, 2.4211E-11, -6.0438E-11, -24420.178,
    2.3587E-11, -5.849E-11, -4144.2998, 2.2963E-11, -5.6543E-11, -3868.7815,
    2.2339E-11, -5.4596E-11, -3593.9832, 2.1715E-11, -5.2648E-11, -3320.2647,
    2.1091E-11, -5.0701E-11, -3047.986, -6037.7762, 4.1639E-11, 4.2715E-11,
    -5902.5279, 4.0973E-11, 4.1349E-11, -5767.2796, 4.0307E-11, 3.9983E-11,
    -5632.0313, 3.9641E-11, 3.8616E-11, -5496.783, 3.8975E-11, 3.725E-11,
    -5361.5346, 3.8309E-11, 3.5884E-11, -5226.2863, 3.7642E-11, 3.4517E-11,
    -5091.038, 3.6976E-11, 3.3151E-11, -4955.7897, 3.631E-11, 3.1785E-11,
    -4820.5414, 3.5644E-11, 3.0419E-11, -4685.2931, 3.4978E-11, 2.9052E-11,
    -4550.0448, 3.4312E-11, 2.7686E-11, -4414.7964, 3.3646E-11, 2.632E-11,
    -4279.5481, 3.298E-11, 2.4954E-11, -4144.2998, 3.2313E-11, 2.3587E-11,
    -24009.0515, 3.1647E-11, 2.2221E-11, -3745.3382, 3.0981E-11, 2.0855E-11,
    -3481.985, 3.0315E-11, 1.9489E-11, -3219.3516, 2.9649E-11, 1.8122E-11,
    -2957.798, 2.8983E-11, 1.6756E-11, 1.0047E-10, -6037.7762, -1.0512E-10,
    9.5885E-11, -5902.5279, -1.0179E-10, 9.1297E-11, -5767.2796, -9.8461E-11,
    8.6708E-11, -5632.0313, -9.513E-11, 8.212E-11, -5496.783, -9.1799E-11,
    7.7531E-11, -5361.5346, -8.8468E-11, 7.2943E-11, -5226.2863, -8.5138E-11,
    6.8355E-11, -5091.038, -8.1807E-11, 6.3766E-11, -4955.7897, -7.8476E-11,
    5.9178E-11, -4820.5414, -7.5145E-11, 5.4589E-11, -4685.2931, -7.1814E-11,
    5.0001E-11, -4550.0448, -6.8483E-11, 4.5413E-11, -4414.7964, -6.5152E-11,
    4.0824E-11, -4279.5481, -6.1821E-11, 3.6236E-11, -4144.2998, -5.849E-11,
    3.1647E-11, -24009.0515, -5.516E-11, 2.7059E-11, -3745.3382, -5.1829E-11,
    2.2471E-11, -3481.985, -4.8498E-11, 1.7882E-11, -3219.3516, -4.5167E-11,
    1.3294E-11, -2957.798, -4.1836E-11, 3.0938E-11, -8.2838E-11, -6037.7762,
    3.0357E-11, -8.0992E-11, -5902.5279, 2.9775E-11, -7.9147E-11, -5767.2796,
    2.9194E-11, -7.7302E-11, -5632.0313, 2.8613E-11, -7.5457E-11, -5496.783,
    2.8032E-11, -7.3612E-11, -5361.5346, 2.7451E-11, -7.1766E-11, -5226.2863,
    2.687E-11, -6.9921E-11, -5091.038, 2.6289E-11, -6.8076E-11, -4955.7897,
    2.5708E-11, -6.6231E-11, -4820.5414, 2.5127E-11, -6.4386E-11, -4685.2931,
    2.4546E-11, -6.254E-11, -4550.0448, 2.3964E-11, -6.0695E-11, -4414.7964,
    2.3383E-11, -5.885E-11, -4279.5481, 2.2802E-11, -5.7005E-11, -4144.2998,
    2.2221E-11, -5.516E-11, -24009.0515, 2.164E-11, -5.3314E-11, -3745.3382,
    2.1059E-11, -5.1469E-11, -3481.985, 2.0478E-11, -4.9624E-11, -3219.3516,
    1.9897E-11, -4.7779E-11, -2957.798, -5596.9877, 3.5392E-11, 4.1489E-11,
    -5473.5444, 3.4837E-11, 4.0166E-11, -5350.1011, 3.4281E-11, 3.8842E-11,
    -5226.6578, 3.3726E-11, 3.7519E-11, -5103.2145, 3.317E-11, 3.6196E-11,
    -4979.7712, 3.2615E-11, 3.4873E-11, -4856.3279, 3.2059E-11, 3.3549E-11,
    -4732.8846, 3.1503E-11, 3.2226E-11, -4609.4413, 3.0948E-11, 3.0903E-11,
    -4485.998, 3.0392E-11, 2.958E-11, -4362.5547, 2.9837E-11, 2.8256E-11,
    -4239.1114, 2.9281E-11, 2.6933E-11, -4115.6681, 2.8726E-11, 2.561E-11,
    -3992.2248, 2.817E-11, 2.4287E-11, -3868.7815, 2.7614E-11, 2.2963E-11,
    -3745.3382, 2.7059E-11, 2.164E-11, -23621.895, 2.6503E-11, 2.0317E-11,
    -3369.9867, 2.5948E-11, 1.8993E-11, -3118.4384, 2.5392E-11, 1.767E-11,
    -2867.6101, 2.4837E-11, 1.6347E-11, 9.8149E-11, -5596.9877, -1.0175E-10,
    9.3671E-11, -5473.5444, -9.8516E-11, 8.9193E-11, -5350.1011, -9.5288E-11,
    8.4715E-11, -5226.6578, -9.2059E-11, 8.0237E-11, -5103.2145, -8.883E-11,
    7.576E-11, -4979.7712, -8.5602E-11, 7.1282E-11, -4856.3279, -8.2373E-11,
    6.6804E-11, -4732.8846, -7.9144E-11, 6.2326E-11, -4609.4413, -7.5915E-11,
    5.7848E-11, -4485.998, -7.2687E-11, 5.337E-11, -4362.5547, -6.9458E-11,
    4.8893E-11, -4239.1114, -6.6229E-11, 4.4415E-11, -4115.6681, -6.3E-11,
    3.9937E-11, -3992.2248, -5.9772E-11, 3.5459E-11, -3868.7815, -5.6543E-11,
    3.0981E-11, -3745.3382, -5.3314E-11, 2.6503E-11, -23621.895, -5.0086E-11,
    2.2026E-11, -3369.9867, -4.6857E-11, 1.7548E-11, -3118.4384, -4.3628E-11,
    1.307E-11, -2867.6101, -4.0399E-11, 2.8926E-11, -7.7974E-11, -5596.9877,
    2.8388E-11, -7.6231E-11, -5473.5444, 2.785E-11, -7.4488E-11, -5350.1011,
    2.7312E-11, -7.2745E-11, -5226.6578, 2.6774E-11, -7.1002E-11, -5103.2145,
    2.6236E-11, -6.9259E-11, -4979.7712, 2.5698E-11, -6.7516E-11, -4856.3279,
    2.516E-11, -6.5773E-11, -4732.8846, 2.4622E-11, -6.403E-11, -4609.4413,
    2.4083E-11, -6.2287E-11, -4485.998, 2.3545E-11, -6.0544E-11, -4362.5547,
    2.3007E-11, -5.8801E-11, -4239.1114, 2.2469E-11, -5.7058E-11, -4115.6681,
    2.1931E-11, -5.5315E-11, -3992.2248, 2.1393E-11, -5.3572E-11, -3868.7815,
    2.0855E-11, -5.1829E-11, -3745.3382, 2.0317E-11, -5.0086E-11, -23621.895,
    1.9779E-11, -4.8343E-11, -3369.9867, 1.9241E-11, -4.6599E-11, -3118.4384,
    1.8702E-11, -4.4856E-11, -2867.6101, -5161.9584, 2.9146E-11, 4.0263E-11,
    -5049.9602, 2.8701E-11, 3.8983E-11, -4937.962, 2.8256E-11, 3.7702E-11,
    -4825.9637, 2.7811E-11, 3.6422E-11, -4713.9655, 2.7366E-11, 3.5142E-11,
    -4601.9673, 2.6921E-11, 3.3861E-11, -4489.969, 2.6476E-11, 3.2581E-11,
    -4377.9708, 2.6031E-11, 3.1301E-11, -4265.9726, 2.5586E-11, 3.0021E-11,
    -4153.9743, 2.5141E-11, 2.874E-11, -4041.9761, 2.4696E-11, 2.746E-11,
    -3929.9779, 2.4251E-11, 2.618E-11, -3817.9796, 2.3806E-11, 2.49E-11,
    -3705.9814, 2.3361E-11, 2.3619E-11, -3593.9832, 2.2916E-11, 2.2339E-11,
    -3481.985, 2.2471E-11, 2.1059E-11, -3369.9867, 2.2026E-11, 1.9779E-11,
    -23257.9885, 2.1581E-11, 1.8498E-11, -3017.5253, 2.1136E-11, 1.7218E-11,
    -2777.4221, 2.0691E-11, 1.5938E-11, 9.5824E-11, -5161.9584, -9.8367E-11,
    9.1457E-11, -5049.9602, -9.5241E-11, 8.709E-11, -4937.962, -9.2114E-11,
    8.2722E-11, -4825.9637, -8.8988E-11, 7.8355E-11, -4713.9655, -8.5861E-11,
    7.3988E-11, -4601.9673, -8.2735E-11, 6.9621E-11, -4489.969, -7.9608E-11,
    6.5253E-11, -4377.9708, -7.6482E-11, 6.0886E-11, -4265.9726, -7.3355E-11,
    5.6519E-11, -4153.9743, -7.0228E-11, 5.2151E-11, -4041.9761, -6.7102E-11,
    4.7784E-11, -3929.9779, -6.3975E-11, 4.3417E-11, -3817.9796, -6.0849E-11,
    3.905E-11, -3705.9814, -5.7722E-11, 3.4682E-11, -3593.9832, -5.4596E-11,
    3.0315E-11, -3481.985, -5.1469E-11, 2.5948E-11, -3369.9867, -4.8343E-11,
    2.1581E-11, -23257.9885, -4.5216E-11, 1.7213E-11, -3017.5253, -4.2089E-11,
    1.2846E-11, -2777.4221, -3.8963E-11, 2.6915E-11, -7.3111E-11, -5161.9584,
    2.642E-11, -7.147E-11, -5049.9602, 2.5925E-11, -6.9829E-11, -4937.962,
    2.543E-11, -6.8188E-11, -4825.9637, 2.4935E-11, -6.6547E-11, -4713.9655,
    2.444E-11, -6.4906E-11, -4601.9673, 2.3944E-11, -6.3266E-11, -4489.969,
    2.3449E-11, -6.1625E-11, -4377.9708, 2.2954E-11, -5.9984E-11, -4265.9726,
    2.2459E-11, -5.8343E-11, -4153.9743, 2.1964E-11, -5.6702E-11, -4041.9761,
    2.1469E-11, -5.5061E-11, -3929.9779, 2.0974E-11, -5.342E-11, -3817.9796,
    2.0479E-11, -5.1779E-11, -3705.9814, 1.9984E-11, -5.0139E-11, -3593.9832,
    1.9489E-11, -4.8498E-11, -3481.985, 1.8993E-11, -4.6857E-11, -3369.9867,
    1.8498E-11, -4.5216E-11, -23257.9885, 1.8003E-11, -4.3575E-11, -3017.5253,
    1.7508E-11, -4.1934E-11, -2777.4221, -4733.0484, 2.2899E-11, 3.9037E-11,
    -4632.1353, 2.2564E-11, 3.7799E-11, -4531.2221, 2.223E-11, 3.6562E-11,
    -4430.309, 2.1895E-11, 3.5325E-11, -4329.3959, 2.1561E-11, 3.4088E-11,
    -4228.4828, 2.1227E-11, 3.285E-11, -4127.5697, 2.0892E-11, 3.1613E-11,
    -4026.6565, 2.0558E-11, 3.0376E-11, -3925.7434, 2.0223E-11, 2.9139E-11,
    -3824.8303, 1.9889E-11, 2.7901E-11, -3723.9172, 1.9554E-11, 2.6664E-11,
    -3623.0041, 1.922E-11, 2.5427E-11, -3522.0909, 1.8885E-11, 2.419E-11,
    -3421.1778, 1.8551E-11, 2.2952E-11, -3320.2647, 1.8217E-11, 2.1715E-11,
    -3219.3516, 1.7882E-11, 2.0478E-11, -3118.4384, 1.7548E-11, 1.9241E-11,
    -3017.5253, 1.7213E-11, 1.8003E-11, -22916.6122, 1.6879E-11, 1.6766E-11,
    -2687.2341, 1.6544E-11, 1.5529E-11, 9.35E-11, -4733.0484, -9.499E-11,
    8.9243E-11, -4632.1353, -9.1965E-11, 8.4986E-11, -4531.2221, -8.8941E-11,
    8.073E-11, -4430.309, -8.5917E-11, 7.6473E-11, -4329.3959, -8.2892E-11,
    7.2216E-11, -4228.4828, -7.9868E-11, 6.7959E-11, -4127.5697, -7.6843E-11,
    6.3703E-11, -4026.6565, -7.3819E-11, 5.9446E-11, -3925.7434, -7.0795E-11,
    5.5189E-11, -3824.8303, -6.777E-11, 5.0933E-11, -3723.9172, -6.4746E-11,
    4.6676E-11, -3623.0041, -6.1721E-11, 4.2419E-11, -3522.0909, -5.8697E-11,
    3.8162E-11, -3421.1778, -5.5673E-11, 3.3906E-11, -3320.2647, -5.2648E-11,
    2.9649E-11, -3219.3516, -4.9624E-11, 2.5392E-11, -3118.4384, -4.6599E-11,
    2.1136E-11, -3017.5253, -4.3575E-11, 1.6879E-11, -22916.6122, -4.0551E-11,
    1.2622E-11, -2687.2341, -3.7526E-11, 2.4904E-11, -6.8247E-11, -4733.0484,
    2.4452E-11, -6.6709E-11, -4632.1353, 2.4E-11, -6.517E-11, -4531.2221,
    2.3547E-11, -6.3631E-11, -4430.309, 2.3095E-11, -6.2093E-11, -4329.3959,
    2.2643E-11, -6.0554E-11, -4228.4828, 2.2191E-11, -5.9015E-11, -4127.5697,
    2.1739E-11, -5.7477E-11, -4026.6565, 2.1287E-11, -5.5938E-11, -3925.7434,
    2.0835E-11, -5.4399E-11, -3824.8303, 2.0383E-11, -5.286E-11, -3723.9172,
    1.9931E-11, -5.1322E-11, -3623.0041, 1.9479E-11, -4.9783E-11, -3522.0909,
    1.9026E-11, -4.8244E-11, -3421.1778, 1.8574E-11, -4.6706E-11, -3320.2647,
    1.8122E-11, -4.5167E-11, -3219.3516, 1.767E-11, -4.3628E-11, -3118.4384,
    1.7218E-11, -4.2089E-11, -3017.5253, 1.6766E-11, -4.0551E-11, -22916.6122,
    1.6314E-11, -3.9012E-11, -2687.2341, -4310.6175, 1.6652E-11, 3.7811E-11,
    -4220.4295, 1.6428E-11, 3.6616E-11, -4130.2416, 1.6204E-11, 3.5422E-11,
    -4040.0536, 1.598E-11, 3.4228E-11, -3949.8656, 1.5756E-11, 3.3033E-11,
    -3859.6777, 1.5532E-11, 3.1839E-11, -3769.4897, 1.5309E-11, 3.0645E-11,
    -3679.3018, 1.5085E-11, 2.9451E-11, -3589.1138, 1.4861E-11, 2.8256E-11,
    -3498.9258, 1.4637E-11, 2.7062E-11, -3408.7379, 1.4413E-11, 2.5868E-11,
    -3318.5499, 1.4189E-11, 2.4674E-11, -3228.3619, 1.3965E-11, 2.3479E-11,
    -3138.174, 1.3741E-11, 2.2285E-11, -3047.986, 1.3518E-11, 2.1091E-11,
    -2957.798, 1.3294E-11, 1.9897E-11, -2867.6101, 1.307E-11, 1.8702E-11,
    -2777.4221, 1.2846E-11, 1.7508E-11, -2687.2341, 1.2622E-11, 1.6314E-11,
    -22597.0462, 1.2398E-11, 1.512E-11, 9.1175E-11, -4310.6175, -9.1612E-11,
    8.7029E-11, -4220.4295, -8.869E-11, 8.2883E-11, -4130.2416, -8.5768E-11,
    7.8737E-11, -4040.0536, -8.2845E-11, 7.4591E-11, -3949.8656, -7.9923E-11,
    7.0444E-11, -3859.6777, -7.7001E-11, 6.6298E-11, -3769.4897, -7.4079E-11,
    6.2152E-11, -3679.3018, -7.1157E-11, 5.8006E-11, -3589.1138, -6.8234E-11,
    5.386E-11, -3498.9258, -6.5312E-11, 4.9714E-11, -3408.7379, -6.239E-11,
    4.5567E-11, -3318.5499, -5.9468E-11, 4.1421E-11, -3228.3619, -5.6545E-11,
    3.7275E-11, -3138.174, -5.3623E-11, 3.3129E-11, -3047.986, -5.0701E-11,
    2.8983E-11, -2957.798, -4.7779E-11, 2.4837E-11, -2867.6101, -4.4856E-11,
    2.0691E-11, -2777.4221, -4.1934E-11, 1.6544E-11, -2687.2341, -3.9012E-11,
    1.2398E-11, -22597.0462, -3.609E-11, 2.2892E-11, -6.3384E-11, -4310.6175,
    2.2483E-11, -6.1948E-11, -4220.4295, 2.2074E-11, -6.0511E-11, -4130.2416,
    2.1665E-11, -5.9074E-11, -4040.0536, 2.1256E-11, -5.7638E-11, -3949.8656,
    2.0847E-11, -5.6201E-11, -3859.6777, 2.0438E-11, -5.4765E-11, -3769.4897,
    2.0029E-11, -5.3328E-11, -3679.3018, 1.962E-11, -5.1892E-11, -3589.1138,
    1.9211E-11, -5.0455E-11, -3498.9258, 1.8802E-11, -4.9019E-11, -3408.7379,
    1.8392E-11, -4.7582E-11, -3318.5499, 1.7983E-11, -4.6146E-11, -3228.3619,
    1.7574E-11, -4.4709E-11, -3138.174, 1.7165E-11, -4.3272E-11, -3047.986,
    1.6756E-11, -4.1836E-11, -2957.798, 1.6347E-11, -4.0399E-11, -2867.6101,
    1.5938E-11, -3.8963E-11, -2777.4221, 1.5529E-11, -3.7526E-11, -2687.2341,
    1.512E-11, -3.609E-11, -22597.0462 };

  static const real_T dv[60] = { 0.070756, -1.8178E-16, 2.6227E-18, 0.070075,
    -1.9839E-16, 4.0625E-18, 0.061681, -2.1628E-16, 5.5283E-18, 0.053784,
    -2.3555E-16, 7.0295E-18, 0.04658, -2.563E-16, 8.5755E-18, 0.040018,
    -2.7867E-16, 1.0176E-17, 0.03405, -3.0278E-16, 1.184E-17, 0.028632,
    -3.2874E-16, 1.3579E-17, 0.023723, -3.5671E-16, 1.5402E-17, 0.019283,
    -3.8683E-16, 1.732E-17, 0.015277, -4.1924E-16, 1.9343E-17, 0.01167,
    -4.5411E-16, 2.1484E-17, 0.0084305, -4.9159E-16, 2.3754E-17, 0.00553,
    -5.3188E-16, 2.6165E-17, 0.0029409, -5.7515E-16, 2.8729E-17, 0.00063751,
    -6.216E-16, 3.146E-17, -0.0014036, -6.7142E-16, 3.4371E-17, -0.0032046,
    -7.2485E-16, 3.7477E-17, -0.0047858, -7.821E-16, 4.0792E-17, -0.0061662,
    -8.434E-16, 4.4332E-17 };

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

      igr2[i] = nlamt[i] * ((d + 0.070757) - mesil[i]);
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

