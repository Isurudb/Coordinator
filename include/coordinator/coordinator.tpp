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

  std::string Estimate_status = "Best";

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
  virtual void RunTest3(ros::NodeHandle *nh) {};
  virtual void RunTest4(ros::NodeHandle *nh) {};
  virtual void RunTest5(ros::NodeHandle *nh) {};
  virtual void RunTest6(ros::NodeHandle *nh) {};

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
double inti_e_x=0;
double inti_e_y=0;
double inti_e_z=0;
float q0_x = 0;
float q0_y = 0;
float q0_z = 0;

double z_nominal[6];
double zp_nextNominal[6];
double v_mpc[3];
bool initial_run = true;
bool initialzation = false;

int count=0;

double kn_tilda[3];
double kN[3];
bool rotation_done = false;
void step_PID();
void step_PID_good();
void step_PID_worst();

// Function Declarations
//void main_MPC_Guidance_v3_sand();
void MPC_Guidance_v3_sand();
void MPC_Guidance_v3_sand_good();
void MPC_Guidance_v3_sand_worst();
void mldivide(double A[3600], double B[60]);
bool rtIsNaN(double value);
void nominal_dynamics();
void nominal_dynamics_good();
void nominal_dynamics_worst();
void tubing_mpc();
void tubing_mpc_good();
void tubing_mpc_worst();


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
  
  uint  speed_ =3;
  flight_mode_.speed=3;
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
      else if(base_status_.test_number  == 3) {
        RunTest3(nh);
      }
      else if(base_status_.test_number  == 4) {
        RunTest4(nh);
      }
      else if(base_status_.test_number  == 5) {
        RunTest5(nh);
      }
      else if(base_status_.test_number  == 6) {
        RunTest6(nh);
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
    else if (test_number_str[2] == '1') {  // tube MPC
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
    double prod_val;
    prod_val=qx*q0_x + qy*q0_y + qz*q0_z;
    double deno;
    deno=sqrt(prod_val*prod_val);
    if ((prod_val==0)){
    deno=1;
    }
    if(prod_val/deno <-0.9 ){
      qx=-qx;
      qy=-qy;
      qz=-qz;
      qw=-qw;
      //ROS_INFO(" Quaternion sign change detected >>>>>>>>>>>>>>>>>>>>> ");
    }
    q0_x=qx;
    q0_y=qy;
    q0_z=qz;
    attitude.x=qx;
    attitude.y=qy;
    attitude.z=qz;
    attitude.w=qw;
    omega.x=wx;
    omega.y=wy;
    omega.z=wz;
   // geometry_msgs::Vector3 torque, axes_rot;
    double r=0, p=0, y=3.14159265/180*90;  // Rotate the previous pose by 45* about Z
   /*  axes_rot.x = 0;
    axes_rot.y = 0;
    axes_rot.z = 1;
    double angle = 45/180*3.14570; */
        q_ref.setRPY(r, p, y);
        //q_ref.setRotation(axes_rot, 45/180*3.14570);
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
    // Initiating PID controller <<<<<<<<<<<<<<<<<<<<<<<ID
    if(Estimate_status=="Best"){
    step_PID(); //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<ID
    }
    else if(Estimate_status=="Good"){
    step_PID_good();
    }
    else if(Estimate_status=="Worst"){
    step_PID_worst();
    }
    // MPC Controller inbound <<<<<<<<<<<<<<<<<<<<<<<ID
    x0[0]=position_error.x;
    x0[1]=position_error.y;
    x0[2]=position_error.z;
    x0[3]=vx;
    x0[4]=vy;
    x0[5]=vz;
    if(Estimate_status=="Best"){
    MPC_Guidance_v3_sand(); //>>>>>>>>>>>>>>>>>>>>>ID
    }
    else if(Estimate_status=="Good"){
     MPC_Guidance_v3_sand_good();
    }
    else if(Estimate_status=="Worst"){
     MPC_Guidance_v3_sand_worst();
    }
    // TRMPC Inbound <<<<<<<<<<<<<<<<<<<<<<<<<<<<ID
    v_mpc[0]=Fx;
    v_mpc[1]=Fy;
    v_mpc[2]=Fz;
    if(Estimate_status=="Best"){
    nominal_dynamics();
    }
    else if(Estimate_status=="Good"){
     nominal_dynamics_good();
    }  
    else if(Estimate_status=="Worst"){
     nominal_dynamics_worst();
    } 
    //sqrt(q_e.getX()*q_e.getX()+q_e.getY()*q_e.getY()+q_e.getZ()*q_e.getZ())>0.005
    if (count==0){
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

    
    kn_tilda[0]=Fx;
    kn_tilda[1]=Fy;
    kn_tilda[2]=Fz;

      /* double sx=x0[0]-zp_nextNominal[0];
      double sy=x0[1]-zp_nextNominal[1];
      double sz=x0[2]-zp_nextNominal[2];
      double svx=x0[3]-zp_nextNominal[3];
      double svy=x0[4]-zp_nextNominal[4];
      double svz=x0[5]-zp_nextNominal[5]; */
     /* -------------------------------------->>>>>>>>>>>>>>>>>>>>>>>>>>>ID   */
     if(Estimate_status=="Best"){
    tubing_mpc();
     }
     else if(Estimate_status=="Good"){
    tubing_mpc_good();
     }
      else if(Estimate_status=="Worst"){
    tubing_mpc_worst();
     }
    count+=1;
    if (count==11){

      count =0;
    } 
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
void CoordinatorBase<T>::tubing_mpc()
{
  //sub-checkpoint --------------------------------------->>>>>
  double a[18] = { 6.46950000, 0.0, 0.0, 0.0, 6.46950000, 0.0, 0.0, 0.0,
    6.46950000, 5.1423000, 0.0, 0.0, 0.0, 5.1423000, 0.0, 0.0, 0.0, 5.1423000 };
    // 2.0202, 0.0, 0.0, 0.0, 2.0202, 0.0, 0.0, 0.0,
    //2.0202, 4.0895, 0.0, 0.0, 0.0, 4.0895, 0.0, 0.0, 0.0, 4.0895 
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

template<typename T>
void CoordinatorBase<T>::tubing_mpc_good()
{
  //sub-checkpoint --------------------------------------->>>>>
  double a[18] = { 4.8997, 0.0, 0.0, 0.0, 4.8997, 0.0, 0.0, 0.0,
    4.8997, 3.8946, 0.0, 0.0, 0.0, 3.8946, 0.0, 0.0, 0.0, 3.8946 };
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

template<typename T>
void CoordinatorBase<T>::tubing_mpc_worst()
{
  //sub-checkpoint --------------------------------------->>>>>
  double a[18] = { 4.1148, 0.0, 0.0, 0.0, 4.1148, 0.0, 0.0, 0.0,
    4.1148, 3.2707, 0.0, 0.0, 0.0, 3.2707, 0.0, 0.0, 0.0, 3.2707 };
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
void CoordinatorBase<T>::nominal_dynamics()
{
  // checkpoint nominal dynamics --------------------------->>>>>>>>
   double b_a[36] = { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.16, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.16, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.16, 0.0, 0.0, 1.0 };

   double a[18] = { 0.10131, 0.0, 0.0, 0.10131, 0.0, 0.0, 0.0,
    0.10131, 0.0, 0.0, 0.10131, 0.0, 0.0, 0.0, 0.10131, 0.0, 0.0,
    0.10131 };

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
void CoordinatorBase<T>::nominal_dynamics_good()
{
  // checkpoint nominal dynamics --------------------------->>>>>>>>
   double b_a[36] = { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 2.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 1.0 };

   double a[18] = { 0.13377, 0.0, 0.0, 0.13377, 0.0, 0.0, 0.0,
    0.13377, 0.0, 0.0, 0.13377, 0.0, 0.0, 0.0, 0.13377, 0.0, 0.0,
    0.10131 };

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
void CoordinatorBase<T>::nominal_dynamics_worst()
{
  // checkpoint nominal dynamics --------------------------->>>>>>>>
   double b_a[36] = { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.16, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.16, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.16, 0.0, 0.0, 1.0 };

   double a[18] = { 0.15929, 0.0, 0.0, 0.15929, 0.0, 0.0, 0.0,
    0.15929, 0.0, 0.0, 0.15929, 0.0, 0.0, 0.0, 0.15929, 0.0, 0.0,
    0.15929 };

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
void CoordinatorBase<T>::step_PID()
{
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


 double qw = q_e.getW();
 double qx = q_e.getX();
 double qy = q_e.getY();
 double qz = q_e.getZ();
 double fx = Fx; 
 double fy = Fy;
 double fz = Fz;
 double Omega1 = arg_omegax;
 double Omega2 = arg_omegay;
 double Omega3 = arg_omegaz;
 //check point attitude controller -------------------------------->>>
 float a1[18] = { -0.0034973916096086594, -0.0, -0.0, -0.0,
    -0.011242698798043299, -0.0, -0.0, -0.0, -0.0094754742284346386,
    -0.17486958048043294, -0.0, -0.0, -0.0, -0.56213493990216479, -0.0, -0.0,
    -0.0, -0.4737737114217318  };

  static const real_T b_a[9] = { 0.0, -0.072060179322222792, 0.0,
    0.072060179322222792, 0.0, 0.14412035864444558, -0.0, -0.14412035864444558,
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
 inti_e_x-=q_e.getX()*0.16;
 inti_e_y-=q_e.getY()*0.16;
 inti_e_z-=q_e.getZ()*0.16;
 double Ki=0.0002;

  arg_tau_x = U[0] + Ki*inti_e_x;
  arg_tau_y = U[1] + Ki*inti_e_y;
  arg_tau_z = U[2] + Ki*inti_e_z;


}

template<typename T>
void CoordinatorBase<T>::step_PID_good()
{
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



 double qw = q_e.getW();
 double qx = q_e.getX();
 double qy = q_e.getY();
 double qz = q_e.getZ();
 double fx = Fx; 
 double fy = Fy;
 double fz = Fz;
 double Omega1 = arg_omegax;
 double Omega2 = arg_omegay;
 double Omega3 = arg_omegaz;
  //check point attitude controller -------------------------------->>>
 float a1[18] = {  -0.07905885393391987, -0.0, -0.0, -0.0,
    -0.23918778841959931, -0.0, -0.0, -0.0, -0.20338311048567945,
    -0.7905885393391987, -0.0, -0.0, -0.0, -2.3918778841959929, -0.0, -0.0, -0.0,
    -2.0338311048567945};

  static const real_T b_a[9] = {  0.0, -0.050293625844425129, 0.0,
    0.050293625844425129, 0.0, 0.10058725168885026, -0.0, -0.10058725168885026,
    0.0};

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
 inti_e_x-=q_e.getX()*0.16;
 inti_e_y-=q_e.getY()*0.16;
 inti_e_z-=q_e.getZ()*0.16;
 double Ki=0.0002;
 
  arg_tau_x = U[0] + Ki*inti_e_x;
  arg_tau_y = U[1] + Ki*inti_e_y;
  arg_tau_z = U[2] + Ki*inti_e_z;


}

template<typename T>
void CoordinatorBase<T>::step_PID_worst()
{
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

 double qw = q_e.getW();
 double qx = q_e.getX();
 double qy = q_e.getY();
 double qz = q_e.getZ();
 double fx = Fx; 
 double fy = Fy;
 double fz = Fz;
 double Omega1 = arg_omegax;
 double Omega2 = arg_omegay;
 double Omega3 = arg_omegaz;

 //check point attitude controller -------------------------------->>>
 float a1[18] = {  -0.079286483982532269, -0.0, -0.0, -0.0,
    -0.2403259386626613, -0.0, -0.0, -0.0, -0.20429363068012904,
    -0.79286483982532263, -0.0, -0.0, -0.0, -2.4032593866266132, -0.0, -0.0,
    -0.0, -2.0429363068012902 };

  static const real_T b_a[9] = { 0.0, -0.033182542210895193, 0.0,
    0.033182542210895193, 0.0, 0.066365084421790385, -0.0, -0.066365084421790385,
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
 inti_e_x-=q_e.getX()*0.16;
 inti_e_y-=q_e.getY()*0.16;
 inti_e_z-=q_e.getZ()*0.16;
 double Ki=0.0002;

  arg_tau_x = U[0] + Ki*inti_e_x;
  arg_tau_y = U[1] + Ki*inti_e_y;
  arg_tau_z = U[2] + Ki*inti_e_z;

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


//template<typename T>
/* void CoordinatorBase<T>::main_MPC_Guidance_v3_sand()
{

  // Initialize function 'MPC_Guidance_v3_sand' input arguments.
  // Initialize function input argument 'x0'.
  // Call the entry-point 'MPC_Guidance_v3_sand'.
  
  MPC_Guidance_v3_sand();
} */


template<typename T>
bool  CoordinatorBase<T>::rtIsNaN(double value){
  return ((value!=value) ? true : false);
}

template<typename T>
void CoordinatorBase<T>::MPC_Guidance_v3_sand()
{
   double b[14400] = { 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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
    100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0,
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
    0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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
    0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0,
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
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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
    100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0,
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
    0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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
    0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0,
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
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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
    100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0,
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
    0.0, 0.0, 0.0, 0.0, 201556.5626, 1.0696E-9, -9.3512E-10, 123391.3808,
    -3.3448E-9, 8.6161E-8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0696E-9,
    201556.5626, -1.0925E-8, 4.2876E-10, 123391.3808, -3.1308E-8, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -9.3512E-10, -1.0925E-8, 201556.5626,
    -1.3582E-8, -1.0697E-7, 123391.3808, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 123391.3808, 4.2876E-10, -1.3582E-8, 4.0282676788E+6, -1.9252E-9,
    4.2611E-8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -3.3448E-9,
    123391.3808, -1.0697E-7, -1.9252E-9, 4.0282676788E+6, -1.6701E-7, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.6161E-8, -3.1308E-8, 123391.3808,
    4.2611E-8, -1.6701E-7, 4.0282676788E+6 };

   double b_b[7200] = { 6.484E-6, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    1.9452E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 3.242E-5, 0.0, 0.0, 0.0008105, 0.0,
    0.0, 4.5388E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 5.8356E-5, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 7.1324E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 8.4292E-5,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 9.726E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    0.00011023, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0001232, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 0.00013616, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00014913, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 0.0001621, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00017507,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00018804, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    0.000201, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00021397, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 0.00022694, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00023991, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 0.00025287, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0,
    6.484E-6, 0.0, 0.0, 0.0008105, 0.0, 0.0, 1.9452E-5, 0.0, 0.0, 0.0008105, 0.0,
    0.0, 3.242E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 4.5388E-5, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 5.8356E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 7.1324E-5, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 8.4292E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 9.726E-5, 0.0,
    0.0, 0.0008105, 0.0, 0.0, 0.00011023, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    0.0001232, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00013616, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 0.00014913, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0001621, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 0.00017507, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00018804,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 0.000201, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    0.00021397, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00022694, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 0.00023991, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00025287, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 0.0, 6.484E-6, 0.0, 0.0, 0.0008105, 0.0, 0.0, 1.9452E-5,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 3.242E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    4.5388E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 5.8356E-5, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 7.1324E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 8.4292E-5, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 9.726E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00011023,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0001232, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    0.00013616, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00014913, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 0.0001621, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00017507, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 0.00018804, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.000201,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00021397, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    0.00022694, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00023991, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 0.00025287, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    6.484E-6, 0.0, 0.0, 0.0008105, 0.0, 0.0, 1.9452E-5, 0.0, 0.0, 0.0008105, 0.0,
    0.0, 3.242E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 4.5388E-5, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 5.8356E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 7.1324E-5, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 8.4292E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 9.726E-5, 0.0,
    0.0, 0.0008105, 0.0, 0.0, 0.00011023, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    0.0001232, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00013616, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 0.00014913, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0001621, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 0.00017507, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00018804,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 0.000201, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    0.00021397, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00022694, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 0.00023991, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 6.484E-6, 0.0, 0.0, 0.0008105, 0.0, 0.0, 1.9452E-5, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 3.242E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 4.5388E-5, 0.0,
    0.0, 0.0008105, 0.0, 0.0, 5.8356E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    7.1324E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 8.4292E-5, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 9.726E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00011023, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 0.0001232, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00013616,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00014913, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    0.0001621, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00017507, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 0.00018804, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.000201, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 0.00021397, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00022694,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00023991, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.484E-6, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    1.9452E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 3.242E-5, 0.0, 0.0, 0.0008105, 0.0,
    0.0, 4.5388E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 5.8356E-5, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 7.1324E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 8.4292E-5,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 9.726E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    0.00011023, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0001232, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 0.00013616, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00014913, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 0.0001621, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00017507,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00018804, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    0.000201, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00021397, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 0.00022694, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00023991, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    6.484E-6, 0.0, 0.0, 0.0008105, 0.0, 0.0, 1.9452E-5, 0.0, 0.0, 0.0008105, 0.0,
    0.0, 3.242E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 4.5388E-5, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 5.8356E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 7.1324E-5, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 8.4292E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 9.726E-5, 0.0,
    0.0, 0.0008105, 0.0, 0.0, 0.00011023, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    0.0001232, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00013616, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 0.00014913, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0001621, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 0.00017507, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00018804,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 0.000201, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    0.00021397, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00022694, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    6.484E-6, 0.0, 0.0, 0.0008105, 0.0, 0.0, 1.9452E-5, 0.0, 0.0, 0.0008105, 0.0,
    0.0, 3.242E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 4.5388E-5, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 5.8356E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 7.1324E-5, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 8.4292E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 9.726E-5, 0.0,
    0.0, 0.0008105, 0.0, 0.0, 0.00011023, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    0.0001232, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00013616, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 0.00014913, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0001621, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 0.00017507, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00018804,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 0.000201, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    0.00021397, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00022694, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    6.484E-6, 0.0, 0.0, 0.0008105, 0.0, 0.0, 1.9452E-5, 0.0, 0.0, 0.0008105, 0.0,
    0.0, 3.242E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 4.5388E-5, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 5.8356E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 7.1324E-5, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 8.4292E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 9.726E-5, 0.0,
    0.0, 0.0008105, 0.0, 0.0, 0.00011023, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    0.0001232, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00013616, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 0.00014913, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0001621, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 0.00017507, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00018804,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 0.000201, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    0.00021397, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00022694, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 6.484E-6, 0.0, 0.0, 0.0008105, 0.0, 0.0, 1.9452E-5, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 3.242E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 4.5388E-5, 0.0,
    0.0, 0.0008105, 0.0, 0.0, 5.8356E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    7.1324E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 8.4292E-5, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 9.726E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00011023, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 0.0001232, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00013616,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00014913, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    0.0001621, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00017507, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 0.00018804, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.000201, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 0.00021397, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 6.484E-6, 0.0, 0.0, 0.0008105, 0.0, 0.0, 1.9452E-5, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 3.242E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 4.5388E-5, 0.0,
    0.0, 0.0008105, 0.0, 0.0, 5.8356E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    7.1324E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 8.4292E-5, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 9.726E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00011023, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 0.0001232, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00013616,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00014913, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    0.0001621, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00017507, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 0.00018804, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.000201, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 0.00021397, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 6.484E-6, 0.0, 0.0, 0.0008105, 0.0, 0.0, 1.9452E-5, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 3.242E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 4.5388E-5, 0.0,
    0.0, 0.0008105, 0.0, 0.0, 5.8356E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    7.1324E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 8.4292E-5, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 9.726E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00011023, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 0.0001232, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00013616,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00014913, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    0.0001621, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00017507, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 0.00018804, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.000201, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 0.00021397, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 6.484E-6, 0.0, 0.0, 0.0008105, 0.0, 0.0, 1.9452E-5,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 3.242E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    4.5388E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 5.8356E-5, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 7.1324E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 8.4292E-5, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 9.726E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00011023,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0001232, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    0.00013616, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00014913, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 0.0001621, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00017507, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 0.00018804, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.000201,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 6.484E-6, 0.0, 0.0, 0.0008105, 0.0, 0.0, 1.9452E-5, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 3.242E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 4.5388E-5, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 5.8356E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 7.1324E-5,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 8.4292E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    9.726E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00011023, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 0.0001232, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00013616, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 0.00014913, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0001621,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00017507, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    0.00018804, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.000201, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.484E-6, 0.0,
    0.0, 0.0008105, 0.0, 0.0, 1.9452E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 3.242E-5,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 4.5388E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    5.8356E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 7.1324E-5, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 8.4292E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 9.726E-5, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 0.00011023, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0001232,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00013616, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    0.00014913, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0001621, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 0.00017507, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00018804, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 0.000201, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.484E-6, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 1.9452E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 3.242E-5, 0.0,
    0.0, 0.0008105, 0.0, 0.0, 4.5388E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    5.8356E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 7.1324E-5, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 8.4292E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 9.726E-5, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 0.00011023, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0001232,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00013616, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    0.00014913, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0001621, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 0.00017507, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00018804, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 6.484E-6, 0.0, 0.0, 0.0008105, 0.0, 0.0, 1.9452E-5,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 3.242E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    4.5388E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 5.8356E-5, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 7.1324E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 8.4292E-5, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 9.726E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00011023,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0001232, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    0.00013616, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00014913, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 0.0001621, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00017507, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 0.00018804, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    6.484E-6, 0.0, 0.0, 0.0008105, 0.0, 0.0, 1.9452E-5, 0.0, 0.0, 0.0008105, 0.0,
    0.0, 3.242E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 4.5388E-5, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 5.8356E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 7.1324E-5, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 8.4292E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 9.726E-5, 0.0,
    0.0, 0.0008105, 0.0, 0.0, 0.00011023, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    0.0001232, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00013616, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 0.00014913, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0001621, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 0.00017507, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00018804,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.484E-6, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 1.9452E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 3.242E-5, 0.0,
    0.0, 0.0008105, 0.0, 0.0, 4.5388E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    5.8356E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 7.1324E-5, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 8.4292E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 9.726E-5, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 0.00011023, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0001232,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00013616, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    0.00014913, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0001621, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 0.00017507, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 6.484E-6, 0.0, 0.0, 0.0008105, 0.0, 0.0, 1.9452E-5, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 3.242E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 4.5388E-5, 0.0,
    0.0, 0.0008105, 0.0, 0.0, 5.8356E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    7.1324E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 8.4292E-5, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 9.726E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00011023, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 0.0001232, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00013616,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00014913, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    0.0001621, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00017507, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.484E-6, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 1.9452E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 3.242E-5, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 4.5388E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 5.8356E-5,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 7.1324E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    8.4292E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 9.726E-5, 0.0, 0.0, 0.0008105, 0.0,
    0.0, 0.00011023, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0001232, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 0.00013616, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00014913,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0001621, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    0.00017507, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 6.484E-6, 0.0, 0.0, 0.0008105, 0.0, 0.0, 1.9452E-5, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 3.242E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 4.5388E-5, 0.0,
    0.0, 0.0008105, 0.0, 0.0, 5.8356E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    7.1324E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 8.4292E-5, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 9.726E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00011023, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 0.0001232, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00013616,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00014913, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    0.0001621, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.484E-6, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    1.9452E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 3.242E-5, 0.0, 0.0, 0.0008105, 0.0,
    0.0, 4.5388E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 5.8356E-5, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 7.1324E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 8.4292E-5,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 9.726E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    0.00011023, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0001232, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 0.00013616, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00014913, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 0.0001621, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.484E-6, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 1.9452E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 3.242E-5, 0.0,
    0.0, 0.0008105, 0.0, 0.0, 4.5388E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    5.8356E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 7.1324E-5, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 8.4292E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 9.726E-5, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 0.00011023, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0001232,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00013616, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    0.00014913, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0001621, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 6.484E-6, 0.0, 0.0, 0.0008105, 0.0, 0.0, 1.9452E-5, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 3.242E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 4.5388E-5, 0.0,
    0.0, 0.0008105, 0.0, 0.0, 5.8356E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    7.1324E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 8.4292E-5, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 9.726E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00011023, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 0.0001232, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00013616,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00014913, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 6.484E-6, 0.0, 0.0, 0.0008105, 0.0, 0.0, 1.9452E-5, 0.0,
    0.0, 0.0008105, 0.0, 0.0, 3.242E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 4.5388E-5,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 5.8356E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    7.1324E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 8.4292E-5, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 9.726E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00011023, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 0.0001232, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00013616,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00014913, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 6.484E-6, 0.0, 0.0, 0.0008105, 0.0, 0.0, 1.9452E-5, 0.0,
    0.0, 0.0008105, 0.0, 0.0, 3.242E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 4.5388E-5,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 5.8356E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    7.1324E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 8.4292E-5, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 9.726E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00011023, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 0.0001232, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00013616,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00014913, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.484E-6, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    1.9452E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 3.242E-5, 0.0, 0.0, 0.0008105, 0.0,
    0.0, 4.5388E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 5.8356E-5, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 7.1324E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 8.4292E-5,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 9.726E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    0.00011023, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0001232, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 0.00013616, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 6.484E-6, 0.0, 0.0, 0.0008105, 0.0, 0.0, 1.9452E-5,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 3.242E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    4.5388E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 5.8356E-5, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 7.1324E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 8.4292E-5, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 9.726E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00011023,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0001232, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    0.00013616, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 6.484E-6, 0.0, 0.0, 0.0008105, 0.0, 0.0, 1.9452E-5, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 3.242E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 4.5388E-5, 0.0,
    0.0, 0.0008105, 0.0, 0.0, 5.8356E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    7.1324E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 8.4292E-5, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 9.726E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00011023, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 0.0001232, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00013616,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 6.484E-6, 0.0, 0.0, 0.0008105, 0.0, 0.0, 1.9452E-5, 0.0,
    0.0, 0.0008105, 0.0, 0.0, 3.242E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 4.5388E-5,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 5.8356E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    7.1324E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 8.4292E-5, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 9.726E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.00011023, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 0.0001232, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.484E-6,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 1.9452E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    3.242E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 4.5388E-5, 0.0, 0.0, 0.0008105, 0.0,
    0.0, 5.8356E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 7.1324E-5, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 8.4292E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 9.726E-5, 0.0,
    0.0, 0.0008105, 0.0, 0.0, 0.00011023, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    0.0001232, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.484E-6, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 1.9452E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 3.242E-5, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 4.5388E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 5.8356E-5,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 7.1324E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    8.4292E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 9.726E-5, 0.0, 0.0, 0.0008105, 0.0,
    0.0, 0.00011023, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0001232, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.484E-6, 0.0, 0.0, 0.0008105, 0.0,
    0.0, 1.9452E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 3.242E-5, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 4.5388E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 5.8356E-5, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 7.1324E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 8.4292E-5,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 9.726E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    0.00011023, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    6.484E-6, 0.0, 0.0, 0.0008105, 0.0, 0.0, 1.9452E-5, 0.0, 0.0, 0.0008105, 0.0,
    0.0, 3.242E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 4.5388E-5, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 5.8356E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 7.1324E-5, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 8.4292E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 9.726E-5, 0.0,
    0.0, 0.0008105, 0.0, 0.0, 0.00011023, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.484E-6, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    1.9452E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 3.242E-5, 0.0, 0.0, 0.0008105, 0.0,
    0.0, 4.5388E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 5.8356E-5, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 7.1324E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 8.4292E-5,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 9.726E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    0.00011023, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 6.484E-6, 0.0, 0.0, 0.0008105, 0.0, 0.0, 1.9452E-5, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 3.242E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 4.5388E-5, 0.0,
    0.0, 0.0008105, 0.0, 0.0, 5.8356E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    7.1324E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 8.4292E-5, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 9.726E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.484E-6, 0.0, 0.0, 0.0008105, 0.0,
    0.0, 1.9452E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 3.242E-5, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 4.5388E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 5.8356E-5, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 7.1324E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 8.4292E-5,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 9.726E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.484E-6, 0.0,
    0.0, 0.0008105, 0.0, 0.0, 1.9452E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 3.242E-5,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 4.5388E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    5.8356E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 7.1324E-5, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 8.4292E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 9.726E-5, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 6.484E-6, 0.0, 0.0, 0.0008105, 0.0, 0.0, 1.9452E-5,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 3.242E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    4.5388E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 5.8356E-5, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 7.1324E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 8.4292E-5, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.484E-6, 0.0, 0.0, 0.0008105, 0.0,
    0.0, 1.9452E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 3.242E-5, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 4.5388E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 5.8356E-5, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 7.1324E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 8.4292E-5,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.484E-6, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 1.9452E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 3.242E-5, 0.0,
    0.0, 0.0008105, 0.0, 0.0, 4.5388E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    5.8356E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 7.1324E-5, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 8.4292E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 6.484E-6, 0.0, 0.0, 0.0008105, 0.0, 0.0, 1.9452E-5, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 3.242E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 4.5388E-5, 0.0,
    0.0, 0.0008105, 0.0, 0.0, 5.8356E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    7.1324E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 6.484E-6, 0.0, 0.0, 0.0008105, 0.0, 0.0, 1.9452E-5, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 3.242E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 4.5388E-5, 0.0,
    0.0, 0.0008105, 0.0, 0.0, 5.8356E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    7.1324E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 6.484E-6, 0.0, 0.0, 0.0008105, 0.0, 0.0, 1.9452E-5, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 3.242E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 4.5388E-5, 0.0,
    0.0, 0.0008105, 0.0, 0.0, 5.8356E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    7.1324E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.484E-6, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    1.9452E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 3.242E-5, 0.0, 0.0, 0.0008105, 0.0,
    0.0, 4.5388E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 5.8356E-5, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 6.484E-6, 0.0, 0.0, 0.0008105, 0.0, 0.0, 1.9452E-5,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 3.242E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    4.5388E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 5.8356E-5, 0.0, 0.0, 0.0008105,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 6.484E-6, 0.0, 0.0, 0.0008105, 0.0, 0.0, 1.9452E-5, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 3.242E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 4.5388E-5, 0.0,
    0.0, 0.0008105, 0.0, 0.0, 5.8356E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 6.484E-6, 0.0, 0.0, 0.0008105, 0.0, 0.0, 1.9452E-5, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 3.242E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 4.5388E-5, 0.0,
    0.0, 0.0008105, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.484E-6, 0.0,
    0.0, 0.0008105, 0.0, 0.0, 1.9452E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 3.242E-5,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 4.5388E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.484E-6, 0.0, 0.0, 0.0008105, 0.0, 0.0,
    1.9452E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 3.242E-5, 0.0, 0.0, 0.0008105, 0.0,
    0.0, 4.5388E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 6.484E-6, 0.0, 0.0, 0.0008105, 0.0, 0.0, 1.9452E-5, 0.0,
    0.0, 0.0008105, 0.0, 0.0, 3.242E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.484E-6, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 1.9452E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 3.242E-5, 0.0,
    0.0, 0.0008105, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 6.484E-6, 0.0, 0.0, 0.0008105, 0.0, 0.0, 1.9452E-5, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 3.242E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.484E-6,
    0.0, 0.0, 0.0008105, 0.0, 0.0, 1.9452E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 6.484E-6, 0.0, 0.0, 0.0008105, 0.0, 0.0, 1.9452E-5, 0.0, 0.0,
    0.0008105, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.484E-6, 0.0, 0.0, 0.0008105, 0.0,
    0.0, 1.9452E-5, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 6.484E-6, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 6.484E-6, 0.0, 0.0, 0.0008105, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.484E-6, 0.0, 0.0, 0.0008105 };

  
   double b_H[3600] = { 207.91568962305743, -3.58787126963432E-15,
    8.5613912980450155E-14, 7.7803882041935317, -3.6038893042877592E-15,
    8.5905374951120982E-14, 7.6450837995329533, -3.6199196985358394E-15,
    8.619706181528774E-14, 7.5097794622411449, -3.6359500927839195E-15,
    8.6488748679454472E-14, 7.3744752259570951, -3.651980487032E-15,
    8.67804355436212E-14, 7.23917414334048, -3.66799852168544E-15,
    8.7071897514292056E-14, 7.103870209457205, -3.68402891593352E-15,
    8.7363584378458788E-14, 6.9685664774727272, -3.7000593101816E-15,
    8.765527124262552E-14, 6.8332629810260386, -3.71608970442968E-15,
    8.7946958106792252E-14, 6.6979597537587221, -3.73212009867776E-15,
    8.8238644970959E-14, 6.5626598483356409, -3.7481381333312E-15,
    8.8530106941629836E-14, 6.42735725981494, -3.7641685275792794E-15,
    8.8821793805796568E-14, 6.292055041354276, -3.7801989218273596E-15,
    8.9113480669963326E-14, 6.1567538304019642, -3.7962268441565124E-15,
    8.9405122555430869E-14, 6.0214530566917084, -3.8122547664856636E-15,
    8.9696764440898438E-14, 5.8861527538625005, -3.8282826888148164E-15,
    8.9988406326366E-14, 5.7508529555403651, -3.8443106111439677E-15,
    9.028004821183355E-14, 5.6155536953513252, -3.86033853347312E-15,
    9.05716900973011E-14, 5.4802550069343736, -3.8763664558022717E-15,
    9.0863331982768662E-14, 5.3449569239285024, -3.8923943781314237E-15,
    9.1154973868236218E-14, -3.58787126963432E-15, 207.91568962305743,
    -2.77498787522325E-13, -3.52461375203816E-15, 7.7803882041935317,
    -2.7517993793200503E-13, -3.4613074246290396E-15, 7.6450837995329533,
    -2.7285929910589E-13, -3.3980010972199203E-15, 7.5097794622411449,
    -2.70538660279775E-13, -3.3346947698108E-15, 7.3744752259570951,
    -2.6821802145366004E-13, -3.27143725221464E-15, 7.23917414334048,
    -2.6589917186334004E-13, -3.20813092480552E-15, 7.103870209457205,
    -2.63578533037225E-13, -3.1448245973963997E-15, 6.9685664774727272,
    -2.6125789421111E-13, -3.08151826998728E-15, 6.8332629810260386,
    -2.58937255384995E-13, -3.01821194257816E-15, 6.6979597537587221,
    -2.5661661655888003E-13, -2.954954424982E-15, 6.5626598483356409,
    -2.5429776696856003E-13, -2.8916480975728796E-15, 6.42735725981494,
    -2.51977128142445E-13, -2.82834177016376E-15, 6.292055041354276,
    -2.4965648931633E-13, -2.7650452047172323E-15, 6.1567538304019642,
    -2.47336208337374E-13, -2.701748639270704E-15, 6.0214530566917084,
    -2.45015927358418E-13, -2.6384520738241761E-15, 5.8861527538625005,
    -2.4269564637946204E-13, -2.5751555083776482E-15, 5.7508529555403651,
    -2.40375365400506E-13, -2.5118589429311203E-15, 5.6155536953513252,
    -2.3805508442155E-13, -2.448562377484592E-15, 5.4802550069343736,
    -2.35734803442594E-13, -2.3852658120380641E-15, 5.3449569239285024,
    -2.33414522463638E-13, 8.5613912980450143E-14, -2.77498787522325E-13,
    207.91568962305743, 8.3809958048241E-14, -2.76769457333925E-13,
    7.7803882041935317, 8.2004611175497722E-14, -2.7603956439075E-13,
    7.6450837995329533, 8.0199264302754473E-14, -2.75309671447575E-13,
    7.5097794622411449, 7.8393917430011212E-14, -2.745797785044E-13,
    7.3744752259570951, 7.6589962497802053E-14, -2.7385044831599997E-13,
    7.23917414334048, 7.4784615625058791E-14, -2.73120555372825E-13,
    7.103870209457205, 7.297926875231553E-14, -2.7239066242965E-13,
    6.9685664774727272, 7.1173921879572256E-14, -2.71660769486475E-13,
    6.8332629810260386, 6.9368575006829E-14, -2.709308765433E-13,
    6.6979597537587221, 6.7564620074619835E-14, -2.7020154635490003E-13,
    6.5626598483356409, 6.5759273201876587E-14, -2.69471653411725E-13,
    6.42735725981494, 6.3953926329133326E-14, -2.6874176046855E-13,
    6.292055041354276, 6.214885784449687E-14, -2.6801198007633E-13,
    6.1567538304019642, 6.0343789359860439E-14, -2.6728219968410996E-13,
    6.0214530566917084, 5.8538720875224E-14, -2.6655241929189003E-13,
    5.8861527538625005, 5.6733652390587558E-14, -2.6582263889967005E-13,
    5.7508529555403651, 5.4928583905951108E-14, -2.6509285850745E-13,
    5.6155536953513252, 5.3123515421314665E-14, -2.6436307811523E-13,
    5.4802550069343736, 5.1318446936678215E-14, -2.6363329772301E-13,
    5.3449569239285024, 7.7803882041935317, -3.52461375203816E-15,
    8.3809958048241E-14, 207.77653715668245, -3.54027248323688E-15,
    8.410110589040905E-14, 7.64130110038759, -3.55594329678992E-15,
    8.43924783836899E-14, 7.5060650778225018, -3.57161411034296E-15,
    8.4683850876970761E-14, 7.3708291226261844, -3.5872849238959995E-15,
    8.4975223370251591E-14, 7.2355962352407879, -3.6029436550947196E-15,
    8.5266371212419671E-14, 7.10036051519322, -3.61861446864776E-15,
    8.5557743705700514E-14, 6.9651249634054349, -3.6342852822008E-15,
    8.5849116198981356E-14, 6.829889613516448, -3.64995609575384E-15,
    8.6140488692262212E-14, 6.6946544991652495, -3.66562690930688E-15,
    8.6431861185543067E-14, 6.5594226207965862, -3.6812856405056E-15,
    8.6723009027711121E-14, 6.424188077942552, -3.69695645405864E-15,
    8.7014381520991977E-14, 6.2889538715173421, -3.71262726761168E-15,
    8.7305754014272832E-14, 6.1537206285128008, -3.728295664693856E-15,
    8.7597081577331112E-14, 6.0184877891113242, -3.7439640617760322E-15,
    8.78884091403894E-14, 5.8832553869519026, -3.7596324588582085E-15,
    8.81797367034477E-14, 5.74802345567353, -3.775300855940384E-15,
    8.8471064266505977E-14, 5.6127920289022279, -3.7909692530225595E-15,
    8.8762391829564282E-14, 5.4775611402640232, -3.8066376501047358E-15,
    8.9053719392622562E-14, 5.3423308233979059, -3.8223060471869121E-15,
    8.9345046955680855E-14, -3.60388930428776E-15, 7.7803882041935317,
    -2.7676945733392495E-13, -3.54027248323688E-15, 207.77653715668245,
    -2.74454277704565E-13, -3.4766065751327197E-15, 7.64130110038759,
    -2.7213731167117E-13, -3.4129406670285602E-15, 7.5060650778225018,
    -2.69820345637775E-13, -3.3492747589244E-15, 7.3708291226261844,
    -2.6750337960438E-13, -3.28565793787352E-15, 7.2355962352407879,
    -2.6518819997502E-13, -3.22199202976936E-15, 7.10036051519322,
    -2.6287123394162504E-13, -3.1583261216652E-15, 6.9651249634054349,
    -2.6055426790823E-13, -3.09466021356104E-15, 6.829889613516448,
    -2.58237301874835E-13, -3.0309943054568802E-15, 6.6946544991652495,
    -2.5592033584144E-13, -2.967377484406E-15, 6.5594226207965862,
    -2.5360515621208E-13, -2.9037115763018395E-15, 6.424188077942552,
    -2.51288190178685E-13, -2.8400456681976797E-15, 6.2889538715173421,
    -2.4897122414529E-13, -2.776389577504176E-15, 6.1537206285128008,
    -2.46654615392702E-13, -2.712733486810672E-15, 6.0184877891113242,
    -2.44338006640114E-13, -2.649077396117168E-15, 5.8832553869519026,
    -2.4202139788752603E-13, -2.5854213054236643E-15, 5.74802345567353,
    -2.39704789134938E-13, -2.5217652147301603E-15, 5.6127920289022279,
    -2.3738818038235E-13, -2.4581091240366563E-15, 5.4775611402640232,
    -2.35071571629762E-13, -2.3944530333431522E-15, 5.3423308233979059,
    -2.32754962877174E-13, 8.5905374951120982E-14, -2.75179937932005E-13,
    7.7803882041935317, 8.4101105890409063E-14, -2.74454277704565E-13,
    207.77653715668245, 8.22954446467799E-14, -2.7372805755411004E-13,
    7.64130110038759, 8.0489783403150762E-14, -2.7300183740365496E-13,
    7.5060650778225018, 7.8684122159521611E-14, -2.722756172532E-13,
    7.3708291226261844, 7.6879853098809667E-14, -2.7154995702575997E-13,
    7.2355962352407879, 7.5074191855180529E-14, -2.70823736875305E-13,
    7.10036051519322, 7.3268530611551365E-14, -2.7009751672485003E-13,
    6.9651249634054349, 7.14628693679222E-14, -2.69371296574395E-13,
    6.829889613516448, 6.9657208124293051E-14, -2.6864507642394E-13,
    6.6946544991652495, 6.785293906358112E-14, -2.679194161965E-13,
    6.5594226207965862, 6.6047277819951969E-14, -2.67193196046045E-13,
    6.424188077942552, 6.4241616576322818E-14, -2.6646697589559E-13,
    6.2889538715173421, 6.2436233769277111E-14, -2.65740867729738E-13,
    6.1537206285128008, 6.06308509622314E-14, -2.65014759563886E-13,
    6.0184877891113242, 5.88254681551857E-14, -2.6428865139803404E-13,
    5.8832553869519026, 5.7020085348139991E-14, -2.6356254323218205E-13,
    5.74802345567353, 5.5214702541094271E-14, -2.6283643506633E-13,
    5.6127920289022279, 5.3409319734048564E-14, -2.62110326900478E-13,
    5.4775611402640232, 5.1603936927002851E-14, -2.6138421873462603E-13,
    5.3423308233979059, 7.6450837995329533, -3.4613074246290396E-15,
    8.2004611175497735E-14, 7.64130110038759, -3.4766065751327197E-15,
    8.2295444646779912E-14, 207.63751548632862, -3.4919175305364797E-15,
    8.2586502526604744E-14, 7.5023478307688318, -3.5072284859402396E-15,
    8.2877560406429576E-14, 7.3671802089388274, -3.5225394413439996E-15,
    8.31686182862544E-14, 7.2320155690229226, -3.5378385918476796E-15,
    8.3459451757536573E-14, 7.0968481150896263, -3.55314954725144E-15,
    8.37505096373614E-14, 6.9616807957771414, -3.5684605026552004E-15,
    8.4041567517186237E-14, 6.82651364472444, -3.58377145805896E-15,
    8.4332625397011069E-14, 6.6913466955705356, -3.59908241346272E-15,
    8.46236832768359E-14, 6.5561828964997506, -3.6143815639664004E-15,
    8.4914516748118078E-14, 6.4210164515831876, -3.62969251937016E-15,
    8.520557462794291E-14, 6.2858503094642364, -3.64500347477392E-15,
    8.5496632507767755E-14, 6.1506850866831755, -3.6603120691976639E-15,
    8.5787645505884048E-14, 6.0155202338661873, -3.6756206636214082E-15,
    8.6078658504000354E-14, 5.8803557846522629, -3.6909292580451517E-15,
    8.6369671502116648E-14, 5.7451917726803936, -3.706237852468896E-15,
    8.6660684500232954E-14, 5.6100282315895731, -3.7215464468926396E-15,
    8.6951697498349247E-14, 5.4748651950058242, -3.7368550413163839E-15,
    8.7242710496465553E-14, 5.3397026965551717, -3.7521636357401274E-15,
    8.7533723494581847E-14, -3.61991969853584E-15, 7.6450837995329533,
    -2.7603956439075E-13, -3.55594329678992E-15, 7.64130110038759,
    -2.7372805755411004E-13, -3.4919175305364797E-15, 207.63751548632862,
    -2.7141476714738E-13, -3.42789176428304E-15, 7.5023478307688318,
    -2.6910147674065E-13, -3.3638659980295997E-15, 7.3671802089388274,
    -2.6678818633392003E-13, -3.29988959628368E-15, 7.2320155690229226,
    -2.6447667949728E-13, -3.23586383003024E-15, 7.0968481150896263,
    -2.6216338909055003E-13, -3.1718380637768E-15, 6.9616807957771414,
    -2.5985009868382E-13, -3.10781229752336E-15, 6.82651364472444,
    -2.5753680827709E-13, -3.04378653126992E-15, 6.6913466955705356,
    -2.5522351787036E-13, -2.9798101295240003E-15, 6.5561828964997506,
    -2.5291201103372E-13, -2.9157843632705596E-15, 6.4210164515831876,
    -2.5059872062699E-13, -2.8517585970171196E-15, 6.2858503094642364,
    -2.4828543022026E-13, -2.7877427036651844E-15, 6.1506850866831755,
    -2.45972496527548E-13, -2.723726810313248E-15, 6.0155202338661873,
    -2.4365956283483603E-13, -2.6597109169613119E-15, 5.8803557846522629,
    -2.4134662914212403E-13, -2.5956950236093763E-15, 5.7451917726803936,
    -2.3903369544941203E-13, -2.5316791302574403E-15, 5.6100282315895731,
    -2.367207617567E-13, -2.4676632369055039E-15, 5.4748651950058242,
    -2.34407828063988E-13, -2.4036473435535683E-15, 5.3397026965551717,
    -2.32094894371276E-13, 8.6197061815287727E-14, -2.7285929910588996E-13,
    7.6450837995329533, 8.43924783836899E-14, -2.7213731167117E-13,
    7.64130110038759, 8.2586502526604732E-14, -2.7141476714738E-13,
    207.63751548632862, 8.078052666951957E-14, -2.7069222262358997E-13,
    7.5023478307688318, 7.8974550812434409E-14, -2.699696780998E-13,
    7.3671802089388274, 7.7169967380836587E-14, -2.6924769066508E-13,
    7.2320155690229226, 7.5363991523751426E-14, -2.6852514614129E-13,
    7.0968481150896263, 7.3558015666666239E-14, -2.678026016175E-13,
    6.9616807957771414, 7.1752039809581065E-14, -2.6708005709371003E-13,
    6.82651364472444, 6.99460639524959E-14, -2.6635751256992E-13,
    6.6913466955705356, 6.8141480520898082E-14, -2.656355251352E-13,
    6.5561828964997506, 6.6335504663812921E-14, -2.6491298061141E-13,
    6.4210164515831876, 6.4529528806727747E-14, -2.6419043608762E-13,
    6.2858503094642364, 6.2723831434740041E-14, -2.63468002981644E-13,
    6.1506850866831755, 6.0918134062752347E-14, -2.62745569875668E-13,
    6.0155202338661873, 5.9112436690764654E-14, -2.6202313676969203E-13,
    5.8803557846522629, 5.730673931877696E-14, -2.61300703663716E-13,
    5.7451917726803936, 5.5501041946789248E-14, -2.6057827055774E-13,
    5.6100282315895731, 5.3695344574801555E-14, -2.59855837451764E-13,
    5.4748651950058242, 5.1889647202813842E-14, -2.59133404345788E-13,
    5.3397026965551717, 7.5097794622411449, -3.3980010972199203E-15,
    8.0199264302754473E-14, 7.5060650778225018, -3.41294066702856E-15,
    8.0489783403150749E-14, 7.5023478307688318, -3.42789176428304E-15,
    8.0780526669519583E-14, 207.49863058380075, -3.44284286153752E-15,
    8.1071269935888391E-14, 7.3635312953318808, -3.4577939587919997E-15,
    8.13620132022572E-14, 7.2284349028802772, -3.4727335286006397E-15,
    8.16525323026535E-14, 7.0933357150560647, -3.48768462585512E-15,
    8.1943275569022309E-14, 6.9582366282136716, -3.5026357231096E-15,
    8.223401883539113E-14, 6.823137675992089, -3.51758682036408E-15,
    8.2524762101759938E-14, 6.6880388920302911, -3.53253791761856E-15,
    8.2815505368128759E-14, 6.5529431722547864, -3.5474774874272002E-15,
    8.3106024468525048E-14, 6.4178448252756946, -3.5624285846816805E-15,
    8.3396767734893869E-14, 6.2827467474500347, -3.57737968193616E-15,
    8.3687511001262677E-14, 6.1476495448794868, -3.5923284737014719E-15,
    8.3978209434436984E-14, 6.0125526786469861, -3.6072772654667842E-15,
    8.4268907867611291E-14, 5.877456182378558, -3.6222260572320965E-15,
    8.4559606300785611E-14, 5.7423600897131939, -3.6371748489974081E-15,
    8.4850304733959917E-14, 5.6072644342898856, -3.65212364076272E-15,
    8.5141003167134237E-14, 5.4721692497476253, -3.667072432528032E-15,
    8.5431701600308544E-14, 5.3370745697124375, -3.6820212242933443E-15,
    8.5722400033482851E-14, -3.63595009278392E-15, 7.5097794622411449,
    -2.75309671447575E-13, -3.57161411034296E-15, 7.5060650778225018,
    -2.73001837403655E-13, -3.5072284859402396E-15, 7.5023478307688318,
    -2.7069222262358997E-13, -3.44284286153752E-15, 207.49863058380075,
    -2.68382607843525E-13, -3.3784572371347998E-15, 7.3635312953318808,
    -2.6607299306346E-13, -3.31412125469384E-15, 7.2284349028802772,
    -2.6376515901954E-13, -3.2497356302911203E-15, 7.0933357150560647,
    -2.61455544239475E-13, -3.1853500058883997E-15, 6.9582366282136716,
    -2.5914592945941E-13, -3.12096438148568E-15, 6.823137675992089,
    -2.56836314679345E-13, -3.0565787570829603E-15, 6.6880388920302911,
    -2.5452669989928E-13, -2.992242774642E-15, 6.5529431722547864,
    -2.5221886585536E-13, -2.9278571502392796E-15, 6.4178448252756946,
    -2.4990925107529503E-13, -2.86347152583656E-15, 6.2827467474500347,
    -2.4759963629523E-13, -2.7990958298261923E-15, 6.1476495448794868,
    -2.45290377662394E-13, -2.7347201338158239E-15, 6.0125526786469861,
    -2.42981119029558E-13, -2.6703444378054559E-15, 5.877456182378558,
    -2.4067186039672203E-13, -2.6059687417950883E-15, 5.7423600897131939,
    -2.38362601763886E-13, -2.5415930457847204E-15, 5.6072644342898856,
    -2.3605334313105E-13, -2.477217349774352E-15, 5.4721692497476253,
    -2.33744084498214E-13, -2.412841653763984E-15, 5.3370745697124375,
    -2.31434825865378E-13, 8.6488748679454459E-14, -2.7053866027977495E-13,
    7.5097794622411449, 8.4683850876970761E-14, -2.69820345637775E-13,
    7.5060650778225018, 8.2877560406429563E-14, -2.6910147674065E-13,
    7.5023478307688318, 8.1071269935888391E-14, -2.68382607843525E-13,
    207.49863058380075, 7.92649794653472E-14, -2.676637389464E-13,
    7.3635312953318808, 7.74600816628635E-14, -2.669454243044E-13,
    7.2284349028802772, 7.5653791192322311E-14, -2.66226555407275E-13,
    7.0933357150560647, 7.3847500721781126E-14, -2.6550768651015003E-13,
    6.9582366282136716, 7.2041210251239928E-14, -2.64788817613025E-13,
    6.823137675992089, 7.0234919780698756E-14, -2.640699487159E-13,
    6.6880388920302911, 6.8430021978215045E-14, -2.633516340739E-13,
    6.5529431722547864, 6.662373150767386E-14, -2.62632765176775E-13,
    6.4178448252756946, 6.4817441037132675E-14, -2.6191389627965E-13,
    6.2827467474500347, 6.3011429100202983E-14, -2.6119513823355003E-13,
    6.1476495448794868, 6.120541716327329E-14, -2.6047638018745E-13,
    6.0125526786469861, 5.939940522634361E-14, -2.5975762214135E-13,
    5.877456182378558, 5.7593393289413918E-14, -2.5903886409525005E-13,
    5.7423600897131939, 5.5787381352484225E-14, -2.5832010604915E-13,
    5.6072644342898856, 5.3981369415554539E-14, -2.5760134800305E-13,
    5.4721692497476253, 5.217535747862484E-14, -2.5688258995695E-13,
    5.3370745697124375, 7.3744752259570951, -3.3346947698108E-15,
    7.8393917430011212E-14, 7.3708291226261844, -3.3492747589244E-15,
    7.86841221595216E-14, 7.3671802089388283, -3.3638659980295997E-15,
    7.8974550812434409E-14, 7.3635312953318808, -3.3784572371348E-15,
    7.92649794653472E-14, 207.35988238181054, -3.39304847624E-15,
    7.955540811826E-14, 7.22485423681804, -3.4076284653535998E-15,
    7.98456128477704E-14, 7.0898233150977212, -3.4222197044588003E-15,
    8.01360415006832E-14, 6.9547924607202329, -3.436810943564E-15,
    8.042647015359601E-14, 6.8197617073245622, -3.4514021826691997E-15,
    8.07168988065088E-14, 6.6847310885497029, -3.4659934217744E-15,
    8.10073274594216E-14, 6.5497034480642924, -3.480573410888E-15,
    8.1297532188932E-14, 6.4146731990200738, -3.4951646499931998E-15,
    8.15879608418448E-14, 6.2796431854877053, -3.5097558890984002E-15,
    8.1878389494757612E-14, 6.1446140031147012, -3.5243448782052802E-15,
    8.216877336298992E-14, 6.0095851234537214, -3.53893386731216E-15,
    8.245915723122224E-14, 5.87455658013079, -3.55352285641904E-15,
    8.2749541099454561E-14, 5.73952840677193, -3.56811184552592E-15,
    8.30399249676869E-14, 5.6045006370161339, -3.5827008346328E-15,
    8.33303088359192E-14, 5.4694733045023947, -3.59728982373968E-15,
    8.3620692704151535E-14, 5.3344464428697025, -3.61187881284656E-15,
    8.3911076572383842E-14, -3.651980487032E-15, 7.3744752259570951,
    -2.745797785044E-13, -3.5872849238959995E-15, 7.3708291226261844,
    -2.7227561725320004E-13, -3.5225394413439996E-15, 7.3671802089388283,
    -2.699696780998E-13, -3.457793958792E-15, 7.3635312953318808,
    -2.676637389464E-13, -3.39304847624E-15, 207.35988238181054,
    -2.65357799793E-13, -3.328352913104E-15, 7.22485423681804,
    -2.630536385418E-13, -3.263607430552E-15, 7.0898233150977212,
    -2.607476993884E-13, -3.198861948E-15, 6.9547924607202329,
    -2.58441760235E-13, -3.134116465448E-15, 6.8197617073245622,
    -2.561358210816E-13, -3.069370982896E-15, 6.6847310885497029,
    -2.538298819282E-13, -3.0046754197600003E-15, 6.5497034480642924,
    -2.5152572067700003E-13, -2.9399299372079996E-15, 6.4146731990200738,
    -2.4921978152360003E-13, -2.8751844546559997E-15, 6.2796431854877053,
    -2.469138423702E-13, -2.8104489559872E-15, 6.1446140031147012,
    -2.4460825879724E-13, -2.7457134573184E-15, 6.0095851234537214,
    -2.4230267522428E-13, -2.6809779586496E-15, 5.87455658013079,
    -2.3999709165132003E-13, -2.6162424599808003E-15, 5.73952840677193,
    -2.3769150807836E-13, -2.5515069613120004E-15, 5.6045006370161339,
    -2.353859245054E-13, -2.4867714626432E-15, 5.4694733045023947,
    -2.3308034093244E-13, -2.4220359639744E-15, 5.3344464428697025,
    -2.3077475735948003E-13, 8.6780435543621191E-14, -2.6821802145366E-13,
    7.3744752259570951, 8.49752233702516E-14, -2.6750337960438E-13,
    7.3708291226261844, 8.31686182862544E-14, -2.6678818633392003E-13,
    7.3671802089388283, 8.13620132022572E-14, -2.6607299306346E-13,
    7.3635312953318808, 7.955540811826E-14, -2.65357799793E-13,
    207.35988238181054, 7.77501959448904E-14, -2.6464315794372E-13,
    7.22485423681804, 7.5943590860893208E-14, -2.6392796467326E-13,
    7.0898233150977212, 7.4136985776896E-14, -2.6321277140280006E-13,
    6.9547924607202329, 7.2330380692898791E-14, -2.6249757813234E-13,
    6.8197617073245622, 7.05237756089016E-14, -2.6178238486188E-13,
    6.6847310885497029, 6.8718563435532E-14, -2.610677430126E-13,
    6.5497034480642924, 6.69119583515348E-14, -2.6035254974214E-13,
    6.4146731990200738, 6.51053532675376E-14, -2.5963735647168E-13,
    6.2796431854877053, 6.3299026765665925E-14, -2.58922273485456E-13,
    6.1446140031147012, 6.1492700263794246E-14, -2.5820719049923196E-13,
    6.0095851234537214, 5.9686373761922567E-14, -2.57492107513008E-13,
    5.87455658013079, 5.7880047260050887E-14, -2.56777024526784E-13,
    5.73952840677193, 5.60737207581792E-14, -2.5606194154056003E-13,
    5.6045006370161339, 5.4267394256307523E-14, -2.5534685855433603E-13,
    5.4694733045023947, 5.2461067754435837E-14, -2.5463177556811204E-13,
    5.3344464428697025, 7.2391741433404793, -3.27143725221464E-15,
    7.6589962497802053E-14, 7.235596235240787, -3.28565793787352E-15,
    7.6879853098809667E-14, 7.2320155690229218, -3.29988959628368E-15,
    7.7169967380836587E-14, 7.2284349028802755, -3.31412125469384E-15,
    7.74600816628635E-14, 7.2248542368180395, -3.3283529131039996E-15,
    7.7750195944890391E-14, 207.22127632865352, -3.34257359876288E-15,
    7.8040086545898018E-14, 7.08631362077408, -3.35680525717304E-15,
    7.8330200827924926E-14, 6.9513509465984713, -3.3710369155832003E-15,
    7.8620315109951847E-14, 6.8163883397656928, -3.38526857399336E-15,
    7.8910429391978755E-14, 6.6814258339147328, -3.39950023240352E-15,
    7.9200543674005676E-14, 6.5464662204967077, -3.4137209180624003E-15,
    7.949043427501329E-14, 6.4115040171243436, -3.42795257647256E-15,
    7.97805485570402E-14, 6.2765420156248357, -3.44218423488272E-15,
    8.0070662839067119E-14, 6.141580801199602, -3.4564136987426242E-15,
    8.0360732384890163E-14, 6.0066198558474007, -3.470643162602528E-15,
    8.0650801930713207E-14, 5.8716592132072236, -3.4848726264624318E-15,
    8.0940871476536264E-14, 5.7366989069050947, -3.4991020903223357E-15,
    8.1230941022359321E-14, 5.6017389705670366, -3.5133315541822395E-15,
    8.1521010568182378E-14, 5.4667794378320433, -3.5275610180421442E-15,
    8.1811080114005422E-14, 5.331820342339106, -3.541790481902048E-15,
    8.2101149659828467E-14, -3.66799852168544E-15, 7.2391741433404793,
    -2.7385044831599997E-13, -3.6029436550947196E-15, 7.235596235240787,
    -2.7154995702576E-13, -3.5378385918476796E-15, 7.2320155690229218,
    -2.6924769066508E-13, -3.47273352860064E-15, 7.2284349028802755,
    -2.669454243044E-13, -3.4076284653535998E-15, 7.2248542368180395,
    -2.6464315794372003E-13, -3.34257359876288E-15, 207.22127632865352,
    -2.6234266665348004E-13, -3.27746853551584E-15, 7.08631362077408,
    -2.600404002928E-13, -3.2123634722687996E-15, 6.9513509465984713,
    -2.5773813393212E-13, -3.14725840902176E-15, 6.8163883397656928,
    -2.5543586757144E-13, -3.08215334577472E-15, 6.6814258339147328,
    -2.5313360121076004E-13, -3.0170984791840002E-15, 6.5464662204967077,
    -2.5083310992052004E-13, -2.9519934159369595E-15, 6.4115040171243436,
    -2.4853084355984003E-13, -2.8868883526899196E-15, 6.2765420156248357,
    -2.4622857719916E-13, -2.8217933287741443E-15, 6.141580801199602,
    -2.43926665852568E-13, -2.7566983048583678E-15, 6.0066198558474007,
    -2.41624754505976E-13, -2.6916032809425921E-15, 5.8716592132072236,
    -2.39322843159384E-13, -2.6265082570268161E-15, 5.7366989069050947,
    -2.37020931812792E-13, -2.5614132331110404E-15, 5.6017389705670366,
    -2.347190204662E-13, -2.4963182091952639E-15, 5.4667794378320433,
    -2.32417109119608E-13, -2.4312231852794882E-15, 5.331820342339106,
    -2.30115197773016E-13, 8.7071897514292043E-14, -2.6589917186334E-13,
    7.2391741433404793, 8.5266371212419671E-14, -2.6518819997502E-13,
    7.235596235240787, 8.3459451757536573E-14, -2.6447667949728E-13,
    7.2320155690229218, 8.1652532302653488E-14, -2.6376515901954E-13,
    7.2284349028802755, 7.98456128477704E-14, -2.630536385418E-13,
    7.2248542368180395, 7.8040086545898018E-14, -2.6234266665348E-13,
    207.22127632865352, 7.6233167091014945E-14, -2.6163114617574E-13,
    7.08631362077408, 7.4426247636131847E-14, -2.6091962569800003E-13,
    6.9513509465984713, 7.261932818124875E-14, -2.6020810522026E-13,
    6.8163883397656928, 7.0812408726365664E-14, -2.5949658474252E-13,
    6.6814258339147328, 6.900688242449328E-14, -2.5878561285420005E-13,
    6.5464662204967077, 6.71999629696102E-14, -2.5807409237646E-13,
    6.4115040171243436, 6.5393043514727109E-14, -2.5736257189872E-13,
    6.2765420156248357, 6.3586402690446154E-14, -2.56651161138864E-13,
    6.141580801199602, 6.1779761866165211E-14, -2.55939750379008E-13,
    6.0066198558474007, 5.9973121041884269E-14, -2.5522833961915203E-13,
    5.8716592132072236, 5.8166480217603326E-14, -2.54516928859296E-13,
    5.7366989069050947, 5.6359839393322371E-14, -2.5380551809944E-13,
    5.6017389705670366, 5.4553198569041428E-14, -2.53094107339584E-13,
    5.4667794378320433, 5.2746557744760467E-14, -2.52382696579728E-13,
    5.331820342339106, 7.1038702094572042, -3.20813092480552E-15,
    7.4784615625058791E-14, 7.100360515193219, -3.2219920297693597E-15,
    7.5074191855180516E-14, 7.0968481150896263, -3.23586383003024E-15,
    7.5363991523751413E-14, 7.0933357150560639, -3.24973563029112E-15,
    7.5653791192322311E-14, 7.0898233150977212, -3.2636074305519998E-15,
    7.59435908608932E-14, 7.08631362077408, -3.27746853551584E-15,
    7.6233167091014933E-14, 207.08280122060566, -3.29134033577672E-15,
    7.6522966759585817E-14, 6.9479067789105216, -3.3052121360376E-15,
    7.6812766428156727E-14, 6.8130123709192159, -3.3190839362984796E-15,
    7.7102566096727612E-14, 6.67811803027074, -3.33295573655936E-15,
    7.7392365765298509E-14, 6.5432264961583746, -3.3468168415232E-15,
    7.7681941995420247E-14, 6.4083323907364491, -3.36068864178408E-15,
    7.7971741663991144E-14, 6.2734384535483878, -3.37456044204496E-15,
    7.8261541332562041E-14, 6.138545259344041, -3.3884301032464321E-15,
    7.85512963134431E-14, 6.003652300576328, -3.402299764447904E-15,
    7.8841051294324157E-14, 5.8687596108816473, -3.4161694256493759E-15,
    7.9130806275205227E-14, 5.733867223898991, -3.4300390868508477E-15,
    7.9420561256086285E-14, 5.5989751732543818, -3.4439087480523196E-15,
    7.9710316236967343E-14, 5.4640834925738444, -3.4577784092537919E-15,
    8.0000071217848413E-14, 5.3291922154963718, -3.4716480704552641E-15,
    8.0289826198729471E-14, -3.68402891593352E-15, 7.1038702094572042,
    -2.7312055537282496E-13, -3.61861446864776E-15, 7.100360515193219,
    -2.7082373687530505E-13, -3.55314954725144E-15, 7.0968481150896263,
    -2.6852514614129E-13, -3.48768462585512E-15, 7.0933357150560639,
    -2.66226555407275E-13, -3.4222197044588E-15, 7.0898233150977212,
    -2.6392796467326E-13, -3.35680525717304E-15, 7.08631362077408,
    -2.6163114617574003E-13, -3.2913403357767202E-15, 207.08280122060566,
    -2.59332555441725E-13, -3.2258754143803998E-15, 6.9479067789105216,
    -2.5703396470771E-13, -3.16041049298408E-15, 6.8130123709192159,
    -2.54735373973695E-13, -3.09494557158776E-15, 6.67811803027074,
    -2.5243678323968E-13, -3.029531124302E-15, 6.5432264961583746,
    -2.5013996474216E-13, -2.9640662029056796E-15, 6.4083323907364491,
    -2.4784137400814503E-13, -2.89860128150936E-15, 6.2734384535483878,
    -2.4554278327413E-13, -2.8331464549351522E-15, 6.138545259344041,
    -2.43244546987414E-13, -2.7676916283609438E-15, 6.003652300576328,
    -2.4094631070069803E-13, -2.7022368017867361E-15, 5.8687596108816473,
    -2.38648074413982E-13, -2.6367819752125281E-15, 5.733867223898991,
    -2.36349838127266E-13, -2.5713271486383204E-15, 5.5989751732543818,
    -2.3405160184055E-13, -2.505872322064112E-15, 5.4640834925738444,
    -2.31753365553834E-13, -2.4404174954899043E-15, 5.3291922154963718,
    -2.29455129267118E-13, 8.7363584378458775E-14, -2.6357853303722497E-13,
    7.1038702094572042, 8.5557743705700526E-14, -2.62871233941625E-13,
    7.100360515193219, 8.37505096373614E-14, -2.6216338909055003E-13,
    7.0968481150896263, 8.1943275569022309E-14, -2.6145554423947497E-13,
    7.0933357150560639, 8.0136041500683213E-14, -2.607476993884E-13,
    7.0898233150977212, 7.8330200827924939E-14, -2.6004040029279997E-13,
    7.08631362077408, 7.6522966759585843E-14, -2.59332555441725E-13,
    207.08280122060566, 7.4715732691246721E-14, -2.5862471059065E-13,
    6.9479067789105216, 7.2908498622907613E-14, -2.57916865739575E-13,
    6.8130123709192159, 7.1101264554568517E-14, -2.572090208885E-13,
    6.67811803027074, 6.9295423881810242E-14, -2.565017217929E-13,
    6.5432264961583746, 6.7488189813471134E-14, -2.55793876941825E-13,
    6.4083323907364491, 6.5680955745132038E-14, -2.5508603209075E-13,
    6.2734384535483878, 6.38740003559091E-14, -2.5437829639077E-13,
    6.138545259344041, 6.2067044966686167E-14, -2.5367056069079E-13,
    6.003652300576328, 6.0260089577463225E-14, -2.5296282499081E-13,
    5.8687596108816473, 5.84531341882403E-14, -2.5225508929083004E-13,
    5.733867223898991, 5.6646178799017348E-14, -2.5154735359085E-13,
    5.5989751732543818, 5.4839223409794412E-14, -2.5083961789087E-13,
    5.4640834925738444, 5.3032268020571464E-14, -2.5013188219089E-13,
    5.3291922154963718, 6.9685664774727263, -3.1448245973963997E-15,
    7.297926875231553E-14, 6.9651249634054349, -3.1583261216652E-15,
    7.3268530611551365E-14, 6.9616807957771414, -3.1718380637768E-15,
    7.3558015666666252E-14, 6.9582366282136707, -3.1853500058884E-15,
    7.3847500721781126E-14, 6.9547924607202321, -3.198861948E-15,
    7.4136985776896E-14, 6.9513509465984713, -3.2123634722688E-15,
    7.4426247636131847E-14, 6.9479067789105216, -3.2258754143804E-15,
    7.4715732691246721E-14, 206.94446261128223, -3.2393873564920003E-15,
    7.5005217746361608E-14, 6.8096364021272082, -3.2528992986036E-15,
    7.5294702801476482E-14, 6.6748102266760263, -3.2664112407152E-15,
    7.5584187856591368E-14, 6.5399867718641325, -3.2799127649840003E-15,
    7.58734497158272E-14, 6.4051607643874586, -3.2934247070956E-15,
    7.616293477094209E-14, 6.270334891505656, -3.3069366492072E-15,
    7.6452419826056964E-14, 6.1355097175170092, -3.32044650775024E-15,
    7.6741860241996035E-14, 6.0006847453311911, -3.33395636629328E-15,
    7.7031300657935106E-14, 5.8658600085820067, -3.3474662248363203E-15,
    7.7320741073874178E-14, 5.7310355409058547, -3.3609760833793598E-15,
    7.7610181489813249E-14, 5.596211375941726, -3.3744859419223997E-15,
    7.7899621905752333E-14, 5.4613875473156464, -3.38799580046544E-15,
    7.81890623216914E-14, 5.3265640886536376, -3.40150565900848E-15,
    7.8478502737630462E-14, -3.7000593101816E-15, 6.9685664774727263,
    -2.7239066242964996E-13, -3.6342852822008E-15, 6.9651249634054349,
    -2.7009751672485003E-13, -3.5684605026552E-15, 6.9616807957771414,
    -2.678026016175E-13, -3.5026357231096E-15, 6.9582366282136707,
    -2.6550768651015E-13, -3.436810943564E-15, 6.9547924607202321,
    -2.632127714028E-13, -3.3710369155832E-15, 6.9513509465984713,
    -2.6091962569800003E-13, -3.3052121360376E-15, 6.9479067789105216,
    -2.5862471059065E-13, -3.239387356492E-15, 206.94446261128223,
    -2.563297954833E-13, -3.1735625769464E-15, 6.8096364021272082,
    -2.5403488037595E-13, -3.1077377974008002E-15, 6.6748102266760263,
    -2.5173996526860004E-13, -3.0419637694200002E-15, 6.5399867718641325,
    -2.494468195638E-13, -2.9761389898743996E-15, 6.4051607643874586,
    -2.4715190445645004E-13, -2.9103142103287998E-15, 6.270334891505656,
    -2.448569893491E-13, -2.84449958109616E-15, 6.1355097175170092,
    -2.4256242812226E-13, -2.77868495186352E-15, 6.0006847453311911,
    -2.4026786689542E-13, -2.71287032263088E-15, 5.8658600085820067,
    -2.3797330566858E-13, -2.64705569339824E-15, 5.7310355409058547,
    -2.3567874444174003E-13, -2.5812410641656004E-15, 5.596211375941726,
    -2.333841832149E-13, -2.51542643493296E-15, 5.4613875473156464,
    -2.3108962198806E-13, -2.4496118057003204E-15, 5.3265640886536376,
    -2.2879506076122003E-13, 8.765527124262552E-14, -2.6125789421110996E-13,
    6.9685664774727263, 8.5849116198981369E-14, -2.6055426790823E-13,
    6.9651249634054349, 8.4041567517186237E-14, -2.5985009868382004E-13,
    6.9616807957771414, 8.2234018835391117E-14, -2.5914592945941E-13,
    6.9582366282136707, 8.042647015359601E-14, -2.58441760235E-13,
    6.9547924607202321, 7.8620315109951847E-14, -2.5773813393212E-13,
    6.9513509465984713, 7.6812766428156727E-14, -2.5703396470771E-13,
    6.9479067789105216, 7.5005217746361608E-14, -2.5632979548330004E-13,
    206.94446261128223, 7.3197669064566476E-14, -2.5562562625889E-13,
    6.8096364021272082, 7.1390120382771356E-14, -2.5492145703448E-13,
    6.6748102266760263, 6.95839653391272E-14, -2.5421783073160005E-13,
    6.5399867718641325, 6.7776416657332086E-14, -2.5351366150719E-13,
    6.4051607643874586, 6.5968867975536966E-14, -2.5280949228278003E-13,
    6.270334891505656, 6.4161598021372025E-14, -2.52105431642676E-13,
    6.1355097175170092, 6.235432806720711E-14, -2.51401371002572E-13,
    6.0006847453311911, 6.0547058113042181E-14, -2.5069731036246805E-13,
    5.8658600085820067, 5.8739788158877253E-14, -2.49993249722364E-13,
    5.7310355409058547, 5.6932518204712325E-14, -2.4928918908226E-13,
    5.596211375941726, 5.5125248250547397E-14, -2.4858512844215603E-13,
    5.4613875473156464, 5.3317978296382456E-14, -2.4788106780205204E-13,
    5.3265640886536376, 6.8332629810260386, -3.08151826998728E-15,
    7.1173921879572269E-14, 6.8298896135164471, -3.0946602135610397E-15,
    7.1462869367922215E-14, 6.8265136447244394, -3.10781229752336E-15,
    7.1752039809581078E-14, 6.823137675992089, -3.12096438148568E-15,
    7.2041210251239941E-14, 6.8197617073245613, -3.1341164654479996E-15,
    7.23303806928988E-14, 6.8163883397656937, -3.14725840902176E-15,
    7.261932818124875E-14, 6.8130123709192159, -3.16041049298408E-15,
    7.2908498622907613E-14, 6.8096364021272091, -3.1735625769464E-15,
    7.3197669064566488E-14, 206.80626043339487, -3.1867146609087197E-15,
    7.3486839506225339E-14, 6.6715024231357818, -3.19986674487104E-15,
    7.3776009947884215E-14, 6.5367470476191691, -3.2130086884448E-15,
    7.4064957436234173E-14, 6.40198913808256, -3.22616077240712E-15,
    7.4354127877893023E-14, 6.2672313295018292, -3.23931285636944E-15,
    7.46432983195519E-14, 6.1324741757236945, -3.252462912254048E-15,
    7.4932424170548971E-14, 5.9977171901145843, -3.265612968138656E-15,
    7.5221550021546056E-14, 5.8629604063083027, -3.2787630240232639E-15,
    7.5510675872543141E-14, 5.7282038579386549, -3.2919130799078718E-15,
    7.5799801723540213E-14, 5.5934475786420386, -3.3050631357924797E-15,
    7.60889275745373E-14, 5.4586916020574474, -3.3182131916770881E-15,
    7.6378053425534382E-14, 5.3239359618109034, -3.331363247561696E-15,
    7.6667179276531467E-14, -3.7160897044296804E-15, 6.8332629810260386,
    -2.71660769486475E-13, -3.64995609575384E-15, 6.8298896135164471,
    -2.69371296574395E-13, -3.58377145805896E-15, 6.8265136447244394,
    -2.6708005709371E-13, -3.51758682036408E-15, 6.823137675992089,
    -2.64788817613025E-13, -3.4514021826691997E-15, 6.8197617073245613,
    -2.6249757813234E-13, -3.38526857399336E-15, 6.8163883397656937,
    -2.6020810522026003E-13, -3.31908393629848E-15, 6.8130123709192159,
    -2.5791686573957505E-13, -3.2528992986035996E-15, 6.8096364021272091,
    -2.5562562625888997E-13, -3.18671466090872E-15, 206.80626043339487,
    -2.53334386778205E-13, -3.12053002321384E-15, 6.6715024231357818,
    -2.5104314729752E-13, -3.054396414538E-15, 6.5367470476191691,
    -2.4875367438544E-13, -2.9882117768431196E-15, 6.40198913808256,
    -2.4646243490475504E-13, -2.9220271391482397E-15, 6.2672313295018292,
    -2.4417119542407E-13, -2.8558527072571681E-15, 6.1324741757236945,
    -2.41880309257106E-13, -2.7896782753660961E-15, 5.9977171901145843,
    -2.39589423090142E-13, -2.7235038434750241E-15, 5.8629604063083027,
    -2.37298536923178E-13, -2.6573294115839521E-15, 5.7282038579386549,
    -2.35007650756214E-13, -2.5911549796928805E-15, 5.5934475786420386,
    -2.3271676458925E-13, -2.5249805478018081E-15, 5.4586916020574474,
    -2.30425878422286E-13, -2.4588061159107361E-15, 5.3239359618109034,
    -2.28134992255322E-13, 8.7946958106792252E-14, -2.58937255384995E-13,
    6.8332629810260386, 8.6140488692262212E-14, -2.58237301874835E-13,
    6.8298896135164471, 8.4332625397011069E-14, -2.5753680827709E-13,
    6.8265136447244394, 8.2524762101759938E-14, -2.56836314679345E-13,
    6.823137675992089, 8.0716898806508808E-14, -2.561358210816E-13,
    6.8197617073245613, 7.8910429391978755E-14, -2.5543586757144E-13,
    6.8163883397656937, 7.7102566096727625E-14, -2.5473537397369503E-13,
    6.8130123709192159, 7.5294702801476482E-14, -2.5403488037595E-13,
    6.8096364021272091, 7.3486839506225339E-14, -2.53334386778205E-13,
    206.80626043339487, 7.1678976210974209E-14, -2.5263389318046003E-13,
    6.6715024231357818, 6.9872506796444156E-14, -2.519339396703E-13,
    6.5367470476191691, 6.8064643501193025E-14, -2.51233446072555E-13,
    6.40198913808256, 6.62567802059419E-14, -2.5053295247481E-13,
    6.2672313295018292, 6.4449195686834967E-14, -2.49832566894582E-13,
    6.1324741757236945, 6.2641611167728053E-14, -2.49132181314354E-13,
    5.9977171901145843, 6.0834026648621138E-14, -2.4843179573412604E-13,
    5.8629604063083027, 5.9026442129514223E-14, -2.4773141015389804E-13,
    5.7282038579386549, 5.72188576104073E-14, -2.4703102457367003E-13,
    5.5934475786420386, 5.5411273091300381E-14, -2.46330638993442E-13,
    5.4586916020574474, 5.3603688572193453E-14, -2.45630253413214E-13,
    5.3239359618109034, 6.6979597537587221, -3.0182119425781598E-15,
    6.9368575006829E-14, 6.6946544991652495, -3.03099430545688E-15,
    6.9657208124293051E-14, 6.6913466955705356, -3.0437865312699197E-15,
    6.9946063952495916E-14, 6.6880388920302911, -3.05657875708296E-15,
    7.0234919780698756E-14, 6.684731088549702, -3.0693709828959997E-15,
    7.05237756089016E-14, 6.6814258339147337, -3.0821533457747198E-15,
    7.0812408726365664E-14, 6.678118030270741, -3.09494557158776E-15,
    7.11012645545685E-14, 6.6748102266760263, -3.1077377974008002E-15,
    7.1390120382771369E-14, 6.6715024231357818, -3.1205300232138397E-15,
    7.1678976210974209E-14, 206.6681946196552, -3.13332224902688E-15,
    7.1967832039177061E-14, 6.533507323428676, -3.1461046119056004E-15,
    7.225646515664113E-14, 6.3988175118269393, -3.1588968377186398E-15,
    7.2545320984843969E-14, 6.2641277675420923, -3.17168906353168E-15,
    7.2834176813046822E-14, 6.1294386339692837, -3.184479316757856E-15,
    7.3122988099101907E-14, 5.9947496349316935, -3.1972695699840319E-15,
    7.3411799385157E-14, 5.8600608040631279, -3.2100598232102079E-15,
    7.3700610671212091E-14, 5.7253721749973909, -3.2228500764363839E-15,
    7.3989421957267189E-14, 5.5906837813682877, -3.23564032966256E-15,
    7.4278233243322275E-14, 5.4559956568122159, -3.2484305828887358E-15,
    7.4567044529377373E-14, 5.3213078349681693, -3.2612208361149117E-15,
    7.4855855815432458E-14, -3.73212009867776E-15, 6.6979597537587221,
    -2.709308765433E-13, -3.66562690930688E-15, 6.6946544991652495,
    -2.6864507642394003E-13, -3.5990824134627196E-15, 6.6913466955705356,
    -2.6635751256992E-13, -3.5325379176185603E-15, 6.6880388920302911,
    -2.640699487159E-13, -3.4659934217743997E-15, 6.684731088549702,
    -2.6178238486188004E-13, -3.39950023240352E-15, 6.6814258339147337,
    -2.5949658474252E-13, -3.33295573655936E-15, 6.678118030270741,
    -2.5720902088850004E-13, -3.2664112407151998E-15, 6.6748102266760263,
    -2.5492145703448E-13, -3.19986674487104E-15, 6.6715024231357818,
    -2.5263389318046003E-13, -3.13332224902688E-15, 206.6681946196552,
    -2.5034632932644E-13, -3.066829059656E-15, 6.533507323428676,
    -2.4806052920708003E-13, -3.0002845638118397E-15, 6.3988175118269393,
    -2.4577296535306005E-13, -2.9337400679676795E-15, 6.2641277675420923,
    -2.4348540149903997E-13, -2.8672058334181764E-15, 6.1294386339692837,
    -2.41198190391952E-13, -2.800671598868672E-15, 5.9947496349316935,
    -2.3891097928486403E-13, -2.734137364319168E-15, 5.8600608040631279,
    -2.36623768177776E-13, -2.6676031297696641E-15, 5.7253721749973909,
    -2.34336557070688E-13, -2.6010688952201605E-15, 5.5906837813682877,
    -2.3204934596359997E-13, -2.5345346606706562E-15, 5.4559956568122159,
    -2.29762134856512E-13, -2.4680004261211522E-15, 5.3213078349681693,
    -2.2747492374942404E-13, 8.8238644970958984E-14, -2.5661661655888E-13,
    6.6979597537587221, 8.6431861185543067E-14, -2.5592033584144E-13,
    6.6946544991652495, 8.46236832768359E-14, -2.5522351787036E-13,
    6.6913466955705356, 8.2815505368128759E-14, -2.5452669989928E-13,
    6.6880388920302911, 8.10073274594216E-14, -2.538298819282E-13,
    6.684731088549702, 7.9200543674005676E-14, -2.5313360121076E-13,
    6.6814258339147337, 7.7392365765298522E-14, -2.5243678323968E-13,
    6.678118030270741, 7.5584187856591368E-14, -2.5173996526860004E-13,
    6.6748102266760263, 7.37760099478842E-14, -2.5104314729752E-13,
    6.6715024231357818, 7.1967832039177061E-14, -2.5034632932644E-13,
    206.6681946196552, 7.0161048253761118E-14, -2.49650048609E-13,
    6.533507323428676, 6.8352870345053977E-14, -2.4895323063792E-13,
    6.3988175118269393, 6.6544692436346824E-14, -2.4825641266684E-13,
    6.2641277675420923, 6.473679335229791E-14, -2.47559702146488E-13,
    6.1294386339692837, 6.2928894268249008E-14, -2.4686299162613596E-13,
    5.9947496349316935, 6.11209951842001E-14, -2.4616628110578403E-13,
    5.8600608040631279, 5.9313096100151193E-14, -2.4546957058543206E-13,
    5.7253721749973909, 5.7505197016102279E-14, -2.4477286006508003E-13,
    5.5906837813682877, 5.5697297932053365E-14, -2.44076149544728E-13,
    5.4559956568122159, 5.3889398848004451E-14, -2.4337943902437603E-13,
    5.3213078349681693, 6.56265984833564, -2.954954424982E-15,
    6.7564620074619848E-14, 6.5594226207965853, -2.9673774844059997E-15,
    6.785293906358112E-14, 6.55618289649975, -2.979810129524E-15,
    6.81414805208981E-14, 6.5529431722547855, -2.992242774642E-15,
    6.8430021978215045E-14, 6.5497034480642924, -3.00467541976E-15,
    6.8718563435532E-14, 6.5464662204967086, -3.0170984791840002E-15,
    6.9006882424493292E-14, 6.5432264961583755, -3.029531124302E-15,
    6.9295423881810242E-14, 6.5399867718641325, -3.0419637694200002E-15,
    6.95839653391272E-14, 6.5367470476191691, -3.054396414538E-15,
    6.9872506796444168E-14, 6.5335073234286751, -3.066829059656E-15,
    7.0161048253761131E-14, 206.53027009586108, -3.07925211908E-15,
    7.0449367242722415E-14, 6.395648329931209, -3.091684764198E-15,
    7.0737908700039365E-14, 6.2610265976792228, -3.104117409316E-15,
    7.1026450157356328E-14, 6.126405432056778, -3.1165481372952003E-15,
    7.131494712100215E-14, 5.9917843673357476, -3.1289788652744E-15,
    7.1603444084647972E-14, 5.8571634371499366, -3.1414095932536E-15,
    7.18919410482938E-14, 5.7225426751331492, -3.1538403212328E-15,
    7.2180438011939629E-14, 5.58792211491919, -3.166271049212E-15,
    7.2468934975585438E-14, 5.4533017901418654, -3.1787017771912E-15,
    7.2757431939231273E-14, 5.3186817344375728, -3.1911325051704E-15,
    7.30459289028771E-14, -3.7481381333312E-15, 6.56265984833564,
    -2.702015463549E-13, -3.6812856405056E-15, 6.5594226207965853,
    -2.679194161965E-13, -3.6143815639664E-15, 6.55618289649975,
    -2.656355251352E-13, -3.5474774874272002E-15, 6.5529431722547855,
    -2.633516340739E-13, -3.4805734108879997E-15, 6.5497034480642924,
    -2.610677430126E-13, -3.4137209180624E-15, 6.5464662204967086,
    -2.587856128542E-13, -3.3468168415232E-15, 6.5432264961583755,
    -2.5650172179290005E-13, -3.279912764984E-15, 6.5399867718641325,
    -2.542178307316E-13, -3.2130086884448E-15, 6.5367470476191691,
    -2.519339396703E-13, -3.1461046119056E-15, 6.5335073234286751,
    -2.49650048609E-13, -3.07925211908E-15, 206.53027009586108,
    -2.4736791845060004E-13, -3.0123480425407996E-15, 6.395648329931209,
    -2.4508402738930004E-13, -2.9454439660015998E-15, 6.2610265976792228,
    -2.42800136328E-13, -2.87855020620512E-15, 6.126405432056778,
    -2.4051659744728E-13, -2.81165644640864E-15, 5.9917843673357476,
    -2.3823305856656003E-13, -2.74476268661216E-15, 5.8571634371499366,
    -2.3594951968584E-13, -2.67786892681568E-15, 5.7225426751331492,
    -2.3366598080512E-13, -2.6109751670192005E-15, 5.58792211491919,
    -2.313824419244E-13, -2.54408140722272E-15, 5.4533017901418654,
    -2.2909890304368E-13, -2.4771876474262403E-15, 5.3186817344375728,
    -2.2681536416296002E-13, 8.8530106941629836E-14, -2.5429776696856E-13,
    6.56265984833564, 8.6723009027711121E-14, -2.5360515621208E-13,
    6.5594226207965853, 8.4914516748118078E-14, -2.5291201103372E-13,
    6.55618289649975, 8.3106024468525048E-14, -2.5221886585536E-13,
    6.5529431722547855, 8.1297532188932E-14, -2.5152572067700003E-13,
    6.5497034480642924, 7.949043427501329E-14, -2.5083310992052E-13,
    6.5464662204967086, 7.7681941995420247E-14, -2.5013996474216E-13,
    6.5432264961583755, 7.58734497158272E-14, -2.494468195638E-13,
    6.5399867718641325, 7.406495743623416E-14, -2.4875367438544E-13,
    6.5367470476191691, 7.2256465156641117E-14, -2.4806052920708003E-13,
    6.5335073234286751, 7.04493672427224E-14, -2.4736791845060004E-13,
    206.53027009586108, 6.8640874963129359E-14, -2.4667477327224E-13,
    6.395648329931209, 6.6832382683536329E-14, -2.4598162809388E-13,
    6.2610265976792228, 6.5024169277078139E-14, -2.45288589799896E-13,
    6.126405432056778, 6.3215955870619974E-14, -2.44595551505912E-13,
    5.9917843673357476, 6.14077424641618E-14, -2.4390251321192805E-13,
    5.8571634371499366, 5.9599529057703632E-14, -2.4320947491794406E-13,
    5.7225426751331492, 5.7791315651245441E-14, -2.4251643662396E-13,
    5.58792211491919, 5.598310224478727E-14, -2.4182339832997603E-13,
    5.4533017901418654, 5.417488883832908E-14, -2.4113036003599204E-13,
    5.3186817344375728, 6.42735725981494, -2.89164809757288E-15,
    6.5759273201876587E-14, 6.424188077942552, -2.9037115763018395E-15,
    6.6047277819951969E-14, 6.4210164515831867, -2.91578436327056E-15,
    6.6335504663812921E-14, 6.4178448252756946, -2.92785715023928E-15,
    6.662373150767386E-14, 6.4146731990200738, -2.9399299372079996E-15,
    6.69119583515348E-14, 6.4115040171243445, -2.95199341593696E-15,
    6.71999629696102E-14, 6.4083323907364491, -2.96406620290568E-15,
    6.7488189813471134E-14, 6.4051607643874586, -2.9761389898744E-15,
    6.7776416657332086E-14, 6.40198913808256, -2.98821177684312E-15,
    6.8064643501193025E-14, 6.3988175118269393, -3.00028456381184E-15,
    6.8352870345053977E-14, 6.395648329931209, -3.0123480425408E-15,
    6.8640874963129372E-14, 206.39247670354331, -3.02442082950952E-15,
    6.8929101806990311E-14, 6.2579230356027757, -3.03649361647824E-15,
    6.9217328650851251E-14, 6.123369890201217, -3.0485645417990083E-15,
    6.9505511049555086E-14, 5.9888168120672685, -3.0606354671197761E-15,
    6.9793693448258922E-14, 5.8542638348347342, -3.072706392440544E-15,
    7.0081875846962757E-14, 5.7197109921374194, -3.0847773177613118E-15,
    7.0370058245666593E-14, 5.5851583176091282, -3.0968482430820797E-15,
    7.0658240644370416E-14, 5.4506058448836665, -3.1089191684028479E-15,
    7.0946423043074251E-14, 5.3160536075948377, -3.1209900937236158E-15,
    7.1234605441778087E-14, -3.76416852757928E-15, 6.42735725981494,
    -2.6947165341172497E-13, -3.6969564540586395E-15, 6.424188077942552,
    -2.6719319604604504E-13, -3.62969251937016E-15, 6.4210164515831867,
    -2.6491298061141E-13, -3.5624285846816805E-15, 6.4178448252756946,
    -2.62632765176775E-13, -3.4951646499931998E-15, 6.4146731990200738,
    -2.6035254974214003E-13, -3.42795257647256E-15, 6.4115040171243445,
    -2.5807409237646E-13, -3.36068864178408E-15, 6.4083323907364491,
    -2.5579387694182505E-13, -3.2934247070955996E-15, 6.4051607643874586,
    -2.5351366150719E-13, -3.22616077240712E-15, 6.40198913808256,
    -2.51233446072555E-13, -3.15889683771864E-15, 6.3988175118269393,
    -2.4895323063792003E-13, -3.0916847641980003E-15, 6.395648329931209,
    -2.4667477327224E-13, -3.0244208295095196E-15, 206.39247670354331,
    -2.4439455783760505E-13, -2.9571568948210397E-15, 6.2579230356027757,
    -2.4211434240297E-13, -2.8899033323661284E-15, 6.123369890201217,
    -2.39834478582126E-13, -2.8226497699112159E-15, 5.9888168120672685,
    -2.37554614761282E-13, -2.7553962074563043E-15, 5.8542638348347342,
    -2.35274750940438E-13, -2.6881426450013922E-15, 5.7197109921374194,
    -2.32994887119594E-13, -2.6208890825464805E-15, 5.5851583176091282,
    -2.3071502329875E-13, -2.5536355200915681E-15, 5.4506058448836665,
    -2.28435159477906E-13, -2.4863819576366564E-15, 5.3160536075948377,
    -2.26155295657062E-13, 8.8821793805796568E-14, -2.5197712814244497E-13,
    6.42735725981494, 8.7014381520991977E-14, -2.51288190178685E-13,
    6.424188077942552, 8.520557462794291E-14, -2.5059872062699E-13,
    6.4210164515831867, 8.3396767734893856E-14, -2.49909251075295E-13,
    6.4178448252756946, 8.1587960841844815E-14, -2.4921978152360003E-13,
    6.4146731990200738, 7.97805485570402E-14, -2.4853084355984E-13,
    6.4115040171243445, 7.7971741663991144E-14, -2.4784137400814503E-13,
    6.4083323907364491, 7.616293477094209E-14, -2.4715190445645004E-13,
    6.4051607643874586, 7.4354127877893023E-14, -2.46462434904755E-13,
    6.40198913808256, 7.2545320984843969E-14, -2.4577296535306E-13,
    6.3988175118269393, 7.0737908700039365E-14, -2.4508402738930004E-13,
    6.395648329931209, 6.8929101806990311E-14, -2.44394557837605E-13,
    206.39247670354331, 6.7120294913941257E-14, -2.4370508828591E-13,
    6.2579230356027757, 6.5311766942541081E-14, -2.43015725051802E-13,
    6.123369890201217, 6.3503238971140917E-14, -2.42326361817694E-13,
    5.9888168120672685, 6.1694710999740753E-14, -2.4163699858358604E-13,
    5.8542638348347342, 5.9886183028340589E-14, -2.4094763534947803E-13,
    5.7197109921374194, 5.8077655056940412E-14, -2.4025827211537E-13,
    5.5851583176091282, 5.6269127085540254E-14, -2.39568908881262E-13,
    5.4506058448836665, 5.4460599114140078E-14, -2.38879545647154E-13,
    5.3160536075948377, 6.292055041354276, -2.82834177016376E-15,
    6.3953926329133326E-14, 6.2889538715173421, -2.8400456681976797E-15,
    6.4241616576322818E-14, 6.2858503094642364, -2.85175859701712E-15,
    6.4529528806727747E-14, 6.2827467474500347, -2.86347152583656E-15,
    6.4817441037132675E-14, 6.2796431854877044, -2.8751844546559997E-15,
    6.51053532675376E-14, 6.2765420156248366, -2.88688835268992E-15,
    6.5393043514727109E-14, 6.2734384535483887, -2.89860128150936E-15,
    6.5680955745132025E-14, 6.2703348915056569, -2.9103142103288E-15,
    6.5968867975536966E-14, 6.2672313295018292, -2.9220271391482397E-15,
    6.6256780205941882E-14, 6.2641277675420923, -2.93374006796768E-15,
    6.6544692436346824E-14, 6.2610265976792236, -2.9454439660016E-15,
    6.6832382683536329E-14, 6.2579230356027757, -2.95715689482104E-15,
    6.7120294913941245E-14, 206.25481947356005, -2.96886982364048E-15,
    6.7408207144346186E-14, 6.1203343483741852, -2.9805809463028162E-15,
    6.7696074978108022E-14, 5.9858492568221315, -2.9922920689651521E-15,
    6.7983942811869871E-14, 5.851364232537688, -3.004003191627488E-15,
    6.8271810645631708E-14, 5.7168793091546579, -3.0157143142898239E-15,
    6.8559678479393557E-14, 5.5823945203068472, -3.0274254369521598E-15,
    6.8847546313155393E-14, 5.447909899628061, -3.039136559614496E-15,
    6.9135414146917242E-14, 5.3134254807521035, -3.0508476822768319E-15,
    6.9423281980679078E-14, -3.78019892182736E-15, 6.292055041354276,
    -2.6874176046854997E-13, -3.71262726761168E-15, 6.2889538715173421,
    -2.6646697589559E-13, -3.6450034747739196E-15, 6.2858503094642364,
    -2.6419043608762E-13, -3.57737968193616E-15, 6.2827467474500347,
    -2.6191389627965E-13, -3.5097558890984E-15, 6.2796431854877044,
    -2.5963735647168E-13, -3.44218423488272E-15, 6.2765420156248366,
    -2.5736257189872E-13, -3.3745604420449603E-15, 6.2734384535483887,
    -2.5508603209075004E-13, -3.3069366492071998E-15, 6.2703348915056569,
    -2.5280949228278E-13, -3.23931285636944E-15, 6.2672313295018292,
    -2.5053295247481E-13, -3.17168906353168E-15, 6.2641277675420923,
    -2.4825641266684E-13, -3.104117409316E-15, 6.2610265976792236,
    -2.4598162809388E-13, -3.0364936164782396E-15, 6.2579230356027757,
    -2.4370508828591006E-13, -2.9688698236404796E-15, 206.25481947356005,
    -2.4142854847794E-13, -2.9012564585271363E-15, 6.1203343483741852,
    -2.39152359716972E-13, -2.8336430934137919E-15, 5.9858492568221315,
    -2.36876170956004E-13, -2.7660297283004482E-15, 5.851364232537688,
    -2.3459998219503605E-13, -2.6984163631871042E-15, 5.7168793091546579,
    -2.3232379343406803E-13, -2.63080299807376E-15, 5.5823945203068472,
    -2.300476046731E-13, -2.5631896329604161E-15, 5.447909899628061,
    -2.27771415912132E-13, -2.4955762678470721E-15, 5.3134254807521035,
    -2.25495227151164E-13, 8.9113480669963313E-14, -2.4965648931632995E-13,
    6.292055041354276, 8.7305754014272819E-14, -2.4897122414529E-13,
    6.2889538715173421, 8.5496632507767742E-14, -2.4828543022026003E-13,
    6.2858503094642364, 8.3687511001262677E-14, -2.4759963629523E-13,
    6.2827467474500347, 8.1878389494757612E-14, -2.4691384237020004E-13,
    6.2796431854877044, 8.00706628390671E-14, -2.4622857719916E-13,
    6.2765420156248366, 7.8261541332562041E-14, -2.4554278327413E-13,
    6.2734384535483887, 7.6452419826056964E-14, -2.448569893491E-13,
    6.2703348915056569, 7.4643298319551887E-14, -2.4417119542407E-13,
    6.2672313295018292, 7.2834176813046822E-14, -2.4348540149904E-13,
    6.2641277675420923, 7.1026450157356316E-14, -2.4280013632800005E-13,
    6.2610265976792236, 6.9217328650851251E-14, -2.4211434240297E-13,
    6.2579230356027757, 6.7408207144346186E-14, -2.4142854847794E-13,
    206.25481947356005, 6.5599364608004023E-14, -2.40742860303708E-13,
    6.1203343483741852, 6.379052207166186E-14, -2.40057172129476E-13,
    5.9858492568221315, 6.1981679535319709E-14, -2.39371483955244E-13,
    5.851364232537688, 6.0172836998977559E-14, -2.3868579578101205E-13,
    5.7168793091546579, 5.83639944626354E-14, -2.3800010760678003E-13,
    5.5823945203068472, 5.6555151926293238E-14, -2.37314419432548E-13,
    5.447909899628061, 5.4746309389951075E-14, -2.3662873125831603E-13,
    5.3134254807521035, 6.1567538304019642, -2.7650452047172319E-15,
    6.2148857844496882E-14, 6.1537206285128008, -2.7763895775041757E-15,
    6.2436233769277111E-14, 6.1506850866831755, -2.787742703665184E-15,
    6.2723831434740053E-14, 6.1476495448794868, -2.7990958298261919E-15,
    6.3011429100202983E-14, 6.1446140031147012, -2.8104489559871998E-15,
    6.3299026765665912E-14, 6.1415808011996029, -2.8217933287741439E-15,
    6.3586402690446154E-14, 6.1385452593440419, -2.8331464549351522E-15,
    6.38740003559091E-14, 6.13550971751701, -2.84449958109616E-15,
    6.4161598021372038E-14, 6.1324741757236945, -2.8558527072571677E-15,
    6.4449195686834967E-14, 6.1294386339692837, -2.867205833418176E-15,
    6.473679335229791E-14, 6.126405432056778, -2.87855020620512E-15,
    6.5024169277078151E-14, 6.123369890201217, -2.889903332366128E-15,
    6.5311766942541081E-14, 6.1203343483741852, -2.9012564585271359E-15,
    6.5599364608004023E-14, 206.11729927454044, -2.9126078340133312E-15,
    6.5886917925330425E-14, 5.9828821591151327, -2.9239592094995266E-15,
    6.6174471242656814E-14, 5.8484650773236364, -2.9353105849857219E-15,
    6.6462024559983216E-14, 5.71404806279975, -2.9466619604719164E-15,
    6.6749577877309618E-14, 5.5796311491772785, -2.9580133359581121E-15,
    6.7037131194636007E-14, 5.4452143700900253, -2.969364711444307E-15,
    6.7324684511962409E-14, 5.3107977591717965, -2.9807160869305023E-15,
    6.76122378292888E-14, -3.7962268441565124E-15, 6.1567538304019642,
    -2.6801198007633E-13, -3.728295664693856E-15, 6.1537206285128008,
    -2.6574086772973803E-13, -3.6603120691976639E-15, 6.1506850866831755,
    -2.63468002981644E-13, -3.5923284737014719E-15, 6.1476495448794868,
    -2.6119513823355003E-13, -3.52434487820528E-15, 6.1446140031147012,
    -2.58922273485456E-13, -3.4564136987426238E-15, 6.1415808011996029,
    -2.5665116113886404E-13, -3.3884301032464321E-15, 6.1385452593440419,
    -2.5437829639077E-13, -3.3204465077502397E-15, 6.13550971751701,
    -2.52105431642676E-13, -3.252462912254048E-15, 6.1324741757236945,
    -2.49832566894582E-13, -3.184479316757856E-15, 6.1294386339692837,
    -2.4755970214648803E-13, -3.1165481372952003E-15, 6.126405432056778,
    -2.45288589799896E-13, -3.0485645417990075E-15, 6.123369890201217,
    -2.4301572505180205E-13, -2.9805809463028159E-15, 6.1203343483741852,
    -2.4074286030370797E-13, -2.9126078340133312E-15, 206.11729927454044,
    -2.3847034603591442E-13, -2.8446347217238462E-15, 5.9828821591151327,
    -2.3619783176812082E-13, -2.7766616094343616E-15, 5.8484650773236364,
    -2.3392531750032721E-13, -2.708688497144877E-15, 5.71404806279975,
    -2.3165280323253361E-13, -2.6407153848553924E-15, 5.5796311491772785,
    -2.2938028896474E-13, -2.5727422725659074E-15, 5.4452143700900253,
    -2.2710777469694641E-13, -2.5047691602764224E-15, 5.3107977591717965,
    -2.248352604291528E-13, 8.9405122555430869E-14, -2.47336208337374E-13,
    6.1567538304019642, 8.7597081577331112E-14, -2.46654615392702E-13,
    6.1537206285128008, 8.5787645505884048E-14, -2.4597249652754803E-13,
    6.1506850866831755, 8.3978209434436984E-14, -2.4529037766239397E-13,
    6.1476495448794868, 8.2168773362989933E-14, -2.4460825879724E-13,
    6.1446140031147012, 8.0360732384890163E-14, -2.43926665852568E-13,
    6.1415808011996029, 7.8551296313443112E-14, -2.43244546987414E-13,
    6.1385452593440419, 7.6741860241996035E-14, -2.4256242812226003E-13,
    6.13550971751701, 7.4932424170548971E-14, -2.41880309257106E-13,
    6.1324741757236945, 7.3122988099101907E-14, -2.41198190391952E-13,
    6.1294386339692837, 7.131494712100215E-14, -2.4051659744728E-13,
    6.126405432056778, 6.9505511049555086E-14, -2.39834478582126E-13,
    6.123369890201217, 6.7696074978108022E-14, -2.39152359716972E-13,
    6.1203343483741852, 6.5886917925330412E-14, -2.3847034603591442E-13,
    206.11729927454044, 6.4077760872552815E-14, -2.3778833235485676E-13,
    5.9828821591151327, 6.2268603819775217E-14, -2.3710631867379925E-13,
    5.8484650773236364, 6.045944676699762E-14, -2.3642430499274164E-13,
    5.71404806279975, 5.8650289714220009E-14, -2.3574229131168403E-13,
    5.5796311491772785, 5.6841132661442406E-14, -2.3506027763062642E-13,
    5.4452143700900253, 5.5031975608664795E-14, -2.3437826394956881E-13,
    5.3107977591717965, 6.0214530566917084, -2.701748639270704E-15,
    6.0343789359860439E-14, 6.0184877891113242, -2.7127334868106716E-15,
    6.06308509622314E-14, 6.0155202338661873, -2.723726810313248E-15,
    6.091813406275236E-14, 6.0125526786469861, -2.7347201338158239E-15,
    6.120541716327329E-14, 6.0095851234537214, -2.7457134573184E-15,
    6.1492700263794233E-14, 6.0066198558474015, -2.7566983048583678E-15,
    6.1779761866165211E-14, 6.0036523005763289, -2.7676916283609442E-15,
    6.2067044966686154E-14, 6.000684745331192, -2.77868495186352E-15,
    6.235432806720711E-14, 5.9977171901145843, -2.7896782753660957E-15,
    6.2641611167728053E-14, 5.9947496349316935, -2.800671598868672E-15,
    6.2928894268249008E-14, 5.9917843673357476, -2.81165644640864E-15,
    6.3215955870619974E-14, 5.9888168120672685, -2.8226497699112159E-15,
    6.3503238971140917E-14, 5.9858492568221315, -2.8336430934137919E-15,
    6.3790522071661872E-14, 5.9828821591151327, -2.8446347217238466E-15,
    6.4077760872552815E-14, 205.97991506140812, -2.855626350033901E-15,
    6.436499967344377E-14, 5.8455659221095857, -2.8666179783439553E-15,
    6.4652238474334725E-14, 5.7112168164448427, -2.8776096066540093E-15,
    6.493947727522568E-14, 5.5768677780477089, -2.8886012349640637E-15,
    6.5226716076116622E-14, 5.44251884055199, -2.8995928632741184E-15,
    6.5513954877007577E-14, 5.30817003759149, -2.9105844915841728E-15,
    6.5801193677898532E-14, -3.8122547664856644E-15, 6.0214530566917084,
    -2.6728219968410996E-13, -3.7439640617760314E-15, 6.0184877891113242,
    -2.6501475956388604E-13, -3.6756206636214074E-15, 6.0155202338661873,
    -2.62745569875668E-13, -3.6072772654667842E-15, 6.0125526786469861,
    -2.6047638018745E-13, -3.5389338673121598E-15, 6.0095851234537214,
    -2.58207190499232E-13, -3.470643162602528E-15, 6.0066198558474015,
    -2.5593975037900803E-13, -3.402299764447904E-15, 6.0036523005763289,
    -2.5367056069079004E-13, -3.33395636629328E-15, 6.000684745331192,
    -2.51401371002572E-13, -3.265612968138656E-15, 5.9977171901145843,
    -2.49132181314354E-13, -3.1972695699840319E-15, 5.9947496349316935,
    -2.46862991626136E-13, -3.1289788652744E-15, 5.9917843673357476,
    -2.4459555150591204E-13, -3.0606354671197754E-15, 5.9888168120672685,
    -2.4232636181769404E-13, -2.9922920689651517E-15, 5.9858492568221315,
    -2.40057172129476E-13, -2.9239592094995266E-15, 5.9828821591151327,
    -2.3778833235485681E-13, -2.8556263500339006E-15, 205.97991506140812,
    -2.3551949258023762E-13, -2.7872934905682754E-15, 5.8455659221095857,
    -2.3325065280561843E-13, -2.71896063110265E-15, 5.7112168164448427,
    -2.3098181303099919E-13, -2.6506277716370243E-15, 5.5768677780477089,
    -2.2871297325638E-13, -2.5822949121713983E-15, 5.44251884055199,
    -2.2644413348176081E-13, -2.5139620527057731E-15, 5.30817003759149,
    -2.2417529370714162E-13, 8.9696764440898425E-14, -2.4501592735841796E-13,
    6.0214530566917084, 8.78884091403894E-14, -2.44338006640114E-13,
    6.0184877891113242, 8.6078658504000341E-14, -2.4365956283483603E-13,
    6.0155202338661873, 8.4268907867611291E-14, -2.42981119029558E-13,
    6.0125526786469861, 8.2459157231222253E-14, -2.4230267522428E-13,
    6.0095851234537214, 8.065080193071322E-14, -2.41624754505976E-13,
    6.0066198558474015, 7.8841051294324169E-14, -2.4094631070069803E-13,
    6.0036523005763289, 7.7031300657935106E-14, -2.4026786689542005E-13,
    6.000684745331192, 7.5221550021546043E-14, -2.39589423090142E-13,
    5.9977171901145843, 7.3411799385157E-14, -2.38910979284864E-13,
    5.9947496349316935, 7.1603444084647972E-14, -2.3823305856656003E-13,
    5.9917843673357476, 6.9793693448258922E-14, -2.37554614761282E-13,
    5.9888168120672685, 6.7983942811869871E-14, -2.36876170956004E-13,
    5.9858492568221315, 6.6174471242656814E-14, -2.3619783176812082E-13,
    5.9828821591151327, 6.436499967344377E-14, -2.3551949258023757E-13,
    205.97991506140812, 6.2555528104230725E-14, -2.3484115339235442E-13,
    5.8455659221095857, 6.0746056535017681E-14, -2.3416281420447123E-13,
    5.7112168164448427, 5.8936584965804624E-14, -2.3348447501658803E-13,
    5.5768677780477089, 5.7127113396591579E-14, -2.3280613582870483E-13,
    5.44251884055199, 5.5317641827378522E-14, -2.3212779664082164E-13,
    5.30817003759149, 5.8861527538625005, -2.6384520738241761E-15,
    5.8538720875224E-14, 5.8832553869519026, -2.649077396117168E-15,
    5.88254681551857E-14, 5.8803557846522621, -2.6597109169613119E-15,
    5.9112436690764654E-14, 5.8774561823785572, -2.6703444378054559E-15,
    5.939940522634361E-14, 5.8745565801307889, -2.6809779586496E-15,
    5.9686373761922554E-14, 5.8716592132072236, -2.6916032809425917E-15,
    5.9973121041884269E-14, 5.8687596108816473, -2.7022368017867361E-15,
    6.0260089577463213E-14, 5.8658600085820067, -2.71287032263088E-15,
    6.0547058113042181E-14, 5.8629604063083027, -2.7235038434750241E-15,
    6.0834026648621138E-14, 5.8600608040631279, -2.734137364319168E-15,
    6.11209951842001E-14, 5.8571634371499366, -2.7447626866121603E-15,
    6.14077424641618E-14, 5.8542638348347342, -2.7553962074563039E-15,
    6.1694710999740753E-14, 5.8513642325376871, -2.7660297283004478E-15,
    6.1981679535319709E-14, 5.8484650773236364, -2.7766616094343616E-15,
    6.2268603819775217E-14, 5.8455659221095857, -2.7872934905682754E-15,
    6.2555528104230725E-14, 205.84266676689555, -2.7979253717021888E-15,
    6.2842452388686233E-14, 5.7083855700899351, -2.8085572528361022E-15,
    6.3129376673141741E-14, 5.5741044069181394, -2.819189133970016E-15,
    6.3416300957597237E-14, 5.4398233110139547, -2.8298210151039298E-15,
    6.3703225242052745E-14, 5.3055423160111834, -2.8404528962378432E-15,
    6.3990149526508253E-14, -3.8282826888148164E-15, 5.8861527538625005,
    -2.6655241929189E-13, -3.7596324588582077E-15, 5.8832553869519026,
    -2.6428865139803404E-13, -3.6909292580451517E-15, 5.8803557846522621,
    -2.62023136769692E-13, -3.6222260572320957E-15, 5.8774561823785572,
    -2.5975762214135E-13, -3.55352285641904E-15, 5.8745565801307889,
    -2.57492107513008E-13, -3.4848726264624318E-15, 5.8716592132072236,
    -2.5522833961915203E-13, -3.4161694256493763E-15, 5.8687596108816473,
    -2.5296282499081E-13, -3.34746622483632E-15, 5.8658600085820067,
    -2.50697310362468E-13, -3.2787630240232639E-15, 5.8629604063083027,
    -2.48431795734126E-13, -3.2100598232102079E-15, 5.8600608040631279,
    -2.4616628110578403E-13, -3.1414095932536E-15, 5.8571634371499366,
    -2.43902513211928E-13, -3.0727063924405436E-15, 5.8542638348347342,
    -2.4163699858358604E-13, -3.0040031916274876E-15, 5.8513642325376871,
    -2.3937148395524397E-13, -2.9353105849857219E-15, 5.8484650773236364,
    -2.371063186737992E-13, -2.8666179783439549E-15, 5.8455659221095857,
    -2.3484115339235442E-13, -2.7979253717021888E-15, 205.84266676689555,
    -2.3257598811090965E-13, -2.7292327650604227E-15, 5.7083855700899351,
    -2.3031082282946482E-13, -2.6605401584186561E-15, 5.5741044069181394,
    -2.2804565754802E-13, -2.5918475517768896E-15, 5.4398233110139547,
    -2.2578049226657519E-13, -2.5231549451351235E-15, 5.3055423160111834,
    -2.2351532698513042E-13, 8.9988406326365981E-14, -2.42695646379462E-13,
    5.8861527538625005, 8.81797367034477E-14, -2.42021397887526E-13,
    5.8832553869519026, 8.6369671502116648E-14, -2.4134662914212403E-13,
    5.8803557846522621, 8.4559606300785611E-14, -2.40671860396722E-13,
    5.8774561823785572, 8.2749541099454573E-14, -2.3999709165132003E-13,
    5.8745565801307889, 8.0940871476536264E-14, -2.39322843159384E-13,
    5.8716592132072236, 7.9130806275205227E-14, -2.38648074413982E-13,
    5.8687596108816473, 7.7320741073874178E-14, -2.3797330566858E-13,
    5.8658600085820067, 7.5510675872543128E-14, -2.37298536923178E-13,
    5.8629604063083027, 7.3700610671212091E-14, -2.36623768177776E-13,
    5.8600608040631279, 7.18919410482938E-14, -2.3594951968584005E-13,
    5.8571634371499366, 7.0081875846962757E-14, -2.35274750940438E-13,
    5.8542638348347342, 6.827181064563172E-14, -2.34599982195036E-13,
    5.8513642325376871, 6.6462024559983216E-14, -2.3392531750032721E-13,
    5.8484650773236364, 6.4652238474334725E-14, -2.3325065280561838E-13,
    5.8455659221095857, 6.2842452388686233E-14, -2.3257598811090965E-13,
    205.84266676689555, 6.1032666303037742E-14, -2.3190132341620081E-13,
    5.7083855700899351, 5.9222880217389238E-14, -2.3122665872149203E-13,
    5.5741044069181394, 5.7413094131740746E-14, -2.305519940267832E-13,
    5.4398233110139547, 5.5603308046092242E-14, -2.2987732933207441E-13,
    5.3055423160111834, 5.7508529555403651, -2.5751555083776478E-15,
    5.6733652390587552E-14, 5.748023455673529, -2.5854213054236639E-15,
    5.7020085348139991E-14, 5.7451917726803927, -2.5956950236093759E-15,
    5.730673931877696E-14, 5.742360089713193, -2.6059687417950879E-15,
    5.7593393289413918E-14, 5.739528406771929, -2.6162424599807996E-15,
    5.7880047260050875E-14, 5.7366989069050947, -2.6265082570268161E-15,
    5.8166480217603314E-14, 5.733867223898991, -2.6367819752125281E-15,
    5.8453134188240283E-14, 5.7310355409058547, -2.64705569339824E-15,
    5.8739788158877253E-14, 5.7282038579386549, -2.6573294115839517E-15,
    5.902644212951421E-14, 5.7253721749973909, -2.6676031297696641E-15,
    5.9313096100151193E-14, 5.7225426751331492, -2.67786892681568E-15,
    5.9599529057703619E-14, 5.7197109921374194, -2.6881426450013918E-15,
    5.9886183028340589E-14, 5.7168793091546579, -2.6984163631871042E-15,
    6.0172836998977559E-14, 5.71404806279975, -2.708688497144877E-15,
    6.0459446766997607E-14, 5.7112168164448418, -2.71896063110265E-15,
    6.0746056535017681E-14, 5.7083855700899342, -2.7292327650604227E-15,
    6.1032666303037729E-14, 205.70555432373502, -2.7395048990181951E-15,
    6.131927607105779E-14, 5.5713410357885707, -2.7497770329759679E-15,
    6.1605885839077851E-14, 5.4371277814759189, -2.7600491669337408E-15,
    6.1892495607097912E-14, 5.3029145944308764, -2.7703213008915136E-15,
    6.2179105375117973E-14, -3.8443106111439684E-15, 5.7508529555403651,
    -2.6582263889967E-13, -3.775300855940384E-15, 5.748023455673529,
    -2.6356254323218205E-13, -3.7062378524688953E-15, 5.7451917726803927,
    -2.6130070366371597E-13, -3.6371748489974081E-15, 5.742360089713193,
    -2.5903886409525E-13, -3.56811184552592E-15, 5.739528406771929,
    -2.56777024526784E-13, -3.4991020903223361E-15, 5.7366989069050947,
    -2.54516928859296E-13, -3.4300390868508481E-15, 5.733867223898991,
    -2.5225508929083004E-13, -3.3609760833793598E-15, 5.7310355409058547,
    -2.49993249722364E-13, -3.2919130799078718E-15, 5.7282038579386549,
    -2.47731410153898E-13, -3.2228500764363842E-15, 5.7253721749973909,
    -2.45469570585432E-13, -3.1538403212328002E-15, 5.7225426751331492,
    -2.43209474917944E-13, -3.0847773177613115E-15, 5.7197109921374194,
    -2.4094763534947803E-13, -3.0157143142898235E-15, 5.7168793091546579,
    -2.38685795781012E-13, -2.9466619604719172E-15, 5.71404806279975,
    -2.3642430499274159E-13, -2.8776096066540097E-15, 5.7112168164448418,
    -2.3416281420447123E-13, -2.8085572528361026E-15, 5.7083855700899342,
    -2.3190132341620081E-13, -2.7395048990181955E-15, 205.70555432373502,
    -2.296398326279304E-13, -2.6704525452002884E-15, 5.5713410357885707,
    -2.2737834183966E-13, -2.6014001913823809E-15, 5.4371277814759189,
    -2.251168510513896E-13, -2.5323478375644738E-15, 5.3029145944308764,
    -2.2285536026311921E-13, 9.0280048211833537E-14, -2.4037536540050597E-13,
    5.7508529555403651, 8.847106426650599E-14, -2.39704789134938E-13,
    5.748023455673529, 8.6660684500232941E-14, -2.3903369544941203E-13,
    5.7451917726803927, 8.4850304733959917E-14, -2.38362601763886E-13,
    5.742360089713193, 8.30399249676869E-14, -2.3769150807836E-13,
    5.739528406771929, 8.1230941022359321E-14, -2.37020931812792E-13,
    5.7366989069050947, 7.94205612560863E-14, -2.36349838127266E-13,
    5.733867223898991, 7.7610181489813249E-14, -2.3567874444174003E-13,
    5.7310355409058547, 7.5799801723540213E-14, -2.35007650756214E-13,
    5.7282038579386549, 7.3989421957267176E-14, -2.34336557070688E-13,
    5.7253721749973909, 7.2180438011939616E-14, -2.3366598080512E-13,
    5.7225426751331492, 7.0370058245666593E-14, -2.32994887119594E-13,
    5.7197109921374194, 6.8559678479393557E-14, -2.32323793434068E-13,
    5.7168793091546579, 6.67495778773096E-14, -2.3165280323253361E-13,
    5.71404806279975, 6.4939477275225667E-14, -2.3098181303099919E-13,
    5.7112168164448418, 6.3129376673141741E-14, -2.3031082282946482E-13,
    5.7083855700899342, 6.13192760710578E-14, -2.2963983262793045E-13,
    205.70555432373502, 5.9509175468973852E-14, -2.2896884242639603E-13,
    5.5713410357885707, 5.7699074866889913E-14, -2.2829785222486161E-13,
    5.4371277814759189, 5.5888974264805962E-14, -2.2762686202332719E-13,
    5.3029145944308764, 5.6155536953513252, -2.51185894293112E-15,
    5.4928583905951108E-14, 5.6127920289022279, -2.52176521473016E-15,
    5.5214702541094277E-14, 5.6100282315895722, -2.53167913025744E-15,
    5.5501041946789254E-14, 5.6072644342898847, -2.54159304578472E-15,
    5.5787381352484225E-14, 5.6045006370161339, -2.5515069613119996E-15,
    5.6073720758179196E-14, 5.6017389705670375, -2.56141323311104E-15,
    5.6359839393322371E-14, 5.5989751732543818, -2.57132714863832E-15,
    5.6646178799017342E-14, 5.596211375941726, -2.5812410641656E-15,
    5.6932518204712325E-14, 5.5934475786420386, -2.5911549796928797E-15,
    5.72188576104073E-14, 5.5906837813682868, -2.60106889522016E-15,
    5.7505197016102279E-14, 5.58792211491919, -2.6109751670192E-15,
    5.7791315651245441E-14, 5.5851583176091282, -2.62088908254648E-15,
    5.8077655056940412E-14, 5.5823945203068472, -2.63080299807376E-15,
    5.83639944626354E-14, 5.5796311491772776, -2.640715384855392E-15,
    5.8650289714220009E-14, 5.5768677780477089, -2.6506277716370243E-15,
    5.8936584965804624E-14, 5.5741044069181394, -2.6605401584186561E-15,
    5.9222880217389238E-14, 5.5713410357885707, -2.670452545200288E-15,
    5.9509175468973852E-14, 205.568577664659, -2.68036493198192E-15,
    5.9795470720558466E-14, 5.4344322519378832, -2.6902773187635521E-15,
    6.008176597214308E-14, 5.3002868728505694, -2.700189705545184E-15,
    6.03680612237277E-14, -3.8603385334731205E-15, 5.6155536953513252,
    -2.6509285850744997E-13, -3.7909692530225595E-15, 5.6127920289022279,
    -2.6283643506633E-13, -3.7215464468926396E-15, 5.6100282315895722,
    -2.6057827055774E-13, -3.6521236407627204E-15, 5.6072644342898847,
    -2.5832010604915E-13, -3.5827008346328E-15, 5.6045006370161339,
    -2.5606194154056003E-13, -3.51333155418224E-15, 5.6017389705670375,
    -2.5380551809944E-13, -3.44390874805232E-15, 5.5989751732543818,
    -2.5154735359085E-13, -3.3744859419223997E-15, 5.596211375941726,
    -2.4928918908225997E-13, -3.30506313579248E-15, 5.5934475786420386,
    -2.4703102457367003E-13, -3.2356403296625602E-15, 5.5906837813682868,
    -2.4477286006508003E-13, -3.166271049212E-15, 5.58792211491919,
    -2.4251643662396E-13, -3.0968482430820797E-15, 5.5851583176091282,
    -2.4025827211537E-13, -3.0274254369521598E-15, 5.5823945203068472,
    -2.3800010760678E-13, -2.9580133359581121E-15, 5.5796311491772776,
    -2.35742291311684E-13, -2.888601234964064E-15, 5.5768677780477089,
    -2.3348447501658803E-13, -2.819189133970016E-15, 5.5741044069181394,
    -2.3122665872149203E-13, -2.7497770329759683E-15, 5.5713410357885707,
    -2.2896884242639603E-13, -2.6803649319819203E-15, 205.568577664659,
    -2.2671102613129998E-13, -2.6109528309878722E-15, 5.4344322519378832,
    -2.24453209836204E-13, -2.5415407299938242E-15, 5.3002868728505694,
    -2.22195393541108E-13, 9.05716900973011E-14, -2.3805508442154995E-13,
    5.6155536953513252, 8.8762391829564282E-14, -2.3738818038235E-13,
    5.6127920289022279, 8.6951697498349247E-14, -2.3672076175670003E-13,
    5.6100282315895722, 8.5141003167134224E-14, -2.3605334313105E-13,
    5.6072644342898847, 8.3330308835919214E-14, -2.353859245054E-13,
    5.6045006370161339, 8.1521010568182378E-14, -2.347190204662E-13,
    5.6017389705670375, 7.9710316236967355E-14, -2.3405160184055003E-13,
    5.5989751732543818, 7.789962190575232E-14, -2.3338418321490005E-13,
    5.596211375941726, 7.60889275745373E-14, -2.3271676458925E-13,
    5.5934475786420386, 7.4278233243322275E-14, -2.320493459636E-13,
    5.5906837813682868, 7.2468934975585438E-14, -2.3138244192440004E-13,
    5.58792211491919, 7.0658240644370416E-14, -2.3071502329875E-13,
    5.5851583176091282, 6.88475463131554E-14, -2.300476046731E-13,
    5.5823945203068472, 6.7037131194636007E-14, -2.2938028896474E-13,
    5.5796311491772776, 6.5226716076116622E-14, -2.2871297325638E-13,
    5.5768677780477089, 6.3416300957597237E-14, -2.2804565754802005E-13,
    5.5741044069181394, 6.1605885839077864E-14, -2.2737834183966004E-13,
    5.5713410357885707, 5.9795470720558466E-14, -2.267110261313E-13,
    205.568577664659, 5.798505560203908E-14, -2.2604371042294E-13,
    5.4344322519378832, 5.6174640483519689E-14, -2.2537639471458002E-13,
    5.3002868728505694, 5.4802550069343736, -2.448562377484592E-15,
    5.3123515421314671E-14, 5.4775611402640223, -2.4581091240366559E-15,
    5.3409319734048564E-14, 5.4748651950058234, -2.4676632369055039E-15,
    5.3695344574801555E-14, 5.4721692497476244, -2.477217349774352E-15,
    5.3981369415554532E-14, 5.4694733045023938, -2.4867714626431997E-15,
    5.4267394256307517E-14, 5.4667794378320433, -2.4963182091952639E-15,
    5.4553198569041422E-14, 5.4640834925738444, -2.505872322064112E-15,
    5.4839223409794406E-14, 5.4613875473156455, -2.51542643493296E-15,
    5.5125248250547397E-14, 5.4586916020574465, -2.5249805478018077E-15,
    5.5411273091300381E-14, 5.4559956568122159, -2.5345346606706562E-15,
    5.5697297932053371E-14, 5.4533017901418646, -2.5440814072227204E-15,
    5.598310224478727E-14, 5.4506058448836656, -2.5536355200915681E-15,
    5.6269127085540254E-14, 5.447909899628061, -2.5631896329604161E-15,
    5.6555151926293245E-14, 5.4452143700900253, -2.5727422725659074E-15,
    5.6841132661442406E-14, 5.4425188405519895, -2.5822949121713987E-15,
    5.7127113396591579E-14, 5.4398233110139538, -2.5918475517768896E-15,
    5.7413094131740746E-14, 5.4371277814759189, -2.6014001913823805E-15,
    5.7699074866889913E-14, 5.4344322519378832, -2.6109528309878718E-15,
    5.798505560203908E-14, 205.43173672239985, -2.6205054705933631E-15,
    5.8271036337188247E-14, 5.2976591512702633, -2.6300581101988544E-15,
    5.8557017072337415E-14, -3.8763664558022725E-15, 5.4802550069343736,
    -2.6436307811523E-13, -3.8066376501047358E-15, 5.4775611402640223,
    -2.62110326900478E-13, -3.7368550413163839E-15, 5.4748651950058234,
    -2.59855837451764E-13, -3.667072432528032E-15, 5.4721692497476244,
    -2.5760134800305E-13, -3.59728982373968E-15, 5.4694733045023938,
    -2.5534685855433603E-13, -3.5275610180421438E-15, 5.4667794378320433,
    -2.53094107339584E-13, -3.4577784092537923E-15, 5.4640834925738444,
    -2.5083961789087004E-13, -3.38799580046544E-15, 5.4613875473156455,
    -2.48585128442156E-13, -3.3182131916770881E-15, 5.4586916020574465,
    -2.46330638993442E-13, -3.2484305828887362E-15, 5.4559956568122159,
    -2.44076149544728E-13, -3.1787017771912002E-15, 5.4533017901418646,
    -2.4182339832997603E-13, -3.1089191684028476E-15, 5.4506058448836656,
    -2.39568908881262E-13, -3.0391365596144957E-15, 5.447909899628061,
    -2.37314419432548E-13, -2.9693647114443074E-15, 5.4452143700900253,
    -2.3506027763062642E-13, -2.8995928632741184E-15, 5.4425188405519895,
    -2.3280613582870483E-13, -2.8298210151039298E-15, 5.4398233110139538,
    -2.305519940267832E-13, -2.7600491669337412E-15, 5.4371277814759189,
    -2.2829785222486161E-13, -2.6902773187635521E-15, 5.4344322519378832,
    -2.2604371042294E-13, -2.6205054705933631E-15, 205.43173672239985,
    -2.2378956862101839E-13, -2.5507336224231745E-15, 5.2976591512702633,
    -2.2153542681909681E-13, 9.0863331982768662E-14, -2.35734803442594E-13,
    5.4802550069343736, 8.9053719392622575E-14, -2.35071571629762E-13,
    5.4775611402640223, 8.7242710496465553E-14, -2.34407828063988E-13,
    5.4748651950058234, 8.5431701600308544E-14, -2.3374408449821397E-13,
    5.4721692497476244, 8.3620692704151535E-14, -2.3308034093244E-13,
    5.4694733045023938, 8.1811080114005422E-14, -2.32417109119608E-13,
    5.4667794378320433, 8.0000071217848413E-14, -2.31753365553834E-13,
    5.4640834925738444, 7.8189062321691391E-14, -2.3108962198806E-13,
    5.4613875473156455, 7.6378053425534369E-14, -2.30425878422286E-13,
    5.4586916020574465, 7.456704452937736E-14, -2.29762134856512E-13,
    5.4559956568122159, 7.275743193923126E-14, -2.2909890304368005E-13,
    5.4533017901418646, 7.0946423043074251E-14, -2.28435159477906E-13,
    5.4506058448836656, 6.9135414146917242E-14, -2.27771415912132E-13,
    5.447909899628061, 6.7324684511962409E-14, -2.2710777469694641E-13,
    5.4452143700900253, 6.5513954877007577E-14, -2.2644413348176079E-13,
    5.4425188405519895, 6.3703225242052745E-14, -2.2578049226657522E-13,
    5.4398233110139538, 6.1892495607097925E-14, -2.2511685105138963E-13,
    5.4371277814759189, 6.008176597214308E-14, -2.24453209836204E-13,
    5.4344322519378832, 5.8271036337188247E-14, -2.2378956862101842E-13,
    205.43173672239985, 5.6460306702233409E-14, -2.2312592740583282E-13,
    5.2976591512702633, 5.3449569239285015, -2.3852658120380641E-15,
    5.1318446936678228E-14, 5.342330823397905, -2.3944530333431518E-15,
    5.1603936927002857E-14, 5.3397026965551708, -2.4036473435535679E-15,
    5.1889647202813855E-14, 5.3370745697124358, -2.412841653763984E-15,
    5.2175357478624846E-14, 5.3344464428697016, -2.4220359639743997E-15,
    5.2461067754435837E-14, 5.331820342339106, -2.4312231852794878E-15,
    5.2746557744760473E-14, 5.3291922154963718, -2.4404174954899039E-15,
    5.3032268020571471E-14, 5.3265640886536376, -2.44961180570032E-15,
    5.3317978296382468E-14, 5.3239359618109026, -2.4588061159107361E-15,
    5.3603688572193459E-14, 5.3213078349681684, -2.4680004261211522E-15,
    5.3889398848004463E-14, 5.3186817344375719, -2.4771876474262403E-15,
    5.4174888838329093E-14, 5.3160536075948377, -2.486381957636656E-15,
    5.4460599114140084E-14, 5.3134254807521026, -2.4955762678470721E-15,
    5.4746309389951082E-14, 5.3107977591717965, -2.5047691602764224E-15,
    5.50319756086648E-14, 5.3081700375914895, -2.5139620527057731E-15,
    5.5317641827378528E-14, 5.3055423160111825, -2.5231549451351235E-15,
    5.5603308046092248E-14, 5.3029145944308764, -2.5323478375644734E-15,
    5.5888974264805975E-14, 5.3002868728505694, -2.5415407299938238E-15,
    5.6174640483519695E-14, 5.2976591512702633, -2.5507336224231745E-15,
    5.6460306702233421E-14, 205.29503142968997, -2.5599265148525248E-15,
    5.6745972920947142E-14, -3.8923943781314245E-15, 5.3449569239285015,
    -2.6363329772300997E-13, -3.8223060471869121E-15, 5.342330823397905,
    -2.6138421873462603E-13, -3.7521636357401274E-15, 5.3397026965551708,
    -2.59133404345788E-13, -3.6820212242933443E-15, 5.3370745697124358,
    -2.5688258995695E-13, -3.61187881284656E-15, 5.3344464428697016,
    -2.5463177556811204E-13, -3.541790481902048E-15, 5.331820342339106,
    -2.52382696579728E-13, -3.4716480704552641E-15, 5.3291922154963718,
    -2.5013188219089E-13, -3.40150565900848E-15, 5.3265640886536376,
    -2.47881067802052E-13, -3.331363247561696E-15, 5.3239359618109026,
    -2.45630253413214E-13, -3.2612208361149121E-15, 5.3213078349681684,
    -2.4337943902437603E-13, -3.1911325051704E-15, 5.3186817344375719,
    -2.4113036003599204E-13, -3.1209900937236154E-15, 5.3160536075948377,
    -2.38879545647154E-13, -3.0508476822768315E-15, 5.3134254807521026,
    -2.36628731258316E-13, -2.9807160869305027E-15, 5.3107977591717965,
    -2.3437826394956881E-13, -2.9105844915841728E-15, 5.3081700375914895,
    -2.3212779664082164E-13, -2.8404528962378432E-15, 5.3055423160111825,
    -2.2987732933207441E-13, -2.770321300891514E-15, 5.3029145944308764,
    -2.2762686202332719E-13, -2.7001897055451844E-15, 5.3002868728505694,
    -2.2537639471458E-13, -2.6300581101988544E-15, 5.2976591512702633,
    -2.231259274058328E-13, -2.5599265148525248E-15, 205.29503142968997,
    -2.208754600970856E-13, 9.1154973868236218E-14, -2.3341452246363796E-13,
    5.3449569239285015, 8.9345046955680868E-14, -2.32754962877174E-13,
    5.342330823397905, 8.7533723494581847E-14, -2.32094894371276E-13,
    5.3397026965551708, 8.5722400033482851E-14, -2.31434825865378E-13,
    5.3370745697124358, 8.3911076572383855E-14, -2.3077475735948003E-13,
    5.3344464428697016, 8.2101149659828479E-14, -2.30115197773016E-13,
    5.331820342339106, 8.0289826198729484E-14, -2.29455129267118E-13,
    5.3291922154963718, 7.8478502737630462E-14, -2.2879506076122003E-13,
    5.3265640886536376, 7.6667179276531454E-14, -2.28134992255322E-13,
    5.3239359618109026, 7.4855855815432458E-14, -2.27474923749424E-13,
    5.3213078349681684, 7.3045928902877083E-14, -2.2681536416296002E-13,
    5.3186817344375719, 7.1234605441778087E-14, -2.26155295657062E-13,
    5.3160536075948377, 6.9423281980679091E-14, -2.25495227151164E-13,
    5.3134254807521026, 6.76122378292888E-14, -2.248352604291528E-13,
    5.3107977591717965, 6.5801193677898532E-14, -2.2417529370714157E-13,
    5.3081700375914895, 6.3990149526508253E-14, -2.2351532698513042E-13,
    5.3055423160111825, 6.2179105375117986E-14, -2.2285536026311924E-13,
    5.3029145944308764, 6.03680612237277E-14, -2.22195393541108E-13,
    5.3002868728505694, 5.8557017072337415E-14, -2.2153542681909681E-13,
    5.2976591512702633, 5.6745972920947135E-14, -2.2087546009708563E-13,
    205.29503142968997 };

   double b_a[3600] = { -207.91568962305743, 3.58787126963432E-15,
    -8.5613912980450155E-14, -7.7803882041935317, 3.6038893042877592E-15,
    -8.5905374951120982E-14, -7.6450837995329533, 3.6199196985358394E-15,
    -8.619706181528774E-14, -7.5097794622411449, 3.6359500927839195E-15,
    -8.6488748679454472E-14, -7.3744752259570951, 3.651980487032E-15,
    -8.67804355436212E-14, -7.23917414334048, 3.66799852168544E-15,
    -8.7071897514292056E-14, -7.103870209457205, 3.68402891593352E-15,
    -8.7363584378458788E-14, -6.9685664774727272, 3.7000593101816E-15,
    -8.765527124262552E-14, -6.8332629810260386, 3.71608970442968E-15,
    -8.7946958106792252E-14, -6.6979597537587221, 3.73212009867776E-15,
    -8.8238644970959E-14, -6.5626598483356409, 3.7481381333312E-15,
    -8.8530106941629836E-14, -6.42735725981494, 3.7641685275792794E-15,
    -8.8821793805796568E-14, -6.292055041354276, 3.7801989218273596E-15,
    -8.9113480669963326E-14, -6.1567538304019642, 3.7962268441565124E-15,
    -8.9405122555430869E-14, -6.0214530566917084, 3.8122547664856636E-15,
    -8.9696764440898438E-14, -5.8861527538625005, 3.8282826888148164E-15,
    -8.9988406326366E-14, -5.7508529555403651, 3.8443106111439677E-15,
    -9.028004821183355E-14, -5.6155536953513252, 3.86033853347312E-15,
    -9.05716900973011E-14, -5.4802550069343736, 3.8763664558022717E-15,
    -9.0863331982768662E-14, -5.3449569239285024, 3.8923943781314237E-15,
    -9.1154973868236218E-14, 3.58787126963432E-15, -207.91568962305743,
    2.77498787522325E-13, 3.52461375203816E-15, -7.7803882041935317,
    2.7517993793200503E-13, 3.4613074246290396E-15, -7.6450837995329533,
    2.7285929910589E-13, 3.3980010972199203E-15, -7.5097794622411449,
    2.70538660279775E-13, 3.3346947698108E-15, -7.3744752259570951,
    2.6821802145366004E-13, 3.27143725221464E-15, -7.23917414334048,
    2.6589917186334004E-13, 3.20813092480552E-15, -7.103870209457205,
    2.63578533037225E-13, 3.1448245973963997E-15, -6.9685664774727272,
    2.6125789421111E-13, 3.08151826998728E-15, -6.8332629810260386,
    2.58937255384995E-13, 3.01821194257816E-15, -6.6979597537587221,
    2.5661661655888003E-13, 2.954954424982E-15, -6.5626598483356409,
    2.5429776696856003E-13, 2.8916480975728796E-15, -6.42735725981494,
    2.51977128142445E-13, 2.82834177016376E-15, -6.292055041354276,
    2.4965648931633E-13, 2.7650452047172323E-15, -6.1567538304019642,
    2.47336208337374E-13, 2.701748639270704E-15, -6.0214530566917084,
    2.45015927358418E-13, 2.6384520738241761E-15, -5.8861527538625005,
    2.4269564637946204E-13, 2.5751555083776482E-15, -5.7508529555403651,
    2.40375365400506E-13, 2.5118589429311203E-15, -5.6155536953513252,
    2.3805508442155E-13, 2.448562377484592E-15, -5.4802550069343736,
    2.35734803442594E-13, 2.3852658120380641E-15, -5.3449569239285024,
    2.33414522463638E-13, -8.5613912980450143E-14, 2.77498787522325E-13,
    -207.91568962305743, -8.3809958048241E-14, 2.76769457333925E-13,
    -7.7803882041935317, -8.2004611175497722E-14, 2.7603956439075E-13,
    -7.6450837995329533, -8.0199264302754473E-14, 2.75309671447575E-13,
    -7.5097794622411449, -7.8393917430011212E-14, 2.745797785044E-13,
    -7.3744752259570951, -7.6589962497802053E-14, 2.7385044831599997E-13,
    -7.23917414334048, -7.4784615625058791E-14, 2.73120555372825E-13,
    -7.103870209457205, -7.297926875231553E-14, 2.7239066242965E-13,
    -6.9685664774727272, -7.1173921879572256E-14, 2.71660769486475E-13,
    -6.8332629810260386, -6.9368575006829E-14, 2.709308765433E-13,
    -6.6979597537587221, -6.7564620074619835E-14, 2.7020154635490003E-13,
    -6.5626598483356409, -6.5759273201876587E-14, 2.69471653411725E-13,
    -6.42735725981494, -6.3953926329133326E-14, 2.6874176046855E-13,
    -6.292055041354276, -6.214885784449687E-14, 2.6801198007633E-13,
    -6.1567538304019642, -6.0343789359860439E-14, 2.6728219968410996E-13,
    -6.0214530566917084, -5.8538720875224E-14, 2.6655241929189003E-13,
    -5.8861527538625005, -5.6733652390587558E-14, 2.6582263889967005E-13,
    -5.7508529555403651, -5.4928583905951108E-14, 2.6509285850745E-13,
    -5.6155536953513252, -5.3123515421314665E-14, 2.6436307811523E-13,
    -5.4802550069343736, -5.1318446936678215E-14, 2.6363329772301E-13,
    -5.3449569239285024, -7.7803882041935317, 3.52461375203816E-15,
    -8.3809958048241E-14, -207.77653715668245, 3.54027248323688E-15,
    -8.410110589040905E-14, -7.64130110038759, 3.55594329678992E-15,
    -8.43924783836899E-14, -7.5060650778225018, 3.57161411034296E-15,
    -8.4683850876970761E-14, -7.3708291226261844, 3.5872849238959995E-15,
    -8.4975223370251591E-14, -7.2355962352407879, 3.6029436550947196E-15,
    -8.5266371212419671E-14, -7.10036051519322, 3.61861446864776E-15,
    -8.5557743705700514E-14, -6.9651249634054349, 3.6342852822008E-15,
    -8.5849116198981356E-14, -6.829889613516448, 3.64995609575384E-15,
    -8.6140488692262212E-14, -6.6946544991652495, 3.66562690930688E-15,
    -8.6431861185543067E-14, -6.5594226207965862, 3.6812856405056E-15,
    -8.6723009027711121E-14, -6.424188077942552, 3.69695645405864E-15,
    -8.7014381520991977E-14, -6.2889538715173421, 3.71262726761168E-15,
    -8.7305754014272832E-14, -6.1537206285128008, 3.728295664693856E-15,
    -8.7597081577331112E-14, -6.0184877891113242, 3.7439640617760322E-15,
    -8.78884091403894E-14, -5.8832553869519026, 3.7596324588582085E-15,
    -8.81797367034477E-14, -5.74802345567353, 3.775300855940384E-15,
    -8.8471064266505977E-14, -5.6127920289022279, 3.7909692530225595E-15,
    -8.8762391829564282E-14, -5.4775611402640232, 3.8066376501047358E-15,
    -8.9053719392622562E-14, -5.3423308233979059, 3.8223060471869121E-15,
    -8.9345046955680855E-14, 3.60388930428776E-15, -7.7803882041935317,
    2.7676945733392495E-13, 3.54027248323688E-15, -207.77653715668245,
    2.74454277704565E-13, 3.4766065751327197E-15, -7.64130110038759,
    2.7213731167117E-13, 3.4129406670285602E-15, -7.5060650778225018,
    2.69820345637775E-13, 3.3492747589244E-15, -7.3708291226261844,
    2.6750337960438E-13, 3.28565793787352E-15, -7.2355962352407879,
    2.6518819997502E-13, 3.22199202976936E-15, -7.10036051519322,
    2.6287123394162504E-13, 3.1583261216652E-15, -6.9651249634054349,
    2.6055426790823E-13, 3.09466021356104E-15, -6.829889613516448,
    2.58237301874835E-13, 3.0309943054568802E-15, -6.6946544991652495,
    2.5592033584144E-13, 2.967377484406E-15, -6.5594226207965862,
    2.5360515621208E-13, 2.9037115763018395E-15, -6.424188077942552,
    2.51288190178685E-13, 2.8400456681976797E-15, -6.2889538715173421,
    2.4897122414529E-13, 2.776389577504176E-15, -6.1537206285128008,
    2.46654615392702E-13, 2.712733486810672E-15, -6.0184877891113242,
    2.44338006640114E-13, 2.649077396117168E-15, -5.8832553869519026,
    2.4202139788752603E-13, 2.5854213054236643E-15, -5.74802345567353,
    2.39704789134938E-13, 2.5217652147301603E-15, -5.6127920289022279,
    2.3738818038235E-13, 2.4581091240366563E-15, -5.4775611402640232,
    2.35071571629762E-13, 2.3944530333431522E-15, -5.3423308233979059,
    2.32754962877174E-13, -8.5905374951120982E-14, 2.75179937932005E-13,
    -7.7803882041935317, -8.4101105890409063E-14, 2.74454277704565E-13,
    -207.77653715668245, -8.22954446467799E-14, 2.7372805755411004E-13,
    -7.64130110038759, -8.0489783403150762E-14, 2.7300183740365496E-13,
    -7.5060650778225018, -7.8684122159521611E-14, 2.722756172532E-13,
    -7.3708291226261844, -7.6879853098809667E-14, 2.7154995702575997E-13,
    -7.2355962352407879, -7.5074191855180529E-14, 2.70823736875305E-13,
    -7.10036051519322, -7.3268530611551365E-14, 2.7009751672485003E-13,
    -6.9651249634054349, -7.14628693679222E-14, 2.69371296574395E-13,
    -6.829889613516448, -6.9657208124293051E-14, 2.6864507642394E-13,
    -6.6946544991652495, -6.785293906358112E-14, 2.679194161965E-13,
    -6.5594226207965862, -6.6047277819951969E-14, 2.67193196046045E-13,
    -6.424188077942552, -6.4241616576322818E-14, 2.6646697589559E-13,
    -6.2889538715173421, -6.2436233769277111E-14, 2.65740867729738E-13,
    -6.1537206285128008, -6.06308509622314E-14, 2.65014759563886E-13,
    -6.0184877891113242, -5.88254681551857E-14, 2.6428865139803404E-13,
    -5.8832553869519026, -5.7020085348139991E-14, 2.6356254323218205E-13,
    -5.74802345567353, -5.5214702541094271E-14, 2.6283643506633E-13,
    -5.6127920289022279, -5.3409319734048564E-14, 2.62110326900478E-13,
    -5.4775611402640232, -5.1603936927002851E-14, 2.6138421873462603E-13,
    -5.3423308233979059, -7.6450837995329533, 3.4613074246290396E-15,
    -8.2004611175497735E-14, -7.64130110038759, 3.4766065751327197E-15,
    -8.2295444646779912E-14, -207.63751548632862, 3.4919175305364797E-15,
    -8.2586502526604744E-14, -7.5023478307688318, 3.5072284859402396E-15,
    -8.2877560406429576E-14, -7.3671802089388274, 3.5225394413439996E-15,
    -8.31686182862544E-14, -7.2320155690229226, 3.5378385918476796E-15,
    -8.3459451757536573E-14, -7.0968481150896263, 3.55314954725144E-15,
    -8.37505096373614E-14, -6.9616807957771414, 3.5684605026552004E-15,
    -8.4041567517186237E-14, -6.82651364472444, 3.58377145805896E-15,
    -8.4332625397011069E-14, -6.6913466955705356, 3.59908241346272E-15,
    -8.46236832768359E-14, -6.5561828964997506, 3.6143815639664004E-15,
    -8.4914516748118078E-14, -6.4210164515831876, 3.62969251937016E-15,
    -8.520557462794291E-14, -6.2858503094642364, 3.64500347477392E-15,
    -8.5496632507767755E-14, -6.1506850866831755, 3.6603120691976639E-15,
    -8.5787645505884048E-14, -6.0155202338661873, 3.6756206636214082E-15,
    -8.6078658504000354E-14, -5.8803557846522629, 3.6909292580451517E-15,
    -8.6369671502116648E-14, -5.7451917726803936, 3.706237852468896E-15,
    -8.6660684500232954E-14, -5.6100282315895731, 3.7215464468926396E-15,
    -8.6951697498349247E-14, -5.4748651950058242, 3.7368550413163839E-15,
    -8.7242710496465553E-14, -5.3397026965551717, 3.7521636357401274E-15,
    -8.7533723494581847E-14, 3.61991969853584E-15, -7.6450837995329533,
    2.7603956439075E-13, 3.55594329678992E-15, -7.64130110038759,
    2.7372805755411004E-13, 3.4919175305364797E-15, -207.63751548632862,
    2.7141476714738E-13, 3.42789176428304E-15, -7.5023478307688318,
    2.6910147674065E-13, 3.3638659980295997E-15, -7.3671802089388274,
    2.6678818633392003E-13, 3.29988959628368E-15, -7.2320155690229226,
    2.6447667949728E-13, 3.23586383003024E-15, -7.0968481150896263,
    2.6216338909055003E-13, 3.1718380637768E-15, -6.9616807957771414,
    2.5985009868382E-13, 3.10781229752336E-15, -6.82651364472444,
    2.5753680827709E-13, 3.04378653126992E-15, -6.6913466955705356,
    2.5522351787036E-13, 2.9798101295240003E-15, -6.5561828964997506,
    2.5291201103372E-13, 2.9157843632705596E-15, -6.4210164515831876,
    2.5059872062699E-13, 2.8517585970171196E-15, -6.2858503094642364,
    2.4828543022026E-13, 2.7877427036651844E-15, -6.1506850866831755,
    2.45972496527548E-13, 2.723726810313248E-15, -6.0155202338661873,
    2.4365956283483603E-13, 2.6597109169613119E-15, -5.8803557846522629,
    2.4134662914212403E-13, 2.5956950236093763E-15, -5.7451917726803936,
    2.3903369544941203E-13, 2.5316791302574403E-15, -5.6100282315895731,
    2.367207617567E-13, 2.4676632369055039E-15, -5.4748651950058242,
    2.34407828063988E-13, 2.4036473435535683E-15, -5.3397026965551717,
    2.32094894371276E-13, -8.6197061815287727E-14, 2.7285929910588996E-13,
    -7.6450837995329533, -8.43924783836899E-14, 2.7213731167117E-13,
    -7.64130110038759, -8.2586502526604732E-14, 2.7141476714738E-13,
    -207.63751548632862, -8.078052666951957E-14, 2.7069222262358997E-13,
    -7.5023478307688318, -7.8974550812434409E-14, 2.699696780998E-13,
    -7.3671802089388274, -7.7169967380836587E-14, 2.6924769066508E-13,
    -7.2320155690229226, -7.5363991523751426E-14, 2.6852514614129E-13,
    -7.0968481150896263, -7.3558015666666239E-14, 2.678026016175E-13,
    -6.9616807957771414, -7.1752039809581065E-14, 2.6708005709371003E-13,
    -6.82651364472444, -6.99460639524959E-14, 2.6635751256992E-13,
    -6.6913466955705356, -6.8141480520898082E-14, 2.656355251352E-13,
    -6.5561828964997506, -6.6335504663812921E-14, 2.6491298061141E-13,
    -6.4210164515831876, -6.4529528806727747E-14, 2.6419043608762E-13,
    -6.2858503094642364, -6.2723831434740041E-14, 2.63468002981644E-13,
    -6.1506850866831755, -6.0918134062752347E-14, 2.62745569875668E-13,
    -6.0155202338661873, -5.9112436690764654E-14, 2.6202313676969203E-13,
    -5.8803557846522629, -5.730673931877696E-14, 2.61300703663716E-13,
    -5.7451917726803936, -5.5501041946789248E-14, 2.6057827055774E-13,
    -5.6100282315895731, -5.3695344574801555E-14, 2.59855837451764E-13,
    -5.4748651950058242, -5.1889647202813842E-14, 2.59133404345788E-13,
    -5.3397026965551717, -7.5097794622411449, 3.3980010972199203E-15,
    -8.0199264302754473E-14, -7.5060650778225018, 3.41294066702856E-15,
    -8.0489783403150749E-14, -7.5023478307688318, 3.42789176428304E-15,
    -8.0780526669519583E-14, -207.49863058380075, 3.44284286153752E-15,
    -8.1071269935888391E-14, -7.3635312953318808, 3.4577939587919997E-15,
    -8.13620132022572E-14, -7.2284349028802772, 3.4727335286006397E-15,
    -8.16525323026535E-14, -7.0933357150560647, 3.48768462585512E-15,
    -8.1943275569022309E-14, -6.9582366282136716, 3.5026357231096E-15,
    -8.223401883539113E-14, -6.823137675992089, 3.51758682036408E-15,
    -8.2524762101759938E-14, -6.6880388920302911, 3.53253791761856E-15,
    -8.2815505368128759E-14, -6.5529431722547864, 3.5474774874272002E-15,
    -8.3106024468525048E-14, -6.4178448252756946, 3.5624285846816805E-15,
    -8.3396767734893869E-14, -6.2827467474500347, 3.57737968193616E-15,
    -8.3687511001262677E-14, -6.1476495448794868, 3.5923284737014719E-15,
    -8.3978209434436984E-14, -6.0125526786469861, 3.6072772654667842E-15,
    -8.4268907867611291E-14, -5.877456182378558, 3.6222260572320965E-15,
    -8.4559606300785611E-14, -5.7423600897131939, 3.6371748489974081E-15,
    -8.4850304733959917E-14, -5.6072644342898856, 3.65212364076272E-15,
    -8.5141003167134237E-14, -5.4721692497476253, 3.667072432528032E-15,
    -8.5431701600308544E-14, -5.3370745697124375, 3.6820212242933443E-15,
    -8.5722400033482851E-14, 3.63595009278392E-15, -7.5097794622411449,
    2.75309671447575E-13, 3.57161411034296E-15, -7.5060650778225018,
    2.73001837403655E-13, 3.5072284859402396E-15, -7.5023478307688318,
    2.7069222262358997E-13, 3.44284286153752E-15, -207.49863058380075,
    2.68382607843525E-13, 3.3784572371347998E-15, -7.3635312953318808,
    2.6607299306346E-13, 3.31412125469384E-15, -7.2284349028802772,
    2.6376515901954E-13, 3.2497356302911203E-15, -7.0933357150560647,
    2.61455544239475E-13, 3.1853500058883997E-15, -6.9582366282136716,
    2.5914592945941E-13, 3.12096438148568E-15, -6.823137675992089,
    2.56836314679345E-13, 3.0565787570829603E-15, -6.6880388920302911,
    2.5452669989928E-13, 2.992242774642E-15, -6.5529431722547864,
    2.5221886585536E-13, 2.9278571502392796E-15, -6.4178448252756946,
    2.4990925107529503E-13, 2.86347152583656E-15, -6.2827467474500347,
    2.4759963629523E-13, 2.7990958298261923E-15, -6.1476495448794868,
    2.45290377662394E-13, 2.7347201338158239E-15, -6.0125526786469861,
    2.42981119029558E-13, 2.6703444378054559E-15, -5.877456182378558,
    2.4067186039672203E-13, 2.6059687417950883E-15, -5.7423600897131939,
    2.38362601763886E-13, 2.5415930457847204E-15, -5.6072644342898856,
    2.3605334313105E-13, 2.477217349774352E-15, -5.4721692497476253,
    2.33744084498214E-13, 2.412841653763984E-15, -5.3370745697124375,
    2.31434825865378E-13, -8.6488748679454459E-14, 2.7053866027977495E-13,
    -7.5097794622411449, -8.4683850876970761E-14, 2.69820345637775E-13,
    -7.5060650778225018, -8.2877560406429563E-14, 2.6910147674065E-13,
    -7.5023478307688318, -8.1071269935888391E-14, 2.68382607843525E-13,
    -207.49863058380075, -7.92649794653472E-14, 2.676637389464E-13,
    -7.3635312953318808, -7.74600816628635E-14, 2.669454243044E-13,
    -7.2284349028802772, -7.5653791192322311E-14, 2.66226555407275E-13,
    -7.0933357150560647, -7.3847500721781126E-14, 2.6550768651015003E-13,
    -6.9582366282136716, -7.2041210251239928E-14, 2.64788817613025E-13,
    -6.823137675992089, -7.0234919780698756E-14, 2.640699487159E-13,
    -6.6880388920302911, -6.8430021978215045E-14, 2.633516340739E-13,
    -6.5529431722547864, -6.662373150767386E-14, 2.62632765176775E-13,
    -6.4178448252756946, -6.4817441037132675E-14, 2.6191389627965E-13,
    -6.2827467474500347, -6.3011429100202983E-14, 2.6119513823355003E-13,
    -6.1476495448794868, -6.120541716327329E-14, 2.6047638018745E-13,
    -6.0125526786469861, -5.939940522634361E-14, 2.5975762214135E-13,
    -5.877456182378558, -5.7593393289413918E-14, 2.5903886409525005E-13,
    -5.7423600897131939, -5.5787381352484225E-14, 2.5832010604915E-13,
    -5.6072644342898856, -5.3981369415554539E-14, 2.5760134800305E-13,
    -5.4721692497476253, -5.217535747862484E-14, 2.5688258995695E-13,
    -5.3370745697124375, -7.3744752259570951, 3.3346947698108E-15,
    -7.8393917430011212E-14, -7.3708291226261844, 3.3492747589244E-15,
    -7.86841221595216E-14, -7.3671802089388283, 3.3638659980295997E-15,
    -7.8974550812434409E-14, -7.3635312953318808, 3.3784572371348E-15,
    -7.92649794653472E-14, -207.35988238181054, 3.39304847624E-15,
    -7.955540811826E-14, -7.22485423681804, 3.4076284653535998E-15,
    -7.98456128477704E-14, -7.0898233150977212, 3.4222197044588003E-15,
    -8.01360415006832E-14, -6.9547924607202329, 3.436810943564E-15,
    -8.042647015359601E-14, -6.8197617073245622, 3.4514021826691997E-15,
    -8.07168988065088E-14, -6.6847310885497029, 3.4659934217744E-15,
    -8.10073274594216E-14, -6.5497034480642924, 3.480573410888E-15,
    -8.1297532188932E-14, -6.4146731990200738, 3.4951646499931998E-15,
    -8.15879608418448E-14, -6.2796431854877053, 3.5097558890984002E-15,
    -8.1878389494757612E-14, -6.1446140031147012, 3.5243448782052802E-15,
    -8.216877336298992E-14, -6.0095851234537214, 3.53893386731216E-15,
    -8.245915723122224E-14, -5.87455658013079, 3.55352285641904E-15,
    -8.2749541099454561E-14, -5.73952840677193, 3.56811184552592E-15,
    -8.30399249676869E-14, -5.6045006370161339, 3.5827008346328E-15,
    -8.33303088359192E-14, -5.4694733045023947, 3.59728982373968E-15,
    -8.3620692704151535E-14, -5.3344464428697025, 3.61187881284656E-15,
    -8.3911076572383842E-14, 3.651980487032E-15, -7.3744752259570951,
    2.745797785044E-13, 3.5872849238959995E-15, -7.3708291226261844,
    2.7227561725320004E-13, 3.5225394413439996E-15, -7.3671802089388283,
    2.699696780998E-13, 3.457793958792E-15, -7.3635312953318808,
    2.676637389464E-13, 3.39304847624E-15, -207.35988238181054,
    2.65357799793E-13, 3.328352913104E-15, -7.22485423681804, 2.630536385418E-13,
    3.263607430552E-15, -7.0898233150977212, 2.607476993884E-13, 3.198861948E-15,
    -6.9547924607202329, 2.58441760235E-13, 3.134116465448E-15,
    -6.8197617073245622, 2.561358210816E-13, 3.069370982896E-15,
    -6.6847310885497029, 2.538298819282E-13, 3.0046754197600003E-15,
    -6.5497034480642924, 2.5152572067700003E-13, 2.9399299372079996E-15,
    -6.4146731990200738, 2.4921978152360003E-13, 2.8751844546559997E-15,
    -6.2796431854877053, 2.469138423702E-13, 2.8104489559872E-15,
    -6.1446140031147012, 2.4460825879724E-13, 2.7457134573184E-15,
    -6.0095851234537214, 2.4230267522428E-13, 2.6809779586496E-15,
    -5.87455658013079, 2.3999709165132003E-13, 2.6162424599808003E-15,
    -5.73952840677193, 2.3769150807836E-13, 2.5515069613120004E-15,
    -5.6045006370161339, 2.353859245054E-13, 2.4867714626432E-15,
    -5.4694733045023947, 2.3308034093244E-13, 2.4220359639744E-15,
    -5.3344464428697025, 2.3077475735948003E-13, -8.6780435543621191E-14,
    2.6821802145366E-13, -7.3744752259570951, -8.49752233702516E-14,
    2.6750337960438E-13, -7.3708291226261844, -8.31686182862544E-14,
    2.6678818633392003E-13, -7.3671802089388283, -8.13620132022572E-14,
    2.6607299306346E-13, -7.3635312953318808, -7.955540811826E-14,
    2.65357799793E-13, -207.35988238181054, -7.77501959448904E-14,
    2.6464315794372E-13, -7.22485423681804, -7.5943590860893208E-14,
    2.6392796467326E-13, -7.0898233150977212, -7.4136985776896E-14,
    2.6321277140280006E-13, -6.9547924607202329, -7.2330380692898791E-14,
    2.6249757813234E-13, -6.8197617073245622, -7.05237756089016E-14,
    2.6178238486188E-13, -6.6847310885497029, -6.8718563435532E-14,
    2.610677430126E-13, -6.5497034480642924, -6.69119583515348E-14,
    2.6035254974214E-13, -6.4146731990200738, -6.51053532675376E-14,
    2.5963735647168E-13, -6.2796431854877053, -6.3299026765665925E-14,
    2.58922273485456E-13, -6.1446140031147012, -6.1492700263794246E-14,
    2.5820719049923196E-13, -6.0095851234537214, -5.9686373761922567E-14,
    2.57492107513008E-13, -5.87455658013079, -5.7880047260050887E-14,
    2.56777024526784E-13, -5.73952840677193, -5.60737207581792E-14,
    2.5606194154056003E-13, -5.6045006370161339, -5.4267394256307523E-14,
    2.5534685855433603E-13, -5.4694733045023947, -5.2461067754435837E-14,
    2.5463177556811204E-13, -5.3344464428697025, -7.2391741433404793,
    3.27143725221464E-15, -7.6589962497802053E-14, -7.235596235240787,
    3.28565793787352E-15, -7.6879853098809667E-14, -7.2320155690229218,
    3.29988959628368E-15, -7.7169967380836587E-14, -7.2284349028802755,
    3.31412125469384E-15, -7.74600816628635E-14, -7.2248542368180395,
    3.3283529131039996E-15, -7.7750195944890391E-14, -207.22127632865352,
    3.34257359876288E-15, -7.8040086545898018E-14, -7.08631362077408,
    3.35680525717304E-15, -7.8330200827924926E-14, -6.9513509465984713,
    3.3710369155832003E-15, -7.8620315109951847E-14, -6.8163883397656928,
    3.38526857399336E-15, -7.8910429391978755E-14, -6.6814258339147328,
    3.39950023240352E-15, -7.9200543674005676E-14, -6.5464662204967077,
    3.4137209180624003E-15, -7.949043427501329E-14, -6.4115040171243436,
    3.42795257647256E-15, -7.97805485570402E-14, -6.2765420156248357,
    3.44218423488272E-15, -8.0070662839067119E-14, -6.141580801199602,
    3.4564136987426242E-15, -8.0360732384890163E-14, -6.0066198558474007,
    3.470643162602528E-15, -8.0650801930713207E-14, -5.8716592132072236,
    3.4848726264624318E-15, -8.0940871476536264E-14, -5.7366989069050947,
    3.4991020903223357E-15, -8.1230941022359321E-14, -5.6017389705670366,
    3.5133315541822395E-15, -8.1521010568182378E-14, -5.4667794378320433,
    3.5275610180421442E-15, -8.1811080114005422E-14, -5.331820342339106,
    3.541790481902048E-15, -8.2101149659828467E-14, 3.66799852168544E-15,
    -7.2391741433404793, 2.7385044831599997E-13, 3.6029436550947196E-15,
    -7.235596235240787, 2.7154995702576E-13, 3.5378385918476796E-15,
    -7.2320155690229218, 2.6924769066508E-13, 3.47273352860064E-15,
    -7.2284349028802755, 2.669454243044E-13, 3.4076284653535998E-15,
    -7.2248542368180395, 2.6464315794372003E-13, 3.34257359876288E-15,
    -207.22127632865352, 2.6234266665348004E-13, 3.27746853551584E-15,
    -7.08631362077408, 2.600404002928E-13, 3.2123634722687996E-15,
    -6.9513509465984713, 2.5773813393212E-13, 3.14725840902176E-15,
    -6.8163883397656928, 2.5543586757144E-13, 3.08215334577472E-15,
    -6.6814258339147328, 2.5313360121076004E-13, 3.0170984791840002E-15,
    -6.5464662204967077, 2.5083310992052004E-13, 2.9519934159369595E-15,
    -6.4115040171243436, 2.4853084355984003E-13, 2.8868883526899196E-15,
    -6.2765420156248357, 2.4622857719916E-13, 2.8217933287741443E-15,
    -6.141580801199602, 2.43926665852568E-13, 2.7566983048583678E-15,
    -6.0066198558474007, 2.41624754505976E-13, 2.6916032809425921E-15,
    -5.8716592132072236, 2.39322843159384E-13, 2.6265082570268161E-15,
    -5.7366989069050947, 2.37020931812792E-13, 2.5614132331110404E-15,
    -5.6017389705670366, 2.347190204662E-13, 2.4963182091952639E-15,
    -5.4667794378320433, 2.32417109119608E-13, 2.4312231852794882E-15,
    -5.331820342339106, 2.30115197773016E-13, -8.7071897514292043E-14,
    2.6589917186334E-13, -7.2391741433404793, -8.5266371212419671E-14,
    2.6518819997502E-13, -7.235596235240787, -8.3459451757536573E-14,
    2.6447667949728E-13, -7.2320155690229218, -8.1652532302653488E-14,
    2.6376515901954E-13, -7.2284349028802755, -7.98456128477704E-14,
    2.630536385418E-13, -7.2248542368180395, -7.8040086545898018E-14,
    2.6234266665348E-13, -207.22127632865352, -7.6233167091014945E-14,
    2.6163114617574E-13, -7.08631362077408, -7.4426247636131847E-14,
    2.6091962569800003E-13, -6.9513509465984713, -7.261932818124875E-14,
    2.6020810522026E-13, -6.8163883397656928, -7.0812408726365664E-14,
    2.5949658474252E-13, -6.6814258339147328, -6.900688242449328E-14,
    2.5878561285420005E-13, -6.5464662204967077, -6.71999629696102E-14,
    2.5807409237646E-13, -6.4115040171243436, -6.5393043514727109E-14,
    2.5736257189872E-13, -6.2765420156248357, -6.3586402690446154E-14,
    2.56651161138864E-13, -6.141580801199602, -6.1779761866165211E-14,
    2.55939750379008E-13, -6.0066198558474007, -5.9973121041884269E-14,
    2.5522833961915203E-13, -5.8716592132072236, -5.8166480217603326E-14,
    2.54516928859296E-13, -5.7366989069050947, -5.6359839393322371E-14,
    2.5380551809944E-13, -5.6017389705670366, -5.4553198569041428E-14,
    2.53094107339584E-13, -5.4667794378320433, -5.2746557744760467E-14,
    2.52382696579728E-13, -5.331820342339106, -7.1038702094572042,
    3.20813092480552E-15, -7.4784615625058791E-14, -7.100360515193219,
    3.2219920297693597E-15, -7.5074191855180516E-14, -7.0968481150896263,
    3.23586383003024E-15, -7.5363991523751413E-14, -7.0933357150560639,
    3.24973563029112E-15, -7.5653791192322311E-14, -7.0898233150977212,
    3.2636074305519998E-15, -7.59435908608932E-14, -7.08631362077408,
    3.27746853551584E-15, -7.6233167091014933E-14, -207.08280122060566,
    3.29134033577672E-15, -7.6522966759585817E-14, -6.9479067789105216,
    3.3052121360376E-15, -7.6812766428156727E-14, -6.8130123709192159,
    3.3190839362984796E-15, -7.7102566096727612E-14, -6.67811803027074,
    3.33295573655936E-15, -7.7392365765298509E-14, -6.5432264961583746,
    3.3468168415232E-15, -7.7681941995420247E-14, -6.4083323907364491,
    3.36068864178408E-15, -7.7971741663991144E-14, -6.2734384535483878,
    3.37456044204496E-15, -7.8261541332562041E-14, -6.138545259344041,
    3.3884301032464321E-15, -7.85512963134431E-14, -6.003652300576328,
    3.402299764447904E-15, -7.8841051294324157E-14, -5.8687596108816473,
    3.4161694256493759E-15, -7.9130806275205227E-14, -5.733867223898991,
    3.4300390868508477E-15, -7.9420561256086285E-14, -5.5989751732543818,
    3.4439087480523196E-15, -7.9710316236967343E-14, -5.4640834925738444,
    3.4577784092537919E-15, -8.0000071217848413E-14, -5.3291922154963718,
    3.4716480704552641E-15, -8.0289826198729471E-14, 3.68402891593352E-15,
    -7.1038702094572042, 2.7312055537282496E-13, 3.61861446864776E-15,
    -7.100360515193219, 2.7082373687530505E-13, 3.55314954725144E-15,
    -7.0968481150896263, 2.6852514614129E-13, 3.48768462585512E-15,
    -7.0933357150560639, 2.66226555407275E-13, 3.4222197044588E-15,
    -7.0898233150977212, 2.6392796467326E-13, 3.35680525717304E-15,
    -7.08631362077408, 2.6163114617574003E-13, 3.2913403357767202E-15,
    -207.08280122060566, 2.59332555441725E-13, 3.2258754143803998E-15,
    -6.9479067789105216, 2.5703396470771E-13, 3.16041049298408E-15,
    -6.8130123709192159, 2.54735373973695E-13, 3.09494557158776E-15,
    -6.67811803027074, 2.5243678323968E-13, 3.029531124302E-15,
    -6.5432264961583746, 2.5013996474216E-13, 2.9640662029056796E-15,
    -6.4083323907364491, 2.4784137400814503E-13, 2.89860128150936E-15,
    -6.2734384535483878, 2.4554278327413E-13, 2.8331464549351522E-15,
    -6.138545259344041, 2.43244546987414E-13, 2.7676916283609438E-15,
    -6.003652300576328, 2.4094631070069803E-13, 2.7022368017867361E-15,
    -5.8687596108816473, 2.38648074413982E-13, 2.6367819752125281E-15,
    -5.733867223898991, 2.36349838127266E-13, 2.5713271486383204E-15,
    -5.5989751732543818, 2.3405160184055E-13, 2.505872322064112E-15,
    -5.4640834925738444, 2.31753365553834E-13, 2.4404174954899043E-15,
    -5.3291922154963718, 2.29455129267118E-13, -8.7363584378458775E-14,
    2.6357853303722497E-13, -7.1038702094572042, -8.5557743705700526E-14,
    2.62871233941625E-13, -7.100360515193219, -8.37505096373614E-14,
    2.6216338909055003E-13, -7.0968481150896263, -8.1943275569022309E-14,
    2.6145554423947497E-13, -7.0933357150560639, -8.0136041500683213E-14,
    2.607476993884E-13, -7.0898233150977212, -7.8330200827924939E-14,
    2.6004040029279997E-13, -7.08631362077408, -7.6522966759585843E-14,
    2.59332555441725E-13, -207.08280122060566, -7.4715732691246721E-14,
    2.5862471059065E-13, -6.9479067789105216, -7.2908498622907613E-14,
    2.57916865739575E-13, -6.8130123709192159, -7.1101264554568517E-14,
    2.572090208885E-13, -6.67811803027074, -6.9295423881810242E-14,
    2.565017217929E-13, -6.5432264961583746, -6.7488189813471134E-14,
    2.55793876941825E-13, -6.4083323907364491, -6.5680955745132038E-14,
    2.5508603209075E-13, -6.2734384535483878, -6.38740003559091E-14,
    2.5437829639077E-13, -6.138545259344041, -6.2067044966686167E-14,
    2.5367056069079E-13, -6.003652300576328, -6.0260089577463225E-14,
    2.5296282499081E-13, -5.8687596108816473, -5.84531341882403E-14,
    2.5225508929083004E-13, -5.733867223898991, -5.6646178799017348E-14,
    2.5154735359085E-13, -5.5989751732543818, -5.4839223409794412E-14,
    2.5083961789087E-13, -5.4640834925738444, -5.3032268020571464E-14,
    2.5013188219089E-13, -5.3291922154963718, -6.9685664774727263,
    3.1448245973963997E-15, -7.297926875231553E-14, -6.9651249634054349,
    3.1583261216652E-15, -7.3268530611551365E-14, -6.9616807957771414,
    3.1718380637768E-15, -7.3558015666666252E-14, -6.9582366282136707,
    3.1853500058884E-15, -7.3847500721781126E-14, -6.9547924607202321,
    3.198861948E-15, -7.4136985776896E-14, -6.9513509465984713,
    3.2123634722688E-15, -7.4426247636131847E-14, -6.9479067789105216,
    3.2258754143804E-15, -7.4715732691246721E-14, -206.94446261128223,
    3.2393873564920003E-15, -7.5005217746361608E-14, -6.8096364021272082,
    3.2528992986036E-15, -7.5294702801476482E-14, -6.6748102266760263,
    3.2664112407152E-15, -7.5584187856591368E-14, -6.5399867718641325,
    3.2799127649840003E-15, -7.58734497158272E-14, -6.4051607643874586,
    3.2934247070956E-15, -7.616293477094209E-14, -6.270334891505656,
    3.3069366492072E-15, -7.6452419826056964E-14, -6.1355097175170092,
    3.32044650775024E-15, -7.6741860241996035E-14, -6.0006847453311911,
    3.33395636629328E-15, -7.7031300657935106E-14, -5.8658600085820067,
    3.3474662248363203E-15, -7.7320741073874178E-14, -5.7310355409058547,
    3.3609760833793598E-15, -7.7610181489813249E-14, -5.596211375941726,
    3.3744859419223997E-15, -7.7899621905752333E-14, -5.4613875473156464,
    3.38799580046544E-15, -7.81890623216914E-14, -5.3265640886536376,
    3.40150565900848E-15, -7.8478502737630462E-14, 3.7000593101816E-15,
    -6.9685664774727263, 2.7239066242964996E-13, 3.6342852822008E-15,
    -6.9651249634054349, 2.7009751672485003E-13, 3.5684605026552E-15,
    -6.9616807957771414, 2.678026016175E-13, 3.5026357231096E-15,
    -6.9582366282136707, 2.6550768651015E-13, 3.436810943564E-15,
    -6.9547924607202321, 2.632127714028E-13, 3.3710369155832E-15,
    -6.9513509465984713, 2.6091962569800003E-13, 3.3052121360376E-15,
    -6.9479067789105216, 2.5862471059065E-13, 3.239387356492E-15,
    -206.94446261128223, 2.563297954833E-13, 3.1735625769464E-15,
    -6.8096364021272082, 2.5403488037595E-13, 3.1077377974008002E-15,
    -6.6748102266760263, 2.5173996526860004E-13, 3.0419637694200002E-15,
    -6.5399867718641325, 2.494468195638E-13, 2.9761389898743996E-15,
    -6.4051607643874586, 2.4715190445645004E-13, 2.9103142103287998E-15,
    -6.270334891505656, 2.448569893491E-13, 2.84449958109616E-15,
    -6.1355097175170092, 2.4256242812226E-13, 2.77868495186352E-15,
    -6.0006847453311911, 2.4026786689542E-13, 2.71287032263088E-15,
    -5.8658600085820067, 2.3797330566858E-13, 2.64705569339824E-15,
    -5.7310355409058547, 2.3567874444174003E-13, 2.5812410641656004E-15,
    -5.596211375941726, 2.333841832149E-13, 2.51542643493296E-15,
    -5.4613875473156464, 2.3108962198806E-13, 2.4496118057003204E-15,
    -5.3265640886536376, 2.2879506076122003E-13, -8.765527124262552E-14,
    2.6125789421110996E-13, -6.9685664774727263, -8.5849116198981369E-14,
    2.6055426790823E-13, -6.9651249634054349, -8.4041567517186237E-14,
    2.5985009868382004E-13, -6.9616807957771414, -8.2234018835391117E-14,
    2.5914592945941E-13, -6.9582366282136707, -8.042647015359601E-14,
    2.58441760235E-13, -6.9547924607202321, -7.8620315109951847E-14,
    2.5773813393212E-13, -6.9513509465984713, -7.6812766428156727E-14,
    2.5703396470771E-13, -6.9479067789105216, -7.5005217746361608E-14,
    2.5632979548330004E-13, -206.94446261128223, -7.3197669064566476E-14,
    2.5562562625889E-13, -6.8096364021272082, -7.1390120382771356E-14,
    2.5492145703448E-13, -6.6748102266760263, -6.95839653391272E-14,
    2.5421783073160005E-13, -6.5399867718641325, -6.7776416657332086E-14,
    2.5351366150719E-13, -6.4051607643874586, -6.5968867975536966E-14,
    2.5280949228278003E-13, -6.270334891505656, -6.4161598021372025E-14,
    2.52105431642676E-13, -6.1355097175170092, -6.235432806720711E-14,
    2.51401371002572E-13, -6.0006847453311911, -6.0547058113042181E-14,
    2.5069731036246805E-13, -5.8658600085820067, -5.8739788158877253E-14,
    2.49993249722364E-13, -5.7310355409058547, -5.6932518204712325E-14,
    2.4928918908226E-13, -5.596211375941726, -5.5125248250547397E-14,
    2.4858512844215603E-13, -5.4613875473156464, -5.3317978296382456E-14,
    2.4788106780205204E-13, -5.3265640886536376, -6.8332629810260386,
    3.08151826998728E-15, -7.1173921879572269E-14, -6.8298896135164471,
    3.0946602135610397E-15, -7.1462869367922215E-14, -6.8265136447244394,
    3.10781229752336E-15, -7.1752039809581078E-14, -6.823137675992089,
    3.12096438148568E-15, -7.2041210251239941E-14, -6.8197617073245613,
    3.1341164654479996E-15, -7.23303806928988E-14, -6.8163883397656937,
    3.14725840902176E-15, -7.261932818124875E-14, -6.8130123709192159,
    3.16041049298408E-15, -7.2908498622907613E-14, -6.8096364021272091,
    3.1735625769464E-15, -7.3197669064566488E-14, -206.80626043339487,
    3.1867146609087197E-15, -7.3486839506225339E-14, -6.6715024231357818,
    3.19986674487104E-15, -7.3776009947884215E-14, -6.5367470476191691,
    3.2130086884448E-15, -7.4064957436234173E-14, -6.40198913808256,
    3.22616077240712E-15, -7.4354127877893023E-14, -6.2672313295018292,
    3.23931285636944E-15, -7.46432983195519E-14, -6.1324741757236945,
    3.252462912254048E-15, -7.4932424170548971E-14, -5.9977171901145843,
    3.265612968138656E-15, -7.5221550021546056E-14, -5.8629604063083027,
    3.2787630240232639E-15, -7.5510675872543141E-14, -5.7282038579386549,
    3.2919130799078718E-15, -7.5799801723540213E-14, -5.5934475786420386,
    3.3050631357924797E-15, -7.60889275745373E-14, -5.4586916020574474,
    3.3182131916770881E-15, -7.6378053425534382E-14, -5.3239359618109034,
    3.331363247561696E-15, -7.6667179276531467E-14, 3.7160897044296804E-15,
    -6.8332629810260386, 2.71660769486475E-13, 3.64995609575384E-15,
    -6.8298896135164471, 2.69371296574395E-13, 3.58377145805896E-15,
    -6.8265136447244394, 2.6708005709371E-13, 3.51758682036408E-15,
    -6.823137675992089, 2.64788817613025E-13, 3.4514021826691997E-15,
    -6.8197617073245613, 2.6249757813234E-13, 3.38526857399336E-15,
    -6.8163883397656937, 2.6020810522026003E-13, 3.31908393629848E-15,
    -6.8130123709192159, 2.5791686573957505E-13, 3.2528992986035996E-15,
    -6.8096364021272091, 2.5562562625888997E-13, 3.18671466090872E-15,
    -206.80626043339487, 2.53334386778205E-13, 3.12053002321384E-15,
    -6.6715024231357818, 2.5104314729752E-13, 3.054396414538E-15,
    -6.5367470476191691, 2.4875367438544E-13, 2.9882117768431196E-15,
    -6.40198913808256, 2.4646243490475504E-13, 2.9220271391482397E-15,
    -6.2672313295018292, 2.4417119542407E-13, 2.8558527072571681E-15,
    -6.1324741757236945, 2.41880309257106E-13, 2.7896782753660961E-15,
    -5.9977171901145843, 2.39589423090142E-13, 2.7235038434750241E-15,
    -5.8629604063083027, 2.37298536923178E-13, 2.6573294115839521E-15,
    -5.7282038579386549, 2.35007650756214E-13, 2.5911549796928805E-15,
    -5.5934475786420386, 2.3271676458925E-13, 2.5249805478018081E-15,
    -5.4586916020574474, 2.30425878422286E-13, 2.4588061159107361E-15,
    -5.3239359618109034, 2.28134992255322E-13, -8.7946958106792252E-14,
    2.58937255384995E-13, -6.8332629810260386, -8.6140488692262212E-14,
    2.58237301874835E-13, -6.8298896135164471, -8.4332625397011069E-14,
    2.5753680827709E-13, -6.8265136447244394, -8.2524762101759938E-14,
    2.56836314679345E-13, -6.823137675992089, -8.0716898806508808E-14,
    2.561358210816E-13, -6.8197617073245613, -7.8910429391978755E-14,
    2.5543586757144E-13, -6.8163883397656937, -7.7102566096727625E-14,
    2.5473537397369503E-13, -6.8130123709192159, -7.5294702801476482E-14,
    2.5403488037595E-13, -6.8096364021272091, -7.3486839506225339E-14,
    2.53334386778205E-13, -206.80626043339487, -7.1678976210974209E-14,
    2.5263389318046003E-13, -6.6715024231357818, -6.9872506796444156E-14,
    2.519339396703E-13, -6.5367470476191691, -6.8064643501193025E-14,
    2.51233446072555E-13, -6.40198913808256, -6.62567802059419E-14,
    2.5053295247481E-13, -6.2672313295018292, -6.4449195686834967E-14,
    2.49832566894582E-13, -6.1324741757236945, -6.2641611167728053E-14,
    2.49132181314354E-13, -5.9977171901145843, -6.0834026648621138E-14,
    2.4843179573412604E-13, -5.8629604063083027, -5.9026442129514223E-14,
    2.4773141015389804E-13, -5.7282038579386549, -5.72188576104073E-14,
    2.4703102457367003E-13, -5.5934475786420386, -5.5411273091300381E-14,
    2.46330638993442E-13, -5.4586916020574474, -5.3603688572193453E-14,
    2.45630253413214E-13, -5.3239359618109034, -6.6979597537587221,
    3.0182119425781598E-15, -6.9368575006829E-14, -6.6946544991652495,
    3.03099430545688E-15, -6.9657208124293051E-14, -6.6913466955705356,
    3.0437865312699197E-15, -6.9946063952495916E-14, -6.6880388920302911,
    3.05657875708296E-15, -7.0234919780698756E-14, -6.684731088549702,
    3.0693709828959997E-15, -7.05237756089016E-14, -6.6814258339147337,
    3.0821533457747198E-15, -7.0812408726365664E-14, -6.678118030270741,
    3.09494557158776E-15, -7.11012645545685E-14, -6.6748102266760263,
    3.1077377974008002E-15, -7.1390120382771369E-14, -6.6715024231357818,
    3.1205300232138397E-15, -7.1678976210974209E-14, -206.6681946196552,
    3.13332224902688E-15, -7.1967832039177061E-14, -6.533507323428676,
    3.1461046119056004E-15, -7.225646515664113E-14, -6.3988175118269393,
    3.1588968377186398E-15, -7.2545320984843969E-14, -6.2641277675420923,
    3.17168906353168E-15, -7.2834176813046822E-14, -6.1294386339692837,
    3.184479316757856E-15, -7.3122988099101907E-14, -5.9947496349316935,
    3.1972695699840319E-15, -7.3411799385157E-14, -5.8600608040631279,
    3.2100598232102079E-15, -7.3700610671212091E-14, -5.7253721749973909,
    3.2228500764363839E-15, -7.3989421957267189E-14, -5.5906837813682877,
    3.23564032966256E-15, -7.4278233243322275E-14, -5.4559956568122159,
    3.2484305828887358E-15, -7.4567044529377373E-14, -5.3213078349681693,
    3.2612208361149117E-15, -7.4855855815432458E-14, 3.73212009867776E-15,
    -6.6979597537587221, 2.709308765433E-13, 3.66562690930688E-15,
    -6.6946544991652495, 2.6864507642394003E-13, 3.5990824134627196E-15,
    -6.6913466955705356, 2.6635751256992E-13, 3.5325379176185603E-15,
    -6.6880388920302911, 2.640699487159E-13, 3.4659934217743997E-15,
    -6.684731088549702, 2.6178238486188004E-13, 3.39950023240352E-15,
    -6.6814258339147337, 2.5949658474252E-13, 3.33295573655936E-15,
    -6.678118030270741, 2.5720902088850004E-13, 3.2664112407151998E-15,
    -6.6748102266760263, 2.5492145703448E-13, 3.19986674487104E-15,
    -6.6715024231357818, 2.5263389318046003E-13, 3.13332224902688E-15,
    -206.6681946196552, 2.5034632932644E-13, 3.066829059656E-15,
    -6.533507323428676, 2.4806052920708003E-13, 3.0002845638118397E-15,
    -6.3988175118269393, 2.4577296535306005E-13, 2.9337400679676795E-15,
    -6.2641277675420923, 2.4348540149903997E-13, 2.8672058334181764E-15,
    -6.1294386339692837, 2.41198190391952E-13, 2.800671598868672E-15,
    -5.9947496349316935, 2.3891097928486403E-13, 2.734137364319168E-15,
    -5.8600608040631279, 2.36623768177776E-13, 2.6676031297696641E-15,
    -5.7253721749973909, 2.34336557070688E-13, 2.6010688952201605E-15,
    -5.5906837813682877, 2.3204934596359997E-13, 2.5345346606706562E-15,
    -5.4559956568122159, 2.29762134856512E-13, 2.4680004261211522E-15,
    -5.3213078349681693, 2.2747492374942404E-13, -8.8238644970958984E-14,
    2.5661661655888E-13, -6.6979597537587221, -8.6431861185543067E-14,
    2.5592033584144E-13, -6.6946544991652495, -8.46236832768359E-14,
    2.5522351787036E-13, -6.6913466955705356, -8.2815505368128759E-14,
    2.5452669989928E-13, -6.6880388920302911, -8.10073274594216E-14,
    2.538298819282E-13, -6.684731088549702, -7.9200543674005676E-14,
    2.5313360121076E-13, -6.6814258339147337, -7.7392365765298522E-14,
    2.5243678323968E-13, -6.678118030270741, -7.5584187856591368E-14,
    2.5173996526860004E-13, -6.6748102266760263, -7.37760099478842E-14,
    2.5104314729752E-13, -6.6715024231357818, -7.1967832039177061E-14,
    2.5034632932644E-13, -206.6681946196552, -7.0161048253761118E-14,
    2.49650048609E-13, -6.533507323428676, -6.8352870345053977E-14,
    2.4895323063792E-13, -6.3988175118269393, -6.6544692436346824E-14,
    2.4825641266684E-13, -6.2641277675420923, -6.473679335229791E-14,
    2.47559702146488E-13, -6.1294386339692837, -6.2928894268249008E-14,
    2.4686299162613596E-13, -5.9947496349316935, -6.11209951842001E-14,
    2.4616628110578403E-13, -5.8600608040631279, -5.9313096100151193E-14,
    2.4546957058543206E-13, -5.7253721749973909, -5.7505197016102279E-14,
    2.4477286006508003E-13, -5.5906837813682877, -5.5697297932053365E-14,
    2.44076149544728E-13, -5.4559956568122159, -5.3889398848004451E-14,
    2.4337943902437603E-13, -5.3213078349681693, -6.56265984833564,
    2.954954424982E-15, -6.7564620074619848E-14, -6.5594226207965853,
    2.9673774844059997E-15, -6.785293906358112E-14, -6.55618289649975,
    2.979810129524E-15, -6.81414805208981E-14, -6.5529431722547855,
    2.992242774642E-15, -6.8430021978215045E-14, -6.5497034480642924,
    3.00467541976E-15, -6.8718563435532E-14, -6.5464662204967086,
    3.0170984791840002E-15, -6.9006882424493292E-14, -6.5432264961583755,
    3.029531124302E-15, -6.9295423881810242E-14, -6.5399867718641325,
    3.0419637694200002E-15, -6.95839653391272E-14, -6.5367470476191691,
    3.054396414538E-15, -6.9872506796444168E-14, -6.5335073234286751,
    3.066829059656E-15, -7.0161048253761131E-14, -206.53027009586108,
    3.07925211908E-15, -7.0449367242722415E-14, -6.395648329931209,
    3.091684764198E-15, -7.0737908700039365E-14, -6.2610265976792228,
    3.104117409316E-15, -7.1026450157356328E-14, -6.126405432056778,
    3.1165481372952003E-15, -7.131494712100215E-14, -5.9917843673357476,
    3.1289788652744E-15, -7.1603444084647972E-14, -5.8571634371499366,
    3.1414095932536E-15, -7.18919410482938E-14, -5.7225426751331492,
    3.1538403212328E-15, -7.2180438011939629E-14, -5.58792211491919,
    3.166271049212E-15, -7.2468934975585438E-14, -5.4533017901418654,
    3.1787017771912E-15, -7.2757431939231273E-14, -5.3186817344375728,
    3.1911325051704E-15, -7.30459289028771E-14, 3.7481381333312E-15,
    -6.56265984833564, 2.702015463549E-13, 3.6812856405056E-15,
    -6.5594226207965853, 2.679194161965E-13, 3.6143815639664E-15,
    -6.55618289649975, 2.656355251352E-13, 3.5474774874272002E-15,
    -6.5529431722547855, 2.633516340739E-13, 3.4805734108879997E-15,
    -6.5497034480642924, 2.610677430126E-13, 3.4137209180624E-15,
    -6.5464662204967086, 2.587856128542E-13, 3.3468168415232E-15,
    -6.5432264961583755, 2.5650172179290005E-13, 3.279912764984E-15,
    -6.5399867718641325, 2.542178307316E-13, 3.2130086884448E-15,
    -6.5367470476191691, 2.519339396703E-13, 3.1461046119056E-15,
    -6.5335073234286751, 2.49650048609E-13, 3.07925211908E-15,
    -206.53027009586108, 2.4736791845060004E-13, 3.0123480425407996E-15,
    -6.395648329931209, 2.4508402738930004E-13, 2.9454439660015998E-15,
    -6.2610265976792228, 2.42800136328E-13, 2.87855020620512E-15,
    -6.126405432056778, 2.4051659744728E-13, 2.81165644640864E-15,
    -5.9917843673357476, 2.3823305856656003E-13, 2.74476268661216E-15,
    -5.8571634371499366, 2.3594951968584E-13, 2.67786892681568E-15,
    -5.7225426751331492, 2.3366598080512E-13, 2.6109751670192005E-15,
    -5.58792211491919, 2.313824419244E-13, 2.54408140722272E-15,
    -5.4533017901418654, 2.2909890304368E-13, 2.4771876474262403E-15,
    -5.3186817344375728, 2.2681536416296002E-13, -8.8530106941629836E-14,
    2.5429776696856E-13, -6.56265984833564, -8.6723009027711121E-14,
    2.5360515621208E-13, -6.5594226207965853, -8.4914516748118078E-14,
    2.5291201103372E-13, -6.55618289649975, -8.3106024468525048E-14,
    2.5221886585536E-13, -6.5529431722547855, -8.1297532188932E-14,
    2.5152572067700003E-13, -6.5497034480642924, -7.949043427501329E-14,
    2.5083310992052E-13, -6.5464662204967086, -7.7681941995420247E-14,
    2.5013996474216E-13, -6.5432264961583755, -7.58734497158272E-14,
    2.494468195638E-13, -6.5399867718641325, -7.406495743623416E-14,
    2.4875367438544E-13, -6.5367470476191691, -7.2256465156641117E-14,
    2.4806052920708003E-13, -6.5335073234286751, -7.04493672427224E-14,
    2.4736791845060004E-13, -206.53027009586108, -6.8640874963129359E-14,
    2.4667477327224E-13, -6.395648329931209, -6.6832382683536329E-14,
    2.4598162809388E-13, -6.2610265976792228, -6.5024169277078139E-14,
    2.45288589799896E-13, -6.126405432056778, -6.3215955870619974E-14,
    2.44595551505912E-13, -5.9917843673357476, -6.14077424641618E-14,
    2.4390251321192805E-13, -5.8571634371499366, -5.9599529057703632E-14,
    2.4320947491794406E-13, -5.7225426751331492, -5.7791315651245441E-14,
    2.4251643662396E-13, -5.58792211491919, -5.598310224478727E-14,
    2.4182339832997603E-13, -5.4533017901418654, -5.417488883832908E-14,
    2.4113036003599204E-13, -5.3186817344375728, -6.42735725981494,
    2.89164809757288E-15, -6.5759273201876587E-14, -6.424188077942552,
    2.9037115763018395E-15, -6.6047277819951969E-14, -6.4210164515831867,
    2.91578436327056E-15, -6.6335504663812921E-14, -6.4178448252756946,
    2.92785715023928E-15, -6.662373150767386E-14, -6.4146731990200738,
    2.9399299372079996E-15, -6.69119583515348E-14, -6.4115040171243445,
    2.95199341593696E-15, -6.71999629696102E-14, -6.4083323907364491,
    2.96406620290568E-15, -6.7488189813471134E-14, -6.4051607643874586,
    2.9761389898744E-15, -6.7776416657332086E-14, -6.40198913808256,
    2.98821177684312E-15, -6.8064643501193025E-14, -6.3988175118269393,
    3.00028456381184E-15, -6.8352870345053977E-14, -6.395648329931209,
    3.0123480425408E-15, -6.8640874963129372E-14, -206.39247670354331,
    3.02442082950952E-15, -6.8929101806990311E-14, -6.2579230356027757,
    3.03649361647824E-15, -6.9217328650851251E-14, -6.123369890201217,
    3.0485645417990083E-15, -6.9505511049555086E-14, -5.9888168120672685,
    3.0606354671197761E-15, -6.9793693448258922E-14, -5.8542638348347342,
    3.072706392440544E-15, -7.0081875846962757E-14, -5.7197109921374194,
    3.0847773177613118E-15, -7.0370058245666593E-14, -5.5851583176091282,
    3.0968482430820797E-15, -7.0658240644370416E-14, -5.4506058448836665,
    3.1089191684028479E-15, -7.0946423043074251E-14, -5.3160536075948377,
    3.1209900937236158E-15, -7.1234605441778087E-14, 3.76416852757928E-15,
    -6.42735725981494, 2.6947165341172497E-13, 3.6969564540586395E-15,
    -6.424188077942552, 2.6719319604604504E-13, 3.62969251937016E-15,
    -6.4210164515831867, 2.6491298061141E-13, 3.5624285846816805E-15,
    -6.4178448252756946, 2.62632765176775E-13, 3.4951646499931998E-15,
    -6.4146731990200738, 2.6035254974214003E-13, 3.42795257647256E-15,
    -6.4115040171243445, 2.5807409237646E-13, 3.36068864178408E-15,
    -6.4083323907364491, 2.5579387694182505E-13, 3.2934247070955996E-15,
    -6.4051607643874586, 2.5351366150719E-13, 3.22616077240712E-15,
    -6.40198913808256, 2.51233446072555E-13, 3.15889683771864E-15,
    -6.3988175118269393, 2.4895323063792003E-13, 3.0916847641980003E-15,
    -6.395648329931209, 2.4667477327224E-13, 3.0244208295095196E-15,
    -206.39247670354331, 2.4439455783760505E-13, 2.9571568948210397E-15,
    -6.2579230356027757, 2.4211434240297E-13, 2.8899033323661284E-15,
    -6.123369890201217, 2.39834478582126E-13, 2.8226497699112159E-15,
    -5.9888168120672685, 2.37554614761282E-13, 2.7553962074563043E-15,
    -5.8542638348347342, 2.35274750940438E-13, 2.6881426450013922E-15,
    -5.7197109921374194, 2.32994887119594E-13, 2.6208890825464805E-15,
    -5.5851583176091282, 2.3071502329875E-13, 2.5536355200915681E-15,
    -5.4506058448836665, 2.28435159477906E-13, 2.4863819576366564E-15,
    -5.3160536075948377, 2.26155295657062E-13, -8.8821793805796568E-14,
    2.5197712814244497E-13, -6.42735725981494, -8.7014381520991977E-14,
    2.51288190178685E-13, -6.424188077942552, -8.520557462794291E-14,
    2.5059872062699E-13, -6.4210164515831867, -8.3396767734893856E-14,
    2.49909251075295E-13, -6.4178448252756946, -8.1587960841844815E-14,
    2.4921978152360003E-13, -6.4146731990200738, -7.97805485570402E-14,
    2.4853084355984E-13, -6.4115040171243445, -7.7971741663991144E-14,
    2.4784137400814503E-13, -6.4083323907364491, -7.616293477094209E-14,
    2.4715190445645004E-13, -6.4051607643874586, -7.4354127877893023E-14,
    2.46462434904755E-13, -6.40198913808256, -7.2545320984843969E-14,
    2.4577296535306E-13, -6.3988175118269393, -7.0737908700039365E-14,
    2.4508402738930004E-13, -6.395648329931209, -6.8929101806990311E-14,
    2.44394557837605E-13, -206.39247670354331, -6.7120294913941257E-14,
    2.4370508828591E-13, -6.2579230356027757, -6.5311766942541081E-14,
    2.43015725051802E-13, -6.123369890201217, -6.3503238971140917E-14,
    2.42326361817694E-13, -5.9888168120672685, -6.1694710999740753E-14,
    2.4163699858358604E-13, -5.8542638348347342, -5.9886183028340589E-14,
    2.4094763534947803E-13, -5.7197109921374194, -5.8077655056940412E-14,
    2.4025827211537E-13, -5.5851583176091282, -5.6269127085540254E-14,
    2.39568908881262E-13, -5.4506058448836665, -5.4460599114140078E-14,
    2.38879545647154E-13, -5.3160536075948377, -6.292055041354276,
    2.82834177016376E-15, -6.3953926329133326E-14, -6.2889538715173421,
    2.8400456681976797E-15, -6.4241616576322818E-14, -6.2858503094642364,
    2.85175859701712E-15, -6.4529528806727747E-14, -6.2827467474500347,
    2.86347152583656E-15, -6.4817441037132675E-14, -6.2796431854877044,
    2.8751844546559997E-15, -6.51053532675376E-14, -6.2765420156248366,
    2.88688835268992E-15, -6.5393043514727109E-14, -6.2734384535483887,
    2.89860128150936E-15, -6.5680955745132025E-14, -6.2703348915056569,
    2.9103142103288E-15, -6.5968867975536966E-14, -6.2672313295018292,
    2.9220271391482397E-15, -6.6256780205941882E-14, -6.2641277675420923,
    2.93374006796768E-15, -6.6544692436346824E-14, -6.2610265976792236,
    2.9454439660016E-15, -6.6832382683536329E-14, -6.2579230356027757,
    2.95715689482104E-15, -6.7120294913941245E-14, -206.25481947356005,
    2.96886982364048E-15, -6.7408207144346186E-14, -6.1203343483741852,
    2.9805809463028162E-15, -6.7696074978108022E-14, -5.9858492568221315,
    2.9922920689651521E-15, -6.7983942811869871E-14, -5.851364232537688,
    3.004003191627488E-15, -6.8271810645631708E-14, -5.7168793091546579,
    3.0157143142898239E-15, -6.8559678479393557E-14, -5.5823945203068472,
    3.0274254369521598E-15, -6.8847546313155393E-14, -5.447909899628061,
    3.039136559614496E-15, -6.9135414146917242E-14, -5.3134254807521035,
    3.0508476822768319E-15, -6.9423281980679078E-14, 3.78019892182736E-15,
    -6.292055041354276, 2.6874176046854997E-13, 3.71262726761168E-15,
    -6.2889538715173421, 2.6646697589559E-13, 3.6450034747739196E-15,
    -6.2858503094642364, 2.6419043608762E-13, 3.57737968193616E-15,
    -6.2827467474500347, 2.6191389627965E-13, 3.5097558890984E-15,
    -6.2796431854877044, 2.5963735647168E-13, 3.44218423488272E-15,
    -6.2765420156248366, 2.5736257189872E-13, 3.3745604420449603E-15,
    -6.2734384535483887, 2.5508603209075004E-13, 3.3069366492071998E-15,
    -6.2703348915056569, 2.5280949228278E-13, 3.23931285636944E-15,
    -6.2672313295018292, 2.5053295247481E-13, 3.17168906353168E-15,
    -6.2641277675420923, 2.4825641266684E-13, 3.104117409316E-15,
    -6.2610265976792236, 2.4598162809388E-13, 3.0364936164782396E-15,
    -6.2579230356027757, 2.4370508828591006E-13, 2.9688698236404796E-15,
    -206.25481947356005, 2.4142854847794E-13, 2.9012564585271363E-15,
    -6.1203343483741852, 2.39152359716972E-13, 2.8336430934137919E-15,
    -5.9858492568221315, 2.36876170956004E-13, 2.7660297283004482E-15,
    -5.851364232537688, 2.3459998219503605E-13, 2.6984163631871042E-15,
    -5.7168793091546579, 2.3232379343406803E-13, 2.63080299807376E-15,
    -5.5823945203068472, 2.300476046731E-13, 2.5631896329604161E-15,
    -5.447909899628061, 2.27771415912132E-13, 2.4955762678470721E-15,
    -5.3134254807521035, 2.25495227151164E-13, -8.9113480669963313E-14,
    2.4965648931632995E-13, -6.292055041354276, -8.7305754014272819E-14,
    2.4897122414529E-13, -6.2889538715173421, -8.5496632507767742E-14,
    2.4828543022026003E-13, -6.2858503094642364, -8.3687511001262677E-14,
    2.4759963629523E-13, -6.2827467474500347, -8.1878389494757612E-14,
    2.4691384237020004E-13, -6.2796431854877044, -8.00706628390671E-14,
    2.4622857719916E-13, -6.2765420156248366, -7.8261541332562041E-14,
    2.4554278327413E-13, -6.2734384535483887, -7.6452419826056964E-14,
    2.448569893491E-13, -6.2703348915056569, -7.4643298319551887E-14,
    2.4417119542407E-13, -6.2672313295018292, -7.2834176813046822E-14,
    2.4348540149904E-13, -6.2641277675420923, -7.1026450157356316E-14,
    2.4280013632800005E-13, -6.2610265976792236, -6.9217328650851251E-14,
    2.4211434240297E-13, -6.2579230356027757, -6.7408207144346186E-14,
    2.4142854847794E-13, -206.25481947356005, -6.5599364608004023E-14,
    2.40742860303708E-13, -6.1203343483741852, -6.379052207166186E-14,
    2.40057172129476E-13, -5.9858492568221315, -6.1981679535319709E-14,
    2.39371483955244E-13, -5.851364232537688, -6.0172836998977559E-14,
    2.3868579578101205E-13, -5.7168793091546579, -5.83639944626354E-14,
    2.3800010760678003E-13, -5.5823945203068472, -5.6555151926293238E-14,
    2.37314419432548E-13, -5.447909899628061, -5.4746309389951075E-14,
    2.3662873125831603E-13, -5.3134254807521035, -6.1567538304019642,
    2.7650452047172319E-15, -6.2148857844496882E-14, -6.1537206285128008,
    2.7763895775041757E-15, -6.2436233769277111E-14, -6.1506850866831755,
    2.787742703665184E-15, -6.2723831434740053E-14, -6.1476495448794868,
    2.7990958298261919E-15, -6.3011429100202983E-14, -6.1446140031147012,
    2.8104489559871998E-15, -6.3299026765665912E-14, -6.1415808011996029,
    2.8217933287741439E-15, -6.3586402690446154E-14, -6.1385452593440419,
    2.8331464549351522E-15, -6.38740003559091E-14, -6.13550971751701,
    2.84449958109616E-15, -6.4161598021372038E-14, -6.1324741757236945,
    2.8558527072571677E-15, -6.4449195686834967E-14, -6.1294386339692837,
    2.867205833418176E-15, -6.473679335229791E-14, -6.126405432056778,
    2.87855020620512E-15, -6.5024169277078151E-14, -6.123369890201217,
    2.889903332366128E-15, -6.5311766942541081E-14, -6.1203343483741852,
    2.9012564585271359E-15, -6.5599364608004023E-14, -206.11729927454044,
    2.9126078340133312E-15, -6.5886917925330425E-14, -5.9828821591151327,
    2.9239592094995266E-15, -6.6174471242656814E-14, -5.8484650773236364,
    2.9353105849857219E-15, -6.6462024559983216E-14, -5.71404806279975,
    2.9466619604719164E-15, -6.6749577877309618E-14, -5.5796311491772785,
    2.9580133359581121E-15, -6.7037131194636007E-14, -5.4452143700900253,
    2.969364711444307E-15, -6.7324684511962409E-14, -5.3107977591717965,
    2.9807160869305023E-15, -6.76122378292888E-14, 3.7962268441565124E-15,
    -6.1567538304019642, 2.6801198007633E-13, 3.728295664693856E-15,
    -6.1537206285128008, 2.6574086772973803E-13, 3.6603120691976639E-15,
    -6.1506850866831755, 2.63468002981644E-13, 3.5923284737014719E-15,
    -6.1476495448794868, 2.6119513823355003E-13, 3.52434487820528E-15,
    -6.1446140031147012, 2.58922273485456E-13, 3.4564136987426238E-15,
    -6.1415808011996029, 2.5665116113886404E-13, 3.3884301032464321E-15,
    -6.1385452593440419, 2.5437829639077E-13, 3.3204465077502397E-15,
    -6.13550971751701, 2.52105431642676E-13, 3.252462912254048E-15,
    -6.1324741757236945, 2.49832566894582E-13, 3.184479316757856E-15,
    -6.1294386339692837, 2.4755970214648803E-13, 3.1165481372952003E-15,
    -6.126405432056778, 2.45288589799896E-13, 3.0485645417990075E-15,
    -6.123369890201217, 2.4301572505180205E-13, 2.9805809463028159E-15,
    -6.1203343483741852, 2.4074286030370797E-13, 2.9126078340133312E-15,
    -206.11729927454044, 2.3847034603591442E-13, 2.8446347217238462E-15,
    -5.9828821591151327, 2.3619783176812082E-13, 2.7766616094343616E-15,
    -5.8484650773236364, 2.3392531750032721E-13, 2.708688497144877E-15,
    -5.71404806279975, 2.3165280323253361E-13, 2.6407153848553924E-15,
    -5.5796311491772785, 2.2938028896474E-13, 2.5727422725659074E-15,
    -5.4452143700900253, 2.2710777469694641E-13, 2.5047691602764224E-15,
    -5.3107977591717965, 2.248352604291528E-13, -8.9405122555430869E-14,
    2.47336208337374E-13, -6.1567538304019642, -8.7597081577331112E-14,
    2.46654615392702E-13, -6.1537206285128008, -8.5787645505884048E-14,
    2.4597249652754803E-13, -6.1506850866831755, -8.3978209434436984E-14,
    2.4529037766239397E-13, -6.1476495448794868, -8.2168773362989933E-14,
    2.4460825879724E-13, -6.1446140031147012, -8.0360732384890163E-14,
    2.43926665852568E-13, -6.1415808011996029, -7.8551296313443112E-14,
    2.43244546987414E-13, -6.1385452593440419, -7.6741860241996035E-14,
    2.4256242812226003E-13, -6.13550971751701, -7.4932424170548971E-14,
    2.41880309257106E-13, -6.1324741757236945, -7.3122988099101907E-14,
    2.41198190391952E-13, -6.1294386339692837, -7.131494712100215E-14,
    2.4051659744728E-13, -6.126405432056778, -6.9505511049555086E-14,
    2.39834478582126E-13, -6.123369890201217, -6.7696074978108022E-14,
    2.39152359716972E-13, -6.1203343483741852, -6.5886917925330412E-14,
    2.3847034603591442E-13, -206.11729927454044, -6.4077760872552815E-14,
    2.3778833235485676E-13, -5.9828821591151327, -6.2268603819775217E-14,
    2.3710631867379925E-13, -5.8484650773236364, -6.045944676699762E-14,
    2.3642430499274164E-13, -5.71404806279975, -5.8650289714220009E-14,
    2.3574229131168403E-13, -5.5796311491772785, -5.6841132661442406E-14,
    2.3506027763062642E-13, -5.4452143700900253, -5.5031975608664795E-14,
    2.3437826394956881E-13, -5.3107977591717965, -6.0214530566917084,
    2.701748639270704E-15, -6.0343789359860439E-14, -6.0184877891113242,
    2.7127334868106716E-15, -6.06308509622314E-14, -6.0155202338661873,
    2.723726810313248E-15, -6.091813406275236E-14, -6.0125526786469861,
    2.7347201338158239E-15, -6.120541716327329E-14, -6.0095851234537214,
    2.7457134573184E-15, -6.1492700263794233E-14, -6.0066198558474015,
    2.7566983048583678E-15, -6.1779761866165211E-14, -6.0036523005763289,
    2.7676916283609442E-15, -6.2067044966686154E-14, -6.000684745331192,
    2.77868495186352E-15, -6.235432806720711E-14, -5.9977171901145843,
    2.7896782753660957E-15, -6.2641611167728053E-14, -5.9947496349316935,
    2.800671598868672E-15, -6.2928894268249008E-14, -5.9917843673357476,
    2.81165644640864E-15, -6.3215955870619974E-14, -5.9888168120672685,
    2.8226497699112159E-15, -6.3503238971140917E-14, -5.9858492568221315,
    2.8336430934137919E-15, -6.3790522071661872E-14, -5.9828821591151327,
    2.8446347217238466E-15, -6.4077760872552815E-14, -205.97991506140812,
    2.855626350033901E-15, -6.436499967344377E-14, -5.8455659221095857,
    2.8666179783439553E-15, -6.4652238474334725E-14, -5.7112168164448427,
    2.8776096066540093E-15, -6.493947727522568E-14, -5.5768677780477089,
    2.8886012349640637E-15, -6.5226716076116622E-14, -5.44251884055199,
    2.8995928632741184E-15, -6.5513954877007577E-14, -5.30817003759149,
    2.9105844915841728E-15, -6.5801193677898532E-14, 3.8122547664856644E-15,
    -6.0214530566917084, 2.6728219968410996E-13, 3.7439640617760314E-15,
    -6.0184877891113242, 2.6501475956388604E-13, 3.6756206636214074E-15,
    -6.0155202338661873, 2.62745569875668E-13, 3.6072772654667842E-15,
    -6.0125526786469861, 2.6047638018745E-13, 3.5389338673121598E-15,
    -6.0095851234537214, 2.58207190499232E-13, 3.470643162602528E-15,
    -6.0066198558474015, 2.5593975037900803E-13, 3.402299764447904E-15,
    -6.0036523005763289, 2.5367056069079004E-13, 3.33395636629328E-15,
    -6.000684745331192, 2.51401371002572E-13, 3.265612968138656E-15,
    -5.9977171901145843, 2.49132181314354E-13, 3.1972695699840319E-15,
    -5.9947496349316935, 2.46862991626136E-13, 3.1289788652744E-15,
    -5.9917843673357476, 2.4459555150591204E-13, 3.0606354671197754E-15,
    -5.9888168120672685, 2.4232636181769404E-13, 2.9922920689651517E-15,
    -5.9858492568221315, 2.40057172129476E-13, 2.9239592094995266E-15,
    -5.9828821591151327, 2.3778833235485681E-13, 2.8556263500339006E-15,
    -205.97991506140812, 2.3551949258023762E-13, 2.7872934905682754E-15,
    -5.8455659221095857, 2.3325065280561843E-13, 2.71896063110265E-15,
    -5.7112168164448427, 2.3098181303099919E-13, 2.6506277716370243E-15,
    -5.5768677780477089, 2.2871297325638E-13, 2.5822949121713983E-15,
    -5.44251884055199, 2.2644413348176081E-13, 2.5139620527057731E-15,
    -5.30817003759149, 2.2417529370714162E-13, -8.9696764440898425E-14,
    2.4501592735841796E-13, -6.0214530566917084, -8.78884091403894E-14,
    2.44338006640114E-13, -6.0184877891113242, -8.6078658504000341E-14,
    2.4365956283483603E-13, -6.0155202338661873, -8.4268907867611291E-14,
    2.42981119029558E-13, -6.0125526786469861, -8.2459157231222253E-14,
    2.4230267522428E-13, -6.0095851234537214, -8.065080193071322E-14,
    2.41624754505976E-13, -6.0066198558474015, -7.8841051294324169E-14,
    2.4094631070069803E-13, -6.0036523005763289, -7.7031300657935106E-14,
    2.4026786689542005E-13, -6.000684745331192, -7.5221550021546043E-14,
    2.39589423090142E-13, -5.9977171901145843, -7.3411799385157E-14,
    2.38910979284864E-13, -5.9947496349316935, -7.1603444084647972E-14,
    2.3823305856656003E-13, -5.9917843673357476, -6.9793693448258922E-14,
    2.37554614761282E-13, -5.9888168120672685, -6.7983942811869871E-14,
    2.36876170956004E-13, -5.9858492568221315, -6.6174471242656814E-14,
    2.3619783176812082E-13, -5.9828821591151327, -6.436499967344377E-14,
    2.3551949258023757E-13, -205.97991506140812, -6.2555528104230725E-14,
    2.3484115339235442E-13, -5.8455659221095857, -6.0746056535017681E-14,
    2.3416281420447123E-13, -5.7112168164448427, -5.8936584965804624E-14,
    2.3348447501658803E-13, -5.5768677780477089, -5.7127113396591579E-14,
    2.3280613582870483E-13, -5.44251884055199, -5.5317641827378522E-14,
    2.3212779664082164E-13, -5.30817003759149, -5.8861527538625005,
    2.6384520738241761E-15, -5.8538720875224E-14, -5.8832553869519026,
    2.649077396117168E-15, -5.88254681551857E-14, -5.8803557846522621,
    2.6597109169613119E-15, -5.9112436690764654E-14, -5.8774561823785572,
    2.6703444378054559E-15, -5.939940522634361E-14, -5.8745565801307889,
    2.6809779586496E-15, -5.9686373761922554E-14, -5.8716592132072236,
    2.6916032809425917E-15, -5.9973121041884269E-14, -5.8687596108816473,
    2.7022368017867361E-15, -6.0260089577463213E-14, -5.8658600085820067,
    2.71287032263088E-15, -6.0547058113042181E-14, -5.8629604063083027,
    2.7235038434750241E-15, -6.0834026648621138E-14, -5.8600608040631279,
    2.734137364319168E-15, -6.11209951842001E-14, -5.8571634371499366,
    2.7447626866121603E-15, -6.14077424641618E-14, -5.8542638348347342,
    2.7553962074563039E-15, -6.1694710999740753E-14, -5.8513642325376871,
    2.7660297283004478E-15, -6.1981679535319709E-14, -5.8484650773236364,
    2.7766616094343616E-15, -6.2268603819775217E-14, -5.8455659221095857,
    2.7872934905682754E-15, -6.2555528104230725E-14, -205.84266676689555,
    2.7979253717021888E-15, -6.2842452388686233E-14, -5.7083855700899351,
    2.8085572528361022E-15, -6.3129376673141741E-14, -5.5741044069181394,
    2.819189133970016E-15, -6.3416300957597237E-14, -5.4398233110139547,
    2.8298210151039298E-15, -6.3703225242052745E-14, -5.3055423160111834,
    2.8404528962378432E-15, -6.3990149526508253E-14, 3.8282826888148164E-15,
    -5.8861527538625005, 2.6655241929189E-13, 3.7596324588582077E-15,
    -5.8832553869519026, 2.6428865139803404E-13, 3.6909292580451517E-15,
    -5.8803557846522621, 2.62023136769692E-13, 3.6222260572320957E-15,
    -5.8774561823785572, 2.5975762214135E-13, 3.55352285641904E-15,
    -5.8745565801307889, 2.57492107513008E-13, 3.4848726264624318E-15,
    -5.8716592132072236, 2.5522833961915203E-13, 3.4161694256493763E-15,
    -5.8687596108816473, 2.5296282499081E-13, 3.34746622483632E-15,
    -5.8658600085820067, 2.50697310362468E-13, 3.2787630240232639E-15,
    -5.8629604063083027, 2.48431795734126E-13, 3.2100598232102079E-15,
    -5.8600608040631279, 2.4616628110578403E-13, 3.1414095932536E-15,
    -5.8571634371499366, 2.43902513211928E-13, 3.0727063924405436E-15,
    -5.8542638348347342, 2.4163699858358604E-13, 3.0040031916274876E-15,
    -5.8513642325376871, 2.3937148395524397E-13, 2.9353105849857219E-15,
    -5.8484650773236364, 2.371063186737992E-13, 2.8666179783439549E-15,
    -5.8455659221095857, 2.3484115339235442E-13, 2.7979253717021888E-15,
    -205.84266676689555, 2.3257598811090965E-13, 2.7292327650604227E-15,
    -5.7083855700899351, 2.3031082282946482E-13, 2.6605401584186561E-15,
    -5.5741044069181394, 2.2804565754802E-13, 2.5918475517768896E-15,
    -5.4398233110139547, 2.2578049226657519E-13, 2.5231549451351235E-15,
    -5.3055423160111834, 2.2351532698513042E-13, -8.9988406326365981E-14,
    2.42695646379462E-13, -5.8861527538625005, -8.81797367034477E-14,
    2.42021397887526E-13, -5.8832553869519026, -8.6369671502116648E-14,
    2.4134662914212403E-13, -5.8803557846522621, -8.4559606300785611E-14,
    2.40671860396722E-13, -5.8774561823785572, -8.2749541099454573E-14,
    2.3999709165132003E-13, -5.8745565801307889, -8.0940871476536264E-14,
    2.39322843159384E-13, -5.8716592132072236, -7.9130806275205227E-14,
    2.38648074413982E-13, -5.8687596108816473, -7.7320741073874178E-14,
    2.3797330566858E-13, -5.8658600085820067, -7.5510675872543128E-14,
    2.37298536923178E-13, -5.8629604063083027, -7.3700610671212091E-14,
    2.36623768177776E-13, -5.8600608040631279, -7.18919410482938E-14,
    2.3594951968584005E-13, -5.8571634371499366, -7.0081875846962757E-14,
    2.35274750940438E-13, -5.8542638348347342, -6.827181064563172E-14,
    2.34599982195036E-13, -5.8513642325376871, -6.6462024559983216E-14,
    2.3392531750032721E-13, -5.8484650773236364, -6.4652238474334725E-14,
    2.3325065280561838E-13, -5.8455659221095857, -6.2842452388686233E-14,
    2.3257598811090965E-13, -205.84266676689555, -6.1032666303037742E-14,
    2.3190132341620081E-13, -5.7083855700899351, -5.9222880217389238E-14,
    2.3122665872149203E-13, -5.5741044069181394, -5.7413094131740746E-14,
    2.305519940267832E-13, -5.4398233110139547, -5.5603308046092242E-14,
    2.2987732933207441E-13, -5.3055423160111834, -5.7508529555403651,
    2.5751555083776478E-15, -5.6733652390587552E-14, -5.748023455673529,
    2.5854213054236639E-15, -5.7020085348139991E-14, -5.7451917726803927,
    2.5956950236093759E-15, -5.730673931877696E-14, -5.742360089713193,
    2.6059687417950879E-15, -5.7593393289413918E-14, -5.739528406771929,
    2.6162424599807996E-15, -5.7880047260050875E-14, -5.7366989069050947,
    2.6265082570268161E-15, -5.8166480217603314E-14, -5.733867223898991,
    2.6367819752125281E-15, -5.8453134188240283E-14, -5.7310355409058547,
    2.64705569339824E-15, -5.8739788158877253E-14, -5.7282038579386549,
    2.6573294115839517E-15, -5.902644212951421E-14, -5.7253721749973909,
    2.6676031297696641E-15, -5.9313096100151193E-14, -5.7225426751331492,
    2.67786892681568E-15, -5.9599529057703619E-14, -5.7197109921374194,
    2.6881426450013918E-15, -5.9886183028340589E-14, -5.7168793091546579,
    2.6984163631871042E-15, -6.0172836998977559E-14, -5.71404806279975,
    2.708688497144877E-15, -6.0459446766997607E-14, -5.7112168164448418,
    2.71896063110265E-15, -6.0746056535017681E-14, -5.7083855700899342,
    2.7292327650604227E-15, -6.1032666303037729E-14, -205.70555432373502,
    2.7395048990181951E-15, -6.131927607105779E-14, -5.5713410357885707,
    2.7497770329759679E-15, -6.1605885839077851E-14, -5.4371277814759189,
    2.7600491669337408E-15, -6.1892495607097912E-14, -5.3029145944308764,
    2.7703213008915136E-15, -6.2179105375117973E-14, 3.8443106111439684E-15,
    -5.7508529555403651, 2.6582263889967E-13, 3.775300855940384E-15,
    -5.748023455673529, 2.6356254323218205E-13, 3.7062378524688953E-15,
    -5.7451917726803927, 2.6130070366371597E-13, 3.6371748489974081E-15,
    -5.742360089713193, 2.5903886409525E-13, 3.56811184552592E-15,
    -5.739528406771929, 2.56777024526784E-13, 3.4991020903223361E-15,
    -5.7366989069050947, 2.54516928859296E-13, 3.4300390868508481E-15,
    -5.733867223898991, 2.5225508929083004E-13, 3.3609760833793598E-15,
    -5.7310355409058547, 2.49993249722364E-13, 3.2919130799078718E-15,
    -5.7282038579386549, 2.47731410153898E-13, 3.2228500764363842E-15,
    -5.7253721749973909, 2.45469570585432E-13, 3.1538403212328002E-15,
    -5.7225426751331492, 2.43209474917944E-13, 3.0847773177613115E-15,
    -5.7197109921374194, 2.4094763534947803E-13, 3.0157143142898235E-15,
    -5.7168793091546579, 2.38685795781012E-13, 2.9466619604719172E-15,
    -5.71404806279975, 2.3642430499274159E-13, 2.8776096066540097E-15,
    -5.7112168164448418, 2.3416281420447123E-13, 2.8085572528361026E-15,
    -5.7083855700899342, 2.3190132341620081E-13, 2.7395048990181955E-15,
    -205.70555432373502, 2.296398326279304E-13, 2.6704525452002884E-15,
    -5.5713410357885707, 2.2737834183966E-13, 2.6014001913823809E-15,
    -5.4371277814759189, 2.251168510513896E-13, 2.5323478375644738E-15,
    -5.3029145944308764, 2.2285536026311921E-13, -9.0280048211833537E-14,
    2.4037536540050597E-13, -5.7508529555403651, -8.847106426650599E-14,
    2.39704789134938E-13, -5.748023455673529, -8.6660684500232941E-14,
    2.3903369544941203E-13, -5.7451917726803927, -8.4850304733959917E-14,
    2.38362601763886E-13, -5.742360089713193, -8.30399249676869E-14,
    2.3769150807836E-13, -5.739528406771929, -8.1230941022359321E-14,
    2.37020931812792E-13, -5.7366989069050947, -7.94205612560863E-14,
    2.36349838127266E-13, -5.733867223898991, -7.7610181489813249E-14,
    2.3567874444174003E-13, -5.7310355409058547, -7.5799801723540213E-14,
    2.35007650756214E-13, -5.7282038579386549, -7.3989421957267176E-14,
    2.34336557070688E-13, -5.7253721749973909, -7.2180438011939616E-14,
    2.3366598080512E-13, -5.7225426751331492, -7.0370058245666593E-14,
    2.32994887119594E-13, -5.7197109921374194, -6.8559678479393557E-14,
    2.32323793434068E-13, -5.7168793091546579, -6.67495778773096E-14,
    2.3165280323253361E-13, -5.71404806279975, -6.4939477275225667E-14,
    2.3098181303099919E-13, -5.7112168164448418, -6.3129376673141741E-14,
    2.3031082282946482E-13, -5.7083855700899342, -6.13192760710578E-14,
    2.2963983262793045E-13, -205.70555432373502, -5.9509175468973852E-14,
    2.2896884242639603E-13, -5.5713410357885707, -5.7699074866889913E-14,
    2.2829785222486161E-13, -5.4371277814759189, -5.5888974264805962E-14,
    2.2762686202332719E-13, -5.3029145944308764, -5.6155536953513252,
    2.51185894293112E-15, -5.4928583905951108E-14, -5.6127920289022279,
    2.52176521473016E-15, -5.5214702541094277E-14, -5.6100282315895722,
    2.53167913025744E-15, -5.5501041946789254E-14, -5.6072644342898847,
    2.54159304578472E-15, -5.5787381352484225E-14, -5.6045006370161339,
    2.5515069613119996E-15, -5.6073720758179196E-14, -5.6017389705670375,
    2.56141323311104E-15, -5.6359839393322371E-14, -5.5989751732543818,
    2.57132714863832E-15, -5.6646178799017342E-14, -5.596211375941726,
    2.5812410641656E-15, -5.6932518204712325E-14, -5.5934475786420386,
    2.5911549796928797E-15, -5.72188576104073E-14, -5.5906837813682868,
    2.60106889522016E-15, -5.7505197016102279E-14, -5.58792211491919,
    2.6109751670192E-15, -5.7791315651245441E-14, -5.5851583176091282,
    2.62088908254648E-15, -5.8077655056940412E-14, -5.5823945203068472,
    2.63080299807376E-15, -5.83639944626354E-14, -5.5796311491772776,
    2.640715384855392E-15, -5.8650289714220009E-14, -5.5768677780477089,
    2.6506277716370243E-15, -5.8936584965804624E-14, -5.5741044069181394,
    2.6605401584186561E-15, -5.9222880217389238E-14, -5.5713410357885707,
    2.670452545200288E-15, -5.9509175468973852E-14, -205.568577664659,
    2.68036493198192E-15, -5.9795470720558466E-14, -5.4344322519378832,
    2.6902773187635521E-15, -6.008176597214308E-14, -5.3002868728505694,
    2.700189705545184E-15, -6.03680612237277E-14, 3.8603385334731205E-15,
    -5.6155536953513252, 2.6509285850744997E-13, 3.7909692530225595E-15,
    -5.6127920289022279, 2.6283643506633E-13, 3.7215464468926396E-15,
    -5.6100282315895722, 2.6057827055774E-13, 3.6521236407627204E-15,
    -5.6072644342898847, 2.5832010604915E-13, 3.5827008346328E-15,
    -5.6045006370161339, 2.5606194154056003E-13, 3.51333155418224E-15,
    -5.6017389705670375, 2.5380551809944E-13, 3.44390874805232E-15,
    -5.5989751732543818, 2.5154735359085E-13, 3.3744859419223997E-15,
    -5.596211375941726, 2.4928918908225997E-13, 3.30506313579248E-15,
    -5.5934475786420386, 2.4703102457367003E-13, 3.2356403296625602E-15,
    -5.5906837813682868, 2.4477286006508003E-13, 3.166271049212E-15,
    -5.58792211491919, 2.4251643662396E-13, 3.0968482430820797E-15,
    -5.5851583176091282, 2.4025827211537E-13, 3.0274254369521598E-15,
    -5.5823945203068472, 2.3800010760678E-13, 2.9580133359581121E-15,
    -5.5796311491772776, 2.35742291311684E-13, 2.888601234964064E-15,
    -5.5768677780477089, 2.3348447501658803E-13, 2.819189133970016E-15,
    -5.5741044069181394, 2.3122665872149203E-13, 2.7497770329759683E-15,
    -5.5713410357885707, 2.2896884242639603E-13, 2.6803649319819203E-15,
    -205.568577664659, 2.2671102613129998E-13, 2.6109528309878722E-15,
    -5.4344322519378832, 2.24453209836204E-13, 2.5415407299938242E-15,
    -5.3002868728505694, 2.22195393541108E-13, -9.05716900973011E-14,
    2.3805508442154995E-13, -5.6155536953513252, -8.8762391829564282E-14,
    2.3738818038235E-13, -5.6127920289022279, -8.6951697498349247E-14,
    2.3672076175670003E-13, -5.6100282315895722, -8.5141003167134224E-14,
    2.3605334313105E-13, -5.6072644342898847, -8.3330308835919214E-14,
    2.353859245054E-13, -5.6045006370161339, -8.1521010568182378E-14,
    2.347190204662E-13, -5.6017389705670375, -7.9710316236967355E-14,
    2.3405160184055003E-13, -5.5989751732543818, -7.789962190575232E-14,
    2.3338418321490005E-13, -5.596211375941726, -7.60889275745373E-14,
    2.3271676458925E-13, -5.5934475786420386, -7.4278233243322275E-14,
    2.320493459636E-13, -5.5906837813682868, -7.2468934975585438E-14,
    2.3138244192440004E-13, -5.58792211491919, -7.0658240644370416E-14,
    2.3071502329875E-13, -5.5851583176091282, -6.88475463131554E-14,
    2.300476046731E-13, -5.5823945203068472, -6.7037131194636007E-14,
    2.2938028896474E-13, -5.5796311491772776, -6.5226716076116622E-14,
    2.2871297325638E-13, -5.5768677780477089, -6.3416300957597237E-14,
    2.2804565754802005E-13, -5.5741044069181394, -6.1605885839077864E-14,
    2.2737834183966004E-13, -5.5713410357885707, -5.9795470720558466E-14,
    2.267110261313E-13, -205.568577664659, -5.798505560203908E-14,
    2.2604371042294E-13, -5.4344322519378832, -5.6174640483519689E-14,
    2.2537639471458002E-13, -5.3002868728505694, -5.4802550069343736,
    2.448562377484592E-15, -5.3123515421314671E-14, -5.4775611402640223,
    2.4581091240366559E-15, -5.3409319734048564E-14, -5.4748651950058234,
    2.4676632369055039E-15, -5.3695344574801555E-14, -5.4721692497476244,
    2.477217349774352E-15, -5.3981369415554532E-14, -5.4694733045023938,
    2.4867714626431997E-15, -5.4267394256307517E-14, -5.4667794378320433,
    2.4963182091952639E-15, -5.4553198569041422E-14, -5.4640834925738444,
    2.505872322064112E-15, -5.4839223409794406E-14, -5.4613875473156455,
    2.51542643493296E-15, -5.5125248250547397E-14, -5.4586916020574465,
    2.5249805478018077E-15, -5.5411273091300381E-14, -5.4559956568122159,
    2.5345346606706562E-15, -5.5697297932053371E-14, -5.4533017901418646,
    2.5440814072227204E-15, -5.598310224478727E-14, -5.4506058448836656,
    2.5536355200915681E-15, -5.6269127085540254E-14, -5.447909899628061,
    2.5631896329604161E-15, -5.6555151926293245E-14, -5.4452143700900253,
    2.5727422725659074E-15, -5.6841132661442406E-14, -5.4425188405519895,
    2.5822949121713987E-15, -5.7127113396591579E-14, -5.4398233110139538,
    2.5918475517768896E-15, -5.7413094131740746E-14, -5.4371277814759189,
    2.6014001913823805E-15, -5.7699074866889913E-14, -5.4344322519378832,
    2.6109528309878718E-15, -5.798505560203908E-14, -205.43173672239985,
    2.6205054705933631E-15, -5.8271036337188247E-14, -5.2976591512702633,
    2.6300581101988544E-15, -5.8557017072337415E-14, 3.8763664558022725E-15,
    -5.4802550069343736, 2.6436307811523E-13, 3.8066376501047358E-15,
    -5.4775611402640223, 2.62110326900478E-13, 3.7368550413163839E-15,
    -5.4748651950058234, 2.59855837451764E-13, 3.667072432528032E-15,
    -5.4721692497476244, 2.5760134800305E-13, 3.59728982373968E-15,
    -5.4694733045023938, 2.5534685855433603E-13, 3.5275610180421438E-15,
    -5.4667794378320433, 2.53094107339584E-13, 3.4577784092537923E-15,
    -5.4640834925738444, 2.5083961789087004E-13, 3.38799580046544E-15,
    -5.4613875473156455, 2.48585128442156E-13, 3.3182131916770881E-15,
    -5.4586916020574465, 2.46330638993442E-13, 3.2484305828887362E-15,
    -5.4559956568122159, 2.44076149544728E-13, 3.1787017771912002E-15,
    -5.4533017901418646, 2.4182339832997603E-13, 3.1089191684028476E-15,
    -5.4506058448836656, 2.39568908881262E-13, 3.0391365596144957E-15,
    -5.447909899628061, 2.37314419432548E-13, 2.9693647114443074E-15,
    -5.4452143700900253, 2.3506027763062642E-13, 2.8995928632741184E-15,
    -5.4425188405519895, 2.3280613582870483E-13, 2.8298210151039298E-15,
    -5.4398233110139538, 2.305519940267832E-13, 2.7600491669337412E-15,
    -5.4371277814759189, 2.2829785222486161E-13, 2.6902773187635521E-15,
    -5.4344322519378832, 2.2604371042294E-13, 2.6205054705933631E-15,
    -205.43173672239985, 2.2378956862101839E-13, 2.5507336224231745E-15,
    -5.2976591512702633, 2.2153542681909681E-13, -9.0863331982768662E-14,
    2.35734803442594E-13, -5.4802550069343736, -8.9053719392622575E-14,
    2.35071571629762E-13, -5.4775611402640223, -8.7242710496465553E-14,
    2.34407828063988E-13, -5.4748651950058234, -8.5431701600308544E-14,
    2.3374408449821397E-13, -5.4721692497476244, -8.3620692704151535E-14,
    2.3308034093244E-13, -5.4694733045023938, -8.1811080114005422E-14,
    2.32417109119608E-13, -5.4667794378320433, -8.0000071217848413E-14,
    2.31753365553834E-13, -5.4640834925738444, -7.8189062321691391E-14,
    2.3108962198806E-13, -5.4613875473156455, -7.6378053425534369E-14,
    2.30425878422286E-13, -5.4586916020574465, -7.456704452937736E-14,
    2.29762134856512E-13, -5.4559956568122159, -7.275743193923126E-14,
    2.2909890304368005E-13, -5.4533017901418646, -7.0946423043074251E-14,
    2.28435159477906E-13, -5.4506058448836656, -6.9135414146917242E-14,
    2.27771415912132E-13, -5.447909899628061, -6.7324684511962409E-14,
    2.2710777469694641E-13, -5.4452143700900253, -6.5513954877007577E-14,
    2.2644413348176079E-13, -5.4425188405519895, -6.3703225242052745E-14,
    2.2578049226657522E-13, -5.4398233110139538, -6.1892495607097925E-14,
    2.2511685105138963E-13, -5.4371277814759189, -6.008176597214308E-14,
    2.24453209836204E-13, -5.4344322519378832, -5.8271036337188247E-14,
    2.2378956862101842E-13, -205.43173672239985, -5.6460306702233409E-14,
    2.2312592740583282E-13, -5.2976591512702633, -5.3449569239285015,
    2.3852658120380641E-15, -5.1318446936678228E-14, -5.342330823397905,
    2.3944530333431518E-15, -5.1603936927002857E-14, -5.3397026965551708,
    2.4036473435535679E-15, -5.1889647202813855E-14, -5.3370745697124358,
    2.412841653763984E-15, -5.2175357478624846E-14, -5.3344464428697016,
    2.4220359639743997E-15, -5.2461067754435837E-14, -5.331820342339106,
    2.4312231852794878E-15, -5.2746557744760473E-14, -5.3291922154963718,
    2.4404174954899039E-15, -5.3032268020571471E-14, -5.3265640886536376,
    2.44961180570032E-15, -5.3317978296382468E-14, -5.3239359618109026,
    2.4588061159107361E-15, -5.3603688572193459E-14, -5.3213078349681684,
    2.4680004261211522E-15, -5.3889398848004463E-14, -5.3186817344375719,
    2.4771876474262403E-15, -5.4174888838329093E-14, -5.3160536075948377,
    2.486381957636656E-15, -5.4460599114140084E-14, -5.3134254807521026,
    2.4955762678470721E-15, -5.4746309389951082E-14, -5.3107977591717965,
    2.5047691602764224E-15, -5.50319756086648E-14, -5.3081700375914895,
    2.5139620527057731E-15, -5.5317641827378528E-14, -5.3055423160111825,
    2.5231549451351235E-15, -5.5603308046092248E-14, -5.3029145944308764,
    2.5323478375644734E-15, -5.5888974264805975E-14, -5.3002868728505694,
    2.5415407299938238E-15, -5.6174640483519695E-14, -5.2976591512702633,
    2.5507336224231745E-15, -5.6460306702233421E-14, -205.29503142968997,
    2.5599265148525248E-15, -5.6745972920947142E-14, 3.8923943781314245E-15,
    -5.3449569239285015, 2.6363329772300997E-13, 3.8223060471869121E-15,
    -5.342330823397905, 2.6138421873462603E-13, 3.7521636357401274E-15,
    -5.3397026965551708, 2.59133404345788E-13, 3.6820212242933443E-15,
    -5.3370745697124358, 2.5688258995695E-13, 3.61187881284656E-15,
    -5.3344464428697016, 2.5463177556811204E-13, 3.541790481902048E-15,
    -5.331820342339106, 2.52382696579728E-13, 3.4716480704552641E-15,
    -5.3291922154963718, 2.5013188219089E-13, 3.40150565900848E-15,
    -5.3265640886536376, 2.47881067802052E-13, 3.331363247561696E-15,
    -5.3239359618109026, 2.45630253413214E-13, 3.2612208361149121E-15,
    -5.3213078349681684, 2.4337943902437603E-13, 3.1911325051704E-15,
    -5.3186817344375719, 2.4113036003599204E-13, 3.1209900937236154E-15,
    -5.3160536075948377, 2.38879545647154E-13, 3.0508476822768315E-15,
    -5.3134254807521026, 2.36628731258316E-13, 2.9807160869305027E-15,
    -5.3107977591717965, 2.3437826394956881E-13, 2.9105844915841728E-15,
    -5.3081700375914895, 2.3212779664082164E-13, 2.8404528962378432E-15,
    -5.3055423160111825, 2.2987732933207441E-13, 2.770321300891514E-15,
    -5.3029145944308764, 2.2762686202332719E-13, 2.7001897055451844E-15,
    -5.3002868728505694, 2.2537639471458E-13, 2.6300581101988544E-15,
    -5.2976591512702633, 2.231259274058328E-13, 2.5599265148525248E-15,
    -205.29503142968997, 2.208754600970856E-13, -9.1154973868236218E-14,
    2.3341452246363796E-13, -5.3449569239285015, -8.9345046955680868E-14,
    2.32754962877174E-13, -5.342330823397905, -8.7533723494581847E-14,
    2.32094894371276E-13, -5.3397026965551708, -8.5722400033482851E-14,
    2.31434825865378E-13, -5.3370745697124358, -8.3911076572383855E-14,
    2.3077475735948003E-13, -5.3344464428697016, -8.2101149659828479E-14,
    2.30115197773016E-13, -5.331820342339106, -8.0289826198729484E-14,
    2.29455129267118E-13, -5.3291922154963718, -7.8478502737630462E-14,
    2.2879506076122003E-13, -5.3265640886536376, -7.6667179276531454E-14,
    2.28134992255322E-13, -5.3239359618109026, -7.4855855815432458E-14,
    2.27474923749424E-13, -5.3213078349681684, -7.3045928902877083E-14,
    2.2681536416296002E-13, -5.3186817344375719, -7.1234605441778087E-14,
    2.26155295657062E-13, -5.3160536075948377, -6.9423281980679091E-14,
    2.25495227151164E-13, -5.3134254807521026, -6.76122378292888E-14,
    2.248352604291528E-13, -5.3107977591717965, -6.5801193677898532E-14,
    2.2417529370714157E-13, -5.3081700375914895, -6.3990149526508253E-14,
    2.2351532698513042E-13, -5.3055423160111825, -6.2179105375117986E-14,
    2.2285536026311924E-13, -5.3029145944308764, -6.03680612237277E-14,
    2.22195393541108E-13, -5.3002868728505694, -5.8557017072337415E-14,
    2.2153542681909681E-13, -5.2976591512702633, -5.6745972920947135E-14,
    2.2087546009708563E-13, -205.29503142968997 };

  
   double dv[720] = { 1.0, 0.0, 0.0, 0.016, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.016, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.016, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0,
    0.0, 0.032, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.032, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.032, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.048, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.048, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.048, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0,
    0.064, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.064, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.064, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.08, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.08,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.08, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.096, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.096, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.096, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 1.0, 0.0, 0.0, 0.112, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.112, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.112, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.128, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.128, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.128, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 1.0, 0.0, 0.0, 0.144, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.144, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.144, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.16, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.16, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.16, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0,
    0.0, 0.0, 0.176, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.176, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.176, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.192, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.192, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.192, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0,
    0.0, 0.208, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.208, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.208, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.224, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.224, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.224, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0,
    0.24, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.24, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.24,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.256, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.256,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.256, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.272, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.272, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.272, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 1.0, 0.0, 0.0, 0.288, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.288, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.288, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.304, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.304, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.304, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 1.0, 0.0, 0.0, 0.32, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.32, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.32, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 };

  
   double B[9] = { 0.6, 0.0, 0.0, 0.0, 0.6, 0.0, 0.0, 0.0,
      0.6 };
  /* double B[9] = { 0.68, 0.0, 0.0, 0.0, 0.332, 0.0, 0.0, 0.0, 0.394
  }; */

   signed char A[400] = { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 };

  double IGA[7200];
  double A_data[3600];
  double a[3600];
  double b_del_lam[240];
  double result[240];
  double varargin_1_data[240];
  double Bineq[120];
  double GAMMA[120];
  double del_lam[120];
  double esig[120];
  double ftol[120];
  double igr2[120];
  double ilam[120];
  double lam[120];
  double mesil[120];
  double H[60];
  double del_z[60];
  double u_max[60];
  double absxk;
  double mu;
  double scale;
  double t;
  int a_tmp;
  int b_i;
  int b_i1;
  int exitflag;
  int i;
  int i1;
  int idx;
  int iter;
  int iy;
  int j2;
  short b_Aineq[7200];
  signed char Aineq[7200];
  signed char At[7200];
  unsigned char ii_data[240];
  boolean_T x[240];

  //  Set Constants - Extract variables from Structure MPCParams
  //  MPCParams=load('MPCParams.mat');
  // MPCParams.X;
  // MPCParams.A;
  // MPCParams.B;
  // MPCParams.P1;
  // MPCParams.PSI;
  // MPCParams.OMEGA;
  // MPCParams.PHI;
  // MPCParams.L1;
  // MPCParams.L2;
  // MPCParams.L3;
  // MPCParams.L4;
  // MPCParams.Fmaxx;
  // MPCParams.Fmaxy;
  // MPCParams.Fmaxz;
  //  MPCParams.ConvThresh;
  // MPCParams.rCollAvoid;
  // MPCParams.maxiter;
  // MPCParams.tol;
  // MPCParams.ObstAvoid;
  // MPCParams.aObst; % Obstacle semi-major axis
  // MPCParams.bObst; % Obstacle semi-minor axis
  // MPCParams.alphaObst; % Orientation of ellipse
  // P_debris = MPCParams.P_debris;
  // MPCParams.pos_of_debris;
  //  MPCParams.AppCone;
  //  MPCParams.phi;
  //  sampT = MPCParams.Ts;
  //  Set Decision Variables - Configure constraints
  //  Compute distance from target
  scale = 3.3121686421112381E-170;
  absxk = std::abs(x0[0]);
  if (absxk > 3.3121686421112381E-170) {
    dr = 1.0;
    scale = absxk;
  } else {
    t = absxk / 3.3121686421112381E-170;
    dr = t * t;
  }

  absxk = std::abs(x0[1]);
  if (absxk > scale) {
    t = scale / absxk;
    dr = dr * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    dr += t * t;
  }

  dr = scale * std::sqrt(dr);

  //  position tolerance to stop control action
  //  Define vectors for Obstacle Avoidance
  //  to current position
  //  to obstacle
  //  to final position
  //  from obstacle to target
  //  from obstacle to FSS
  //  Angle between unit vectors
  //  Approach Cone variables
  //  point where cone is placed
  //  orientation of cone (inertial)
  //  Unit vector of approach axis
  //  Angle between FSS and approach axis
  dock_flag = 0.0;
  CollAvoid_flag = 0.0;

  //  Set obstacle avoidance flag
  //  Define target point for current iteration
  //      xfinal = xfinal;
  //  to final position
  for (i = 0; i < 6; i++) {
    target_state[i] = 0.0;
  }

  //  Define and Solve QP
  //  dock_complete = 0;
  pt_sel = 0.0;
  dock_complete = 0.0;

  //  Compute QP Matrices
  // Dimension of GAMMA is n*horizon x n
  // for i = 1:horizon
  for (b_i = 0; b_i < 120; b_i++) {
    t = 0.0;
    for (i1 = 0; i1 < 6; i1++) {
      t += 2.0 * x0[i1] * dv[i1 + 6 * b_i];
    }

    GAMMA[b_i] = t;
  }

  for (b_i = 0; b_i < 120; b_i++) {
    t = 0.0;
    for (i1 = 0; i1 < 120; i1++) {
      t += GAMMA[i1] * b[i1 + 120 * b_i];
    }

    Bineq[b_i] = t;
  }

  for (b_i = 0; b_i < 60; b_i++) {
    t = 0.0;
    for (i1 = 0; i1 < 120; i1++) {
      t += Bineq[i1] * b_b[i1 + 120 * b_i];
    }

    H[b_i] = t;
  }

  //  Control Consraints
  i = -1;
  for (idx = 0; idx < 20; idx++) {
    for (j2 = 0; j2 < 3; j2++) {
      t = B[3 * j2];
      scale = B[3 * j2 + 1];
      absxk = B[3 * j2 + 2];
      for (b_i1 = 0; b_i1 < 20; b_i1++) {
        i++;
        a_tmp = A[b_i1 + 20 * idx];
        a[i] = static_cast<double>(a_tmp) * t;
        i++;
        a[i] = static_cast<double>(a_tmp) * scale;
        i++;
        a[i] = static_cast<double>(a_tmp) * absxk;
      }
    }
  }

  for (b_i = 0; b_i < 60; b_i++) {
    t = 0.0;
    for (i1 = 0; i1 < 60; i1++) {
      t += a[b_i + 60 * i1];
    }

    u_max[b_i] = t;
  }

  std::memset(&Aineq[0], 0, 7200U * sizeof(signed char));
  std::memset(&a[0], 0, 3600U * sizeof(double));
  for (i = 0; i < 60; i++) {
    a[i + 60 * i] = 1.0;
  }

  for (b_i = 0; b_i < 60; b_i++) {
    for (i1 = 0; i1 < 60; i1++) {
      idx = static_cast<int>(a[i1 + 60 * b_i]);
      i = i1 + 120 * b_i;
      Aineq[i] = static_cast<signed char>(-idx);
      Aineq[i + 60] = static_cast<signed char>(idx);
    }
  }

  std::memset(&Bineq[0], 0, 120U * sizeof(double));

  //  Control Constraints only
  //      Aineq = Aineq;
  //      Bineq = Bineq; % Joh and me comment this on spet 23
  //  Control Constraints only
  //      Aineq = Aineq;
  //      Bineq = Bineq;
  //  Call QP Solver
  for (i = 0; i < 60; i++) {
    t = u_max[i];
    Bineq[i] = t;
    Bineq[i + 60] = t;
  }
   /* if (X_QP[0].empty()){
   for (i = 0; i < 60; i++){
      X_QP[i] = 0.0;

    }
  } */
  //  Solve quadratic programming problem using Wright's (1997) Method
  //  Minimise J(x) = 1/2x'Hx + f'x
  //  Subject to: Ax <= b
  //  Supporting Functions
  //  Reference: S. J. Wright, "Applying New Optimization Algorithms to Model
  //  Predictive Control," in Chemical Process Control-V, CACHE, AIChE
  //  Symposium, 1997, pp. 147-155.
  // Number of decision variables
  //  p = 0;
  // Test for Cold Start
  // Warm Start
  // to tune
  // to tune
  // Default Values
  mu = 10000.0;
  for (i = 0; i < 120; i++) {
    lam[i] = 100.0;
    ftol[i] = 100.0;
    esig[i] = 0.001;
    for (b_i = 0; b_i < 60; b_i++) {
      At[b_i + 60 * i] = Aineq[i + 120 * b_i];
    }
  }

  //  %Linsolve options
  //  opU.UT = true;
  //  opUT.UT = true;
  //  opUT.TRANSA = true;
  // Begin Searching
  //  for iter = 1:maxiter
  iter = 0;
  exitflag = 0;
  while ((iter <= 100) && (exitflag != 1)) {
    boolean_T exitg1;
    boolean_T y;

    // Create common matrices
    for (i = 0; i < 120; i++) {
      t = lam[i];
      scale = 1.0 / t;
      ilam[i] = scale;
      GAMMA[i] = -t / ftol[i];
      mesil[i] = mu * esig[i] * scale;
    }

    // RHS
    for (i = 0; i < 60; i++) {
      for (b_i = 0; b_i < 120; b_i++) {
        idx = b_i + 120 * i;
        IGA[idx] = GAMMA[b_i] * static_cast<double>(Aineq[idx]);
      }

      t = 0.0;
      for (b_i = 0; b_i < 60; b_i++) {
        t += b_a[i + 60 * b_i] * X_QP[b_i];
      }

      scale = 0.0;
      for (b_i = 0; b_i < 120; b_i++) {
        scale += static_cast<double>(At[i + 60 * b_i]) * lam[b_i];
      }

      u_max[i] = (t - scale) - H[i];
    }

    for (b_i = 0; b_i < 7200; b_i++) {
      b_Aineq[b_i] = static_cast<short>(-Aineq[b_i]);
    }

    for (b_i = 0; b_i < 120; b_i++) {
      t = 0.0;
      for (i1 = 0; i1 < 60; i1++) {
        t += static_cast<double>(b_Aineq[b_i + 120 * i1]) * X_QP[i1];
      }

      igr2[b_i] = GAMMA[b_i] * ((t + Bineq[b_i]) - mesil[b_i]);
    }

    // Solve
    for (b_i = 0; b_i < 60; b_i++) {
      for (i1 = 0; i1 < 60; i1++) {
        t = 0.0;
        for (idx = 0; idx < 120; idx++) {
          t += static_cast<double>(At[b_i + 60 * idx]) * IGA[idx + 120 * i1];
        }

        a[b_i + 60 * i1] = t;
      }
    }

    for (b_i = 0; b_i < 3600; b_i++) {
      A_data[b_i] = b_H[b_i] - a[b_i];
    }

    a_tmp = -1;
    idx = 0;
    exitg1 = false;
    while ((!exitg1) && (idx < 60)) {
      int idxA1j;
      int idxAjj;
      int ix;
      idxA1j = idx * 60;
      idxAjj = idxA1j + idx;
      scale = 0.0;
      if (idx >= 1) {
        ix = idxA1j;
        iy = idxA1j;
        for (i = 0; i < idx; i++) {
          scale += A_data[ix] * A_data[iy];
          ix++;
          iy++;
        }
      }

      scale = A_data[idxAjj] - scale;
      if (scale > 0.0) {
        scale = std::sqrt(scale);
        A_data[idxAjj] = scale;
        if (idx + 1 < 60) {
          int idxAjjp1;
          i = idxA1j + 61;
          idxAjjp1 = idxAjj + 61;
          if (idx != 0) {
            iy = idxAjj + 60;
            b_i = (idxA1j + 60 * (58 - idx)) + 61;
            for (j2 = i; j2 <= b_i; j2 += 60) {
              ix = idxA1j;
              absxk = 0.0;
              i1 = (j2 + idx) - 1;
              for (b_i1 = j2; b_i1 <= i1; b_i1++) {
                absxk += A_data[b_i1 - 1] * A_data[ix];
                ix++;
              }

              A_data[iy] += -absxk;
              iy += 60;
            }
          }

          scale = 1.0 / scale;
          b_i = (idxAjj + 60 * (58 - idx)) + 61;
          for (i = idxAjjp1; i <= b_i; i += 60) {
            A_data[i - 1] *= scale;
          }
        }

        idx++;
      } else {
        a_tmp = idx;
        exitg1 = true;
      }
    }

    // [R] = chol(H-At*IGA);
    if (a_tmp + 1 == 0) {
      for (b_i = 0; b_i < 60; b_i++) {
        t = 0.0;
        for (i1 = 0; i1 < 120; i1++) {
          t += static_cast<double>(At[b_i + 60 * i1]) * igr2[i1];
        }

        del_z[b_i] = u_max[b_i] - t;
      }

      for (b_i = 0; b_i < 3600; b_i++) {
        a[b_i] = b_H[b_i] - a[b_i];
      }

      mldivide(a, del_z);

      // old method (LU?)
      //   del_z = linsolve (R, linsolve (R, (r1-At*igr2), opUT), opU); %exploit matrix properties for solving 
    } else {
      // Not Positive Definite (problem? eg infeasible)
      for (b_i = 0; b_i < 60; b_i++) {
        t = 0.0;
        for (i1 = 0; i1 < 120; i1++) {
          t += static_cast<double>(At[b_i + 60 * i1]) * igr2[i1];
        }

        del_z[b_i] = u_max[b_i] - t;
      }

      for (b_i = 0; b_i < 3600; b_i++) {
        a[b_i] = b_H[b_i] - a[b_i];
      }

      mldivide(a, del_z);

      // old method (LU?)
    }

    // Decide on suitable alpha (from Wright's paper)
    // Try Max Increment (alpha = 1)
    // Check lam and ftol > 0
    for (i = 0; i < 120; i++) {
      t = 0.0;
      for (b_i = 0; b_i < 60; b_i++) {
        t += IGA[i + 120 * b_i] * del_z[b_i];
      }

      t = igr2[i] - t;
      del_lam[i] = t;
      scale = ftol[i];
      absxk = (-scale + mesil[i]) - ilam[i] * scale * t;
      mesil[i] = absxk;
      t += lam[i];
      GAMMA[i] = t;
      scale += absxk;
      ilam[i] = scale;
      x[i] = (t < 2.2204460492503131E-16);
      x[i + 120] = (scale < 2.2204460492503131E-16);
    }

    y = false;
    i = 0;
    exitg1 = false;
    while ((!exitg1) && (i < 240)) {
      if (!x[i]) {
        i++;
      } else {
        y = true;
        exitg1 = true;
      }
    }

    if (!y) {
      // KKT met
      std::memcpy(&lam[0], &GAMMA[0], 120U * sizeof(double));
      std::memcpy(&ftol[0], &ilam[0], 120U * sizeof(double));
      for (b_i = 0; b_i < 60; b_i++) {
        X_QP[b_i] += del_z[b_i];
      }
    } else {
      // KKT failed - solve by finding minimum ratio
      for (b_i = 0; b_i < 120; b_i++) {
        result[b_i] = GAMMA[b_i];
        result[b_i + 120] = ilam[b_i];
      }

      idx = 0;
      i = 0;
      exitg1 = false;
      while ((!exitg1) && (i < 240)) {
        if (result[i] < 2.2204460492503131E-16) {
          idx++;
          ii_data[idx - 1] = static_cast<unsigned char>(i + 1);
          if (idx >= 240) {
            exitg1 = true;
          } else {
            i++;
          }
        } else {
          i++;
        }
      }

      if (1 > idx) {
        j2 = 0;
      } else {
        j2 = idx;
      }

      // detects elements breaking KKT condition
      for (b_i = 0; b_i < 120; b_i++) {
        b_del_lam[b_i] = del_lam[b_i];
        b_del_lam[b_i + 120] = mesil[b_i];
      }

      for (b_i = 0; b_i < j2; b_i++) {
        i = ii_data[b_i] - 1;
        varargin_1_data[b_i] = 1.0 - result[i] / b_del_lam[i];
      }

      if (j2 <= 2) {
        if (j2 == 1) {
          scale = varargin_1_data[0];
        } else if ((varargin_1_data[0] > varargin_1_data[1]) || (
                    (varargin_1_data[0]) && (!rtIsNaN(varargin_1_data[1])))) {
          scale = varargin_1_data[1];
        } else {
          scale = varargin_1_data[0];
        }
      } else {
        if (!rtIsNaN(varargin_1_data[0])) {
          idx = 1;
        } else {
          idx = 0;
          i = 2;
          exitg1 = false;
          while ((!exitg1) && (i <= j2)) {
            if (!rtIsNaN(varargin_1_data[i - 1])) {
              idx = i;
              exitg1 = true;
            } else {
              i++;
            }
          }
        }

        if (idx == 0) {
          scale = varargin_1_data[0];
        } else {
          scale = varargin_1_data[idx - 1];
          b_i = idx + 1;
          for (i = b_i; i <= j2; i++) {
            t = varargin_1_data[i - 1];
            if (scale > t) {
              scale = t;
            }
          }
        }
      }

      scale *= 0.995;

      // solves for min ratio (max value of alpha allowed)
      // Increment
      for (b_i = 0; b_i < 120; b_i++) {
        lam[b_i] += scale * del_lam[b_i];
        ftol[b_i] += scale * mesil[b_i];
      }

      for (b_i = 0; b_i < 60; b_i++) {
        X_QP[b_i] += scale * del_z[b_i];
      }
    }

    // Complimentary Gap
    absxk = mu;
    scale = 0.0;
    for (b_i = 0; b_i < 120; b_i++) {
      scale += ftol[b_i] * lam[b_i];
    }

    mu = scale / 120.0;

    //      if(mu < tol)
    //          exitflag = 1;
    //          return
    //      end
    //      %Solve for new Sigma
    //      sigma = mu/mu_old;
    //      if(sigma > 0.1) %to tune
    //          sigma = 0.1;
    //      end
    //      esig = sigma*ones(mc,1);
    if (mu < 0.001) {
      exitflag = 1;
    } else {
      // Solve for new Sigma
      scale = mu / absxk;
      if (scale > 0.1) {
        // to tune
        scale = 0.1;
      }

      for (i = 0; i < 120; i++) {
        esig[i] = scale;
      }
    }

    iter++;
  }

  // Check for failure
  num_iter = iter;

  //  solution to warm-start next iteration
  //  Extract first control
  Fx = X_QP[0];
  Fy = X_QP[1];
  Fz = X_QP[2];
}

template<typename T>
void CoordinatorBase<T>::MPC_Guidance_v3_sand_good()
{
   double b[14400] = { 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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
    100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0,
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
    0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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
    0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0,
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
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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
    100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0,
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
    0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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
    0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0,
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
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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
    100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0,
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
    0.0, 0.0, 0.0, 0.0, 200625.9621, 5.0006E-8, -8.434E-11, 93457.1261,
    1.5777E-8, -8.046E-12, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0006E-8,
    200625.9621, -6.0596E-11, 2.3592E-8, 93457.1261, -2.0259E-11, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -8.434E-11, -6.0596E-11, 200625.9621,
    -3.7164E-11, 1.0179E-10, 93457.1261, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 93457.1261, 2.3592E-8, -3.7164E-11, 3.0492404765E+6, 1.1459E-8,
    -2.5253E-12, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.5777E-8, 93457.1261,
    1.0179E-10, 1.1459E-8, 3.0492404765E+6, -3.3164E-10, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, -8.046E-12, -2.0259E-11, 93457.1261, -2.5253E-12,
    -3.3164E-10, 3.0492404765E+6 };

  
   double b_b[7200] = { 8.5613E-6, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    2.5684E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 4.2807E-5, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 5.9929E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 7.7052E-5, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 9.4174E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0001113,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00012842, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    0.00014554, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00016266, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 0.00017979, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00019691, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 0.00021403, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00023116,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00024828, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    0.0002654, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00028252, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 0.00029965, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00031677, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 0.00033389, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0,
    8.5613E-6, 0.0, 0.0, 0.0010702, 0.0, 0.0, 2.5684E-5, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 4.2807E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 5.9929E-5, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 7.7052E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 9.4174E-5,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0001113, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    0.00012842, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00014554, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 0.00016266, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00017979, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 0.00019691, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00021403,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00023116, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    0.00024828, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0002654, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 0.00028252, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00029965, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 0.00031677, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00033389,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0, 8.5613E-6, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    2.5684E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 4.2807E-5, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 5.9929E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 7.7052E-5, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 9.4174E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0001113,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00012842, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    0.00014554, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00016266, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 0.00017979, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00019691, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 0.00021403, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00023116,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00024828, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    0.0002654, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00028252, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 0.00029965, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00031677, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 0.00033389, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 8.5613E-6, 0.0, 0.0, 0.0010702, 0.0, 0.0, 2.5684E-5, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 4.2807E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 5.9929E-5,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 7.7052E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    9.4174E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0001113, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 0.00012842, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00014554, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 0.00016266, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00017979,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00019691, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    0.00021403, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00023116, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 0.00024828, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0002654, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 0.00028252, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00029965,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00031677, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.5613E-6, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    2.5684E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 4.2807E-5, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 5.9929E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 7.7052E-5, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 9.4174E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0001113,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00012842, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    0.00014554, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00016266, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 0.00017979, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00019691, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 0.00021403, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00023116,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00024828, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    0.0002654, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00028252, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 0.00029965, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00031677, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.5613E-6, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 2.5684E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 4.2807E-5,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 5.9929E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    7.7052E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 9.4174E-5, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 0.0001113, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00012842, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 0.00014554, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00016266,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00017979, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    0.00019691, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00021403, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 0.00023116, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00024828, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 0.0002654, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00028252,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00029965, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    0.00031677, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 8.5613E-6, 0.0, 0.0, 0.0010702, 0.0, 0.0, 2.5684E-5, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 4.2807E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 5.9929E-5,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 7.7052E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    9.4174E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0001113, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 0.00012842, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00014554, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 0.00016266, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00017979,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00019691, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    0.00021403, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00023116, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 0.00024828, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0002654, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 0.00028252, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00029965,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 8.5613E-6, 0.0, 0.0, 0.0010702, 0.0, 0.0, 2.5684E-5, 0.0,
    0.0, 0.0010702, 0.0, 0.0, 4.2807E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    5.9929E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 7.7052E-5, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 9.4174E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0001113, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 0.00012842, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00014554,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00016266, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    0.00017979, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00019691, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 0.00021403, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00023116, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 0.00024828, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0002654,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00028252, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    0.00029965, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.5613E-6, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    2.5684E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 4.2807E-5, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 5.9929E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 7.7052E-5, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 9.4174E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0001113,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00012842, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    0.00014554, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00016266, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 0.00017979, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00019691, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 0.00021403, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00023116,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00024828, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    0.0002654, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00028252, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 0.00029965, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.5613E-6, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 2.5684E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 4.2807E-5,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 5.9929E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    7.7052E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 9.4174E-5, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 0.0001113, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00012842, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 0.00014554, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00016266,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00017979, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    0.00019691, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00021403, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 0.00023116, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00024828, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 0.0002654, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00028252,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.5613E-6, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 2.5684E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 4.2807E-5,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 5.9929E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    7.7052E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 9.4174E-5, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 0.0001113, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00012842, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 0.00014554, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00016266,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00017979, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    0.00019691, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00021403, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 0.00023116, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00024828, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 0.0002654, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00028252,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.5613E-6, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 2.5684E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 4.2807E-5,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 5.9929E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    7.7052E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 9.4174E-5, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 0.0001113, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00012842, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 0.00014554, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00016266,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00017979, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    0.00019691, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00021403, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 0.00023116, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00024828, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 0.0002654, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00028252,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.5613E-6,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 2.5684E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    4.2807E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 5.9929E-5, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 7.7052E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 9.4174E-5, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 0.0001113, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00012842,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00014554, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    0.00016266, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00017979, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 0.00019691, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00021403, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 0.00023116, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00024828,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0002654, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.5613E-6, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 2.5684E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 4.2807E-5, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 5.9929E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 7.7052E-5,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 9.4174E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    0.0001113, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00012842, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 0.00014554, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00016266, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 0.00017979, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00019691,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00021403, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    0.00023116, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00024828, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 0.0002654, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 8.5613E-6, 0.0, 0.0, 0.0010702, 0.0, 0.0, 2.5684E-5,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 4.2807E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    5.9929E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 7.7052E-5, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 9.4174E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0001113, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 0.00012842, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00014554,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00016266, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    0.00017979, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00019691, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 0.00021403, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00023116, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 0.00024828, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0002654,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 8.5613E-6, 0.0, 0.0, 0.0010702, 0.0, 0.0, 2.5684E-5, 0.0,
    0.0, 0.0010702, 0.0, 0.0, 4.2807E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    5.9929E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 7.7052E-5, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 9.4174E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0001113, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 0.00012842, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00014554,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00016266, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    0.00017979, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00019691, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 0.00021403, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00023116, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 0.00024828, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    8.5613E-6, 0.0, 0.0, 0.0010702, 0.0, 0.0, 2.5684E-5, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 4.2807E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 5.9929E-5, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 7.7052E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 9.4174E-5,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0001113, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    0.00012842, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00014554, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 0.00016266, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00017979, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 0.00019691, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00021403,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00023116, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    0.00024828, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.5613E-6, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 2.5684E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 4.2807E-5, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 5.9929E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 7.7052E-5,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 9.4174E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    0.0001113, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00012842, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 0.00014554, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00016266, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 0.00017979, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00019691,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00021403, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    0.00023116, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00024828, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.5613E-6, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    2.5684E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 4.2807E-5, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 5.9929E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 7.7052E-5, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 9.4174E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0001113,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00012842, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    0.00014554, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00016266, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 0.00017979, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00019691, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 0.00021403, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00023116,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.5613E-6,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 2.5684E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    4.2807E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 5.9929E-5, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 7.7052E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 9.4174E-5, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 0.0001113, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00012842,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00014554, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    0.00016266, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00017979, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 0.00019691, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00021403, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 0.00023116, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 8.5613E-6, 0.0, 0.0, 0.0010702, 0.0, 0.0, 2.5684E-5,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 4.2807E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    5.9929E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 7.7052E-5, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 9.4174E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0001113, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 0.00012842, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00014554,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00016266, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    0.00017979, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00019691, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 0.00021403, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00023116, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    8.5613E-6, 0.0, 0.0, 0.0010702, 0.0, 0.0, 2.5684E-5, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 4.2807E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 5.9929E-5, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 7.7052E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 9.4174E-5,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0001113, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    0.00012842, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00014554, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 0.00016266, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00017979, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 0.00019691, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00021403,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 8.5613E-6, 0.0, 0.0, 0.0010702, 0.0, 0.0, 2.5684E-5, 0.0,
    0.0, 0.0010702, 0.0, 0.0, 4.2807E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    5.9929E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 7.7052E-5, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 9.4174E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0001113, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 0.00012842, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00014554,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00016266, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    0.00017979, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00019691, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 0.00021403, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.5613E-6, 0.0, 0.0, 0.0010702, 0.0,
    0.0, 2.5684E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 4.2807E-5, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 5.9929E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 7.7052E-5,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 9.4174E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    0.0001113, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00012842, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 0.00014554, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00016266, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 0.00017979, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00019691,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00021403, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 8.5613E-6, 0.0, 0.0, 0.0010702, 0.0, 0.0, 2.5684E-5, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 4.2807E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 5.9929E-5,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 7.7052E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    9.4174E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0001113, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 0.00012842, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00014554, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 0.00016266, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00017979,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00019691, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 8.5613E-6, 0.0, 0.0, 0.0010702, 0.0, 0.0, 2.5684E-5, 0.0,
    0.0, 0.0010702, 0.0, 0.0, 4.2807E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    5.9929E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 7.7052E-5, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 9.4174E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0001113, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 0.00012842, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00014554,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00016266, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    0.00017979, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00019691, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.5613E-6, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    2.5684E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 4.2807E-5, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 5.9929E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 7.7052E-5, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 9.4174E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0001113,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00012842, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    0.00014554, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00016266, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 0.00017979, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00019691, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.5613E-6, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 2.5684E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 4.2807E-5,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 5.9929E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    7.7052E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 9.4174E-5, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 0.0001113, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00012842, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 0.00014554, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00016266,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00017979, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.5613E-6, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 2.5684E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 4.2807E-5,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 5.9929E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    7.7052E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 9.4174E-5, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 0.0001113, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00012842, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 0.00014554, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00016266,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00017979, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.5613E-6, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 2.5684E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 4.2807E-5,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 5.9929E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    7.7052E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 9.4174E-5, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 0.0001113, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00012842, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 0.00014554, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00016266,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00017979, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.5613E-6,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 2.5684E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    4.2807E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 5.9929E-5, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 7.7052E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 9.4174E-5, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 0.0001113, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00012842,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00014554, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    0.00016266, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.5613E-6, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 2.5684E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 4.2807E-5, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 5.9929E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 7.7052E-5,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 9.4174E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    0.0001113, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00012842, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 0.00014554, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00016266, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 8.5613E-6, 0.0, 0.0, 0.0010702, 0.0, 0.0, 2.5684E-5,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 4.2807E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    5.9929E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 7.7052E-5, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 9.4174E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0001113, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 0.00012842, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00014554,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00016266, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 8.5613E-6, 0.0, 0.0, 0.0010702, 0.0, 0.0, 2.5684E-5, 0.0,
    0.0, 0.0010702, 0.0, 0.0, 4.2807E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    5.9929E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 7.7052E-5, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 9.4174E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0001113, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 0.00012842, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00014554,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.5613E-6,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 2.5684E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    4.2807E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 5.9929E-5, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 7.7052E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 9.4174E-5, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 0.0001113, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00012842,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00014554, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.5613E-6, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    2.5684E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 4.2807E-5, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 5.9929E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 7.7052E-5, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 9.4174E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0001113,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00012842, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    0.00014554, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 8.5613E-6, 0.0, 0.0, 0.0010702, 0.0, 0.0, 2.5684E-5, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 4.2807E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 5.9929E-5,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 7.7052E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    9.4174E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0001113, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 0.00012842, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.5613E-6, 0.0, 0.0, 0.0010702, 0.0,
    0.0, 2.5684E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 4.2807E-5, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 5.9929E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 7.7052E-5,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 9.4174E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    0.0001113, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.00012842, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    8.5613E-6, 0.0, 0.0, 0.0010702, 0.0, 0.0, 2.5684E-5, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 4.2807E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 5.9929E-5, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 7.7052E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 9.4174E-5,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0001113, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    0.00012842, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.5613E-6, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 2.5684E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 4.2807E-5, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 5.9929E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 7.7052E-5,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 9.4174E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    0.0001113, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.5613E-6, 0.0,
    0.0, 0.0010702, 0.0, 0.0, 2.5684E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    4.2807E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 5.9929E-5, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 7.7052E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 9.4174E-5, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 0.0001113, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 8.5613E-6, 0.0, 0.0, 0.0010702, 0.0, 0.0, 2.5684E-5, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 4.2807E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 5.9929E-5,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 7.7052E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    9.4174E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0001113, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.5613E-6, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 2.5684E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 4.2807E-5, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 5.9929E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 7.7052E-5,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 9.4174E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.5613E-6, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 2.5684E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 4.2807E-5, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 5.9929E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 7.7052E-5,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 9.4174E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.5613E-6, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 2.5684E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 4.2807E-5, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 5.9929E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 7.7052E-5,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 9.4174E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.5613E-6, 0.0,
    0.0, 0.0010702, 0.0, 0.0, 2.5684E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    4.2807E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 5.9929E-5, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 7.7052E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.5613E-6, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 2.5684E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 4.2807E-5,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 5.9929E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    7.7052E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.5613E-6, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 2.5684E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 4.2807E-5, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 5.9929E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 7.7052E-5,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.5613E-6, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 2.5684E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 4.2807E-5,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 5.9929E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.5613E-6, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    2.5684E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 4.2807E-5, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 5.9929E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 8.5613E-6, 0.0, 0.0, 0.0010702, 0.0, 0.0, 2.5684E-5, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 4.2807E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 5.9929E-5,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 8.5613E-6, 0.0, 0.0, 0.0010702, 0.0, 0.0, 2.5684E-5, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 4.2807E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.5613E-6, 0.0, 0.0,
    0.0010702, 0.0, 0.0, 2.5684E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 4.2807E-5,
    0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 8.5613E-6, 0.0, 0.0, 0.0010702, 0.0, 0.0, 2.5684E-5, 0.0,
    0.0, 0.0010702, 0.0, 0.0, 4.2807E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    8.5613E-6, 0.0, 0.0, 0.0010702, 0.0, 0.0, 2.5684E-5, 0.0, 0.0, 0.0010702,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.5613E-6, 0.0, 0.0, 0.0010702, 0.0, 0.0,
    2.5684E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.5613E-6, 0.0,
    0.0, 0.0010702, 0.0, 0.0, 2.5684E-5, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 8.5613E-6, 0.0, 0.0, 0.0010702, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.5613E-6, 0.0, 0.0, 0.0010702, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.5613E-6, 0.0, 0.0,
    0.0010702 };

  
   double b_H[3600] = { 211.5154541352035, 6.5533595964669213E-14,
    -5.6899162546212E-17, 11.280659712816917, 6.40974107466916E-14,
    -5.4573129935716004E-17, 11.045865348887073, 6.2661225528714E-14,
    -5.224709732522E-17, 10.811067762614019, 6.12242014175016E-14,
    -4.9919706050816008E-17, 10.576273692184389, 5.9788016199524E-14,
    -4.7593673440320003E-17, 10.341479856129229, 5.83518309815464E-14,
    -4.5267640829824E-17, 10.106686313076017, 5.6915645763568807E-14,
    -4.2941608219328E-17, 9.8718897822287275, 5.5478621652356411E-14,
    -4.0614216944924E-17, 9.6370970017757, 5.404243643437881E-14,
    -3.8288184334428005E-17, 9.4023046902394611, 5.2606251216401209E-14,
    -3.5962151723932006E-17, 9.1675095668102937, 5.1169227105188813E-14,
    -3.3634760449527996E-17, 8.9327183696834176, 4.9733041887211206E-14,
    -3.1308727839032004E-17, 8.6979278173984422, 4.8296856669233612E-14,
    -2.8982695228536E-17, 8.4631379685742569, 4.6860671451256005E-14,
    -2.665666261804E-17, 8.2283468781897078, 4.5423982897337528E-14,
    -2.43298148091992E-17, 7.993557944730604, 4.3987629900712964E-14,
    -2.20035104659216E-17, 7.7587695568412611, 4.2551193014764926E-14,
    -1.9677070256253204E-17, 7.5239824411033318, 4.1114840018140368E-14,
    -1.73507659129756E-17, 7.289195988207168, 3.9678403132192323E-14,
    -1.5024325703307202E-17, 7.0544106909915394, 3.8241991413041333E-14,
    -1.2697926253556042E-17, 6.5533595964669213E-14, 211.5154541352035,
    -7.1491717928462731E-16, 6.4383780855811606E-14, 11.280659712816917,
    -7.1795437889196179E-16, 6.3233965746954E-14, 11.045865348887073,
    -7.2099157849929608E-16, 6.2083479017121611E-14, 10.811067762614019,
    -7.2403055217182159E-16, 6.0933663908264E-14, 10.576273692184389,
    -7.2706775177915607E-16, 5.97838487994064E-14, 10.341479856129229,
    -7.3010495138649046E-16, 5.86340336905488E-14, 10.106686313076017,
    -7.3314215099382494E-16, 5.74835469607164E-14, 9.8718897822287275,
    -7.3618112466635055E-16, 5.6333731851858809E-14, 9.6370970017757,
    -7.3921832427368493E-16, 5.5183916743001208E-14, 9.4023046902394611,
    -7.4225552388101942E-16, 5.40334300131688E-14, 9.1675095668102937,
    -7.4529449755354493E-16, 5.2883614904311194E-14, 8.9327183696834176,
    -7.4833169716087921E-16, 5.17337997954536E-14, 8.6979278173984422,
    -7.5136889676821369E-16, 5.0583984686596006E-14, 8.4631379685742569,
    -7.5440609637554808E-16, 4.9433766605153528E-14, 8.2283468781897078,
    -7.5744436042199714E-16, 4.8283817172100973E-14, 7.993557944730604,
    -7.6048191484236985E-16, 4.7133800576950924E-14, 7.7587695568412611,
    -7.6351964666926162E-16, 4.5983851143898369E-14, 7.5239824411033318,
    -7.6655720108963423E-16, 4.4833834548748326E-14, 7.289195988207168,
    -7.695949329165261E-16, 4.3683838102227525E-14, 7.0544106909915394,
    -7.72632611521462E-16, -5.6899162546212E-17, -7.1491717928462731E-16,
    211.5154541352035, -5.5640119678180009E-17, -7.1348205829720955E-16,
    11.280659712816917, -5.4381076810148E-17, -7.1204693730979209E-16,
    11.045865348887073, -5.3121298519880006E-17, -7.1061097805076971E-16,
    10.811067762614019, -5.1862255651848E-17, -7.0917585706335215E-16,
    10.576273692184389, -5.0603212783816006E-17, -7.0774073607593449E-16,
    10.341479856129229, -4.9344169915784E-17, -7.0630561508851693E-16,
    10.106686313076017, -4.8084391625516004E-17, -7.0486965582949455E-16,
    9.8718897822287275, -4.682534875748401E-17, -7.034345348420768E-16,
    9.6370970017757, -4.5566305889452004E-17, -7.0199941385465934E-16,
    9.4023046902394611, -4.4306527599184E-17, -7.00563454595637E-16,
    9.1675095668102937, -4.3047484731152E-17, -6.991283336082193E-16,
    8.9327183696834176, -4.1788441863120008E-17, -6.9769321262080174E-16,
    8.6979278173984422, -4.0529398995088E-17, -6.9625809163338418E-16,
    8.4631379685742569, -3.92699148737144E-17, -6.9482246768300367E-16,
    8.2283468781897078, -3.80107249212352E-17, -6.93387179041265E-16,
    7.993557944730604, -3.6751461426532404E-17, -6.9195180657236614E-16,
    7.7587695568412611, -3.5492271474053207E-17, -6.905165179306275E-16,
    7.5239824411033318, -3.4233007979350409E-17, -6.8908114546172831E-16,
    7.289195988207168, -3.2973766547314683E-17, -6.876457981409776E-16,
    7.0544106909915394, 11.280659712816917, 6.4383780855811606E-14,
    -5.5640119678180009E-17, 211.27504957325488, 6.297690859496681E-14,
    -5.3363526231076009E-17, 11.040373840180749, 6.1570036334122E-14,
    -5.1086932783972E-17, 10.805694894794851, 6.0162342302096809E-14,
    -4.8809009550976007E-17, 10.571019337930627, 5.8755470041252E-14,
    -4.6532416103872008E-17, 10.336343956813446, 5.73485977804072E-14,
    -4.4255822656768E-17, 10.101668810070741, 5.59417255195624E-14,
    -4.1979229209664008E-17, 9.866990685561472, 5.4534031487537215E-14,
    -3.9701305976668006E-17, 9.6323161841230132, 5.3127159226692406E-14,
    -3.7424712529564006E-17, 9.3976420929824371, 5.172028696584761E-14,
    -3.514811908246E-17, 9.1629651999901434, 5.0312592933822411E-14,
    -3.2870195849464004E-17, 8.9282921059698044, 4.8905720672977608E-14,
    -3.0593602402360005E-17, 8.6936195981553812, 4.7498848412132812E-14,
    -2.8317008955256005E-17, 8.458947735182857, 4.6091976151288E-14,
    -2.6040415508152002E-17, 8.22427461321002, 4.4684610827734969E-14,
    -2.3763024189512803E-17, 7.9896035620420482, 4.3277574212654087E-14,
    -2.1486164785230402E-17, 7.7549330046841192, 4.1870455420455166E-14,
    -1.9209172402358806E-17, 7.5202636539824406, 4.0463418805374284E-14,
    -1.6932312998076405E-17, 7.2855949143456868, 3.9056300013175362E-14,
    -1.4655320615204802E-17, 7.0509272696833909, 3.764920587411186E-14,
    -1.2378368125909962E-17, 6.40974107466916E-14, 11.280659712816917,
    -7.1348205829720975E-16, 6.297690859496681E-14, 211.27504957325488,
    -7.16554778601069E-16, 6.1856406443242014E-14, 11.040373840180749,
    -7.19627498904928E-16, 6.0735249792596817E-14, 10.805694894794851,
    -7.2270201402204886E-16, 5.9614747640872008E-14, 10.571019337930627,
    -7.2577473432590809E-16, 5.84942454891472E-14, 10.336343956813446,
    -7.2884745462976733E-16, 5.73737433374224E-14, 10.101668810070741,
    -7.3192017493362656E-16, 5.6252586686777206E-14, 9.866990685561472,
    -7.3499469005074731E-16, 5.5132084535052411E-14, 9.6323161841230132,
    -7.3806741035460654E-16, 5.4011582383327615E-14, 9.3976420929824371,
    -7.4114013065846577E-16, 5.2890425732682405E-14, 9.1629651999901434,
    -7.4421464577558652E-16, 5.1769923580957596E-14, 8.9282921059698044,
    -7.4728736607944565E-16, 5.06494214292328E-14, 8.6936195981553812,
    -7.5036008638330488E-16, 4.9528919277508004E-14, 8.458947735182857,
    -7.5343280668716412E-16, 4.8408024426430969E-14, 8.22427461321002,
    -7.565066038789802E-16, 4.7287391374922089E-14, 7.9896035620420482,
    -7.5957968314549181E-16, 4.6166692873521164E-14, 7.7549330046841192,
    -7.6265294189332942E-16, 4.504605982201229E-14, 7.5202636539824406,
    -7.6572602115984094E-16, 4.3925361320611365E-14, 7.2855949143456868,
    -7.6879927990767874E-16, 4.2804682454178059E-14, 7.0509272696833909,
    -7.718724848111185E-16, -5.4573129935716004E-17, -7.1795437889196169E-16,
    11.280659712816917, -5.3363526231076015E-17, -7.1655477860106879E-16,
    211.27504957325488, -5.2153922526436008E-17, -7.1515517831017617E-16,
    11.040373840180749, -5.0943612277576019E-17, -7.1375476049574893E-16,
    10.805694894794851, -4.9734008572936006E-17, -7.1235516020485612E-16,
    10.571019337930627, -4.8524404868296004E-17, -7.1095555991396331E-16,
    10.336343956813446, -4.7314801163656003E-17, -7.095559596230705E-16,
    10.101668810070741, -4.6104490914796008E-17, -7.0815554180864326E-16,
    9.866990685561472, -4.4894887210156007E-17, -7.0675594151775045E-16,
    9.6323161841230132, -4.3685283505516006E-17, -7.0535634122685774E-16,
    9.3976420929824371, -4.2474973256656005E-17, -7.0395592341243049E-16,
    9.1629651999901434, -4.1265369552016004E-17, -7.0255632312153768E-16,
    8.9282921059698044, -4.0055765847376009E-17, -7.0115672283064487E-16,
    8.6936195981553812, -3.8846162142736007E-17, -6.9975712253975216E-16,
    8.458947735182857, -3.7636134511564006E-17, -6.9835703173473867E-16,
    8.22427461321002, -3.6426389498080005E-17, -6.9695726793913894E-16,
    7.9896035620420482, -3.5216573830174007E-17, -6.95557422391186E-16,
    7.7549330046841192, -3.4006828816690006E-17, -6.9415765859558615E-16,
    7.5202636539824406, -3.2797013148784008E-17, -6.92757813047633E-16,
    7.2855949143456868, -3.1587218677204607E-17, -6.91357992025386E-16,
    7.0509272696833909, 11.045865348887073, 6.3233965746954E-14,
    -5.4381076810148E-17, 11.040373840180749, 6.1856406443242014E-14,
    -5.2153922526436008E-17, 211.0348823313123, 6.047884713953E-14,
    -4.9926768242724E-17, 10.800322026823309, 5.9100483186692E-14,
    -4.7698313051136007E-17, 10.565764983533031, 5.772292388298E-14,
    -4.5471158767424005E-17, 10.331208057362378, 5.6345364579268006E-14,
    -4.3244004483712004E-17, 10.096651306938771, 5.4967805275556E-14,
    -4.10168502E-17, 9.8620915887812188, 5.3589441322718012E-14,
    -3.8788395008412009E-17, 9.6275353663675709, 5.2211882019006008E-14,
    -3.6561240724700008E-17, 9.3929794956226829, 5.083432271529401E-14,
    -3.4334086440988007E-17, 9.158420833067261, 4.9455958762456008E-14,
    -3.21056312494E-17, 8.9238658421705868, 4.8078399458744004E-14,
    -2.9878476965688E-17, 8.68931137884382, 4.6700840155032012E-14,
    -2.7651322681976E-17, 8.4547575017229679, 4.5323280851320008E-14,
    -2.5424168398264003E-17, 8.220202348178967, 4.3945238758132404E-14,
    -2.31962335698264E-17, 7.9856491793192488, 4.2567518524595204E-14,
    -2.09688191045392E-17, 7.75109645249273, 4.1189717826145406E-14,
    -1.87412745484644E-17, 7.5165448668273047, 3.9811997592608206E-14,
    -1.6513860083177202E-17, 7.2819938404670834, 3.8434196894158408E-14,
    -1.4286315527102403E-17, 7.0474438483752415, 3.7056420335182388E-14,
    -1.2058809998263881E-17, 6.2661225528714012E-14, 11.045865348887073,
    -7.1204693730979209E-16, 6.1570036334122E-14, 11.040373840180749,
    -7.1515517831017617E-16, 6.047884713953E-14, 211.0348823313123,
    -7.1826341931056006E-16, 5.938702056807201E-14, 10.800322026823309,
    -7.21373475872276E-16, 5.8295831373480012E-14, 10.565764983533031,
    -7.2448171687266012E-16, 5.7204642178888E-14, 10.331208057362378,
    -7.275899578730441E-16, 5.6113452984296004E-14, 10.096651306938771,
    -7.3069819887342818E-16, 5.5021626412838004E-14, 9.8620915887812188,
    -7.3380825543514416E-16, 5.3930437218246012E-14, 9.6275353663675709,
    -7.3691649643552814E-16, 5.2839248023654008E-14, 9.3929794956226829,
    -7.4002473743591213E-16, 5.1747421452196008E-14, 9.158420833067261,
    -7.431347939976281E-16, 5.0656232257604E-14, 8.9238658421705868,
    -7.46243034998012E-16, 4.9565043063012E-14, 8.68931137884382,
    -7.4935127599839607E-16, 4.847385386842E-14, 8.4547575017229679,
    -7.5245951699878006E-16, 4.738228224770841E-14, 8.220202348178967,
    -7.5556884733596326E-16, 4.6290965577743205E-14, 7.9856491793192488,
    -7.5867745144861368E-16, 4.5199585170091404E-14, 7.75109645249273,
    -7.6178623711739732E-16, 4.4108268500126205E-14, 7.5165448668273047,
    -7.6489484123004764E-16, 4.301688809247441E-14, 7.2819938404670834,
    -7.6800362689883138E-16, 4.1925526806128586E-14, 7.0474438483752415,
    -7.71112358100775E-16, -5.2247097325220006E-17, -7.2099157849929608E-16,
    11.045865348887073, -5.1086932783972009E-17, -7.19627498904928E-16,
    11.040373840180749, -4.9926768242724E-17, -7.1826341931056016E-16,
    211.0348823313123, -4.8765926035272008E-17, -7.1689854294072815E-16,
    10.800322026823309, -4.7605761494024E-17, -7.1553446334636019E-16,
    10.565764983533031, -4.6445596952776E-17, -7.1417038375199213E-16,
    10.331208057362378, -4.5285432411528E-17, -7.1280630415762407E-16,
    10.096651306938771, -4.4124590204076E-17, -7.1144142778779216E-16,
    9.8620915887812188, -4.2964425662828004E-17, -7.10077348193424E-16,
    9.6275353663675709, -4.1804261121579995E-17, -7.0871326859905614E-16,
    9.3929794956226829, -4.0643418914128E-17, -7.0734839222922413E-16,
    9.158420833067261, -3.948325437288E-17, -7.0598431263485607E-16,
    8.9238658421705868, -3.832308983163201E-17, -7.0462023304048811E-16,
    8.68931137884382, -3.7162925290384007E-17, -7.0325615344612014E-16,
    8.4547575017229679, -3.60023541494136E-17, -7.0189159578647377E-16,
    8.220202348178967, -3.48420540749248E-17, -7.0052735683701285E-16,
    7.9856491793192488, -3.3681686233815604E-17, -6.9916303821000583E-16,
    7.75109645249273, -3.2521386159326805E-17, -6.977987992605449E-16,
    7.5165448668273047, -3.1361018318217607E-17, -6.9643448063353768E-16,
    7.2819938404670834, -3.0200670807094524E-17, -6.950701859097945E-16,
    7.0474438483752415, 10.811067762614019, 6.2083479017121611E-14,
    -5.3121298519880006E-17, 10.805694894794851, 6.0735249792596817E-14,
    -5.0943612277576007E-17, 10.800322026823309, 5.938702056807201E-14,
    -4.8765926035272E-17, 210.7949460252814, 5.8038003826476808E-14,
    -4.6586967779976004E-17, 10.560507564234927, 5.6689774601952E-14,
    -4.4409281537672005E-17, 10.326069161679428, 5.5341545377427206E-14,
    -4.2231595295368006E-17, 10.091630876243558, 5.39933161529024E-14,
    -4.0053909053064007E-17, 9.8571896331465076, 5.264429941130721E-14,
    -3.7874950797768003E-17, 9.6227517584314572, 5.1296070186782415E-14,
    -3.5697264555464004E-17, 9.38831417675259, 4.9947840962257608E-14,
    -3.3519578313160005E-17, 9.1538738133342612, 4.8598824220662412E-14,
    -3.1340620057864E-17, 8.91943699422134, 4.7250594996137605E-14,
    -2.916293381556E-17, 8.6850006440594463, 4.5902365771612811E-14,
    -2.6985247573256004E-17, 8.4505648214674611, 4.455413654708801E-14,
    -2.4807561330952004E-17, 8.2161277050364561, 4.3205434812320967E-14,
    -2.2629111880852803E-17, 7.9816924871702248, 4.185704808438209E-14,
    -2.0451171235950403E-17, 7.7472576595643634, 4.0508582604736168E-14,
    -1.8273103389748802E-17, 7.5128239076032974, 3.9160195876797285E-14,
    -1.6095162744846402E-17, 7.2783906631917308, 3.7811730397151363E-14,
    -1.3917094898644802E-17, 7.0439583923583848, 3.6463288543017558E-14,
    -1.1739065212832962E-17, 6.1224201417501615E-14, 10.811067762614019,
    -7.1061097805076971E-16, 6.0162342302096809E-14, 10.805694894794851,
    -7.13754760495749E-16, 5.9100483186692014E-14, 10.800322026823309,
    -7.16898542940728E-16, 5.8038003826476808E-14, 210.7949460252814,
    -7.2004416170722881E-16, 5.6976144711072014E-14, 10.560507564234927,
    -7.2318794415220813E-16, 5.59142855956672E-14, 10.326069161679428,
    -7.2633172659718725E-16, 5.4852426480262406E-14, 10.091630876243558,
    -7.2947550904216656E-16, 5.3789947120047206E-14, 9.8571896331465076,
    -7.3262112780866732E-16, 5.2728088004642412E-14, 9.6227517584314572,
    -7.3576491025364654E-16, 5.1666228889237611E-14, 9.38831417675259,
    -7.3890869269862576E-16, 5.0603749529022411E-14, 9.1538738133342612,
    -7.4205431146512652E-16, 4.9541890413617604E-14, 8.91943699422134,
    -7.4519809391010564E-16, 4.8480031298212804E-14, 8.6850006440594463,
    -7.4834187635508486E-16, 4.741817218280801E-14, 8.4505648214674611,
    -7.5148565880006407E-16, 4.6355940920516972E-14, 8.2161277050364561,
    -7.546305430379562E-16, 4.529395775615009E-14, 7.9816924871702248,
    -7.5777469274723978E-16, 4.4231912567302164E-14, 7.7472576595643634,
    -7.609190260886755E-16, 4.3169929402935289E-14, 7.5128239076032974,
    -7.64063175797959E-16, 4.210788421408737E-14, 7.2783906631917308,
    -7.6720750913939471E-16, 4.1045857632583755E-14, 7.0439583923583848,
    -7.7035178739118468E-16, -4.9919706050816008E-17, -7.2403055217182169E-16,
    10.811067762614019, -4.8809009550976013E-17, -7.2270201402204876E-16,
    10.805694894794851, -4.7698313051136007E-17, -7.2137347587227613E-16,
    10.800322026823309, -4.6586967779976017E-17, -7.2004416170722891E-16,
    210.7949460252814, -4.5476271280136004E-17, -7.1871562355745618E-16,
    10.560507564234927, -4.4365574780296E-17, -7.1738708540768326E-16,
    10.326069161679428, -4.3254878280456E-17, -7.1605854725791053E-16,
    10.091630876243558, -4.2143533009296006E-17, -7.1472923309286331E-16,
    9.8571896331465076, -4.1032836509456012E-17, -7.1340069494309038E-16,
    9.6227517584314572, -3.9922140009616005E-17, -7.1207215679331775E-16,
    9.38831417675259, -3.8810794738456E-17, -7.1074284262827053E-16,
    9.1538738133342612, -3.7700098238616E-17, -7.094143044784977E-16,
    8.91943699422134, -3.6589401738776014E-17, -7.0808576632872487E-16,
    8.6850006440594463, -3.5478705238936007E-17, -7.0675722817895215E-16,
    8.4505648214674611, -3.4367619476304009E-17, -7.054282244200147E-16,
    8.2161277050364561, -3.3256793222200003E-17, -7.0409953106718694E-16,
    7.9816924871702248, -3.2145902090964004E-17, -7.027707601128319E-16,
    7.7472576595643634, -3.1035075836860004E-17, -7.0144206676000413E-16,
    7.5128239076032974, -2.9924184705624005E-17, -7.00113295805649E-16,
    7.2783906631917308, -2.8813313037527604E-17, -6.9878454813175221E-16,
    7.0439583923583848, 10.576273692184389, 6.0933663908264E-14,
    -5.1862255651848006E-17, 10.571019337930625, 5.9614747640872008E-14,
    -4.9734008572936006E-17, 10.565764983533031, 5.829583137348E-14,
    -4.7605761494024E-17, 10.560507564234925, 5.6976144711072E-14,
    -4.5476271280136004E-17, 210.55525321008787, 5.5657228443680005E-14,
    -4.3348024201224003E-17, 10.320933262462949, 5.4338312176288008E-14,
    -4.1219777122312E-17, 10.086613373329026, 5.3019395908896E-14,
    -3.90915300434E-17, 9.8522905365665938, 5.1699709246488014E-14,
    -3.6962039829512007E-17, 9.61797094086438, 5.0380792979096011E-14,
    -3.4833792750600006E-17, 9.3836515795709019, 4.9061876711704008E-14,
    -3.2705545671688E-17, 9.1493294465654724, 4.774219004929601E-14,
    -3.05760554578E-17, 8.9150107305419777, 4.6423273781904E-14,
    -2.8447808378888004E-17, 8.6806924248506281, 4.5104357514512005E-14,
    -2.6319561299976003E-17, 8.4463745881103076, 4.3785441247120008E-14,
    -2.4191314221064002E-17, 8.2120554400910155, 4.2466062742718408E-14,
    -2.2062321261166402E-17, 7.9777381045159146, 4.1146992396323207E-14,
    -1.99338255552592E-17, 7.7434211074414661, 3.9827845010426408E-14,
    -1.7805205535854404E-17, 7.5091051204995294, 3.8508774664031207E-14,
    -1.5676709829947203E-17, 7.2747895893302488, 3.7189627278134408E-14,
    -1.3548089810542402E-17, 7.0404749710502355, 3.5870503004088086E-14,
    -1.1419507085186881E-17, 5.9788016199524008E-14, 10.576273692184389,
    -7.0917585706335215E-16, 5.8755470041252012E-14, 10.571019337930625,
    -7.1235516020485622E-16, 5.772292388298E-14, 10.565764983533031,
    -7.1553446334636009E-16, 5.6689774601952007E-14, 10.560507564234925,
    -7.18715623557456E-16, 5.5657228443680011E-14, 210.55525321008787,
    -7.2189492669896015E-16, 5.4624682285407996E-14, 10.320933262462949,
    -7.2507422984046412E-16, 5.3592136127136E-14, 10.086613373329026,
    -7.2825353298196818E-16, 5.2558986846108009E-14, 9.8522905365665938,
    -7.3143469319306408E-16, 5.1526440687836014E-14, 9.61797094086438,
    -7.3461399633456815E-16, 5.0493894529564011E-14, 9.3836515795709019,
    -7.3779329947607221E-16, 4.9460745248536E-14, 9.1493294465654724,
    -7.4097445968716811E-16, 4.8428199090264E-14, 8.9150107305419777,
    -7.44153762828672E-16, 4.7395652931992004E-14, 8.6806924248506281,
    -7.4733306597017614E-16, 4.6363106773720008E-14, 8.4463745881103076,
    -7.5051236911168011E-16, 4.5330198741794406E-14, 8.2120554400910155,
    -7.5369278649493926E-16, 4.4297531958971206E-14, 7.9777381045159146,
    -7.5687246105036175E-16, 4.3264804863872404E-14, 7.7434211074414661,
    -7.6005232131274331E-16, 4.2232138081049204E-14, 7.5091051204995294,
    -7.632319958681657E-16, 4.11994109859504E-14, 7.2747895893302488,
    -7.6641185613054735E-16, 4.0166701984534282E-14, 7.0404749710502355,
    -7.6959166068084114E-16, -4.759367344032001E-17, -7.2706775177915607E-16,
    10.576273692184389, -4.6532416103872008E-17, -7.25774734325908E-16,
    10.571019337930625, -4.5471158767424005E-17, -7.2448171687266012E-16,
    10.565764983533031, -4.4409281537672005E-17, -7.2318794415220813E-16,
    10.560507564234925, -4.3348024201224003E-17, -7.2189492669896015E-16,
    210.55525321008787, -4.2286766864776E-17, -7.2060190924571207E-16,
    10.320933262462949, -4.1225509528328E-17, -7.193088917924641E-16,
    10.086613373329026, -4.0163632298576005E-17, -7.1801511907201211E-16,
    9.8522905365665938, -3.9102374962128009E-17, -7.16722101618764E-16,
    9.61797094086438, -3.8041117625680007E-17, -7.1542908416551615E-16,
    9.3836515795709019, -3.6979240395928E-17, -7.1413531144506416E-16,
    9.1493294465654724, -3.5917983059480004E-17, -7.1284229399181608E-16,
    8.9150107305419777, -3.4856725723032E-17, -7.1154927653856811E-16,
    8.6806924248506281, -3.3795468386584006E-17, -7.1025625908532013E-16,
    8.4463745881103076, -3.2733839114153604E-17, -7.0896278847174971E-16,
    8.2120554400910155, -3.16724577990448E-17, -7.0766961996506085E-16,
    7.9777381045159146, -3.06110144946056E-17, -7.0637637593165174E-16,
    7.7434211074414661, -2.95496331794968E-17, -7.0508320742496289E-16,
    7.5091051204995294, -2.848818987505761E-17, -7.0378996339155359E-16,
    7.2747895893302488, -2.7426765167417521E-17, -7.0249674201616071E-16,
    7.0404749710502355, 10.341479856129229, 5.9783848799406411E-14,
    -5.0603212783816006E-17, 10.336343956813447, 5.8494245489147212E-14,
    -4.8524404868296004E-17, 10.331208057362378, 5.7204642178888E-14,
    -4.6445596952776E-17, 10.326069161679429, 5.5914285595667207E-14,
    -4.4365574780296E-17, 10.320933262462951, 5.4624682285407996E-14,
    -4.2286766864776E-17, 210.31579736311858, 5.3335078975148804E-14,
    -4.0207958949256E-17, 10.08159587029636, 5.20454756648896E-14,
    -3.8129151033736E-17, 9.8473914398770823, 5.0755119081668811E-14,
    -3.6049128861256004E-17, 9.6131901231962775, 4.9465515771409607E-14,
    -3.3970320945736E-17, 9.3789889822967556, 4.8175912461150409E-14,
    -3.1891513030216E-17, 9.1447850797179253, 4.6885555877929607E-14,
    -2.9811490857736E-17, 8.9105844667941234, 4.55959525676704E-14,
    -2.7732682942216E-17, 8.67638420557332, 4.4306349257411205E-14,
    -2.5653875026696002E-17, 8.4421843546846613, 4.3016745947152007E-14,
    -2.3575067111176E-17, 8.2079831750942063, 4.1726690673115843E-14,
    -2.1495530641480002E-17, 7.97378372182736, 4.0436936708264323E-14,
    -1.9416479874568E-17, 7.7395845552843232, 3.9147107416116648E-14,
    -1.7337307681960003E-17, 7.5053863333786381, 3.7857353451265122E-14,
    -1.5258256915048003E-17, 7.2711885154687677, 3.6567524159117441E-14,
    -1.3179084722440003E-17, 7.036991549742087, 3.5277717465158614E-14,
    -1.1099948957540802E-17, 5.8351830981546414E-14, 10.341479856129229,
    -7.0774073607593449E-16, 5.73485977804072E-14, 10.336343956813447,
    -7.1095555991396341E-16, 5.6345364579268006E-14, 10.331208057362378,
    -7.14170383751992E-16, 5.5341545377427206E-14, 10.326069161679429,
    -7.1738708540768326E-16, 5.4338312176288008E-14, 10.320933262462951,
    -7.2060190924571207E-16, 5.33350789751488E-14, 210.31579736311858,
    -7.2381673308374089E-16, 5.23318457740096E-14, 10.08159587029636,
    -7.2703155692176971E-16, 5.13280265721688E-14, 9.8473914398770823,
    -7.3024825857746094E-16, 5.03247933710296E-14, 9.6131901231962775,
    -7.3346308241548965E-16, 4.9321560169890405E-14, 9.3789889822967556,
    -7.3667790625351857E-16, 4.8317740968049605E-14, 9.1447850797179253,
    -7.398946079092098E-16, 4.7314507766910395E-14, 8.9105844667941234,
    -7.4310943174723842E-16, 4.63112745657712E-14, 8.67638420557332,
    -7.4632425558526733E-16, 4.5308041364632006E-14, 8.4421843546846613,
    -7.4953907942329605E-16, 4.4304456563071847E-14, 8.2079831750942063,
    -7.5275502995192232E-16, 4.3301106161792329E-14, 7.97378372182736,
    -7.5597022935348362E-16, 4.2297697160442644E-14, 7.7395845552843232,
    -7.5918561653681111E-16, 4.1294346759163126E-14, 7.5053863333786381,
    -7.6240081593837241E-16, 4.0290937757813447E-14, 7.2711885154687677,
    -7.656162031217E-16, 3.928754633648481E-14, 7.036991549742087,
    -7.688315339704976E-16, -4.5267640829824E-17, -7.3010495138649046E-16,
    10.341479856129229, -4.4255822656768008E-17, -7.2884745462976723E-16,
    10.336343956813447, -4.3244004483712E-17, -7.275899578730441E-16,
    10.331208057362378, -4.2231595295368006E-17, -7.2633172659718734E-16,
    10.326069161679429, -4.1219777122311996E-17, -7.2507422984046412E-16,
    10.320933262462951, -4.0207958949256E-17, -7.2381673308374089E-16,
    210.31579736311858, -3.9196140776199995E-17, -7.2255923632701766E-16,
    10.08159587029636, -3.8183731587856003E-17, -7.2130100505116091E-16,
    9.8473914398770823, -3.7171913414800006E-17, -7.2004350829443768E-16,
    9.6131901231962775, -3.6160095241743996E-17, -7.1878601153771455E-16,
    9.3789889822967556, -3.51476860534E-17, -7.1752778026185779E-16,
    9.1447850797179253, -3.4135867880344E-17, -7.1627028350513447E-16,
    8.9105844667941234, -3.3124049707288003E-17, -7.1501278674841134E-16,
    8.67638420557332, -3.2112231534232E-17, -7.1375528999168811E-16,
    8.4421843546846613, -3.11000587520032E-17, -7.1249735252348481E-16,
    8.2079831750942063, -3.00881223758896E-17, -7.1123970886293476E-16,
    7.97378372182736, -2.90761268982472E-17, -7.0998199175047159E-16,
    7.7395845552843232, -2.80641905221336E-17, -7.0872434808992164E-16,
    7.5053863333786381, -2.7052195044491203E-17, -7.0746663097745827E-16,
    7.2711885154687677, -2.6040217297307439E-17, -7.0620893590056912E-16,
    7.036991549742087, 10.106686313076018, 5.86340336905488E-14,
    -4.9344169915784007E-17, 10.101668810070741, 5.7373743337422417E-14,
    -4.7314801163656009E-17, 10.096651306938771, 5.6113452984296004E-14,
    -4.5285432411528E-17, 10.091630876243558, 5.48524264802624E-14,
    -4.3254878280456009E-17, 10.086613373329026, 5.3592136127136E-14,
    -4.1225509528328005E-17, 10.08159587029636, 5.2331845774009607E-14,
    -3.91961407762E-17, 210.0765783671358, 5.10715554208832E-14,
    -3.7166772024072004E-17, 9.84249234306944, 4.9810528916849609E-14,
    -3.5136217893000007E-17, 9.60840930541858, 4.8550238563723209E-14,
    -3.3106849140872003E-17, 9.3743263849215861, 4.7289948210596809E-14,
    -3.1077480388744E-17, 9.14024071277792, 4.6028921706563211E-14,
    -2.9046926257672E-17, 8.9061582029675073, 4.4768631353436805E-14,
    -2.7017557505544005E-17, 8.672075986227517, 4.3508341000310405E-14,
    -2.4988188753416002E-17, 8.4379941211905258, 4.2248050647184005E-14,
    -2.2958820001288004E-17, 8.2039109100289078, 4.098731860351329E-14,
    -2.0928740021793602E-17, 7.9698293390874371, 3.9726881020205446E-14,
    -1.88991341938768E-17, 7.7357480030929349, 3.8466369821806889E-14,
    -1.68694098280656E-17, 7.5016675462235023, 3.7205932238499044E-14,
    -1.4839804000148803E-17, 7.2675874415901642, 3.5945421040100487E-14,
    -1.2810079634337603E-17, 7.0335081284339376, 3.4684931926229141E-14,
    -1.0780390829894723E-17, 5.6915645763568807E-14, 10.106686313076018,
    -7.0630561508851693E-16, 5.5941725519562407E-14, 10.101668810070741,
    -7.095559596230706E-16, 5.4967805275556008E-14, 10.096651306938771,
    -7.1280630415762407E-16, 5.3993316152902412E-14, 10.091630876243558,
    -7.1605854725791043E-16, 5.3019395908896012E-14, 10.086613373329026,
    -7.193088917924641E-16, 5.20454756648896E-14, 10.08159587029636,
    -7.2255923632701766E-16, 5.10715554208832E-14, 210.0765783671358,
    -7.2580958086157133E-16, 5.0097066298229604E-14, 9.84249234306944,
    -7.2906182396185769E-16, 4.9123146054223205E-14, 9.60840930541858,
    -7.3231216849641126E-16, 4.8149225810216811E-14, 9.3743263849215861,
    -7.3556251303096493E-16, 4.7174736687563209E-14, 9.14024071277792,
    -7.3881475613125139E-16, 4.6200816443556797E-14, 8.9061582029675073,
    -7.4206510066580486E-16, 4.5226896199550403E-14, 8.672075986227517,
    -7.4531544520035852E-16, 4.425297595554401E-14, 8.4379941211905258,
    -7.4856578973491209E-16, 4.3278714384349288E-14, 8.2039109100289078,
    -7.5181727340890537E-16, 4.2304680364613451E-14, 7.9698293390874371,
    -7.5506799765660558E-16, 4.1330589457012884E-14, 7.7357480030929349,
    -7.58318911760879E-16, 4.0356555437277047E-14, 7.5016675462235023,
    -7.6156963600857911E-16, 3.9382464529676486E-14, 7.2675874415901642,
    -7.6482055011285264E-16, 3.8408390688435343E-14, 7.0335081284339376,
    -7.6807140726015406E-16, -4.2941608219328E-17, -7.3314215099382484E-16,
    10.106686313076018, -4.1979229209664008E-17, -7.3192017493362636E-16,
    10.101668810070741, -4.10168502E-17, -7.3069819887342809E-16,
    10.096651306938771, -4.0053909053064007E-17, -7.2947550904216656E-16,
    10.091630876243558, -3.90915300434E-17, -7.2825353298196818E-16,
    10.086613373329026, -3.8129151033736009E-17, -7.2703155692176971E-16,
    10.08159587029636, -3.7166772024072E-17, -7.2580958086157133E-16,
    210.0765783671358, -3.6203830877136008E-17, -7.2458689103030971E-16,
    9.84249234306944, -3.5241451867472E-17, -7.2336491497011123E-16,
    9.60840930541858, -3.4279072857808E-17, -7.22142938909913E-16,
    9.3743263849215861, -3.3316131710872E-17, -7.2092024907865133E-16,
    9.14024071277792, -3.2353752701208E-17, -7.1969827301845285E-16,
    8.9061582029675073, -3.1391373691544004E-17, -7.1847629695825447E-16,
    8.672075986227517, -3.0428994681880005E-17, -7.172543208980561E-16,
    8.4379941211905258, -2.94662783898528E-17, -7.1603191657521981E-16,
    8.2039109100289078, -2.8503786952734406E-17, -7.1480979776080867E-16,
    7.9698293390874371, -2.7541239301888804E-17, -7.1358760756929143E-16,
    7.7357480030929349, -2.6578747864770404E-17, -7.1236548875488039E-16,
    7.5016675462235023, -2.5616200213924805E-17, -7.11143298563363E-16,
    7.2675874415901642, -2.4653669427197363E-17, -7.0992112978497752E-16,
    7.0335081284339376, 9.8718897822287275, 5.74835469607164E-14,
    -4.8084391625516004E-17, 9.866990685561472, 5.6252586686777213E-14,
    -4.6104490914796008E-17, 9.8620915887812188, 5.5021626412838004E-14,
    -4.4124590204076E-17, 9.85718963314651, 5.3789947120047206E-14,
    -4.2143533009296006E-17, 9.8522905365665956, 5.2558986846108003E-14,
    -4.0163632298576005E-17, 9.8473914398770823, 5.1328026572168807E-14,
    -3.8183731587856003E-17, 9.84249234306944, 5.0097066298229604E-14,
    -3.6203830877136008E-17, 209.83759038740612, 4.8865387005438813E-14,
    -3.4222773682356007E-17, 9.603625697455076, 4.763442673149961E-14,
    -3.2242872971636006E-17, 9.36966106602924, 4.6403466457560413E-14,
    -3.0262972260916004E-17, 9.1356936930295056, 4.5171787164769609E-14,
    -2.8281915066136004E-17, 8.9017293550011374, 4.3940826890830407E-14,
    -2.6302014355416002E-17, 8.6677652514157444, 4.270986661689121E-14,
    -2.4322113644696E-17, 8.4338014409007762, 4.1478906342952007E-14,
    -2.2342212933976003E-17, 8.19983626685215, 4.0247514657701848E-14,
    -2.0361618332820005E-17, 7.9658726469041685, 3.9016410579992326E-14,
    -1.8381486325288E-17, 7.7319092101474451, 3.7785234600397651E-14,
    -1.6401238669350002E-17, 7.497946586999495, 3.6554130522688123E-14,
    -1.4421106661818004E-17, 7.2639842643148116, 3.5322954543093442E-14,
    -1.2440859005880003E-17, 7.0300226724170809, 3.4091800134064311E-14,
    -1.0460646044463803E-17, 5.5478621652356411E-14, 9.8718897822287275,
    -7.0486965582949455E-16, 5.45340314875372E-14, 9.866990685561472,
    -7.0815554180864336E-16, 5.3589441322718006E-14, 9.8620915887812188,
    -7.1144142778779206E-16, 5.264429941130721E-14, 9.85718963314651,
    -7.1472923309286321E-16, 5.1699709246488014E-14, 9.8522905365665956,
    -7.1801511907201211E-16, 5.0755119081668805E-14, 9.8473914398770823,
    -7.2130100505116091E-16, 4.9810528916849609E-14, 9.84249234306944,
    -7.2458689103030971E-16, 4.8865387005438806E-14, 209.83759038740612,
    -7.27874696335381E-16, 4.792079684061961E-14, 9.603625697455076,
    -7.3116058231452966E-16, 4.6976206675800408E-14, 9.36966106602924,
    -7.3444646829367856E-16, 4.6031064764389605E-14, 9.1356936930295056,
    -7.3773427359874971E-16, 4.50864745995704E-14, 8.9017293550011374,
    -7.4102015957789841E-16, 4.4141884434751207E-14, 8.6677652514157444,
    -7.4430604555704731E-16, 4.3197294269932011E-14, 8.4338014409007762,
    -7.4759193153619611E-16, 4.2252373057157849E-14, 8.19983626685215,
    -7.5087896911089832E-16, 4.130767254302033E-14, 7.9658726469041685,
    -7.5416523895523169E-16, 4.0362916854223645E-14, 7.7319092101474451,
    -7.5745170073215719E-16, 3.9418216340086125E-14, 7.497946586999495,
    -7.6073797057649036E-16, 3.8473460651289446E-14, 7.2639842643148116,
    -7.64024432353416E-16, 3.7528721514890512E-14, 7.0300226724170809,
    -7.6731083655056378E-16, -4.0614216944924E-17, -7.3618112466635045E-16,
    9.8718897822287275, -3.9701305976668012E-17, -7.3499469005074721E-16,
    9.866990685561472, -3.8788395008412E-17, -7.3380825543514416E-16,
    9.8620915887812188, -3.787495079776801E-17, -7.3262112780866732E-16,
    9.85718963314651, -3.6962039829512E-17, -7.3143469319306418E-16,
    9.8522905365665956, -3.6049128861256004E-17, -7.3024825857746084E-16,
    9.8473914398770823, -3.5136217893E-17, -7.2906182396185769E-16,
    9.84249234306944, -3.4222773682356007E-17, -7.2787469633538086E-16,
    209.83759038740612, -3.3309862714100004E-17, -7.2668826171977761E-16,
    9.603625697455076, -3.2396951745844E-17, -7.2550182710417457E-16,
    9.36966106602924, -3.14835075352E-17, -7.2431469947769773E-16,
    9.1356936930295056, -3.0570596566944005E-17, -7.2312826486209449E-16,
    8.9017293550011374, -2.9657685598688008E-17, -7.2194183024649134E-16,
    8.6677652514157444, -2.8744774630432005E-17, -7.207553956308881E-16,
    8.4338014409007762, -2.78315437167432E-17, -7.1956854520876074E-16,
    8.19983626685215, -2.6918526100009604E-17, -7.1838197199098276E-16,
    7.9658726469041685, -2.6005455159037204E-17, -7.171953294721176E-16,
    7.7319092101474451, -2.5092437542303603E-17, -7.1600875625433962E-16,
    7.497946586999495, -2.4179366601331206E-17, -7.1482211373547427E-16,
    7.2639842643148116, -2.3266311657630442E-17, -7.1363549200693533E-16,
    7.0300226724170809, 9.6370970017757, 5.6333731851858809E-14,
    -4.682534875748401E-17, 9.6323161841230132, 5.5132084535052411E-14,
    -4.4894887210156007E-17, 9.6275353663675709, 5.3930437218246E-14,
    -4.2964425662828004E-17, 9.6227517584314572, 5.2728088004642405E-14,
    -4.1032836509456006E-17, 9.6179709408643816, 5.1526440687836E-14,
    -3.9102374962128E-17, 9.6131901231962757, 5.03247933710296E-14,
    -3.7171913414800006E-17, 9.60840930541858, 4.9123146054223205E-14,
    -3.5241451867472E-17, 9.603625697455076, 4.792079684061961E-14,
    -3.3309862714100004E-17, 209.5988448798594, 4.6719149523813206E-14,
    -3.1379401166772008E-17, 9.3649984688201577, 4.5517502207006814E-14,
    -2.9448939619444005E-17, 9.13114932623846, 4.4315152993403207E-14,
    -2.7517350466072E-17, 8.89730309130636, 4.3113505676596809E-14,
    -2.5586888918744003E-17, 8.6634570321898057, 4.191185835979041E-14,
    -2.3656427371416E-17, 8.4296112075162242, 4.0710211042984006E-14,
    -2.1725965824088003E-17, 8.1957640018724653, 3.9508142588099289E-14,
    -1.9794827713133604E-17, 7.9619182642156119, 3.8306354891933449E-14,
    -1.7864140644596802E-17, 7.7280726579903023, 3.7104497006087891E-14,
    -1.5933340815455604E-17, 7.4942277998786047, 3.5902709309922045E-14,
    -1.4002653746918804E-17, 7.26038319045333, 3.4700851424076487E-14,
    -1.2071853917777603E-17, 7.0265392511089324, 3.3499014595134839E-14,
    -1.0141087916817723E-17, 5.4042436434378816E-14, 9.6370970017757,
    -7.034345348420769E-16, 5.3127159226692406E-14, 9.6323161841230132,
    -7.0675594151775055E-16, 5.2211882019006008E-14, 9.6275353663675709,
    -7.100773481934241E-16, 5.1296070186782409E-14, 9.6227517584314572,
    -7.1340069494309048E-16, 5.0380792979096011E-14, 9.6179709408643816,
    -7.1672210161876413E-16, 4.94655157714096E-14, 9.6131901231962757,
    -7.2004350829443768E-16, 4.8550238563723209E-14, 9.60840930541858,
    -7.2336491497011133E-16, 4.7634426731499604E-14, 9.603625697455076,
    -7.2668826171977771E-16, 4.6719149523813206E-14, 209.5988448798594,
    -7.3000966839545126E-16, 4.5803872316126808E-14, 9.3649984688201577,
    -7.33331075071125E-16, 4.48880604839032E-14, 9.13114932623846,
    -7.3665442182079139E-16, 4.39727832762168E-14, 8.89730309130636,
    -7.3997582849646485E-16, 4.3057506068530407E-14, 8.6634570321898057,
    -7.432972351721385E-16, 4.2142228860844009E-14, 8.4296112075162242,
    -7.46618641847812E-16, 4.122663087843529E-14, 8.1957640018724653,
    -7.4994121256788138E-16, 4.0311246745841452E-14, 7.9619182642156119,
    -7.5326300725835355E-16, 3.9395809150793884E-14, 7.7280726579903023,
    -7.56584995956225E-16, 3.8480425018200046E-14, 7.4942277998786047,
    -7.5990679064669707E-16, 3.7564987423152485E-14, 7.26038319045333,
    -7.6322877934456861E-16, 3.6649565866841039E-14, 7.0265392511089324,
    -7.6655070984022024E-16, -3.8288184334428005E-17, -7.3921832427368483E-16,
    9.6370970017757, -3.7424712529564006E-17, -7.3806741035460644E-16,
    9.6323161841230132, -3.6561240724700008E-17, -7.3691649643552814E-16,
    9.6275353663675709, -3.5697264555464011E-17, -7.3576491025364654E-16,
    9.6227517584314572, -3.4833792750600006E-17, -7.3461399633456815E-16,
    9.6179709408643816, -3.3970320945736008E-17, -7.3346308241548965E-16,
    9.6131901231962757, -3.3106849140872003E-17, -7.3231216849641126E-16,
    9.60840930541858, -3.2242872971636006E-17, -7.3116058231452966E-16,
    9.603625697455076, -3.1379401166772008E-17, -7.3000966839545126E-16,
    209.5988448798594, -3.0515929361908003E-17, -7.28858754476373E-16,
    9.3649984688201577, -2.9651953192672006E-17, -7.2770716829449136E-16,
    9.13114932623846, -2.8788481387808E-17, -7.2655625437541287E-16,
    8.89730309130636, -2.7925009582944009E-17, -7.2540534045633448E-16,
    8.6634570321898057, -2.7061537778080005E-17, -7.2425442653725618E-16,
    8.4296112075162242, -2.6197763354592803E-17, -7.2310310926049584E-16,
    8.1957640018724653, -2.5334190676854404E-17, -7.2195206088885667E-16,
    7.9619182642156119, -2.4470567562678804E-17, -7.2080094529093745E-16,
    7.7280726579903023, -2.3606994884940405E-17, -7.1964989691929837E-16,
    7.4942277998786047, -2.2743371770764805E-17, -7.1849878132137886E-16,
    7.26038319045333, -2.1879763787520366E-17, -7.1734768589134373E-16,
    7.0265392511089324, 9.4023046902394611, 5.5183916743001208E-14,
    -4.5566305889452004E-17, 9.3976420929824371, 5.4011582383327615E-14,
    -4.3685283505516006E-17, 9.3929794956226829, 5.2839248023654E-14,
    -4.180426112158E-17, 9.38831417675259, 5.1666228889237605E-14,
    -3.9922140009616005E-17, 9.3836515795709019, 5.0493894529564005E-14,
    -3.8041117625680007E-17, 9.3789889822967538, 4.9321560169890405E-14,
    -3.6160095241744E-17, 9.3743263849215843, 4.8149225810216805E-14,
    -3.4279072857808004E-17, 9.3696610660292379, 4.6976206675800414E-14,
    -3.2396951745844008E-17, 9.3649984688201577, 4.5803872316126808E-14,
    -3.0515929361908003E-17, 209.36033587151744, 4.4631537956453208E-14,
    -2.8634906977972005E-17, 9.1266049593635223, 4.3458518822036811E-14,
    -2.6752785866008002E-17, 8.89287682753625, 4.2286184462363204E-14,
    -2.4871763482072004E-17, 8.6591488128970813, 4.1113850102689611E-14,
    -2.2990741098136003E-17, 8.4254209740734574, 3.9941515743016004E-14,
    -2.11097187142E-17, 8.19169173684826, 3.876877051849673E-14,
    -1.9228037093447204E-17, 7.9579638814928124, 3.7596299203874565E-14,
    -1.73467949639056E-17, 7.724236105798914, 3.6423759411778131E-14,
    -1.5465442961561202E-17, 7.4905090127234679, 3.5251288097155967E-14,
    -1.3584200832019603E-17, 7.2567821165747262, 3.4078748305059526E-14,
    -1.1702848829675204E-17, 7.023055829800783, 3.2906229056205367E-14,
    -9.8215297891716419E-18, 5.2606251216401209E-14, 9.4023046902394611,
    -7.0199941385465934E-16, 5.172028696584761E-14, 9.3976420929824371,
    -7.0535634122685774E-16, 5.083432271529401E-14, 9.3929794956226829,
    -7.08713268599056E-16, 4.9947840962257608E-14, 9.38831417675259,
    -7.1207215679331765E-16, 4.9061876711704008E-14, 9.3836515795709019,
    -7.1542908416551615E-16, 4.81759124611504E-14, 9.3789889822967538,
    -7.1878601153771445E-16, 4.72899482105968E-14, 9.3743263849215843,
    -7.22142938909913E-16, 4.6403466457560407E-14, 9.3696610660292379,
    -7.2550182710417457E-16, 4.5517502207006807E-14, 9.3649984688201577,
    -7.2885875447637287E-16, 4.4631537956453208E-14, 209.36033587151744,
    -7.3221568184857137E-16, 4.3745056203416806E-14, 9.1266049593635223,
    -7.35574570042833E-16, 4.28590919528632E-14, 8.89287682753625,
    -7.3893149741503119E-16, 4.1973127702309607E-14, 8.6591488128970813,
    -7.4228842478722968E-16, 4.1087163451756007E-14, 8.4254209740734574,
    -7.4564535215942808E-16, 4.0200888699712731E-14, 8.19169173684826,
    -7.4900345602486434E-16, 3.9314820948662568E-14, 7.9579638814928124,
    -7.5236077556147552E-16, 3.8428701447364124E-14, 7.724236105798914,
    -7.5571829118029279E-16, 3.7542633696313968E-14, 7.4905090127234679,
    -7.5907561071690388E-16, 3.6656514195015524E-14, 7.2567821165747262,
    -7.6243312633572125E-16, 3.5770410218791567E-14, 7.023055829800783,
    -7.657905831298767E-16, -3.5962151723932006E-17, -7.4225552388101932E-16,
    9.4023046902394611, -3.5148119082460007E-17, -7.4114013065846558E-16,
    9.3976420929824371, -3.4334086440988007E-17, -7.4002473743591213E-16,
    9.3929794956226829, -3.3519578313160011E-17, -7.3890869269862576E-16,
    9.38831417675259, -3.2705545671688006E-17, -7.3779329947607211E-16,
    9.3836515795709019, -3.1891513030216006E-17, -7.3667790625351847E-16,
    9.3789889822967538, -3.1077480388744E-17, -7.3556251303096493E-16,
    9.3743263849215843, -3.0262972260916004E-17, -7.3444646829367846E-16,
    9.3696610660292379, -2.9448939619444005E-17, -7.3333107507112481E-16,
    9.3649984688201577, -2.8634906977972E-17, -7.3221568184857137E-16,
    209.36033587151744, -2.7820398850144003E-17, -7.310996371112849E-16,
    9.1266049593635223, -2.7006366208672004E-17, -7.2998424388873125E-16,
    8.89287682753625, -2.6192333567200004E-17, -7.2886885066617771E-16,
    8.6591488128970813, -2.5378300925728004E-17, -7.2775345744362417E-16,
    8.4254209740734574, -2.4563982992442405E-17, -7.2663767331223085E-16,
    8.19169173684826, -2.3749855253699204E-17, -7.2552214978673068E-16,
    7.9579638814928124, -2.2935679966320404E-17, -7.2440656110975729E-16,
    7.724236105798914, -2.2121552227577203E-17, -7.23291037584257E-16,
    7.4905090127234679, -2.1307376940198407E-17, -7.2217544890728354E-16,
    7.2567821165747262, -2.0493215917410284E-17, -7.2105987977575223E-16,
    7.023055829800783, 9.1675095668102937, 5.4033430013168807E-14,
    -4.4306527599184007E-17, 9.1629651999901416, 5.2890425732682411E-14,
    -4.2474973256656005E-17, 9.158420833067261, 5.1747421452196E-14,
    -4.0643418914128E-17, 9.1538738133342612, 5.0603749529022405E-14,
    -3.8810794738456E-17, 9.1493294465654724, 4.9460745248536E-14,
    -3.6979240395928E-17, 9.1447850797179235, 4.8317740968049605E-14,
    -3.51476860534E-17, 9.14024071277792, 4.71747366875632E-14,
    -3.3316131710872E-17, 9.1356936930295038, 4.6031064764389611E-14,
    -3.14835075352E-17, 9.13114932623846, 4.4888060483903209E-14,
    -2.9651953192672006E-17, 9.1266049593635241, 4.3745056203416812E-14,
    -2.7820398850144003E-17, 209.12205793968755, 4.2601384280243209E-14,
    -2.5987774674472E-17, 8.8884479796349467, 4.1458379999756806E-14,
    -2.4156220331944E-17, 8.6548380781418146, 4.0315375719270409E-14,
    -2.2324665989416002E-17, 8.4212282938367871, 3.9172371438784007E-14,
    -2.0493111646888003E-17, 8.1876170937228725, 3.8028966572685287E-14,
    -1.86609154044736E-17, 7.9540071893506381, 3.6885828763661445E-14,
    -1.6829147095316802E-17, 7.7203973128705465, 3.5742624190368887E-14,
    -1.4997271802845603E-17, 7.4867880534994615, 3.4599486381345046E-14,
    -1.3165503493688802E-17, 7.2531789392993735, 3.3456281808052488E-14,
    -1.1333628201217603E-17, 7.0195703737839263, 3.2313097264040543E-14,
    -9.5017850037407223E-18, 5.1169227105188813E-14, 9.1675095668102937,
    -7.0056345459563686E-16, 5.0312592933822404E-14, 9.1629651999901416,
    -7.0395592341243059E-16, 4.9455958762456008E-14, 9.158420833067261,
    -7.07348392229224E-16, 4.8598824220662406E-14, 9.1538738133342612,
    -7.1074284262827043E-16, 4.774219004929601E-14, 9.1493294465654724,
    -7.1413531144506406E-16, 4.68855558779296E-14, 9.1447850797179235,
    -7.1752778026185769E-16, 4.6028921706563205E-14, 9.14024071277792,
    -7.2092024907865133E-16, 4.51717871647696E-14, 9.1356936930295038,
    -7.2431469947769773E-16, 4.4315152993403207E-14, 9.13114932623846,
    -7.2770716829449126E-16, 4.3458518822036811E-14, 9.1266049593635241,
    -7.31099637111285E-16, 4.26013842802432E-14, 209.12205793968755,
    -7.344940875103314E-16, 4.17447501088768E-14, 8.8884479796349467,
    -7.3788655632712483E-16, 4.0888115937510404E-14, 8.6548380781418146,
    -7.4127902514391847E-16, 4.0031481766144008E-14, 8.4212282938367871,
    -7.446714939607121E-16, 3.9174547372521287E-14, 8.1876170937228725,
    -7.4806515172685738E-16, 3.8317813127069447E-14, 7.9540071893506381,
    -7.5145801686010152E-16, 3.7461028844574885E-14, 7.7203973128705465,
    -7.54851080151571E-16, 3.6604294599123045E-14, 7.4867880534994615,
    -7.5824394528481512E-16, 3.5747510316628484E-14, 7.2531789392993735,
    -7.6163700857628468E-16, 3.4890741045246742E-14, 7.0195703737839263,
    -7.6503001242028642E-16, -3.3634760449528008E-17, -7.4529449755354483E-16,
    9.1675095668102937, -3.2870195849464004E-17, -7.4421464577558642E-16,
    9.1629651999901416, -3.21056312494E-17, -7.431347939976281E-16,
    9.158420833067261, -3.1340620057864008E-17, -7.4205431146512652E-16,
    9.1538738133342612, -3.0576055457800004E-17, -7.4097445968716811E-16,
    9.1493294465654724, -2.9811490857736006E-17, -7.398946079092097E-16,
    9.1447850797179235, -2.9046926257671996E-17, -7.3881475613125129E-16,
    9.14024071277792, -2.8281915066136004E-17, -7.3773427359874971E-16,
    9.1356936930295038, -2.7517350466072003E-17, -7.366544218207912E-16,
    9.13114932623846, -2.6752785866008E-17, -7.35574570042833E-16,
    9.1266049593635241, -2.5987774674472E-17, -7.344940875103314E-16,
    209.12205793968755, -2.5223210074408E-17, -7.3341423573237289E-16,
    8.8884479796349467, -2.4458645474344002E-17, -7.3233438395441448E-16,
    8.6548380781418146, -2.3694080874280004E-17, -7.3125453217645617E-16,
    8.4212282938367871, -2.29292483193328E-17, -7.3017430194577187E-16,
    8.1876170937228725, -2.2164594400974402E-17, -7.2909432401690467E-16,
    7.9540071893506381, -2.1399895823468804E-17, -7.2801428301258346E-16,
    7.7203973128705465, -2.0635241905110403E-17, -7.2693430508371636E-16,
    7.4867880534994615, -1.9870543327604802E-17, -7.25854264079395E-16,
    7.2531789392993735, -1.9105858147843363E-17, -7.2477424199770994E-16,
    7.0195703737839263, 8.9327183696834176, 5.2883614904311207E-14,
    -4.3047484731152008E-17, 8.9282921059698044, 5.1769923580957609E-14,
    -4.126536955201601E-17, 8.923865842170585, 5.0656232257604E-14,
    -3.948325437288E-17, 8.91943699422134, 4.9541890413617604E-14,
    -3.7700098238616E-17, 8.9150107305419777, 4.8428199090264E-14,
    -3.5917983059480004E-17, 8.9105844667941216, 4.7314507766910407E-14,
    -3.4135867880344E-17, 8.9061582029675055, 4.6200816443556797E-14,
    -3.2353752701208E-17, 8.9017293550011374, 4.5086474599570409E-14,
    -3.0570596566944005E-17, 8.89730309130636, 4.3972783276216805E-14,
    -2.8788481387808E-17, 8.89287682753625, 4.2859091952863213E-14,
    -2.7006366208672004E-17, 8.8884479796349467, 4.1744750108876806E-14,
    -2.5223210074408E-17, 208.8840217160126, 4.06310587855232E-14,
    -2.3441094895272002E-17, 8.65052985898094, 3.951736746216961E-14,
    -2.1658979716136E-17, 8.4170380605087409, 3.8403676138816005E-14,
    -1.9876864537E-17, 8.1835448287962667, 3.7289594503082728E-14,
    -1.8094124784787203E-17, 7.950052806713451, 3.6175773075602562E-14,
    -1.63118014146256E-17, 7.7165607607544988, 3.5061886596059127E-14,
    -1.4529373948951202E-17, 7.4830692663956926, 3.3948065168578968E-14,
    -1.2747050578789602E-17, 7.2495778654378915, 3.2834178689035527E-14,
    -1.0964623113115202E-17, 7.016086952475777, 3.1720311725111071E-14,
    -9.1822268760946417E-18, 4.9733041887211212E-14, 8.9327183696834176,
    -6.991283336082193E-16, 4.89057206729776E-14, 8.9282921059698044,
    -7.0255632312153778E-16, 4.8078399458744004E-14, 8.923865842170585,
    -7.0598431263485607E-16, 4.7250594996137605E-14, 8.91943699422134,
    -7.094143044784976E-16, 4.6423273781904014E-14, 8.9150107305419777,
    -7.1284229399181608E-16, 4.55959525676704E-14, 8.9105844667941216,
    -7.1627028350513447E-16, 4.4768631353436805E-14, 8.9061582029675055,
    -7.1969827301845295E-16, 4.3940826890830407E-14, 8.9017293550011374,
    -7.2312826486209449E-16, 4.3113505676596809E-14, 8.89730309130636,
    -7.2655625437541287E-16, 4.2286184462363211E-14, 8.89287682753625,
    -7.2998424388873135E-16, 4.14583799997568E-14, 8.8884479796349467,
    -7.33414235732373E-16, 4.06310587855232E-14, 208.8840217160126,
    -7.3684222524569127E-16, 3.9803737571289604E-14, 8.65052985898094,
    -7.4027021475900966E-16, 3.8976416357056006E-14, 8.4170380605087409,
    -7.43698204272328E-16, 3.8148805193798727E-14, 8.1835448287962667,
    -7.4712739518384034E-16, 3.7321387329890569E-14, 7.950052806713451,
    -7.5055578516322349E-16, 3.6493921141145125E-14, 7.7165607607544988,
    -7.5398437537563878E-16, 3.5666503277236967E-14, 7.4830692663956926,
    -7.5741276535502183E-16, 3.4839037088491522E-14, 7.2495778654378915,
    -7.6084135556743732E-16, 3.4011585397197269E-14, 7.016086952475777,
    -7.6426988570994288E-16, -3.1308727839032004E-17, -7.4833169716087931E-16,
    8.9327183696834176, -3.0593602402360005E-17, -7.4728736607944565E-16,
    8.9282921059698044, -2.9878476965688E-17, -7.4624303499801209E-16,
    8.923865842170585, -2.916293381556E-17, -7.4519809391010574E-16,
    8.91943699422134, -2.8447808378888004E-17, -7.4415376282867218E-16,
    8.9150107305419777, -2.7732682942216E-17, -7.4310943174723852E-16,
    8.9105844667941216, -2.7017557505544E-17, -7.42065100665805E-16,
    8.9061582029675055, -2.6302014355416002E-17, -7.4102015957789851E-16,
    8.9017293550011374, -2.5586888918744003E-17, -7.3997582849646485E-16,
    8.89730309130636, -2.4871763482072E-17, -7.3893149741503138E-16,
    8.89287682753625, -2.4156220331944E-17, -7.3788655632712493E-16,
    8.8884479796349467, -2.3441094895272002E-17, -7.3684222524569127E-16,
    208.8840217160126, -2.2725969458600003E-17, -7.3579789416425771E-16,
    8.65052985898094, -2.2010844021928004E-17, -7.3475356308282415E-16,
    8.4170380605087409, -2.12954679571824E-17, -7.3370886599750688E-16,
    8.1835448287962667, -2.0580258977819202E-17, -7.3266441291477868E-16,
    7.950052806713451, -1.98650082271104E-17, -7.3161989883140331E-16,
    7.7165607607544988, -1.9149799247747202E-17, -7.3057544574867511E-16,
    7.4830692663956926, -1.8434548497038404E-17, -7.2953093166529954E-16,
    7.2495778654378915, -1.7719310277733281E-17, -7.2848643588211845E-16,
    7.016086952475777, 8.6979278173984422, 5.1733799795453607E-14,
    -4.1788441863120008E-17, 8.69361959815538, 5.0649421429232813E-14,
    -4.0055765847376009E-17, 8.68931137884382, 4.9565043063012E-14,
    -3.8323089831632E-17, 8.6850006440594463, 4.8480031298212804E-14,
    -3.6589401738776E-17, 8.6806924248506281, 4.7395652931992E-14,
    -3.4856725723032E-17, 8.6763842055733189, 4.6311274565771203E-14,
    -3.3124049707288003E-17, 8.6720759862275152, 4.52268961995504E-14,
    -3.1391373691544004E-17, 8.6677652514157444, 4.4141884434751213E-14,
    -2.9657685598688008E-17, 8.6634570321898057, 4.3057506068530407E-14,
    -2.7925009582944E-17, 8.6591488128970813, 4.1973127702309613E-14,
    -2.6192333567200004E-17, 8.6548380781418146, 4.088811593751041E-14,
    -2.4458645474344002E-17, 8.65052985898094, 3.9803737571289604E-14,
    -2.2725969458600003E-17, 208.64622163976065, 3.871935920506881E-14,
    -2.0993293442856E-17, 8.4128478271310385, 3.7634980838848003E-14,
    -1.9260617427112002E-17, 8.1794725638285666, 3.6550222433480169E-14,
    -1.7527334165100803E-17, 7.9460984240437309, 3.5465717387543685E-14,
    -1.57944557339344E-17, 7.7127242086144783, 3.4381149001749368E-14,
    -1.4061476095056803E-17, 7.479350479281651, 3.3296643955812883E-14,
    -1.2328597663890403E-17, 7.24597679157641, 3.2212075570018566E-14,
    -1.0595618025012802E-17, 7.0126035311676285, 3.11275261861816E-14,
    -8.862668748448561E-18, 4.8296856669233612E-14, 8.6979278173984422,
    -6.9769321262080174E-16, 4.7498848412132806E-14, 8.69361959815538,
    -7.01156722830645E-16, 4.6700840155032006E-14, 8.68931137884382,
    -7.04620233040488E-16, 4.5902365771612804E-14, 8.6850006440594463,
    -7.0808576632872487E-16, 4.5104357514512011E-14, 8.6806924248506281,
    -7.1154927653856811E-16, 4.43063492574112E-14, 8.6763842055733189,
    -7.1501278674841124E-16, 4.3508341000310405E-14, 8.6720759862275152,
    -7.1847629695825457E-16, 4.2709866616891204E-14, 8.6677652514157444,
    -7.2194183024649134E-16, 4.191185835979041E-14, 8.6634570321898057,
    -7.2540534045633448E-16, 4.1113850102689611E-14, 8.6591488128970813,
    -7.2886885066617781E-16, 4.03153757192704E-14, 8.6548380781418146,
    -7.3233438395441458E-16, 3.95173674621696E-14, 8.65052985898094,
    -7.3579789416425761E-16, 3.8719359205068804E-14, 208.64622163976065,
    -7.3926140437410094E-16, 3.7921350947968004E-14, 8.4128478271310385,
    -7.4272491458394408E-16, 3.7123063015076168E-14, 8.1794725638285666,
    -7.4618963864082339E-16, 3.6324961532711686E-14, 7.9460984240437309,
    -7.4965355346634545E-16, 3.5526813437715365E-14, 7.7127242086144783,
    -7.5311767059970668E-16, 3.4728711955350888E-14, 7.479350479281651,
    -7.5658158542522854E-16, 3.3930563860354568E-14, 7.24597679157641,
    -7.6004570255858986E-16, 3.3132429749147796E-14, 7.0126035311676285,
    -7.6350975899959934E-16, -2.8982695228536006E-17, -7.5136889676821369E-16,
    8.6979278173984422, -2.8317008955256005E-17, -7.5036008638330479E-16,
    8.69361959815538, -2.7651322681976004E-17, -7.4935127599839617E-16,
    8.68931137884382, -2.6985247573256007E-17, -7.48341876355085E-16,
    8.6850006440594463, -2.6319561299976003E-17, -7.4733306597017614E-16,
    8.6806924248506281, -2.5653875026696002E-17, -7.4632425558526724E-16,
    8.6763842055733189, -2.4988188753416002E-17, -7.4531544520035852E-16,
    8.6720759862275152, -2.4322113644696E-17, -7.4430604555704731E-16,
    8.6677652514157444, -2.3656427371416003E-17, -7.432972351721384E-16,
    8.6634570321898057, -2.2990741098136E-17, -7.4228842478722978E-16,
    8.6591488128970813, -2.2324665989416E-17, -7.4127902514391857E-16,
    8.6548380781418146, -2.1658979716136E-17, -7.4027021475900966E-16,
    8.65052985898094, -2.0993293442856004E-17, -7.3926140437410094E-16,
    208.64622163976065, -2.0327607169576003E-17, -7.3825259398919213E-16,
    8.4128478271310385, -1.9661687595032002E-17, -7.3724343004924188E-16,
    8.1794725638285666, -1.8995923554664002E-17, -7.3623450181265259E-16,
    7.9460984240437309, -1.8330120630752E-17, -7.3522551465022315E-16,
    7.7127242086144783, -1.7664356590384E-17, -7.3421658641363376E-16,
    7.479350479281651, -1.6998553666472003E-17, -7.3320759925120423E-16,
    7.24597679157641, -1.6332762407623202E-17, -7.3219862976652685E-16,
    7.0126035311676285, 8.4631379685742587, 5.0583984686596006E-14,
    -4.0529398995088E-17, 8.458947735182857, 4.9528919277508011E-14,
    -3.8846162142736007E-17, 8.4547575017229679, 4.847385386842E-14,
    -3.7162925290384E-17, 8.4505648214674611, 4.7418172182808003E-14,
    -3.5478705238936007E-17, 8.4463745881103076, 4.636310677372E-14,
    -3.3795468386584E-17, 8.4421843546846613, 4.5308041364632006E-14,
    -3.2112231534232E-17, 8.4379941211905258, 4.4252975955544E-14,
    -3.0428994681880005E-17, 8.4338014409007762, 4.3197294269932011E-14,
    -2.8744774630432005E-17, 8.4296112075162242, 4.2142228860844009E-14,
    -2.7061537778080005E-17, 8.4254209740734574, 4.1087163451756013E-14,
    -2.5378300925728004E-17, 8.4212282938367871, 4.0031481766144008E-14,
    -2.369408087428E-17, 8.4170380605087409, 3.8976416357056006E-14,
    -2.2010844021928004E-17, 8.4128478271310385, 3.792135094796801E-14,
    -2.0327607169576003E-17, 208.40865759369393, 3.6866285538880008E-14,
    -1.8644370317224003E-17, 8.17540029881121, 3.5810850363877604E-14,
    -1.6960543545414403E-17, 7.9421440413329165, 3.4755661699484808E-14,
    -1.52771100532432E-17, 7.7088876564419255, 3.3700411407439608E-14,
    -1.3593578241162402E-17, 7.4756316921436365, 3.2645222743046805E-14,
    -1.1910144748991203E-17, 7.2423757177046557, 3.1589972451001605E-14,
    -1.0226612936910401E-17, 7.0091201098594791, 3.0534740647252126E-14,
    -8.5431106208024819E-18, 4.6860671451256011E-14, 8.4631379685742587,
    -6.9625809163338408E-16, 4.6091976151288E-14, 8.458947735182857,
    -6.9975712253975216E-16, 4.5323280851320008E-14, 8.4547575017229679,
    -7.0325615344612E-16, 4.455413654708801E-14, 8.4505648214674611,
    -7.06757228178952E-16, 4.3785441247120008E-14, 8.4463745881103076,
    -7.1025625908532013E-16, 4.3016745947152E-14, 8.4421843546846613,
    -7.1375528999168811E-16, 4.2248050647184005E-14, 8.4379941211905258,
    -7.172543208980562E-16, 4.1478906342952E-14, 8.4338014409007762,
    -7.207553956308881E-16, 4.0710211042984006E-14, 8.4296112075162242,
    -7.2425442653725608E-16, 3.9941515743016011E-14, 8.4254209740734574,
    -7.2775345744362417E-16, 3.9172371438784007E-14, 8.4212282938367871,
    -7.3125453217645617E-16, 3.8403676138816E-14, 8.4170380605087409,
    -7.3475356308282405E-16, 3.7634980838848003E-14, 8.4128478271310385,
    -7.3825259398919213E-16, 3.6866285538880008E-14, 208.40865759369393,
    -7.4175162489556012E-16, 3.6097320836353609E-14, 8.17540029881121,
    -7.4525188209780645E-16, 3.5328535735532808E-14, 7.9421440413329165,
    -7.4875132176946732E-16, 3.4559705734285604E-14, 7.7088876564419255,
    -7.5225096582377448E-16, 3.3790920633464804E-14, 7.4756316921436365,
    -7.5575040549543525E-16, 3.3022090632217606E-14, 7.2423757177046557,
    -7.592500495497425E-16, 3.2253274101098324E-14, 7.0091201098594791,
    -7.627496322892558E-16, -2.6656662618040004E-17, -7.5440609637554808E-16,
    8.4631379685742587, -2.6040415508152005E-17, -7.53432806687164E-16,
    8.458947735182857, -2.5424168398264003E-17, -7.5245951699878016E-16,
    8.4547575017229679, -2.4807561330952008E-17, -7.5148565880006417E-16,
    8.4505648214674611, -2.4191314221064002E-17, -7.5051236911168011E-16,
    8.4463745881103076, -2.3575067111176E-17, -7.4953907942329605E-16,
    8.4421843546846613, -2.2958820001287998E-17, -7.4856578973491209E-16,
    8.4379941211905258, -2.2342212933976003E-17, -7.4759193153619611E-16,
    8.4338014409007762, -2.1725965824088003E-17, -7.46618641847812E-16,
    8.4296112075162242, -2.11097187142E-17, -7.4564535215942818E-16,
    8.4254209740734574, -2.0493111646888003E-17, -7.446714939607122E-16,
    8.4212282938367871, -1.9876864537E-17, -7.43698204272328E-16,
    8.4170380605087409, -1.9260617427112005E-17, -7.4272491458394408E-16,
    8.4128478271310385, -1.8644370317224003E-17, -7.4175162489556012E-16,
    208.40865759369393, -1.8027907232881604E-17, -7.40777994100977E-16,
    8.17540029881121, -1.7411588131508802E-17, -7.398045907105265E-16,
    7.9421440413329165, -1.67952330343936E-17, -7.38831130469043E-16,
    7.7088876564419255, -1.6178913933020802E-17, -7.3785772707859251E-16,
    7.4756316921436365, -1.5562558835905605E-17, -7.3688426683710882E-16,
    7.2423757177046557, -1.4946214537513123E-17, -7.3591082365093535E-16,
    7.0091201098594791, 8.2283468781897078, 4.9433766605153528E-14,
    -3.9269914873714405E-17, 8.22427461321002, 4.8408024426430969E-14,
    -3.7636134511564006E-17, 8.220202348178967, 4.73822822477084E-14,
    -3.60023541494136E-17, 8.2161277050364578, 4.6355940920516965E-14,
    -3.4367619476304E-17, 8.2120554400910155, 4.53301987417944E-14,
    -3.2733839114153604E-17, 8.2079831750942063, 4.4304456563071847E-14,
    -3.1100058752003206E-17, 8.2039109100289078, 4.3278714384349281E-14,
    -2.94662783898528E-17, 8.19983626685215, 4.2252373057157849E-14,
    -2.7831543716743208E-17, 8.1957640018724653, 4.122663087843529E-14,
    -2.6197763354592803E-17, 8.1916917368482611, 4.0200888699712731E-14,
    -2.45639829924424E-17, 8.1876170937228743, 3.9174547372521287E-14,
    -2.29292483193328E-17, 8.1835448287962667, 3.8148805193798727E-14,
    -2.12954679571824E-17, 8.1794725638285666, 3.7123063015076168E-14,
    -1.9661687595032002E-17, 8.17540029881121, 3.6097320836353609E-14,
    -1.80279072328816E-17, 208.17132660690865, 3.5071219168549716E-14,
    -1.6393554284155682E-17, 7.9381882729475333, 3.4045357160133381E-14,
    -1.475958305981344E-17, 7.7050497598066379, 3.3019435236870161E-14,
    -1.3125516404375281E-17, 7.4719116017523151, 3.199357322845382E-14,
    -1.1491545180033043E-17, 7.238773381791427, 3.0967651305190594E-14,
    -9.8574785245948825E-18, 7.0056354677261057, 2.9941747356381442E-14,
    -8.223440498485498E-18, 4.5423982897337534E-14, 8.2283468781897078,
    -6.9482246768300367E-16, 4.4684610827734963E-14, 8.22427461321002,
    -6.9835703173473877E-16, 4.3945238758132404E-14, 8.220202348178967,
    -7.0189159578647368E-16, 4.3205434812320967E-14, 8.2161277050364578,
    -7.054282244200146E-16, 4.2466062742718408E-14, 8.2120554400910155,
    -7.0896278847174971E-16, 4.1726690673115843E-14, 8.2079831750942063,
    -7.1249735252348471E-16, 4.0987318603513284E-14, 8.2039109100289078,
    -7.1603191657521981E-16, 4.0247514657701841E-14, 8.19983626685215,
    -7.1956854520876074E-16, 3.9508142588099289E-14, 8.1957640018724653,
    -7.2310310926049574E-16, 3.876877051849673E-14, 8.1916917368482611,
    -7.2663767331223094E-16, 3.8028966572685281E-14, 8.1876170937228743,
    -7.3017430194577187E-16, 3.7289594503082722E-14, 8.1835448287962667,
    -7.3370886599750678E-16, 3.6550222433480163E-14, 8.1794725638285666,
    -7.3724343004924188E-16, 3.5810850363877604E-14, 8.17540029881121,
    -7.4077799410097688E-16, 3.5071219168549722E-14, 208.17132660690865,
    -7.443137969017954E-16, 3.4331760723705383E-14, 7.9381882729475333,
    -7.4784877386989171E-16, 3.3592259091240159E-14, 7.7050497598066379,
    -7.5138395729616857E-16, 3.285280064639582E-14, 7.4719116017523151,
    -7.5491893426426478E-16, 3.2113299013930596E-14, 7.238773381791427,
    -7.5845411769054164E-16, 3.1373810337751641E-14, 7.0056354677261057,
    -7.6198923917936419E-16, -2.4329814809199204E-17, -7.5744436042199724E-16,
    8.2283468781897078, -2.3763024189512807E-17, -7.565066038789802E-16,
    8.22427461321002, -2.3196233569826403E-17, -7.5556884733596336E-16,
    8.220202348178967, -2.2629111880852806E-17, -7.546305430379563E-16,
    8.2161277050364578, -2.2062321261166402E-17, -7.5369278649493936E-16,
    8.2120554400910155, -2.1495530641480002E-17, -7.5275502995192232E-16,
    8.2079831750942063, -2.0928740021793602E-17, -7.5181727340890537E-16,
    8.2039109100289078, -2.0361618332820005E-17, -7.5087896911089832E-16,
    8.19983626685215, -1.9794827713133604E-17, -7.4994121256788128E-16,
    8.1957640018724653, -1.92280370934472E-17, -7.4900345602486453E-16,
    8.1916917368482611, -1.86609154044736E-17, -7.4806515172685747E-16,
    8.1876170937228743, -1.80941247847872E-17, -7.4712739518384043E-16,
    8.1835448287962667, -1.7527334165100803E-17, -7.4618963864082349E-16,
    8.1794725638285666, -1.69605435454144E-17, -7.4525188209780655E-16,
    8.17540029881121, -1.6393554284155682E-17, -7.443137969017955E-16,
    208.17132660690865, -1.582669745061184E-17, -7.4337593080778046E-16,
    7.9381882729475333, -1.5259807510139283E-17, -7.4243800993826666E-16,
    7.7050497598066379, -1.4692950676595441E-17, -7.4150014384425161E-16,
    7.4719116017523151, -1.4126060736122884E-17, -7.4056222297473752E-16,
    7.238773381791427, -1.3559180727728938E-17, -7.3962431853787332E-16,
    7.0056354677261057, 7.9935579447306031, 4.8283817172100966E-14,
    -3.80107249212352E-17, 7.9896035620420474, 4.7287391374922089E-14,
    -3.6426389498080005E-17, 7.9856491793192479, 4.62909655777432E-14,
    -3.48420540749248E-17, 7.9816924871702239, 4.5293957756150084E-14,
    -3.3256793222200003E-17, 7.9777381045159146, 4.42975319589712E-14,
    -3.16724577990448E-17, 7.973783721827358, 4.3301106161792322E-14,
    -3.00881223758896E-17, 7.9698293390874362, 4.2304680364613439E-14,
    -2.8503786952734406E-17, 7.9658726469041667, 4.130767254302033E-14,
    -2.6918526100009604E-17, 7.9619182642156119, 4.0311246745841446E-14,
    -2.5334190676854404E-17, 7.9579638814928124, 3.9314820948662568E-14,
    -2.37498552536992E-17, 7.9540071893506372, 3.8317813127069447E-14,
    -2.2164594400974402E-17, 7.950052806713451, 3.7321387329890563E-14,
    -2.0580258977819202E-17, 7.94609842404373, 3.6324961532711686E-14,
    -1.8995923554664002E-17, 7.9421440413329156, 3.5328535735532808E-14,
    -1.7411588131508802E-17, 7.9381882729475315, 3.4331760723705383E-14,
    -1.5826697450611843E-17, 207.93423342836567, 3.3335218521643654E-14,
    -1.424217694154272E-17, 7.70121275949525, 3.2338618117140506E-14,
    -1.2657563889516642E-17, 7.4681923802056644, 3.1342075915078771E-14,
    -1.1073043380447522E-17, 7.2351718872420347, 3.0345475510575623E-14,
    -9.4884303284214412E-18, 7.0021516394762147, 2.93488925668049E-14,
    -7.90384503928245E-18, 4.3987629900712971E-14, 7.9935579447306031,
    -6.93387179041265E-16, 4.3277574212654081E-14, 7.9896035620420474,
    -6.96957267939139E-16, 4.2567518524595204E-14, 7.9856491793192479,
    -7.0052735683701285E-16, 4.1857048084382084E-14, 7.9816924871702239,
    -7.0409953106718694E-16, 4.1146992396323207E-14, 7.9777381045159146,
    -7.0766961996506095E-16, 4.0436936708264323E-14, 7.973783721827358,
    -7.1123970886293476E-16, 3.9726881020205446E-14, 7.9698293390874362,
    -7.1480979776080877E-16, 3.901641057999232E-14, 7.9658726469041667,
    -7.1838197199098286E-16, 3.8306354891933449E-14, 7.9619182642156119,
    -7.2195206088885677E-16, 3.7596299203874572E-14, 7.9579638814928124,
    -7.2552214978673078E-16, 3.6885828763661445E-14, 7.9540071893506372,
    -7.2909432401690477E-16, 3.6175773075602562E-14, 7.950052806713451,
    -7.3266441291477858E-16, 3.5465717387543685E-14, 7.94609842404373,
    -7.3623450181265259E-16, 3.4755661699484808E-14, 7.9421440413329156,
    -7.398045907105265E-16, 3.4045357160133381E-14, 7.9381882729475315,
    -7.4337593080778046E-16, 3.3335218521643654E-14, 207.93423342836567,
    -7.4694643677211448E-16, 3.26250384079385E-14, 7.70121275949525,
    -7.5051715126967847E-16, 3.1914899769448775E-14, 7.4681923802056644,
    -7.5408765723401239E-16, 3.1204719655743617E-14, 7.2351718872420347,
    -7.5765837173157638E-16, 3.0494551984603094E-14, 7.0021516394762147,
    -7.6122902366917136E-16, -2.2003510465921604E-17, -7.6048191484236985E-16,
    7.9935579447306031, -2.1486164785230405E-17, -7.5957968314549171E-16,
    7.9896035620420474, -2.0968819104539204E-17, -7.5867745144861378E-16,
    7.9856491793192479, -2.0451171235950406E-17, -7.5777469274723978E-16,
    7.9816924871702239, -1.99338255552592E-17, -7.5687246105036175E-16,
    7.9777381045159146, -1.9416479874568E-17, -7.5597022935348362E-16,
    7.973783721827358, -1.88991341938768E-17, -7.5506799765660558E-16,
    7.9698293390874362, -1.8381486325288003E-17, -7.5416523895523159E-16,
    7.9658726469041667, -1.7864140644596802E-17, -7.5326300725835345E-16,
    7.9619182642156119, -1.73467949639056E-17, -7.5236077556147552E-16,
    7.9579638814928124, -1.6829147095316802E-17, -7.5145801686010162E-16,
    7.9540071893506372, -1.63118014146256E-17, -7.5055578516322349E-16,
    7.950052806713451, -1.5794455733934402E-17, -7.4965355346634535E-16,
    7.94609842404373, -1.52771100532432E-17, -7.4875132176946732E-16,
    7.9421440413329156, -1.475958305981344E-17, -7.4784877386989181E-16,
    7.9381882729475315, -1.4242176941542723E-17, -7.4694643677211448E-16,
    207.93423342836567, -1.372474060448224E-17, -7.4604404697388771E-16,
    7.70121275949525, -1.3207334486211522E-17, -7.4514170987611038E-16,
    7.4681923802056644, -1.2689898149151042E-17, -7.4423932007788351E-16,
    7.2351718872420347, -1.217247087772749E-17, -7.4333694608979164E-16,
    7.0021516394762147, 7.7587695568412611, 4.7133800576950924E-14,
    -3.6751461426532404E-17, 7.7549330046841183, 4.616669287352117E-14,
    -3.5216573830174007E-17, 7.75109645249273, 4.5199585170091404E-14,
    -3.36816862338156E-17, 7.7472576595643634, 4.4231912567302164E-14,
    -3.2145902090964004E-17, 7.7434211074414661, 4.3264804863872404E-14,
    -3.06110144946056E-17, 7.7395845552843223, 4.2297697160442644E-14,
    -2.9076126898247204E-17, 7.735748003092934, 4.1330589457012884E-14,
    -2.7541239301888804E-17, 7.7319092101474443, 4.0362916854223651E-14,
    -2.6005455159037207E-17, 7.7280726579903023, 3.9395809150793891E-14,
    -2.4470567562678804E-17, 7.724236105798914, 3.8428701447364131E-14,
    -2.29356799663204E-17, 7.7203973128705465, 3.7461028844574891E-14,
    -2.13998958234688E-17, 7.7165607607544988, 3.6493921141145125E-14,
    -1.98650082271104E-17, 7.7127242086144792, 3.5526813437715371E-14,
    -1.8330120630752E-17, 7.7088876564419255, 3.4559705734285604E-14,
    -1.67952330343936E-17, 7.705049759806637, 3.3592259091240159E-14,
    -1.5259807510139283E-17, 7.7012127594952506, 3.26250384079385E-14,
    -1.372474060448224E-17, 207.69737553510555, 3.16577612347009E-14,
    -1.2189584044175883E-17, 7.4644729414504143, 3.0690540551399242E-14,
    -1.0654517138518842E-17, 7.231570182352967, 2.972326337816164E-14,
    -9.1193605782124829E-18, 6.9986676077554533, 2.8756003151904819E-14,
    -7.5842309143009172E-18, 4.2551193014764932E-14, 7.7587695568412611,
    -6.91951806572366E-16, 4.1870455420455166E-14, 7.7549330046841183,
    -6.95557422391186E-16, 4.1189717826145406E-14, 7.75109645249273,
    -6.9916303821000563E-16, 4.0508582604736168E-14, 7.7472576595643634,
    -7.027707601128318E-16, 3.9827845010426408E-14, 7.7434211074414661,
    -7.0637637593165174E-16, 3.9147107416116642E-14, 7.7395845552843223,
    -7.0998199175047149E-16, 3.8466369821806882E-14, 7.735748003092934,
    -7.1358760756929143E-16, 3.7785234600397645E-14, 7.7319092101474443,
    -7.171953294721176E-16, 3.7104497006087885E-14, 7.7280726579903023,
    -7.2080094529093735E-16, 3.6423759411778131E-14, 7.724236105798914,
    -7.2440656110975729E-16, 3.5742624190368881E-14, 7.7203973128705465,
    -7.2801428301258346E-16, 3.5061886596059121E-14, 7.7165607607544988,
    -7.3161989883140311E-16, 3.4381149001749361E-14, 7.7127242086144792,
    -7.3522551465022306E-16, 3.3700411407439608E-14, 7.7088876564419255,
    -7.388311304690429E-16, 3.3019435236870161E-14, 7.705049759806637,
    -7.4243800993826646E-16, 3.2338618117140506E-14, 7.7012127594952506,
    -7.4604404697388771E-16, 3.16577612347009E-14, 207.69737553510555,
    -7.4965029461790941E-16, 3.0976944114971244E-14, 7.4644729414504143,
    -7.5325633165353046E-16, 3.0296087232531638E-14, 7.231570182352967,
    -7.5686257929755226E-16, 2.9615242278905019E-14, 6.9986676077554533,
    -7.6046876375905379E-16, -1.9677070256253204E-17, -7.6351964666926162E-16,
    7.7587695568412611, -1.9209172402358806E-17, -7.6265294189332932E-16,
    7.7549330046841183, -1.8741274548464405E-17, -7.6178623711739732E-16,
    7.75109645249273, -1.8273103389748805E-17, -7.609190260886755E-16,
    7.7472576595643634, -1.7805205535854404E-17, -7.6005232131274331E-16,
    7.7434211074414661, -1.7337307681960003E-17, -7.5918561653681111E-16,
    7.7395845552843223, -1.68694098280656E-17, -7.58318911760879E-16,
    7.735748003092934, -1.6401238669350002E-17, -7.5745170073215719E-16,
    7.7319092101474443, -1.5933340815455604E-17, -7.5658499595622489E-16,
    7.7280726579903023, -1.54654429615612E-17, -7.5571829118029289E-16,
    7.724236105798914, -1.49972718028456E-17, -7.5485108015157108E-16,
    7.7203973128705465, -1.4529373948951202E-17, -7.5398437537563878E-16,
    7.7165607607544988, -1.4061476095056803E-17, -7.5311767059970668E-16,
    7.7127242086144792, -1.3593578241162402E-17, -7.5225096582377458E-16,
    7.7088876564419255, -1.3125516404375281E-17, -7.5138395729616857E-16,
    7.705049759806637, -1.2657563889516642E-17, -7.5051715126967837E-16,
    7.7012127594952506, -1.2189584044175881E-17, -7.4965029461790951E-16,
    207.69737553510555, -1.1721631529317242E-17, -7.4878348859141931E-16,
    7.4644729414504143, -1.1253651683976483E-17, -7.4791663193965016E-16,
    7.231570182352967, -1.0785680037780358E-17, -7.4704979047546493E-16,
    6.9986676077554533, 7.5239824411033318, 4.5983851143898369E-14,
    -3.5492271474053207E-17, 7.5202636539824406, 4.504605982201229E-14,
    -3.4006828816690006E-17, 7.5165448668273047, 4.41082685001262E-14,
    -3.25213861593268E-17, 7.5128239076032974, 4.3169929402935283E-14,
    -3.1035075836860004E-17, 7.5091051204995294, 4.22321380810492E-14,
    -2.95496331794968E-17, 7.5053863333786373, 4.1294346759163126E-14,
    -2.80641905221336E-17, 7.5016675462235014, 4.0356555437277041E-14,
    -2.6578747864770404E-17, 7.497946586999495, 3.9418216340086131E-14,
    -2.5092437542303606E-17, 7.4942277998786038, 3.8480425018200046E-14,
    -2.3606994884940402E-17, 7.4905090127234679, 3.7542633696313974E-14,
    -2.2121552227577203E-17, 7.4867880534994615, 3.6604294599123045E-14,
    -2.0635241905110403E-17, 7.4830692663956935, 3.5666503277236967E-14,
    -1.9149799247747202E-17, 7.479350479281651, 3.4728711955350888E-14,
    -1.7664356590384E-17, 7.4756316921436365, 3.3790920633464804E-14,
    -1.6178913933020802E-17, 7.4719116017523151, 3.285280064639582E-14,
    -1.4692950676595441E-17, 7.4681923802056653, 3.1914899769448775E-14,
    -1.320733448621152E-17, 7.4644729414504143, 3.0976944114971244E-14,
    -1.1721631529317242E-17, 207.46075371990256, 3.0039043238024193E-14,
    -1.0236015338933323E-17, 7.2279686878035747, 2.9101087583546662E-14,
    -8.7503123820390416E-18, 6.9951837795055631, 2.8163148362328274E-14,
    -7.2646354550978688E-18, 4.1114840018140375E-14, 7.5239824411033318,
    -6.905165179306275E-16, 4.0463418805374284E-14, 7.5202636539824406,
    -6.9415765859558625E-16, 3.9811997592608206E-14, 7.5165448668273047,
    -6.977987992605449E-16, 3.9160195876797285E-14, 7.5128239076032974,
    -7.0144206676000413E-16, 3.8508774664031213E-14, 7.5091051204995294,
    -7.0508320742496289E-16, 3.7857353451265122E-14, 7.5053863333786373,
    -7.0872434808992164E-16, 3.7205932238499044E-14, 7.5016675462235014,
    -7.1236548875488039E-16, 3.6554130522688123E-14, 7.497946586999495,
    -7.1600875625433962E-16, 3.5902709309922045E-14, 7.4942277998786038,
    -7.1964989691929837E-16, 3.5251288097155967E-14, 7.4905090127234679,
    -7.2329103758425712E-16, 3.4599486381345046E-14, 7.4867880534994615,
    -7.2693430508371636E-16, 3.3948065168578961E-14, 7.4830692663956935,
    -7.30575445748675E-16, 3.3296643955812883E-14, 7.479350479281651,
    -7.3421658641363376E-16, 3.2645222743046805E-14, 7.4756316921436365,
    -7.3785772707859251E-16, 3.199357322845382E-14, 7.4719116017523151,
    -7.4150014384425152E-16, 3.1342075915078777E-14, 7.4681923802056653,
    -7.4514170987611038E-16, 3.0690540551399242E-14, 7.4644729414504143,
    -7.4878348859141931E-16, 3.0039043238024193E-14, 207.46075371990256,
    -7.5242505462327808E-16, 2.9387507874344658E-14, 7.2279686878035747,
    -7.56066833338587E-16, 2.8735983925756472E-14, 6.9951837795055631,
    -7.5970854824886086E-16, -1.7350765912975604E-17, -7.6655720108963423E-16,
    7.5239824411033318, -1.6932312998076405E-17, -7.6572602115984094E-16,
    7.5202636539824406, -1.6513860083177205E-17, -7.6489484123004774E-16,
    7.5165448668273047, -1.6095162744846405E-17, -7.64063175797959E-16,
    7.5128239076032974, -1.5676709829947203E-17, -7.632319958681658E-16,
    7.5091051204995294, -1.5258256915048003E-17, -7.6240081593837241E-16,
    7.5053863333786373, -1.48398040001488E-17, -7.6156963600857911E-16,
    7.5016675462235014, -1.4421106661818004E-17, -7.6073797057649046E-16,
    7.497946586999495, -1.4002653746918804E-17, -7.5990679064669707E-16,
    7.4942277998786038, -1.3584200832019601E-17, -7.5907561071690388E-16,
    7.4905090127234679, -1.3165503493688802E-17, -7.5824394528481522E-16,
    7.4867880534994615, -1.2747050578789602E-17, -7.5741276535502183E-16,
    7.4830692663956935, -1.2328597663890403E-17, -7.5658158542522864E-16,
    7.479350479281651, -1.1910144748991203E-17, -7.5575040549543534E-16,
    7.4756316921436365, -1.1491545180033043E-17, -7.5491893426426478E-16,
    7.4719116017523151, -1.1073043380447522E-17, -7.5408765723401239E-16,
    7.4681923802056653, -1.0654517138518842E-17, -7.5325633165353056E-16,
    7.4644729414504143, -1.0236015338933321E-17, -7.5242505462327808E-16,
    207.46075371990256, -9.8174890970046427E-18, -7.5159372904279615E-16,
    7.2279686878035747, -9.39897018777891E-18, -7.5076241802738315E-16,
    6.9951837795055631, 7.289195988207168, 4.4833834548748326E-14,
    -3.4233007979350409E-17, 7.2855949143456868, 4.3925361320611371E-14,
    -3.2797013148784008E-17, 7.2819938404670825, 4.30168880924744E-14,
    -3.13610183182176E-17, 7.27839066319173, 4.2107884214087363E-14,
    -2.9924184705624005E-17, 7.27478958933025, 4.11994109859504E-14,
    -2.8488189875057604E-17, 7.2711885154687668, 4.0290937757813441E-14,
    -2.7052195044491203E-17, 7.2675874415901633, 3.938246452967648E-14,
    -2.5616200213924802E-17, 7.2639842643148107, 3.8473460651289452E-14,
    -2.4179366601331206E-17, 7.26038319045333, 3.7564987423152485E-14,
    -2.2743371770764802E-17, 7.2567821165747262, 3.665651419501553E-14,
    -2.1307376940198404E-17, 7.2531789392993735, 3.574751031662849E-14,
    -1.9870543327604802E-17, 7.2495778654378924, 3.4839037088491522E-14,
    -1.84345484970384E-17, 7.24597679157641, 3.3930563860354568E-14,
    -1.6998553666472003E-17, 7.2423757177046557, 3.3022090632217606E-14,
    -1.5562558835905602E-17, 7.238773381791427, 3.2113299013930596E-14,
    -1.4126060736122884E-17, 7.2351718872420347, 3.1204719655743623E-14,
    -1.2689898149151041E-17, 7.231570182352967, 3.0296087232531638E-14,
    -1.1253651683976483E-17, 7.2279686878035747, 2.9387507874344658E-14,
    -9.8174890970046427E-18, 207.2243669829133, 2.8478875451132679E-14,
    -8.3812426318300817E-18, 6.9916997477848017, 2.7570258947428195E-14,
    -6.945021330116338E-18, 3.9678403132192329E-14, 7.289195988207168,
    -6.8908114546172841E-16, 3.9056300013175362E-14, 7.2855949143456868,
    -6.927578130476331E-16, 3.8434196894158408E-14, 7.2819938404670825,
    -6.9643448063353768E-16, 3.7811730397151363E-14, 7.27839066319173,
    -7.00113295805649E-16, 3.7189627278134408E-14, 7.27478958933025,
    -7.0378996339155368E-16, 3.6567524159117441E-14, 7.2711885154687668,
    -7.0746663097745827E-16, 3.5945421040100487E-14, 7.2675874415901633,
    -7.1114329856336306E-16, 3.5322954543093442E-14, 7.2639842643148107,
    -7.1482211373547437E-16, 3.4700851424076487E-14, 7.26038319045333,
    -7.18498781321379E-16, 3.4078748305059533E-14, 7.2567821165747262,
    -7.2217544890728374E-16, 3.3456281808052481E-14, 7.2531789392993735,
    -7.2585426407939506E-16, 3.2834178689035521E-14, 7.2495778654378924,
    -7.2953093166529954E-16, 3.2212075570018566E-14, 7.24597679157641,
    -7.3320759925120423E-16, 3.1589972451001605E-14, 7.2423757177046557,
    -7.3688426683710892E-16, 3.09676513051906E-14, 7.238773381791427,
    -7.4056222297473752E-16, 3.0345475510575623E-14, 7.2351718872420347,
    -7.4423932007788361E-16, 2.972326337816164E-14, 7.231570182352967,
    -7.4791663193965026E-16, 2.9101087583546662E-14, 7.2279686878035747,
    -7.5159372904279615E-16, 2.8478875451132679E-14, 207.2243669829133,
    -7.5527104090456289E-16, 2.7856674220058395E-14, 6.9916997477848017,
    -7.5894828833874329E-16, -1.5024325703307205E-17, -7.69594932916526E-16,
    7.289195988207168, -1.4655320615204805E-17, -7.6879927990767854E-16,
    7.2855949143456868, -1.4286315527102403E-17, -7.6800362689883138E-16,
    7.2819938404670825, -1.3917094898644805E-17, -7.6720750913939471E-16,
    7.27839066319173, -1.3548089810542402E-17, -7.6641185613054735E-16,
    7.27478958933025, -1.3179084722440001E-17, -7.656162031216999E-16,
    7.2711885154687668, -1.2810079634337602E-17, -7.6482055011285264E-16,
    7.2675874415901633, -1.2440859005880003E-17, -7.64024432353416E-16,
    7.2639842643148107, -1.2071853917777603E-17, -7.6322877934456851E-16,
    7.26038319045333, -1.17028488296752E-17, -7.6243312633572125E-16,
    7.2567821165747262, -1.1333628201217601E-17, -7.6163700857628468E-16,
    7.2531789392993735, -1.0964623113115202E-17, -7.6084135556743722E-16,
    7.2495778654378924, -1.0595618025012802E-17, -7.6004570255858986E-16,
    7.24597679157641, -1.0226612936910401E-17, -7.592500495497425E-16,
    7.2423757177046557, -9.857478524594881E-18, -7.5845411769054164E-16,
    7.238773381791427, -9.4884303284214412E-18, -7.5765837173157638E-16,
    7.2351718872420347, -9.1193605782124813E-18, -7.5686257929755226E-16,
    7.231570182352967, -8.7503123820390416E-18, -7.56066833338587E-16,
    7.2279686878035747, -8.3812426318300817E-18, -7.5527104090456279E-16,
    207.2243669829133, -8.0121793478317786E-18, -7.5447526241305643E-16,
    6.9916997477848017, 7.05441069099154, 4.3683838102227532E-14,
    -3.2973766547314689E-17, 7.0509272696833909, 4.2804682454178059E-14,
    -3.1587218677204607E-17, 7.0474438483752424, 4.192552680612858E-14,
    -3.0200670807094518E-17, 7.0439583923583857, 4.1045857632583755E-14,
    -2.8813313037527604E-17, 7.0404749710502372, 4.0166701984534282E-14,
    -2.7426765167417521E-17, 7.036991549742087, 3.928754633648481E-14,
    -2.6040217297307442E-17, 7.0335081284339385, 3.8408390688435337E-14,
    -2.4653669427197363E-17, 7.0300226724170809, 3.7528721514890518E-14,
    -2.3266311657630445E-17, 7.0265392511089324, 3.6649565866841046E-14,
    -2.1879763787520363E-17, 7.0230558298007839, 3.5770410218791573E-14,
    -2.0493215917410281E-17, 7.0195703737839272, 3.4890741045246742E-14,
    -1.910585814784336E-17, 7.0160869524757787, 3.4011585397197269E-14,
    -1.7719310277733281E-17, 7.0126035311676294, 3.31324297491478E-14,
    -1.6332762407623202E-17, 7.00912010985948, 3.2253274101098324E-14,
    -1.4946214537513123E-17, 7.0056354677261057, 3.1373810337751641E-14,
    -1.3559180727728938E-17, 7.0021516394762164, 3.0494551984603094E-14,
    -1.217247087772749E-17, 6.9986676077554542, 2.9615242278905019E-14,
    -1.0785680037780358E-17, 6.9951837795055631, 2.8735983925756472E-14,
    -9.39897018777891E-18, 6.9916997477848026, 2.7856674220058395E-14,
    -8.0121793478317786E-18, 206.9882157771053, 2.6977379920125177E-14,
    -6.625412804868351E-18, 3.8241991413041339E-14, 7.05441069099154,
    -6.876457981409775E-16, 3.7649205874111854E-14, 7.0509272696833909,
    -6.91357992025386E-16, 3.7056420335182388E-14, 7.0474438483752424,
    -6.950701859097944E-16, 3.6463288543017558E-14, 7.0439583923583857,
    -6.9878454813175211E-16, 3.5870503004088092E-14, 7.0404749710502372,
    -7.0249674201616062E-16, 3.5277717465158607E-14, 7.036991549742087,
    -7.06208935900569E-16, 3.4684931926229141E-14, 7.0335081284339385,
    -7.0992112978497752E-16, 3.4091800134064311E-14, 7.0300226724170809,
    -7.1363549200693533E-16, 3.3499014595134845E-14, 7.0265392511089324,
    -7.1734768589134373E-16, 3.2906229056205373E-14, 7.0230558298007839,
    -7.2105987977575223E-16, 3.2313097264040537E-14, 7.0195703737839272,
    -7.2477424199770994E-16, 3.1720311725111065E-14, 7.0160869524757787,
    -7.2848643588211825E-16, 3.1127526186181592E-14, 7.0126035311676294,
    -7.3219862976652685E-16, 3.0534740647252126E-14, 7.00912010985948,
    -7.3591082365093525E-16, 2.9941747356381442E-14, 7.0056354677261057,
    -7.3962431853787322E-16, 2.93488925668049E-14, 7.0021516394762164,
    -7.4333694608979164E-16, 2.8756003151904819E-14, 6.9986676077554542,
    -7.4704979047546483E-16, 2.8163148362328274E-14, 6.9951837795055631,
    -7.5076241802738315E-16, 2.7570258947428195E-14, 6.9916997477848026,
    -7.5447526241305643E-16, 2.6977379920125177E-14, 206.9882157771053,
    -7.5818804174860314E-16, -1.2697926253556044E-17, -7.72632611521462E-16,
    7.05441069099154, -1.2378368125909965E-17, -7.718724848111184E-16,
    7.0509272696833909, -1.2058809998263884E-17, -7.71112358100775E-16,
    7.0474438483752424, -1.1739065212832965E-17, -7.7035178739118468E-16,
    7.0439583923583857, -1.1419507085186881E-17, -7.6959166068084124E-16,
    7.0404749710502372, -1.1099948957540802E-17, -7.688315339704976E-16,
    7.036991549742087, -1.0780390829894721E-17, -7.6807140726015406E-16,
    7.0335081284339385, -1.0460646044463803E-17, -7.6731083655056378E-16,
    7.0300226724170809, -1.0141087916817723E-17, -7.6655070984022014E-16,
    7.0265392511089324, -9.82152978917164E-18, -7.657905831298768E-16,
    7.0230558298007839, -9.5017850037407208E-18, -7.6503001242028652E-16,
    7.0195703737839272, -9.1822268760946417E-18, -7.6426988570994288E-16,
    7.0160869524757787, -8.862668748448561E-18, -7.6350975899959934E-16,
    7.0126035311676294, -8.5431106208024819E-18, -7.627496322892558E-16,
    7.00912010985948, -8.223440498485498E-18, -7.6198923917936429E-16,
    7.0056354677261057, -7.90384503928245E-18, -7.6122902366917126E-16,
    7.0021516394762164, -7.5842309143009172E-18, -7.6046876375905389E-16,
    6.9986676077554542, -7.2646354550978688E-18, -7.5970854824886086E-16,
    6.9951837795055631, -6.945021330116338E-18, -7.5894828833874319E-16,
    6.9916997477848026, -6.6254128048683518E-18, -7.5818804174860324E-16,
    206.9882157771053 };

  
   double b_a[3600] = { -211.5154541352035, -6.5533595964669213E-14,
    5.6899162546212E-17, -11.280659712816917, -6.40974107466916E-14,
    5.4573129935716004E-17, -11.045865348887073, -6.2661225528714E-14,
    5.224709732522E-17, -10.811067762614019, -6.12242014175016E-14,
    4.9919706050816008E-17, -10.576273692184389, -5.9788016199524E-14,
    4.7593673440320003E-17, -10.341479856129229, -5.83518309815464E-14,
    4.5267640829824E-17, -10.106686313076017, -5.6915645763568807E-14,
    4.2941608219328E-17, -9.8718897822287275, -5.5478621652356411E-14,
    4.0614216944924E-17, -9.6370970017757, -5.404243643437881E-14,
    3.8288184334428005E-17, -9.4023046902394611, -5.2606251216401209E-14,
    3.5962151723932006E-17, -9.1675095668102937, -5.1169227105188813E-14,
    3.3634760449527996E-17, -8.9327183696834176, -4.9733041887211206E-14,
    3.1308727839032004E-17, -8.6979278173984422, -4.8296856669233612E-14,
    2.8982695228536E-17, -8.4631379685742569, -4.6860671451256005E-14,
    2.665666261804E-17, -8.2283468781897078, -4.5423982897337528E-14,
    2.43298148091992E-17, -7.993557944730604, -4.3987629900712964E-14,
    2.20035104659216E-17, -7.7587695568412611, -4.2551193014764926E-14,
    1.9677070256253204E-17, -7.5239824411033318, -4.1114840018140368E-14,
    1.73507659129756E-17, -7.289195988207168, -3.9678403132192323E-14,
    1.5024325703307202E-17, -7.0544106909915394, -3.8241991413041333E-14,
    1.2697926253556042E-17, -6.5533595964669213E-14, -211.5154541352035,
    7.1491717928462731E-16, -6.4383780855811606E-14, -11.280659712816917,
    7.1795437889196179E-16, -6.3233965746954E-14, -11.045865348887073,
    7.2099157849929608E-16, -6.2083479017121611E-14, -10.811067762614019,
    7.2403055217182159E-16, -6.0933663908264E-14, -10.576273692184389,
    7.2706775177915607E-16, -5.97838487994064E-14, -10.341479856129229,
    7.3010495138649046E-16, -5.86340336905488E-14, -10.106686313076017,
    7.3314215099382494E-16, -5.74835469607164E-14, -9.8718897822287275,
    7.3618112466635055E-16, -5.6333731851858809E-14, -9.6370970017757,
    7.3921832427368493E-16, -5.5183916743001208E-14, -9.4023046902394611,
    7.4225552388101942E-16, -5.40334300131688E-14, -9.1675095668102937,
    7.4529449755354493E-16, -5.2883614904311194E-14, -8.9327183696834176,
    7.4833169716087921E-16, -5.17337997954536E-14, -8.6979278173984422,
    7.5136889676821369E-16, -5.0583984686596006E-14, -8.4631379685742569,
    7.5440609637554808E-16, -4.9433766605153528E-14, -8.2283468781897078,
    7.5744436042199714E-16, -4.8283817172100973E-14, -7.993557944730604,
    7.6048191484236985E-16, -4.7133800576950924E-14, -7.7587695568412611,
    7.6351964666926162E-16, -4.5983851143898369E-14, -7.5239824411033318,
    7.6655720108963423E-16, -4.4833834548748326E-14, -7.289195988207168,
    7.695949329165261E-16, -4.3683838102227525E-14, -7.0544106909915394,
    7.72632611521462E-16, 5.6899162546212E-17, 7.1491717928462731E-16,
    -211.5154541352035, 5.5640119678180009E-17, 7.1348205829720955E-16,
    -11.280659712816917, 5.4381076810148E-17, 7.1204693730979209E-16,
    -11.045865348887073, 5.3121298519880006E-17, 7.1061097805076971E-16,
    -10.811067762614019, 5.1862255651848E-17, 7.0917585706335215E-16,
    -10.576273692184389, 5.0603212783816006E-17, 7.0774073607593449E-16,
    -10.341479856129229, 4.9344169915784E-17, 7.0630561508851693E-16,
    -10.106686313076017, 4.8084391625516004E-17, 7.0486965582949455E-16,
    -9.8718897822287275, 4.682534875748401E-17, 7.034345348420768E-16,
    -9.6370970017757, 4.5566305889452004E-17, 7.0199941385465934E-16,
    -9.4023046902394611, 4.4306527599184E-17, 7.00563454595637E-16,
    -9.1675095668102937, 4.3047484731152E-17, 6.991283336082193E-16,
    -8.9327183696834176, 4.1788441863120008E-17, 6.9769321262080174E-16,
    -8.6979278173984422, 4.0529398995088E-17, 6.9625809163338418E-16,
    -8.4631379685742569, 3.92699148737144E-17, 6.9482246768300367E-16,
    -8.2283468781897078, 3.80107249212352E-17, 6.93387179041265E-16,
    -7.993557944730604, 3.6751461426532404E-17, 6.9195180657236614E-16,
    -7.7587695568412611, 3.5492271474053207E-17, 6.905165179306275E-16,
    -7.5239824411033318, 3.4233007979350409E-17, 6.8908114546172831E-16,
    -7.289195988207168, 3.2973766547314683E-17, 6.876457981409776E-16,
    -7.0544106909915394, -11.280659712816917, -6.4383780855811606E-14,
    5.5640119678180009E-17, -211.27504957325488, -6.297690859496681E-14,
    5.3363526231076009E-17, -11.040373840180749, -6.1570036334122E-14,
    5.1086932783972E-17, -10.805694894794851, -6.0162342302096809E-14,
    4.8809009550976007E-17, -10.571019337930627, -5.8755470041252E-14,
    4.6532416103872008E-17, -10.336343956813446, -5.73485977804072E-14,
    4.4255822656768E-17, -10.101668810070741, -5.59417255195624E-14,
    4.1979229209664008E-17, -9.866990685561472, -5.4534031487537215E-14,
    3.9701305976668006E-17, -9.6323161841230132, -5.3127159226692406E-14,
    3.7424712529564006E-17, -9.3976420929824371, -5.172028696584761E-14,
    3.514811908246E-17, -9.1629651999901434, -5.0312592933822411E-14,
    3.2870195849464004E-17, -8.9282921059698044, -4.8905720672977608E-14,
    3.0593602402360005E-17, -8.6936195981553812, -4.7498848412132812E-14,
    2.8317008955256005E-17, -8.458947735182857, -4.6091976151288E-14,
    2.6040415508152002E-17, -8.22427461321002, -4.4684610827734969E-14,
    2.3763024189512803E-17, -7.9896035620420482, -4.3277574212654087E-14,
    2.1486164785230402E-17, -7.7549330046841192, -4.1870455420455166E-14,
    1.9209172402358806E-17, -7.5202636539824406, -4.0463418805374284E-14,
    1.6932312998076405E-17, -7.2855949143456868, -3.9056300013175362E-14,
    1.4655320615204802E-17, -7.0509272696833909, -3.764920587411186E-14,
    1.2378368125909962E-17, -6.40974107466916E-14, -11.280659712816917,
    7.1348205829720975E-16, -6.297690859496681E-14, -211.27504957325488,
    7.16554778601069E-16, -6.1856406443242014E-14, -11.040373840180749,
    7.19627498904928E-16, -6.0735249792596817E-14, -10.805694894794851,
    7.2270201402204886E-16, -5.9614747640872008E-14, -10.571019337930627,
    7.2577473432590809E-16, -5.84942454891472E-14, -10.336343956813446,
    7.2884745462976733E-16, -5.73737433374224E-14, -10.101668810070741,
    7.3192017493362656E-16, -5.6252586686777206E-14, -9.866990685561472,
    7.3499469005074731E-16, -5.5132084535052411E-14, -9.6323161841230132,
    7.3806741035460654E-16, -5.4011582383327615E-14, -9.3976420929824371,
    7.4114013065846577E-16, -5.2890425732682405E-14, -9.1629651999901434,
    7.4421464577558652E-16, -5.1769923580957596E-14, -8.9282921059698044,
    7.4728736607944565E-16, -5.06494214292328E-14, -8.6936195981553812,
    7.5036008638330488E-16, -4.9528919277508004E-14, -8.458947735182857,
    7.5343280668716412E-16, -4.8408024426430969E-14, -8.22427461321002,
    7.565066038789802E-16, -4.7287391374922089E-14, -7.9896035620420482,
    7.5957968314549181E-16, -4.6166692873521164E-14, -7.7549330046841192,
    7.6265294189332942E-16, -4.504605982201229E-14, -7.5202636539824406,
    7.6572602115984094E-16, -4.3925361320611365E-14, -7.2855949143456868,
    7.6879927990767874E-16, -4.2804682454178059E-14, -7.0509272696833909,
    7.718724848111185E-16, 5.4573129935716004E-17, 7.1795437889196169E-16,
    -11.280659712816917, 5.3363526231076015E-17, 7.1655477860106879E-16,
    -211.27504957325488, 5.2153922526436008E-17, 7.1515517831017617E-16,
    -11.040373840180749, 5.0943612277576019E-17, 7.1375476049574893E-16,
    -10.805694894794851, 4.9734008572936006E-17, 7.1235516020485612E-16,
    -10.571019337930627, 4.8524404868296004E-17, 7.1095555991396331E-16,
    -10.336343956813446, 4.7314801163656003E-17, 7.095559596230705E-16,
    -10.101668810070741, 4.6104490914796008E-17, 7.0815554180864326E-16,
    -9.866990685561472, 4.4894887210156007E-17, 7.0675594151775045E-16,
    -9.6323161841230132, 4.3685283505516006E-17, 7.0535634122685774E-16,
    -9.3976420929824371, 4.2474973256656005E-17, 7.0395592341243049E-16,
    -9.1629651999901434, 4.1265369552016004E-17, 7.0255632312153768E-16,
    -8.9282921059698044, 4.0055765847376009E-17, 7.0115672283064487E-16,
    -8.6936195981553812, 3.8846162142736007E-17, 6.9975712253975216E-16,
    -8.458947735182857, 3.7636134511564006E-17, 6.9835703173473867E-16,
    -8.22427461321002, 3.6426389498080005E-17, 6.9695726793913894E-16,
    -7.9896035620420482, 3.5216573830174007E-17, 6.95557422391186E-16,
    -7.7549330046841192, 3.4006828816690006E-17, 6.9415765859558615E-16,
    -7.5202636539824406, 3.2797013148784008E-17, 6.92757813047633E-16,
    -7.2855949143456868, 3.1587218677204607E-17, 6.91357992025386E-16,
    -7.0509272696833909, -11.045865348887073, -6.3233965746954E-14,
    5.4381076810148E-17, -11.040373840180749, -6.1856406443242014E-14,
    5.2153922526436008E-17, -211.0348823313123, -6.047884713953E-14,
    4.9926768242724E-17, -10.800322026823309, -5.9100483186692E-14,
    4.7698313051136007E-17, -10.565764983533031, -5.772292388298E-14,
    4.5471158767424005E-17, -10.331208057362378, -5.6345364579268006E-14,
    4.3244004483712004E-17, -10.096651306938771, -5.4967805275556E-14,
    4.10168502E-17, -9.8620915887812188, -5.3589441322718012E-14,
    3.8788395008412009E-17, -9.6275353663675709, -5.2211882019006008E-14,
    3.6561240724700008E-17, -9.3929794956226829, -5.083432271529401E-14,
    3.4334086440988007E-17, -9.158420833067261, -4.9455958762456008E-14,
    3.21056312494E-17, -8.9238658421705868, -4.8078399458744004E-14,
    2.9878476965688E-17, -8.68931137884382, -4.6700840155032012E-14,
    2.7651322681976E-17, -8.4547575017229679, -4.5323280851320008E-14,
    2.5424168398264003E-17, -8.220202348178967, -4.3945238758132404E-14,
    2.31962335698264E-17, -7.9856491793192488, -4.2567518524595204E-14,
    2.09688191045392E-17, -7.75109645249273, -4.1189717826145406E-14,
    1.87412745484644E-17, -7.5165448668273047, -3.9811997592608206E-14,
    1.6513860083177202E-17, -7.2819938404670834, -3.8434196894158408E-14,
    1.4286315527102403E-17, -7.0474438483752415, -3.7056420335182388E-14,
    1.2058809998263881E-17, -6.2661225528714012E-14, -11.045865348887073,
    7.1204693730979209E-16, -6.1570036334122E-14, -11.040373840180749,
    7.1515517831017617E-16, -6.047884713953E-14, -211.0348823313123,
    7.1826341931056006E-16, -5.938702056807201E-14, -10.800322026823309,
    7.21373475872276E-16, -5.8295831373480012E-14, -10.565764983533031,
    7.2448171687266012E-16, -5.7204642178888E-14, -10.331208057362378,
    7.275899578730441E-16, -5.6113452984296004E-14, -10.096651306938771,
    7.3069819887342818E-16, -5.5021626412838004E-14, -9.8620915887812188,
    7.3380825543514416E-16, -5.3930437218246012E-14, -9.6275353663675709,
    7.3691649643552814E-16, -5.2839248023654008E-14, -9.3929794956226829,
    7.4002473743591213E-16, -5.1747421452196008E-14, -9.158420833067261,
    7.431347939976281E-16, -5.0656232257604E-14, -8.9238658421705868,
    7.46243034998012E-16, -4.9565043063012E-14, -8.68931137884382,
    7.4935127599839607E-16, -4.847385386842E-14, -8.4547575017229679,
    7.5245951699878006E-16, -4.738228224770841E-14, -8.220202348178967,
    7.5556884733596326E-16, -4.6290965577743205E-14, -7.9856491793192488,
    7.5867745144861368E-16, -4.5199585170091404E-14, -7.75109645249273,
    7.6178623711739732E-16, -4.4108268500126205E-14, -7.5165448668273047,
    7.6489484123004764E-16, -4.301688809247441E-14, -7.2819938404670834,
    7.6800362689883138E-16, -4.1925526806128586E-14, -7.0474438483752415,
    7.71112358100775E-16, 5.2247097325220006E-17, 7.2099157849929608E-16,
    -11.045865348887073, 5.1086932783972009E-17, 7.19627498904928E-16,
    -11.040373840180749, 4.9926768242724E-17, 7.1826341931056016E-16,
    -211.0348823313123, 4.8765926035272008E-17, 7.1689854294072815E-16,
    -10.800322026823309, 4.7605761494024E-17, 7.1553446334636019E-16,
    -10.565764983533031, 4.6445596952776E-17, 7.1417038375199213E-16,
    -10.331208057362378, 4.5285432411528E-17, 7.1280630415762407E-16,
    -10.096651306938771, 4.4124590204076E-17, 7.1144142778779216E-16,
    -9.8620915887812188, 4.2964425662828004E-17, 7.10077348193424E-16,
    -9.6275353663675709, 4.1804261121579995E-17, 7.0871326859905614E-16,
    -9.3929794956226829, 4.0643418914128E-17, 7.0734839222922413E-16,
    -9.158420833067261, 3.948325437288E-17, 7.0598431263485607E-16,
    -8.9238658421705868, 3.832308983163201E-17, 7.0462023304048811E-16,
    -8.68931137884382, 3.7162925290384007E-17, 7.0325615344612014E-16,
    -8.4547575017229679, 3.60023541494136E-17, 7.0189159578647377E-16,
    -8.220202348178967, 3.48420540749248E-17, 7.0052735683701285E-16,
    -7.9856491793192488, 3.3681686233815604E-17, 6.9916303821000583E-16,
    -7.75109645249273, 3.2521386159326805E-17, 6.977987992605449E-16,
    -7.5165448668273047, 3.1361018318217607E-17, 6.9643448063353768E-16,
    -7.2819938404670834, 3.0200670807094524E-17, 6.950701859097945E-16,
    -7.0474438483752415, -10.811067762614019, -6.2083479017121611E-14,
    5.3121298519880006E-17, -10.805694894794851, -6.0735249792596817E-14,
    5.0943612277576007E-17, -10.800322026823309, -5.938702056807201E-14,
    4.8765926035272E-17, -210.7949460252814, -5.8038003826476808E-14,
    4.6586967779976004E-17, -10.560507564234927, -5.6689774601952E-14,
    4.4409281537672005E-17, -10.326069161679428, -5.5341545377427206E-14,
    4.2231595295368006E-17, -10.091630876243558, -5.39933161529024E-14,
    4.0053909053064007E-17, -9.8571896331465076, -5.264429941130721E-14,
    3.7874950797768003E-17, -9.6227517584314572, -5.1296070186782415E-14,
    3.5697264555464004E-17, -9.38831417675259, -4.9947840962257608E-14,
    3.3519578313160005E-17, -9.1538738133342612, -4.8598824220662412E-14,
    3.1340620057864E-17, -8.91943699422134, -4.7250594996137605E-14,
    2.916293381556E-17, -8.6850006440594463, -4.5902365771612811E-14,
    2.6985247573256004E-17, -8.4505648214674611, -4.455413654708801E-14,
    2.4807561330952004E-17, -8.2161277050364561, -4.3205434812320967E-14,
    2.2629111880852803E-17, -7.9816924871702248, -4.185704808438209E-14,
    2.0451171235950403E-17, -7.7472576595643634, -4.0508582604736168E-14,
    1.8273103389748802E-17, -7.5128239076032974, -3.9160195876797285E-14,
    1.6095162744846402E-17, -7.2783906631917308, -3.7811730397151363E-14,
    1.3917094898644802E-17, -7.0439583923583848, -3.6463288543017558E-14,
    1.1739065212832962E-17, -6.1224201417501615E-14, -10.811067762614019,
    7.1061097805076971E-16, -6.0162342302096809E-14, -10.805694894794851,
    7.13754760495749E-16, -5.9100483186692014E-14, -10.800322026823309,
    7.16898542940728E-16, -5.8038003826476808E-14, -210.7949460252814,
    7.2004416170722881E-16, -5.6976144711072014E-14, -10.560507564234927,
    7.2318794415220813E-16, -5.59142855956672E-14, -10.326069161679428,
    7.2633172659718725E-16, -5.4852426480262406E-14, -10.091630876243558,
    7.2947550904216656E-16, -5.3789947120047206E-14, -9.8571896331465076,
    7.3262112780866732E-16, -5.2728088004642412E-14, -9.6227517584314572,
    7.3576491025364654E-16, -5.1666228889237611E-14, -9.38831417675259,
    7.3890869269862576E-16, -5.0603749529022411E-14, -9.1538738133342612,
    7.4205431146512652E-16, -4.9541890413617604E-14, -8.91943699422134,
    7.4519809391010564E-16, -4.8480031298212804E-14, -8.6850006440594463,
    7.4834187635508486E-16, -4.741817218280801E-14, -8.4505648214674611,
    7.5148565880006407E-16, -4.6355940920516972E-14, -8.2161277050364561,
    7.546305430379562E-16, -4.529395775615009E-14, -7.9816924871702248,
    7.5777469274723978E-16, -4.4231912567302164E-14, -7.7472576595643634,
    7.609190260886755E-16, -4.3169929402935289E-14, -7.5128239076032974,
    7.64063175797959E-16, -4.210788421408737E-14, -7.2783906631917308,
    7.6720750913939471E-16, -4.1045857632583755E-14, -7.0439583923583848,
    7.7035178739118468E-16, 4.9919706050816008E-17, 7.2403055217182169E-16,
    -10.811067762614019, 4.8809009550976013E-17, 7.2270201402204876E-16,
    -10.805694894794851, 4.7698313051136007E-17, 7.2137347587227613E-16,
    -10.800322026823309, 4.6586967779976017E-17, 7.2004416170722891E-16,
    -210.7949460252814, 4.5476271280136004E-17, 7.1871562355745618E-16,
    -10.560507564234927, 4.4365574780296E-17, 7.1738708540768326E-16,
    -10.326069161679428, 4.3254878280456E-17, 7.1605854725791053E-16,
    -10.091630876243558, 4.2143533009296006E-17, 7.1472923309286331E-16,
    -9.8571896331465076, 4.1032836509456012E-17, 7.1340069494309038E-16,
    -9.6227517584314572, 3.9922140009616005E-17, 7.1207215679331775E-16,
    -9.38831417675259, 3.8810794738456E-17, 7.1074284262827053E-16,
    -9.1538738133342612, 3.7700098238616E-17, 7.094143044784977E-16,
    -8.91943699422134, 3.6589401738776014E-17, 7.0808576632872487E-16,
    -8.6850006440594463, 3.5478705238936007E-17, 7.0675722817895215E-16,
    -8.4505648214674611, 3.4367619476304009E-17, 7.054282244200147E-16,
    -8.2161277050364561, 3.3256793222200003E-17, 7.0409953106718694E-16,
    -7.9816924871702248, 3.2145902090964004E-17, 7.027707601128319E-16,
    -7.7472576595643634, 3.1035075836860004E-17, 7.0144206676000413E-16,
    -7.5128239076032974, 2.9924184705624005E-17, 7.00113295805649E-16,
    -7.2783906631917308, 2.8813313037527604E-17, 6.9878454813175221E-16,
    -7.0439583923583848, -10.576273692184389, -6.0933663908264E-14,
    5.1862255651848006E-17, -10.571019337930625, -5.9614747640872008E-14,
    4.9734008572936006E-17, -10.565764983533031, -5.829583137348E-14,
    4.7605761494024E-17, -10.560507564234925, -5.6976144711072E-14,
    4.5476271280136004E-17, -210.55525321008787, -5.5657228443680005E-14,
    4.3348024201224003E-17, -10.320933262462949, -5.4338312176288008E-14,
    4.1219777122312E-17, -10.086613373329026, -5.3019395908896E-14,
    3.90915300434E-17, -9.8522905365665938, -5.1699709246488014E-14,
    3.6962039829512007E-17, -9.61797094086438, -5.0380792979096011E-14,
    3.4833792750600006E-17, -9.3836515795709019, -4.9061876711704008E-14,
    3.2705545671688E-17, -9.1493294465654724, -4.774219004929601E-14,
    3.05760554578E-17, -8.9150107305419777, -4.6423273781904E-14,
    2.8447808378888004E-17, -8.6806924248506281, -4.5104357514512005E-14,
    2.6319561299976003E-17, -8.4463745881103076, -4.3785441247120008E-14,
    2.4191314221064002E-17, -8.2120554400910155, -4.2466062742718408E-14,
    2.2062321261166402E-17, -7.9777381045159146, -4.1146992396323207E-14,
    1.99338255552592E-17, -7.7434211074414661, -3.9827845010426408E-14,
    1.7805205535854404E-17, -7.5091051204995294, -3.8508774664031207E-14,
    1.5676709829947203E-17, -7.2747895893302488, -3.7189627278134408E-14,
    1.3548089810542402E-17, -7.0404749710502355, -3.5870503004088086E-14,
    1.1419507085186881E-17, -5.9788016199524008E-14, -10.576273692184389,
    7.0917585706335215E-16, -5.8755470041252012E-14, -10.571019337930625,
    7.1235516020485622E-16, -5.772292388298E-14, -10.565764983533031,
    7.1553446334636009E-16, -5.6689774601952007E-14, -10.560507564234925,
    7.18715623557456E-16, -5.5657228443680011E-14, -210.55525321008787,
    7.2189492669896015E-16, -5.4624682285407996E-14, -10.320933262462949,
    7.2507422984046412E-16, -5.3592136127136E-14, -10.086613373329026,
    7.2825353298196818E-16, -5.2558986846108009E-14, -9.8522905365665938,
    7.3143469319306408E-16, -5.1526440687836014E-14, -9.61797094086438,
    7.3461399633456815E-16, -5.0493894529564011E-14, -9.3836515795709019,
    7.3779329947607221E-16, -4.9460745248536E-14, -9.1493294465654724,
    7.4097445968716811E-16, -4.8428199090264E-14, -8.9150107305419777,
    7.44153762828672E-16, -4.7395652931992004E-14, -8.6806924248506281,
    7.4733306597017614E-16, -4.6363106773720008E-14, -8.4463745881103076,
    7.5051236911168011E-16, -4.5330198741794406E-14, -8.2120554400910155,
    7.5369278649493926E-16, -4.4297531958971206E-14, -7.9777381045159146,
    7.5687246105036175E-16, -4.3264804863872404E-14, -7.7434211074414661,
    7.6005232131274331E-16, -4.2232138081049204E-14, -7.5091051204995294,
    7.632319958681657E-16, -4.11994109859504E-14, -7.2747895893302488,
    7.6641185613054735E-16, -4.0166701984534282E-14, -7.0404749710502355,
    7.6959166068084114E-16, 4.759367344032001E-17, 7.2706775177915607E-16,
    -10.576273692184389, 4.6532416103872008E-17, 7.25774734325908E-16,
    -10.571019337930625, 4.5471158767424005E-17, 7.2448171687266012E-16,
    -10.565764983533031, 4.4409281537672005E-17, 7.2318794415220813E-16,
    -10.560507564234925, 4.3348024201224003E-17, 7.2189492669896015E-16,
    -210.55525321008787, 4.2286766864776E-17, 7.2060190924571207E-16,
    -10.320933262462949, 4.1225509528328E-17, 7.193088917924641E-16,
    -10.086613373329026, 4.0163632298576005E-17, 7.1801511907201211E-16,
    -9.8522905365665938, 3.9102374962128009E-17, 7.16722101618764E-16,
    -9.61797094086438, 3.8041117625680007E-17, 7.1542908416551615E-16,
    -9.3836515795709019, 3.6979240395928E-17, 7.1413531144506416E-16,
    -9.1493294465654724, 3.5917983059480004E-17, 7.1284229399181608E-16,
    -8.9150107305419777, 3.4856725723032E-17, 7.1154927653856811E-16,
    -8.6806924248506281, 3.3795468386584006E-17, 7.1025625908532013E-16,
    -8.4463745881103076, 3.2733839114153604E-17, 7.0896278847174971E-16,
    -8.2120554400910155, 3.16724577990448E-17, 7.0766961996506085E-16,
    -7.9777381045159146, 3.06110144946056E-17, 7.0637637593165174E-16,
    -7.7434211074414661, 2.95496331794968E-17, 7.0508320742496289E-16,
    -7.5091051204995294, 2.848818987505761E-17, 7.0378996339155359E-16,
    -7.2747895893302488, 2.7426765167417521E-17, 7.0249674201616071E-16,
    -7.0404749710502355, -10.341479856129229, -5.9783848799406411E-14,
    5.0603212783816006E-17, -10.336343956813447, -5.8494245489147212E-14,
    4.8524404868296004E-17, -10.331208057362378, -5.7204642178888E-14,
    4.6445596952776E-17, -10.326069161679429, -5.5914285595667207E-14,
    4.4365574780296E-17, -10.320933262462951, -5.4624682285407996E-14,
    4.2286766864776E-17, -210.31579736311858, -5.3335078975148804E-14,
    4.0207958949256E-17, -10.08159587029636, -5.20454756648896E-14,
    3.8129151033736E-17, -9.8473914398770823, -5.0755119081668811E-14,
    3.6049128861256004E-17, -9.6131901231962775, -4.9465515771409607E-14,
    3.3970320945736E-17, -9.3789889822967556, -4.8175912461150409E-14,
    3.1891513030216E-17, -9.1447850797179253, -4.6885555877929607E-14,
    2.9811490857736E-17, -8.9105844667941234, -4.55959525676704E-14,
    2.7732682942216E-17, -8.67638420557332, -4.4306349257411205E-14,
    2.5653875026696002E-17, -8.4421843546846613, -4.3016745947152007E-14,
    2.3575067111176E-17, -8.2079831750942063, -4.1726690673115843E-14,
    2.1495530641480002E-17, -7.97378372182736, -4.0436936708264323E-14,
    1.9416479874568E-17, -7.7395845552843232, -3.9147107416116648E-14,
    1.7337307681960003E-17, -7.5053863333786381, -3.7857353451265122E-14,
    1.5258256915048003E-17, -7.2711885154687677, -3.6567524159117441E-14,
    1.3179084722440003E-17, -7.036991549742087, -3.5277717465158614E-14,
    1.1099948957540802E-17, -5.8351830981546414E-14, -10.341479856129229,
    7.0774073607593449E-16, -5.73485977804072E-14, -10.336343956813447,
    7.1095555991396341E-16, -5.6345364579268006E-14, -10.331208057362378,
    7.14170383751992E-16, -5.5341545377427206E-14, -10.326069161679429,
    7.1738708540768326E-16, -5.4338312176288008E-14, -10.320933262462951,
    7.2060190924571207E-16, -5.33350789751488E-14, -210.31579736311858,
    7.2381673308374089E-16, -5.23318457740096E-14, -10.08159587029636,
    7.2703155692176971E-16, -5.13280265721688E-14, -9.8473914398770823,
    7.3024825857746094E-16, -5.03247933710296E-14, -9.6131901231962775,
    7.3346308241548965E-16, -4.9321560169890405E-14, -9.3789889822967556,
    7.3667790625351857E-16, -4.8317740968049605E-14, -9.1447850797179253,
    7.398946079092098E-16, -4.7314507766910395E-14, -8.9105844667941234,
    7.4310943174723842E-16, -4.63112745657712E-14, -8.67638420557332,
    7.4632425558526733E-16, -4.5308041364632006E-14, -8.4421843546846613,
    7.4953907942329605E-16, -4.4304456563071847E-14, -8.2079831750942063,
    7.5275502995192232E-16, -4.3301106161792329E-14, -7.97378372182736,
    7.5597022935348362E-16, -4.2297697160442644E-14, -7.7395845552843232,
    7.5918561653681111E-16, -4.1294346759163126E-14, -7.5053863333786381,
    7.6240081593837241E-16, -4.0290937757813447E-14, -7.2711885154687677,
    7.656162031217E-16, -3.928754633648481E-14, -7.036991549742087,
    7.688315339704976E-16, 4.5267640829824E-17, 7.3010495138649046E-16,
    -10.341479856129229, 4.4255822656768008E-17, 7.2884745462976723E-16,
    -10.336343956813447, 4.3244004483712E-17, 7.275899578730441E-16,
    -10.331208057362378, 4.2231595295368006E-17, 7.2633172659718734E-16,
    -10.326069161679429, 4.1219777122311996E-17, 7.2507422984046412E-16,
    -10.320933262462951, 4.0207958949256E-17, 7.2381673308374089E-16,
    -210.31579736311858, 3.9196140776199995E-17, 7.2255923632701766E-16,
    -10.08159587029636, 3.8183731587856003E-17, 7.2130100505116091E-16,
    -9.8473914398770823, 3.7171913414800006E-17, 7.2004350829443768E-16,
    -9.6131901231962775, 3.6160095241743996E-17, 7.1878601153771455E-16,
    -9.3789889822967556, 3.51476860534E-17, 7.1752778026185779E-16,
    -9.1447850797179253, 3.4135867880344E-17, 7.1627028350513447E-16,
    -8.9105844667941234, 3.3124049707288003E-17, 7.1501278674841134E-16,
    -8.67638420557332, 3.2112231534232E-17, 7.1375528999168811E-16,
    -8.4421843546846613, 3.11000587520032E-17, 7.1249735252348481E-16,
    -8.2079831750942063, 3.00881223758896E-17, 7.1123970886293476E-16,
    -7.97378372182736, 2.90761268982472E-17, 7.0998199175047159E-16,
    -7.7395845552843232, 2.80641905221336E-17, 7.0872434808992164E-16,
    -7.5053863333786381, 2.7052195044491203E-17, 7.0746663097745827E-16,
    -7.2711885154687677, 2.6040217297307439E-17, 7.0620893590056912E-16,
    -7.036991549742087, -10.106686313076018, -5.86340336905488E-14,
    4.9344169915784007E-17, -10.101668810070741, -5.7373743337422417E-14,
    4.7314801163656009E-17, -10.096651306938771, -5.6113452984296004E-14,
    4.5285432411528E-17, -10.091630876243558, -5.48524264802624E-14,
    4.3254878280456009E-17, -10.086613373329026, -5.3592136127136E-14,
    4.1225509528328005E-17, -10.08159587029636, -5.2331845774009607E-14,
    3.91961407762E-17, -210.0765783671358, -5.10715554208832E-14,
    3.7166772024072004E-17, -9.84249234306944, -4.9810528916849609E-14,
    3.5136217893000007E-17, -9.60840930541858, -4.8550238563723209E-14,
    3.3106849140872003E-17, -9.3743263849215861, -4.7289948210596809E-14,
    3.1077480388744E-17, -9.14024071277792, -4.6028921706563211E-14,
    2.9046926257672E-17, -8.9061582029675073, -4.4768631353436805E-14,
    2.7017557505544005E-17, -8.672075986227517, -4.3508341000310405E-14,
    2.4988188753416002E-17, -8.4379941211905258, -4.2248050647184005E-14,
    2.2958820001288004E-17, -8.2039109100289078, -4.098731860351329E-14,
    2.0928740021793602E-17, -7.9698293390874371, -3.9726881020205446E-14,
    1.88991341938768E-17, -7.7357480030929349, -3.8466369821806889E-14,
    1.68694098280656E-17, -7.5016675462235023, -3.7205932238499044E-14,
    1.4839804000148803E-17, -7.2675874415901642, -3.5945421040100487E-14,
    1.2810079634337603E-17, -7.0335081284339376, -3.4684931926229141E-14,
    1.0780390829894723E-17, -5.6915645763568807E-14, -10.106686313076018,
    7.0630561508851693E-16, -5.5941725519562407E-14, -10.101668810070741,
    7.095559596230706E-16, -5.4967805275556008E-14, -10.096651306938771,
    7.1280630415762407E-16, -5.3993316152902412E-14, -10.091630876243558,
    7.1605854725791043E-16, -5.3019395908896012E-14, -10.086613373329026,
    7.193088917924641E-16, -5.20454756648896E-14, -10.08159587029636,
    7.2255923632701766E-16, -5.10715554208832E-14, -210.0765783671358,
    7.2580958086157133E-16, -5.0097066298229604E-14, -9.84249234306944,
    7.2906182396185769E-16, -4.9123146054223205E-14, -9.60840930541858,
    7.3231216849641126E-16, -4.8149225810216811E-14, -9.3743263849215861,
    7.3556251303096493E-16, -4.7174736687563209E-14, -9.14024071277792,
    7.3881475613125139E-16, -4.6200816443556797E-14, -8.9061582029675073,
    7.4206510066580486E-16, -4.5226896199550403E-14, -8.672075986227517,
    7.4531544520035852E-16, -4.425297595554401E-14, -8.4379941211905258,
    7.4856578973491209E-16, -4.3278714384349288E-14, -8.2039109100289078,
    7.5181727340890537E-16, -4.2304680364613451E-14, -7.9698293390874371,
    7.5506799765660558E-16, -4.1330589457012884E-14, -7.7357480030929349,
    7.58318911760879E-16, -4.0356555437277047E-14, -7.5016675462235023,
    7.6156963600857911E-16, -3.9382464529676486E-14, -7.2675874415901642,
    7.6482055011285264E-16, -3.8408390688435343E-14, -7.0335081284339376,
    7.6807140726015406E-16, 4.2941608219328E-17, 7.3314215099382484E-16,
    -10.106686313076018, 4.1979229209664008E-17, 7.3192017493362636E-16,
    -10.101668810070741, 4.10168502E-17, 7.3069819887342809E-16,
    -10.096651306938771, 4.0053909053064007E-17, 7.2947550904216656E-16,
    -10.091630876243558, 3.90915300434E-17, 7.2825353298196818E-16,
    -10.086613373329026, 3.8129151033736009E-17, 7.2703155692176971E-16,
    -10.08159587029636, 3.7166772024072E-17, 7.2580958086157133E-16,
    -210.0765783671358, 3.6203830877136008E-17, 7.2458689103030971E-16,
    -9.84249234306944, 3.5241451867472E-17, 7.2336491497011123E-16,
    -9.60840930541858, 3.4279072857808E-17, 7.22142938909913E-16,
    -9.3743263849215861, 3.3316131710872E-17, 7.2092024907865133E-16,
    -9.14024071277792, 3.2353752701208E-17, 7.1969827301845285E-16,
    -8.9061582029675073, 3.1391373691544004E-17, 7.1847629695825447E-16,
    -8.672075986227517, 3.0428994681880005E-17, 7.172543208980561E-16,
    -8.4379941211905258, 2.94662783898528E-17, 7.1603191657521981E-16,
    -8.2039109100289078, 2.8503786952734406E-17, 7.1480979776080867E-16,
    -7.9698293390874371, 2.7541239301888804E-17, 7.1358760756929143E-16,
    -7.7357480030929349, 2.6578747864770404E-17, 7.1236548875488039E-16,
    -7.5016675462235023, 2.5616200213924805E-17, 7.11143298563363E-16,
    -7.2675874415901642, 2.4653669427197363E-17, 7.0992112978497752E-16,
    -7.0335081284339376, -9.8718897822287275, -5.74835469607164E-14,
    4.8084391625516004E-17, -9.866990685561472, -5.6252586686777213E-14,
    4.6104490914796008E-17, -9.8620915887812188, -5.5021626412838004E-14,
    4.4124590204076E-17, -9.85718963314651, -5.3789947120047206E-14,
    4.2143533009296006E-17, -9.8522905365665956, -5.2558986846108003E-14,
    4.0163632298576005E-17, -9.8473914398770823, -5.1328026572168807E-14,
    3.8183731587856003E-17, -9.84249234306944, -5.0097066298229604E-14,
    3.6203830877136008E-17, -209.83759038740612, -4.8865387005438813E-14,
    3.4222773682356007E-17, -9.603625697455076, -4.763442673149961E-14,
    3.2242872971636006E-17, -9.36966106602924, -4.6403466457560413E-14,
    3.0262972260916004E-17, -9.1356936930295056, -4.5171787164769609E-14,
    2.8281915066136004E-17, -8.9017293550011374, -4.3940826890830407E-14,
    2.6302014355416002E-17, -8.6677652514157444, -4.270986661689121E-14,
    2.4322113644696E-17, -8.4338014409007762, -4.1478906342952007E-14,
    2.2342212933976003E-17, -8.19983626685215, -4.0247514657701848E-14,
    2.0361618332820005E-17, -7.9658726469041685, -3.9016410579992326E-14,
    1.8381486325288E-17, -7.7319092101474451, -3.7785234600397651E-14,
    1.6401238669350002E-17, -7.497946586999495, -3.6554130522688123E-14,
    1.4421106661818004E-17, -7.2639842643148116, -3.5322954543093442E-14,
    1.2440859005880003E-17, -7.0300226724170809, -3.4091800134064311E-14,
    1.0460646044463803E-17, -5.5478621652356411E-14, -9.8718897822287275,
    7.0486965582949455E-16, -5.45340314875372E-14, -9.866990685561472,
    7.0815554180864336E-16, -5.3589441322718006E-14, -9.8620915887812188,
    7.1144142778779206E-16, -5.264429941130721E-14, -9.85718963314651,
    7.1472923309286321E-16, -5.1699709246488014E-14, -9.8522905365665956,
    7.1801511907201211E-16, -5.0755119081668805E-14, -9.8473914398770823,
    7.2130100505116091E-16, -4.9810528916849609E-14, -9.84249234306944,
    7.2458689103030971E-16, -4.8865387005438806E-14, -209.83759038740612,
    7.27874696335381E-16, -4.792079684061961E-14, -9.603625697455076,
    7.3116058231452966E-16, -4.6976206675800408E-14, -9.36966106602924,
    7.3444646829367856E-16, -4.6031064764389605E-14, -9.1356936930295056,
    7.3773427359874971E-16, -4.50864745995704E-14, -8.9017293550011374,
    7.4102015957789841E-16, -4.4141884434751207E-14, -8.6677652514157444,
    7.4430604555704731E-16, -4.3197294269932011E-14, -8.4338014409007762,
    7.4759193153619611E-16, -4.2252373057157849E-14, -8.19983626685215,
    7.5087896911089832E-16, -4.130767254302033E-14, -7.9658726469041685,
    7.5416523895523169E-16, -4.0362916854223645E-14, -7.7319092101474451,
    7.5745170073215719E-16, -3.9418216340086125E-14, -7.497946586999495,
    7.6073797057649036E-16, -3.8473460651289446E-14, -7.2639842643148116,
    7.64024432353416E-16, -3.7528721514890512E-14, -7.0300226724170809,
    7.6731083655056378E-16, 4.0614216944924E-17, 7.3618112466635045E-16,
    -9.8718897822287275, 3.9701305976668012E-17, 7.3499469005074721E-16,
    -9.866990685561472, 3.8788395008412E-17, 7.3380825543514416E-16,
    -9.8620915887812188, 3.787495079776801E-17, 7.3262112780866732E-16,
    -9.85718963314651, 3.6962039829512E-17, 7.3143469319306418E-16,
    -9.8522905365665956, 3.6049128861256004E-17, 7.3024825857746084E-16,
    -9.8473914398770823, 3.5136217893E-17, 7.2906182396185769E-16,
    -9.84249234306944, 3.4222773682356007E-17, 7.2787469633538086E-16,
    -209.83759038740612, 3.3309862714100004E-17, 7.2668826171977761E-16,
    -9.603625697455076, 3.2396951745844E-17, 7.2550182710417457E-16,
    -9.36966106602924, 3.14835075352E-17, 7.2431469947769773E-16,
    -9.1356936930295056, 3.0570596566944005E-17, 7.2312826486209449E-16,
    -8.9017293550011374, 2.9657685598688008E-17, 7.2194183024649134E-16,
    -8.6677652514157444, 2.8744774630432005E-17, 7.207553956308881E-16,
    -8.4338014409007762, 2.78315437167432E-17, 7.1956854520876074E-16,
    -8.19983626685215, 2.6918526100009604E-17, 7.1838197199098276E-16,
    -7.9658726469041685, 2.6005455159037204E-17, 7.171953294721176E-16,
    -7.7319092101474451, 2.5092437542303603E-17, 7.1600875625433962E-16,
    -7.497946586999495, 2.4179366601331206E-17, 7.1482211373547427E-16,
    -7.2639842643148116, 2.3266311657630442E-17, 7.1363549200693533E-16,
    -7.0300226724170809, -9.6370970017757, -5.6333731851858809E-14,
    4.682534875748401E-17, -9.6323161841230132, -5.5132084535052411E-14,
    4.4894887210156007E-17, -9.6275353663675709, -5.3930437218246E-14,
    4.2964425662828004E-17, -9.6227517584314572, -5.2728088004642405E-14,
    4.1032836509456006E-17, -9.6179709408643816, -5.1526440687836E-14,
    3.9102374962128E-17, -9.6131901231962757, -5.03247933710296E-14,
    3.7171913414800006E-17, -9.60840930541858, -4.9123146054223205E-14,
    3.5241451867472E-17, -9.603625697455076, -4.792079684061961E-14,
    3.3309862714100004E-17, -209.5988448798594, -4.6719149523813206E-14,
    3.1379401166772008E-17, -9.3649984688201577, -4.5517502207006814E-14,
    2.9448939619444005E-17, -9.13114932623846, -4.4315152993403207E-14,
    2.7517350466072E-17, -8.89730309130636, -4.3113505676596809E-14,
    2.5586888918744003E-17, -8.6634570321898057, -4.191185835979041E-14,
    2.3656427371416E-17, -8.4296112075162242, -4.0710211042984006E-14,
    2.1725965824088003E-17, -8.1957640018724653, -3.9508142588099289E-14,
    1.9794827713133604E-17, -7.9619182642156119, -3.8306354891933449E-14,
    1.7864140644596802E-17, -7.7280726579903023, -3.7104497006087891E-14,
    1.5933340815455604E-17, -7.4942277998786047, -3.5902709309922045E-14,
    1.4002653746918804E-17, -7.26038319045333, -3.4700851424076487E-14,
    1.2071853917777603E-17, -7.0265392511089324, -3.3499014595134839E-14,
    1.0141087916817723E-17, -5.4042436434378816E-14, -9.6370970017757,
    7.034345348420769E-16, -5.3127159226692406E-14, -9.6323161841230132,
    7.0675594151775055E-16, -5.2211882019006008E-14, -9.6275353663675709,
    7.100773481934241E-16, -5.1296070186782409E-14, -9.6227517584314572,
    7.1340069494309048E-16, -5.0380792979096011E-14, -9.6179709408643816,
    7.1672210161876413E-16, -4.94655157714096E-14, -9.6131901231962757,
    7.2004350829443768E-16, -4.8550238563723209E-14, -9.60840930541858,
    7.2336491497011133E-16, -4.7634426731499604E-14, -9.603625697455076,
    7.2668826171977771E-16, -4.6719149523813206E-14, -209.5988448798594,
    7.3000966839545126E-16, -4.5803872316126808E-14, -9.3649984688201577,
    7.33331075071125E-16, -4.48880604839032E-14, -9.13114932623846,
    7.3665442182079139E-16, -4.39727832762168E-14, -8.89730309130636,
    7.3997582849646485E-16, -4.3057506068530407E-14, -8.6634570321898057,
    7.432972351721385E-16, -4.2142228860844009E-14, -8.4296112075162242,
    7.46618641847812E-16, -4.122663087843529E-14, -8.1957640018724653,
    7.4994121256788138E-16, -4.0311246745841452E-14, -7.9619182642156119,
    7.5326300725835355E-16, -3.9395809150793884E-14, -7.7280726579903023,
    7.56584995956225E-16, -3.8480425018200046E-14, -7.4942277998786047,
    7.5990679064669707E-16, -3.7564987423152485E-14, -7.26038319045333,
    7.6322877934456861E-16, -3.6649565866841039E-14, -7.0265392511089324,
    7.6655070984022024E-16, 3.8288184334428005E-17, 7.3921832427368483E-16,
    -9.6370970017757, 3.7424712529564006E-17, 7.3806741035460644E-16,
    -9.6323161841230132, 3.6561240724700008E-17, 7.3691649643552814E-16,
    -9.6275353663675709, 3.5697264555464011E-17, 7.3576491025364654E-16,
    -9.6227517584314572, 3.4833792750600006E-17, 7.3461399633456815E-16,
    -9.6179709408643816, 3.3970320945736008E-17, 7.3346308241548965E-16,
    -9.6131901231962757, 3.3106849140872003E-17, 7.3231216849641126E-16,
    -9.60840930541858, 3.2242872971636006E-17, 7.3116058231452966E-16,
    -9.603625697455076, 3.1379401166772008E-17, 7.3000966839545126E-16,
    -209.5988448798594, 3.0515929361908003E-17, 7.28858754476373E-16,
    -9.3649984688201577, 2.9651953192672006E-17, 7.2770716829449136E-16,
    -9.13114932623846, 2.8788481387808E-17, 7.2655625437541287E-16,
    -8.89730309130636, 2.7925009582944009E-17, 7.2540534045633448E-16,
    -8.6634570321898057, 2.7061537778080005E-17, 7.2425442653725618E-16,
    -8.4296112075162242, 2.6197763354592803E-17, 7.2310310926049584E-16,
    -8.1957640018724653, 2.5334190676854404E-17, 7.2195206088885667E-16,
    -7.9619182642156119, 2.4470567562678804E-17, 7.2080094529093745E-16,
    -7.7280726579903023, 2.3606994884940405E-17, 7.1964989691929837E-16,
    -7.4942277998786047, 2.2743371770764805E-17, 7.1849878132137886E-16,
    -7.26038319045333, 2.1879763787520366E-17, 7.1734768589134373E-16,
    -7.0265392511089324, -9.4023046902394611, -5.5183916743001208E-14,
    4.5566305889452004E-17, -9.3976420929824371, -5.4011582383327615E-14,
    4.3685283505516006E-17, -9.3929794956226829, -5.2839248023654E-14,
    4.180426112158E-17, -9.38831417675259, -5.1666228889237605E-14,
    3.9922140009616005E-17, -9.3836515795709019, -5.0493894529564005E-14,
    3.8041117625680007E-17, -9.3789889822967538, -4.9321560169890405E-14,
    3.6160095241744E-17, -9.3743263849215843, -4.8149225810216805E-14,
    3.4279072857808004E-17, -9.3696610660292379, -4.6976206675800414E-14,
    3.2396951745844008E-17, -9.3649984688201577, -4.5803872316126808E-14,
    3.0515929361908003E-17, -209.36033587151744, -4.4631537956453208E-14,
    2.8634906977972005E-17, -9.1266049593635223, -4.3458518822036811E-14,
    2.6752785866008002E-17, -8.89287682753625, -4.2286184462363204E-14,
    2.4871763482072004E-17, -8.6591488128970813, -4.1113850102689611E-14,
    2.2990741098136003E-17, -8.4254209740734574, -3.9941515743016004E-14,
    2.11097187142E-17, -8.19169173684826, -3.876877051849673E-14,
    1.9228037093447204E-17, -7.9579638814928124, -3.7596299203874565E-14,
    1.73467949639056E-17, -7.724236105798914, -3.6423759411778131E-14,
    1.5465442961561202E-17, -7.4905090127234679, -3.5251288097155967E-14,
    1.3584200832019603E-17, -7.2567821165747262, -3.4078748305059526E-14,
    1.1702848829675204E-17, -7.023055829800783, -3.2906229056205367E-14,
    9.8215297891716419E-18, -5.2606251216401209E-14, -9.4023046902394611,
    7.0199941385465934E-16, -5.172028696584761E-14, -9.3976420929824371,
    7.0535634122685774E-16, -5.083432271529401E-14, -9.3929794956226829,
    7.08713268599056E-16, -4.9947840962257608E-14, -9.38831417675259,
    7.1207215679331765E-16, -4.9061876711704008E-14, -9.3836515795709019,
    7.1542908416551615E-16, -4.81759124611504E-14, -9.3789889822967538,
    7.1878601153771445E-16, -4.72899482105968E-14, -9.3743263849215843,
    7.22142938909913E-16, -4.6403466457560407E-14, -9.3696610660292379,
    7.2550182710417457E-16, -4.5517502207006807E-14, -9.3649984688201577,
    7.2885875447637287E-16, -4.4631537956453208E-14, -209.36033587151744,
    7.3221568184857137E-16, -4.3745056203416806E-14, -9.1266049593635223,
    7.35574570042833E-16, -4.28590919528632E-14, -8.89287682753625,
    7.3893149741503119E-16, -4.1973127702309607E-14, -8.6591488128970813,
    7.4228842478722968E-16, -4.1087163451756007E-14, -8.4254209740734574,
    7.4564535215942808E-16, -4.0200888699712731E-14, -8.19169173684826,
    7.4900345602486434E-16, -3.9314820948662568E-14, -7.9579638814928124,
    7.5236077556147552E-16, -3.8428701447364124E-14, -7.724236105798914,
    7.5571829118029279E-16, -3.7542633696313968E-14, -7.4905090127234679,
    7.5907561071690388E-16, -3.6656514195015524E-14, -7.2567821165747262,
    7.6243312633572125E-16, -3.5770410218791567E-14, -7.023055829800783,
    7.657905831298767E-16, 3.5962151723932006E-17, 7.4225552388101932E-16,
    -9.4023046902394611, 3.5148119082460007E-17, 7.4114013065846558E-16,
    -9.3976420929824371, 3.4334086440988007E-17, 7.4002473743591213E-16,
    -9.3929794956226829, 3.3519578313160011E-17, 7.3890869269862576E-16,
    -9.38831417675259, 3.2705545671688006E-17, 7.3779329947607211E-16,
    -9.3836515795709019, 3.1891513030216006E-17, 7.3667790625351847E-16,
    -9.3789889822967538, 3.1077480388744E-17, 7.3556251303096493E-16,
    -9.3743263849215843, 3.0262972260916004E-17, 7.3444646829367846E-16,
    -9.3696610660292379, 2.9448939619444005E-17, 7.3333107507112481E-16,
    -9.3649984688201577, 2.8634906977972E-17, 7.3221568184857137E-16,
    -209.36033587151744, 2.7820398850144003E-17, 7.310996371112849E-16,
    -9.1266049593635223, 2.7006366208672004E-17, 7.2998424388873125E-16,
    -8.89287682753625, 2.6192333567200004E-17, 7.2886885066617771E-16,
    -8.6591488128970813, 2.5378300925728004E-17, 7.2775345744362417E-16,
    -8.4254209740734574, 2.4563982992442405E-17, 7.2663767331223085E-16,
    -8.19169173684826, 2.3749855253699204E-17, 7.2552214978673068E-16,
    -7.9579638814928124, 2.2935679966320404E-17, 7.2440656110975729E-16,
    -7.724236105798914, 2.2121552227577203E-17, 7.23291037584257E-16,
    -7.4905090127234679, 2.1307376940198407E-17, 7.2217544890728354E-16,
    -7.2567821165747262, 2.0493215917410284E-17, 7.2105987977575223E-16,
    -7.023055829800783, -9.1675095668102937, -5.4033430013168807E-14,
    4.4306527599184007E-17, -9.1629651999901416, -5.2890425732682411E-14,
    4.2474973256656005E-17, -9.158420833067261, -5.1747421452196E-14,
    4.0643418914128E-17, -9.1538738133342612, -5.0603749529022405E-14,
    3.8810794738456E-17, -9.1493294465654724, -4.9460745248536E-14,
    3.6979240395928E-17, -9.1447850797179235, -4.8317740968049605E-14,
    3.51476860534E-17, -9.14024071277792, -4.71747366875632E-14,
    3.3316131710872E-17, -9.1356936930295038, -4.6031064764389611E-14,
    3.14835075352E-17, -9.13114932623846, -4.4888060483903209E-14,
    2.9651953192672006E-17, -9.1266049593635241, -4.3745056203416812E-14,
    2.7820398850144003E-17, -209.12205793968755, -4.2601384280243209E-14,
    2.5987774674472E-17, -8.8884479796349467, -4.1458379999756806E-14,
    2.4156220331944E-17, -8.6548380781418146, -4.0315375719270409E-14,
    2.2324665989416002E-17, -8.4212282938367871, -3.9172371438784007E-14,
    2.0493111646888003E-17, -8.1876170937228725, -3.8028966572685287E-14,
    1.86609154044736E-17, -7.9540071893506381, -3.6885828763661445E-14,
    1.6829147095316802E-17, -7.7203973128705465, -3.5742624190368887E-14,
    1.4997271802845603E-17, -7.4867880534994615, -3.4599486381345046E-14,
    1.3165503493688802E-17, -7.2531789392993735, -3.3456281808052488E-14,
    1.1333628201217603E-17, -7.0195703737839263, -3.2313097264040543E-14,
    9.5017850037407223E-18, -5.1169227105188813E-14, -9.1675095668102937,
    7.0056345459563686E-16, -5.0312592933822404E-14, -9.1629651999901416,
    7.0395592341243059E-16, -4.9455958762456008E-14, -9.158420833067261,
    7.07348392229224E-16, -4.8598824220662406E-14, -9.1538738133342612,
    7.1074284262827043E-16, -4.774219004929601E-14, -9.1493294465654724,
    7.1413531144506406E-16, -4.68855558779296E-14, -9.1447850797179235,
    7.1752778026185769E-16, -4.6028921706563205E-14, -9.14024071277792,
    7.2092024907865133E-16, -4.51717871647696E-14, -9.1356936930295038,
    7.2431469947769773E-16, -4.4315152993403207E-14, -9.13114932623846,
    7.2770716829449126E-16, -4.3458518822036811E-14, -9.1266049593635241,
    7.31099637111285E-16, -4.26013842802432E-14, -209.12205793968755,
    7.344940875103314E-16, -4.17447501088768E-14, -8.8884479796349467,
    7.3788655632712483E-16, -4.0888115937510404E-14, -8.6548380781418146,
    7.4127902514391847E-16, -4.0031481766144008E-14, -8.4212282938367871,
    7.446714939607121E-16, -3.9174547372521287E-14, -8.1876170937228725,
    7.4806515172685738E-16, -3.8317813127069447E-14, -7.9540071893506381,
    7.5145801686010152E-16, -3.7461028844574885E-14, -7.7203973128705465,
    7.54851080151571E-16, -3.6604294599123045E-14, -7.4867880534994615,
    7.5824394528481512E-16, -3.5747510316628484E-14, -7.2531789392993735,
    7.6163700857628468E-16, -3.4890741045246742E-14, -7.0195703737839263,
    7.6503001242028642E-16, 3.3634760449528008E-17, 7.4529449755354483E-16,
    -9.1675095668102937, 3.2870195849464004E-17, 7.4421464577558642E-16,
    -9.1629651999901416, 3.21056312494E-17, 7.431347939976281E-16,
    -9.158420833067261, 3.1340620057864008E-17, 7.4205431146512652E-16,
    -9.1538738133342612, 3.0576055457800004E-17, 7.4097445968716811E-16,
    -9.1493294465654724, 2.9811490857736006E-17, 7.398946079092097E-16,
    -9.1447850797179235, 2.9046926257671996E-17, 7.3881475613125129E-16,
    -9.14024071277792, 2.8281915066136004E-17, 7.3773427359874971E-16,
    -9.1356936930295038, 2.7517350466072003E-17, 7.366544218207912E-16,
    -9.13114932623846, 2.6752785866008E-17, 7.35574570042833E-16,
    -9.1266049593635241, 2.5987774674472E-17, 7.344940875103314E-16,
    -209.12205793968755, 2.5223210074408E-17, 7.3341423573237289E-16,
    -8.8884479796349467, 2.4458645474344002E-17, 7.3233438395441448E-16,
    -8.6548380781418146, 2.3694080874280004E-17, 7.3125453217645617E-16,
    -8.4212282938367871, 2.29292483193328E-17, 7.3017430194577187E-16,
    -8.1876170937228725, 2.2164594400974402E-17, 7.2909432401690467E-16,
    -7.9540071893506381, 2.1399895823468804E-17, 7.2801428301258346E-16,
    -7.7203973128705465, 2.0635241905110403E-17, 7.2693430508371636E-16,
    -7.4867880534994615, 1.9870543327604802E-17, 7.25854264079395E-16,
    -7.2531789392993735, 1.9105858147843363E-17, 7.2477424199770994E-16,
    -7.0195703737839263, -8.9327183696834176, -5.2883614904311207E-14,
    4.3047484731152008E-17, -8.9282921059698044, -5.1769923580957609E-14,
    4.126536955201601E-17, -8.923865842170585, -5.0656232257604E-14,
    3.948325437288E-17, -8.91943699422134, -4.9541890413617604E-14,
    3.7700098238616E-17, -8.9150107305419777, -4.8428199090264E-14,
    3.5917983059480004E-17, -8.9105844667941216, -4.7314507766910407E-14,
    3.4135867880344E-17, -8.9061582029675055, -4.6200816443556797E-14,
    3.2353752701208E-17, -8.9017293550011374, -4.5086474599570409E-14,
    3.0570596566944005E-17, -8.89730309130636, -4.3972783276216805E-14,
    2.8788481387808E-17, -8.89287682753625, -4.2859091952863213E-14,
    2.7006366208672004E-17, -8.8884479796349467, -4.1744750108876806E-14,
    2.5223210074408E-17, -208.8840217160126, -4.06310587855232E-14,
    2.3441094895272002E-17, -8.65052985898094, -3.951736746216961E-14,
    2.1658979716136E-17, -8.4170380605087409, -3.8403676138816005E-14,
    1.9876864537E-17, -8.1835448287962667, -3.7289594503082728E-14,
    1.8094124784787203E-17, -7.950052806713451, -3.6175773075602562E-14,
    1.63118014146256E-17, -7.7165607607544988, -3.5061886596059127E-14,
    1.4529373948951202E-17, -7.4830692663956926, -3.3948065168578968E-14,
    1.2747050578789602E-17, -7.2495778654378915, -3.2834178689035527E-14,
    1.0964623113115202E-17, -7.016086952475777, -3.1720311725111071E-14,
    9.1822268760946417E-18, -4.9733041887211212E-14, -8.9327183696834176,
    6.991283336082193E-16, -4.89057206729776E-14, -8.9282921059698044,
    7.0255632312153778E-16, -4.8078399458744004E-14, -8.923865842170585,
    7.0598431263485607E-16, -4.7250594996137605E-14, -8.91943699422134,
    7.094143044784976E-16, -4.6423273781904014E-14, -8.9150107305419777,
    7.1284229399181608E-16, -4.55959525676704E-14, -8.9105844667941216,
    7.1627028350513447E-16, -4.4768631353436805E-14, -8.9061582029675055,
    7.1969827301845295E-16, -4.3940826890830407E-14, -8.9017293550011374,
    7.2312826486209449E-16, -4.3113505676596809E-14, -8.89730309130636,
    7.2655625437541287E-16, -4.2286184462363211E-14, -8.89287682753625,
    7.2998424388873135E-16, -4.14583799997568E-14, -8.8884479796349467,
    7.33414235732373E-16, -4.06310587855232E-14, -208.8840217160126,
    7.3684222524569127E-16, -3.9803737571289604E-14, -8.65052985898094,
    7.4027021475900966E-16, -3.8976416357056006E-14, -8.4170380605087409,
    7.43698204272328E-16, -3.8148805193798727E-14, -8.1835448287962667,
    7.4712739518384034E-16, -3.7321387329890569E-14, -7.950052806713451,
    7.5055578516322349E-16, -3.6493921141145125E-14, -7.7165607607544988,
    7.5398437537563878E-16, -3.5666503277236967E-14, -7.4830692663956926,
    7.5741276535502183E-16, -3.4839037088491522E-14, -7.2495778654378915,
    7.6084135556743732E-16, -3.4011585397197269E-14, -7.016086952475777,
    7.6426988570994288E-16, 3.1308727839032004E-17, 7.4833169716087931E-16,
    -8.9327183696834176, 3.0593602402360005E-17, 7.4728736607944565E-16,
    -8.9282921059698044, 2.9878476965688E-17, 7.4624303499801209E-16,
    -8.923865842170585, 2.916293381556E-17, 7.4519809391010574E-16,
    -8.91943699422134, 2.8447808378888004E-17, 7.4415376282867218E-16,
    -8.9150107305419777, 2.7732682942216E-17, 7.4310943174723852E-16,
    -8.9105844667941216, 2.7017557505544E-17, 7.42065100665805E-16,
    -8.9061582029675055, 2.6302014355416002E-17, 7.4102015957789851E-16,
    -8.9017293550011374, 2.5586888918744003E-17, 7.3997582849646485E-16,
    -8.89730309130636, 2.4871763482072E-17, 7.3893149741503138E-16,
    -8.89287682753625, 2.4156220331944E-17, 7.3788655632712493E-16,
    -8.8884479796349467, 2.3441094895272002E-17, 7.3684222524569127E-16,
    -208.8840217160126, 2.2725969458600003E-17, 7.3579789416425771E-16,
    -8.65052985898094, 2.2010844021928004E-17, 7.3475356308282415E-16,
    -8.4170380605087409, 2.12954679571824E-17, 7.3370886599750688E-16,
    -8.1835448287962667, 2.0580258977819202E-17, 7.3266441291477868E-16,
    -7.950052806713451, 1.98650082271104E-17, 7.3161989883140331E-16,
    -7.7165607607544988, 1.9149799247747202E-17, 7.3057544574867511E-16,
    -7.4830692663956926, 1.8434548497038404E-17, 7.2953093166529954E-16,
    -7.2495778654378915, 1.7719310277733281E-17, 7.2848643588211845E-16,
    -7.016086952475777, -8.6979278173984422, -5.1733799795453607E-14,
    4.1788441863120008E-17, -8.69361959815538, -5.0649421429232813E-14,
    4.0055765847376009E-17, -8.68931137884382, -4.9565043063012E-14,
    3.8323089831632E-17, -8.6850006440594463, -4.8480031298212804E-14,
    3.6589401738776E-17, -8.6806924248506281, -4.7395652931992E-14,
    3.4856725723032E-17, -8.6763842055733189, -4.6311274565771203E-14,
    3.3124049707288003E-17, -8.6720759862275152, -4.52268961995504E-14,
    3.1391373691544004E-17, -8.6677652514157444, -4.4141884434751213E-14,
    2.9657685598688008E-17, -8.6634570321898057, -4.3057506068530407E-14,
    2.7925009582944E-17, -8.6591488128970813, -4.1973127702309613E-14,
    2.6192333567200004E-17, -8.6548380781418146, -4.088811593751041E-14,
    2.4458645474344002E-17, -8.65052985898094, -3.9803737571289604E-14,
    2.2725969458600003E-17, -208.64622163976065, -3.871935920506881E-14,
    2.0993293442856E-17, -8.4128478271310385, -3.7634980838848003E-14,
    1.9260617427112002E-17, -8.1794725638285666, -3.6550222433480169E-14,
    1.7527334165100803E-17, -7.9460984240437309, -3.5465717387543685E-14,
    1.57944557339344E-17, -7.7127242086144783, -3.4381149001749368E-14,
    1.4061476095056803E-17, -7.479350479281651, -3.3296643955812883E-14,
    1.2328597663890403E-17, -7.24597679157641, -3.2212075570018566E-14,
    1.0595618025012802E-17, -7.0126035311676285, -3.11275261861816E-14,
    8.862668748448561E-18, -4.8296856669233612E-14, -8.6979278173984422,
    6.9769321262080174E-16, -4.7498848412132806E-14, -8.69361959815538,
    7.01156722830645E-16, -4.6700840155032006E-14, -8.68931137884382,
    7.04620233040488E-16, -4.5902365771612804E-14, -8.6850006440594463,
    7.0808576632872487E-16, -4.5104357514512011E-14, -8.6806924248506281,
    7.1154927653856811E-16, -4.43063492574112E-14, -8.6763842055733189,
    7.1501278674841124E-16, -4.3508341000310405E-14, -8.6720759862275152,
    7.1847629695825457E-16, -4.2709866616891204E-14, -8.6677652514157444,
    7.2194183024649134E-16, -4.191185835979041E-14, -8.6634570321898057,
    7.2540534045633448E-16, -4.1113850102689611E-14, -8.6591488128970813,
    7.2886885066617781E-16, -4.03153757192704E-14, -8.6548380781418146,
    7.3233438395441458E-16, -3.95173674621696E-14, -8.65052985898094,
    7.3579789416425761E-16, -3.8719359205068804E-14, -208.64622163976065,
    7.3926140437410094E-16, -3.7921350947968004E-14, -8.4128478271310385,
    7.4272491458394408E-16, -3.7123063015076168E-14, -8.1794725638285666,
    7.4618963864082339E-16, -3.6324961532711686E-14, -7.9460984240437309,
    7.4965355346634545E-16, -3.5526813437715365E-14, -7.7127242086144783,
    7.5311767059970668E-16, -3.4728711955350888E-14, -7.479350479281651,
    7.5658158542522854E-16, -3.3930563860354568E-14, -7.24597679157641,
    7.6004570255858986E-16, -3.3132429749147796E-14, -7.0126035311676285,
    7.6350975899959934E-16, 2.8982695228536006E-17, 7.5136889676821369E-16,
    -8.6979278173984422, 2.8317008955256005E-17, 7.5036008638330479E-16,
    -8.69361959815538, 2.7651322681976004E-17, 7.4935127599839617E-16,
    -8.68931137884382, 2.6985247573256007E-17, 7.48341876355085E-16,
    -8.6850006440594463, 2.6319561299976003E-17, 7.4733306597017614E-16,
    -8.6806924248506281, 2.5653875026696002E-17, 7.4632425558526724E-16,
    -8.6763842055733189, 2.4988188753416002E-17, 7.4531544520035852E-16,
    -8.6720759862275152, 2.4322113644696E-17, 7.4430604555704731E-16,
    -8.6677652514157444, 2.3656427371416003E-17, 7.432972351721384E-16,
    -8.6634570321898057, 2.2990741098136E-17, 7.4228842478722978E-16,
    -8.6591488128970813, 2.2324665989416E-17, 7.4127902514391857E-16,
    -8.6548380781418146, 2.1658979716136E-17, 7.4027021475900966E-16,
    -8.65052985898094, 2.0993293442856004E-17, 7.3926140437410094E-16,
    -208.64622163976065, 2.0327607169576003E-17, 7.3825259398919213E-16,
    -8.4128478271310385, 1.9661687595032002E-17, 7.3724343004924188E-16,
    -8.1794725638285666, 1.8995923554664002E-17, 7.3623450181265259E-16,
    -7.9460984240437309, 1.8330120630752E-17, 7.3522551465022315E-16,
    -7.7127242086144783, 1.7664356590384E-17, 7.3421658641363376E-16,
    -7.479350479281651, 1.6998553666472003E-17, 7.3320759925120423E-16,
    -7.24597679157641, 1.6332762407623202E-17, 7.3219862976652685E-16,
    -7.0126035311676285, -8.4631379685742587, -5.0583984686596006E-14,
    4.0529398995088E-17, -8.458947735182857, -4.9528919277508011E-14,
    3.8846162142736007E-17, -8.4547575017229679, -4.847385386842E-14,
    3.7162925290384E-17, -8.4505648214674611, -4.7418172182808003E-14,
    3.5478705238936007E-17, -8.4463745881103076, -4.636310677372E-14,
    3.3795468386584E-17, -8.4421843546846613, -4.5308041364632006E-14,
    3.2112231534232E-17, -8.4379941211905258, -4.4252975955544E-14,
    3.0428994681880005E-17, -8.4338014409007762, -4.3197294269932011E-14,
    2.8744774630432005E-17, -8.4296112075162242, -4.2142228860844009E-14,
    2.7061537778080005E-17, -8.4254209740734574, -4.1087163451756013E-14,
    2.5378300925728004E-17, -8.4212282938367871, -4.0031481766144008E-14,
    2.369408087428E-17, -8.4170380605087409, -3.8976416357056006E-14,
    2.2010844021928004E-17, -8.4128478271310385, -3.792135094796801E-14,
    2.0327607169576003E-17, -208.40865759369393, -3.6866285538880008E-14,
    1.8644370317224003E-17, -8.17540029881121, -3.5810850363877604E-14,
    1.6960543545414403E-17, -7.9421440413329165, -3.4755661699484808E-14,
    1.52771100532432E-17, -7.7088876564419255, -3.3700411407439608E-14,
    1.3593578241162402E-17, -7.4756316921436365, -3.2645222743046805E-14,
    1.1910144748991203E-17, -7.2423757177046557, -3.1589972451001605E-14,
    1.0226612936910401E-17, -7.0091201098594791, -3.0534740647252126E-14,
    8.5431106208024819E-18, -4.6860671451256011E-14, -8.4631379685742587,
    6.9625809163338408E-16, -4.6091976151288E-14, -8.458947735182857,
    6.9975712253975216E-16, -4.5323280851320008E-14, -8.4547575017229679,
    7.0325615344612E-16, -4.455413654708801E-14, -8.4505648214674611,
    7.06757228178952E-16, -4.3785441247120008E-14, -8.4463745881103076,
    7.1025625908532013E-16, -4.3016745947152E-14, -8.4421843546846613,
    7.1375528999168811E-16, -4.2248050647184005E-14, -8.4379941211905258,
    7.172543208980562E-16, -4.1478906342952E-14, -8.4338014409007762,
    7.207553956308881E-16, -4.0710211042984006E-14, -8.4296112075162242,
    7.2425442653725608E-16, -3.9941515743016011E-14, -8.4254209740734574,
    7.2775345744362417E-16, -3.9172371438784007E-14, -8.4212282938367871,
    7.3125453217645617E-16, -3.8403676138816E-14, -8.4170380605087409,
    7.3475356308282405E-16, -3.7634980838848003E-14, -8.4128478271310385,
    7.3825259398919213E-16, -3.6866285538880008E-14, -208.40865759369393,
    7.4175162489556012E-16, -3.6097320836353609E-14, -8.17540029881121,
    7.4525188209780645E-16, -3.5328535735532808E-14, -7.9421440413329165,
    7.4875132176946732E-16, -3.4559705734285604E-14, -7.7088876564419255,
    7.5225096582377448E-16, -3.3790920633464804E-14, -7.4756316921436365,
    7.5575040549543525E-16, -3.3022090632217606E-14, -7.2423757177046557,
    7.592500495497425E-16, -3.2253274101098324E-14, -7.0091201098594791,
    7.627496322892558E-16, 2.6656662618040004E-17, 7.5440609637554808E-16,
    -8.4631379685742587, 2.6040415508152005E-17, 7.53432806687164E-16,
    -8.458947735182857, 2.5424168398264003E-17, 7.5245951699878016E-16,
    -8.4547575017229679, 2.4807561330952008E-17, 7.5148565880006417E-16,
    -8.4505648214674611, 2.4191314221064002E-17, 7.5051236911168011E-16,
    -8.4463745881103076, 2.3575067111176E-17, 7.4953907942329605E-16,
    -8.4421843546846613, 2.2958820001287998E-17, 7.4856578973491209E-16,
    -8.4379941211905258, 2.2342212933976003E-17, 7.4759193153619611E-16,
    -8.4338014409007762, 2.1725965824088003E-17, 7.46618641847812E-16,
    -8.4296112075162242, 2.11097187142E-17, 7.4564535215942818E-16,
    -8.4254209740734574, 2.0493111646888003E-17, 7.446714939607122E-16,
    -8.4212282938367871, 1.9876864537E-17, 7.43698204272328E-16,
    -8.4170380605087409, 1.9260617427112005E-17, 7.4272491458394408E-16,
    -8.4128478271310385, 1.8644370317224003E-17, 7.4175162489556012E-16,
    -208.40865759369393, 1.8027907232881604E-17, 7.40777994100977E-16,
    -8.17540029881121, 1.7411588131508802E-17, 7.398045907105265E-16,
    -7.9421440413329165, 1.67952330343936E-17, 7.38831130469043E-16,
    -7.7088876564419255, 1.6178913933020802E-17, 7.3785772707859251E-16,
    -7.4756316921436365, 1.5562558835905605E-17, 7.3688426683710882E-16,
    -7.2423757177046557, 1.4946214537513123E-17, 7.3591082365093535E-16,
    -7.0091201098594791, -8.2283468781897078, -4.9433766605153528E-14,
    3.9269914873714405E-17, -8.22427461321002, -4.8408024426430969E-14,
    3.7636134511564006E-17, -8.220202348178967, -4.73822822477084E-14,
    3.60023541494136E-17, -8.2161277050364578, -4.6355940920516965E-14,
    3.4367619476304E-17, -8.2120554400910155, -4.53301987417944E-14,
    3.2733839114153604E-17, -8.2079831750942063, -4.4304456563071847E-14,
    3.1100058752003206E-17, -8.2039109100289078, -4.3278714384349281E-14,
    2.94662783898528E-17, -8.19983626685215, -4.2252373057157849E-14,
    2.7831543716743208E-17, -8.1957640018724653, -4.122663087843529E-14,
    2.6197763354592803E-17, -8.1916917368482611, -4.0200888699712731E-14,
    2.45639829924424E-17, -8.1876170937228743, -3.9174547372521287E-14,
    2.29292483193328E-17, -8.1835448287962667, -3.8148805193798727E-14,
    2.12954679571824E-17, -8.1794725638285666, -3.7123063015076168E-14,
    1.9661687595032002E-17, -8.17540029881121, -3.6097320836353609E-14,
    1.80279072328816E-17, -208.17132660690865, -3.5071219168549716E-14,
    1.6393554284155682E-17, -7.9381882729475333, -3.4045357160133381E-14,
    1.475958305981344E-17, -7.7050497598066379, -3.3019435236870161E-14,
    1.3125516404375281E-17, -7.4719116017523151, -3.199357322845382E-14,
    1.1491545180033043E-17, -7.238773381791427, -3.0967651305190594E-14,
    9.8574785245948825E-18, -7.0056354677261057, -2.9941747356381442E-14,
    8.223440498485498E-18, -4.5423982897337534E-14, -8.2283468781897078,
    6.9482246768300367E-16, -4.4684610827734963E-14, -8.22427461321002,
    6.9835703173473877E-16, -4.3945238758132404E-14, -8.220202348178967,
    7.0189159578647368E-16, -4.3205434812320967E-14, -8.2161277050364578,
    7.054282244200146E-16, -4.2466062742718408E-14, -8.2120554400910155,
    7.0896278847174971E-16, -4.1726690673115843E-14, -8.2079831750942063,
    7.1249735252348471E-16, -4.0987318603513284E-14, -8.2039109100289078,
    7.1603191657521981E-16, -4.0247514657701841E-14, -8.19983626685215,
    7.1956854520876074E-16, -3.9508142588099289E-14, -8.1957640018724653,
    7.2310310926049574E-16, -3.876877051849673E-14, -8.1916917368482611,
    7.2663767331223094E-16, -3.8028966572685281E-14, -8.1876170937228743,
    7.3017430194577187E-16, -3.7289594503082722E-14, -8.1835448287962667,
    7.3370886599750678E-16, -3.6550222433480163E-14, -8.1794725638285666,
    7.3724343004924188E-16, -3.5810850363877604E-14, -8.17540029881121,
    7.4077799410097688E-16, -3.5071219168549722E-14, -208.17132660690865,
    7.443137969017954E-16, -3.4331760723705383E-14, -7.9381882729475333,
    7.4784877386989171E-16, -3.3592259091240159E-14, -7.7050497598066379,
    7.5138395729616857E-16, -3.285280064639582E-14, -7.4719116017523151,
    7.5491893426426478E-16, -3.2113299013930596E-14, -7.238773381791427,
    7.5845411769054164E-16, -3.1373810337751641E-14, -7.0056354677261057,
    7.6198923917936419E-16, 2.4329814809199204E-17, 7.5744436042199724E-16,
    -8.2283468781897078, 2.3763024189512807E-17, 7.565066038789802E-16,
    -8.22427461321002, 2.3196233569826403E-17, 7.5556884733596336E-16,
    -8.220202348178967, 2.2629111880852806E-17, 7.546305430379563E-16,
    -8.2161277050364578, 2.2062321261166402E-17, 7.5369278649493936E-16,
    -8.2120554400910155, 2.1495530641480002E-17, 7.5275502995192232E-16,
    -8.2079831750942063, 2.0928740021793602E-17, 7.5181727340890537E-16,
    -8.2039109100289078, 2.0361618332820005E-17, 7.5087896911089832E-16,
    -8.19983626685215, 1.9794827713133604E-17, 7.4994121256788128E-16,
    -8.1957640018724653, 1.92280370934472E-17, 7.4900345602486453E-16,
    -8.1916917368482611, 1.86609154044736E-17, 7.4806515172685747E-16,
    -8.1876170937228743, 1.80941247847872E-17, 7.4712739518384043E-16,
    -8.1835448287962667, 1.7527334165100803E-17, 7.4618963864082349E-16,
    -8.1794725638285666, 1.69605435454144E-17, 7.4525188209780655E-16,
    -8.17540029881121, 1.6393554284155682E-17, 7.443137969017955E-16,
    -208.17132660690865, 1.582669745061184E-17, 7.4337593080778046E-16,
    -7.9381882729475333, 1.5259807510139283E-17, 7.4243800993826666E-16,
    -7.7050497598066379, 1.4692950676595441E-17, 7.4150014384425161E-16,
    -7.4719116017523151, 1.4126060736122884E-17, 7.4056222297473752E-16,
    -7.238773381791427, 1.3559180727728938E-17, 7.3962431853787332E-16,
    -7.0056354677261057, -7.9935579447306031, -4.8283817172100966E-14,
    3.80107249212352E-17, -7.9896035620420474, -4.7287391374922089E-14,
    3.6426389498080005E-17, -7.9856491793192479, -4.62909655777432E-14,
    3.48420540749248E-17, -7.9816924871702239, -4.5293957756150084E-14,
    3.3256793222200003E-17, -7.9777381045159146, -4.42975319589712E-14,
    3.16724577990448E-17, -7.973783721827358, -4.3301106161792322E-14,
    3.00881223758896E-17, -7.9698293390874362, -4.2304680364613439E-14,
    2.8503786952734406E-17, -7.9658726469041667, -4.130767254302033E-14,
    2.6918526100009604E-17, -7.9619182642156119, -4.0311246745841446E-14,
    2.5334190676854404E-17, -7.9579638814928124, -3.9314820948662568E-14,
    2.37498552536992E-17, -7.9540071893506372, -3.8317813127069447E-14,
    2.2164594400974402E-17, -7.950052806713451, -3.7321387329890563E-14,
    2.0580258977819202E-17, -7.94609842404373, -3.6324961532711686E-14,
    1.8995923554664002E-17, -7.9421440413329156, -3.5328535735532808E-14,
    1.7411588131508802E-17, -7.9381882729475315, -3.4331760723705383E-14,
    1.5826697450611843E-17, -207.93423342836567, -3.3335218521643654E-14,
    1.424217694154272E-17, -7.70121275949525, -3.2338618117140506E-14,
    1.2657563889516642E-17, -7.4681923802056644, -3.1342075915078771E-14,
    1.1073043380447522E-17, -7.2351718872420347, -3.0345475510575623E-14,
    9.4884303284214412E-18, -7.0021516394762147, -2.93488925668049E-14,
    7.90384503928245E-18, -4.3987629900712971E-14, -7.9935579447306031,
    6.93387179041265E-16, -4.3277574212654081E-14, -7.9896035620420474,
    6.96957267939139E-16, -4.2567518524595204E-14, -7.9856491793192479,
    7.0052735683701285E-16, -4.1857048084382084E-14, -7.9816924871702239,
    7.0409953106718694E-16, -4.1146992396323207E-14, -7.9777381045159146,
    7.0766961996506095E-16, -4.0436936708264323E-14, -7.973783721827358,
    7.1123970886293476E-16, -3.9726881020205446E-14, -7.9698293390874362,
    7.1480979776080877E-16, -3.901641057999232E-14, -7.9658726469041667,
    7.1838197199098286E-16, -3.8306354891933449E-14, -7.9619182642156119,
    7.2195206088885677E-16, -3.7596299203874572E-14, -7.9579638814928124,
    7.2552214978673078E-16, -3.6885828763661445E-14, -7.9540071893506372,
    7.2909432401690477E-16, -3.6175773075602562E-14, -7.950052806713451,
    7.3266441291477858E-16, -3.5465717387543685E-14, -7.94609842404373,
    7.3623450181265259E-16, -3.4755661699484808E-14, -7.9421440413329156,
    7.398045907105265E-16, -3.4045357160133381E-14, -7.9381882729475315,
    7.4337593080778046E-16, -3.3335218521643654E-14, -207.93423342836567,
    7.4694643677211448E-16, -3.26250384079385E-14, -7.70121275949525,
    7.5051715126967847E-16, -3.1914899769448775E-14, -7.4681923802056644,
    7.5408765723401239E-16, -3.1204719655743617E-14, -7.2351718872420347,
    7.5765837173157638E-16, -3.0494551984603094E-14, -7.0021516394762147,
    7.6122902366917136E-16, 2.2003510465921604E-17, 7.6048191484236985E-16,
    -7.9935579447306031, 2.1486164785230405E-17, 7.5957968314549171E-16,
    -7.9896035620420474, 2.0968819104539204E-17, 7.5867745144861378E-16,
    -7.9856491793192479, 2.0451171235950406E-17, 7.5777469274723978E-16,
    -7.9816924871702239, 1.99338255552592E-17, 7.5687246105036175E-16,
    -7.9777381045159146, 1.9416479874568E-17, 7.5597022935348362E-16,
    -7.973783721827358, 1.88991341938768E-17, 7.5506799765660558E-16,
    -7.9698293390874362, 1.8381486325288003E-17, 7.5416523895523159E-16,
    -7.9658726469041667, 1.7864140644596802E-17, 7.5326300725835345E-16,
    -7.9619182642156119, 1.73467949639056E-17, 7.5236077556147552E-16,
    -7.9579638814928124, 1.6829147095316802E-17, 7.5145801686010162E-16,
    -7.9540071893506372, 1.63118014146256E-17, 7.5055578516322349E-16,
    -7.950052806713451, 1.5794455733934402E-17, 7.4965355346634535E-16,
    -7.94609842404373, 1.52771100532432E-17, 7.4875132176946732E-16,
    -7.9421440413329156, 1.475958305981344E-17, 7.4784877386989181E-16,
    -7.9381882729475315, 1.4242176941542723E-17, 7.4694643677211448E-16,
    -207.93423342836567, 1.372474060448224E-17, 7.4604404697388771E-16,
    -7.70121275949525, 1.3207334486211522E-17, 7.4514170987611038E-16,
    -7.4681923802056644, 1.2689898149151042E-17, 7.4423932007788351E-16,
    -7.2351718872420347, 1.217247087772749E-17, 7.4333694608979164E-16,
    -7.0021516394762147, -7.7587695568412611, -4.7133800576950924E-14,
    3.6751461426532404E-17, -7.7549330046841183, -4.616669287352117E-14,
    3.5216573830174007E-17, -7.75109645249273, -4.5199585170091404E-14,
    3.36816862338156E-17, -7.7472576595643634, -4.4231912567302164E-14,
    3.2145902090964004E-17, -7.7434211074414661, -4.3264804863872404E-14,
    3.06110144946056E-17, -7.7395845552843223, -4.2297697160442644E-14,
    2.9076126898247204E-17, -7.735748003092934, -4.1330589457012884E-14,
    2.7541239301888804E-17, -7.7319092101474443, -4.0362916854223651E-14,
    2.6005455159037207E-17, -7.7280726579903023, -3.9395809150793891E-14,
    2.4470567562678804E-17, -7.724236105798914, -3.8428701447364131E-14,
    2.29356799663204E-17, -7.7203973128705465, -3.7461028844574891E-14,
    2.13998958234688E-17, -7.7165607607544988, -3.6493921141145125E-14,
    1.98650082271104E-17, -7.7127242086144792, -3.5526813437715371E-14,
    1.8330120630752E-17, -7.7088876564419255, -3.4559705734285604E-14,
    1.67952330343936E-17, -7.705049759806637, -3.3592259091240159E-14,
    1.5259807510139283E-17, -7.7012127594952506, -3.26250384079385E-14,
    1.372474060448224E-17, -207.69737553510555, -3.16577612347009E-14,
    1.2189584044175883E-17, -7.4644729414504143, -3.0690540551399242E-14,
    1.0654517138518842E-17, -7.231570182352967, -2.972326337816164E-14,
    9.1193605782124829E-18, -6.9986676077554533, -2.8756003151904819E-14,
    7.5842309143009172E-18, -4.2551193014764932E-14, -7.7587695568412611,
    6.91951806572366E-16, -4.1870455420455166E-14, -7.7549330046841183,
    6.95557422391186E-16, -4.1189717826145406E-14, -7.75109645249273,
    6.9916303821000563E-16, -4.0508582604736168E-14, -7.7472576595643634,
    7.027707601128318E-16, -3.9827845010426408E-14, -7.7434211074414661,
    7.0637637593165174E-16, -3.9147107416116642E-14, -7.7395845552843223,
    7.0998199175047149E-16, -3.8466369821806882E-14, -7.735748003092934,
    7.1358760756929143E-16, -3.7785234600397645E-14, -7.7319092101474443,
    7.171953294721176E-16, -3.7104497006087885E-14, -7.7280726579903023,
    7.2080094529093735E-16, -3.6423759411778131E-14, -7.724236105798914,
    7.2440656110975729E-16, -3.5742624190368881E-14, -7.7203973128705465,
    7.2801428301258346E-16, -3.5061886596059121E-14, -7.7165607607544988,
    7.3161989883140311E-16, -3.4381149001749361E-14, -7.7127242086144792,
    7.3522551465022306E-16, -3.3700411407439608E-14, -7.7088876564419255,
    7.388311304690429E-16, -3.3019435236870161E-14, -7.705049759806637,
    7.4243800993826646E-16, -3.2338618117140506E-14, -7.7012127594952506,
    7.4604404697388771E-16, -3.16577612347009E-14, -207.69737553510555,
    7.4965029461790941E-16, -3.0976944114971244E-14, -7.4644729414504143,
    7.5325633165353046E-16, -3.0296087232531638E-14, -7.231570182352967,
    7.5686257929755226E-16, -2.9615242278905019E-14, -6.9986676077554533,
    7.6046876375905379E-16, 1.9677070256253204E-17, 7.6351964666926162E-16,
    -7.7587695568412611, 1.9209172402358806E-17, 7.6265294189332932E-16,
    -7.7549330046841183, 1.8741274548464405E-17, 7.6178623711739732E-16,
    -7.75109645249273, 1.8273103389748805E-17, 7.609190260886755E-16,
    -7.7472576595643634, 1.7805205535854404E-17, 7.6005232131274331E-16,
    -7.7434211074414661, 1.7337307681960003E-17, 7.5918561653681111E-16,
    -7.7395845552843223, 1.68694098280656E-17, 7.58318911760879E-16,
    -7.735748003092934, 1.6401238669350002E-17, 7.5745170073215719E-16,
    -7.7319092101474443, 1.5933340815455604E-17, 7.5658499595622489E-16,
    -7.7280726579903023, 1.54654429615612E-17, 7.5571829118029289E-16,
    -7.724236105798914, 1.49972718028456E-17, 7.5485108015157108E-16,
    -7.7203973128705465, 1.4529373948951202E-17, 7.5398437537563878E-16,
    -7.7165607607544988, 1.4061476095056803E-17, 7.5311767059970668E-16,
    -7.7127242086144792, 1.3593578241162402E-17, 7.5225096582377458E-16,
    -7.7088876564419255, 1.3125516404375281E-17, 7.5138395729616857E-16,
    -7.705049759806637, 1.2657563889516642E-17, 7.5051715126967837E-16,
    -7.7012127594952506, 1.2189584044175881E-17, 7.4965029461790951E-16,
    -207.69737553510555, 1.1721631529317242E-17, 7.4878348859141931E-16,
    -7.4644729414504143, 1.1253651683976483E-17, 7.4791663193965016E-16,
    -7.231570182352967, 1.0785680037780358E-17, 7.4704979047546493E-16,
    -6.9986676077554533, -7.5239824411033318, -4.5983851143898369E-14,
    3.5492271474053207E-17, -7.5202636539824406, -4.504605982201229E-14,
    3.4006828816690006E-17, -7.5165448668273047, -4.41082685001262E-14,
    3.25213861593268E-17, -7.5128239076032974, -4.3169929402935283E-14,
    3.1035075836860004E-17, -7.5091051204995294, -4.22321380810492E-14,
    2.95496331794968E-17, -7.5053863333786373, -4.1294346759163126E-14,
    2.80641905221336E-17, -7.5016675462235014, -4.0356555437277041E-14,
    2.6578747864770404E-17, -7.497946586999495, -3.9418216340086131E-14,
    2.5092437542303606E-17, -7.4942277998786038, -3.8480425018200046E-14,
    2.3606994884940402E-17, -7.4905090127234679, -3.7542633696313974E-14,
    2.2121552227577203E-17, -7.4867880534994615, -3.6604294599123045E-14,
    2.0635241905110403E-17, -7.4830692663956935, -3.5666503277236967E-14,
    1.9149799247747202E-17, -7.479350479281651, -3.4728711955350888E-14,
    1.7664356590384E-17, -7.4756316921436365, -3.3790920633464804E-14,
    1.6178913933020802E-17, -7.4719116017523151, -3.285280064639582E-14,
    1.4692950676595441E-17, -7.4681923802056653, -3.1914899769448775E-14,
    1.320733448621152E-17, -7.4644729414504143, -3.0976944114971244E-14,
    1.1721631529317242E-17, -207.46075371990256, -3.0039043238024193E-14,
    1.0236015338933323E-17, -7.2279686878035747, -2.9101087583546662E-14,
    8.7503123820390416E-18, -6.9951837795055631, -2.8163148362328274E-14,
    7.2646354550978688E-18, -4.1114840018140375E-14, -7.5239824411033318,
    6.905165179306275E-16, -4.0463418805374284E-14, -7.5202636539824406,
    6.9415765859558625E-16, -3.9811997592608206E-14, -7.5165448668273047,
    6.977987992605449E-16, -3.9160195876797285E-14, -7.5128239076032974,
    7.0144206676000413E-16, -3.8508774664031213E-14, -7.5091051204995294,
    7.0508320742496289E-16, -3.7857353451265122E-14, -7.5053863333786373,
    7.0872434808992164E-16, -3.7205932238499044E-14, -7.5016675462235014,
    7.1236548875488039E-16, -3.6554130522688123E-14, -7.497946586999495,
    7.1600875625433962E-16, -3.5902709309922045E-14, -7.4942277998786038,
    7.1964989691929837E-16, -3.5251288097155967E-14, -7.4905090127234679,
    7.2329103758425712E-16, -3.4599486381345046E-14, -7.4867880534994615,
    7.2693430508371636E-16, -3.3948065168578961E-14, -7.4830692663956935,
    7.30575445748675E-16, -3.3296643955812883E-14, -7.479350479281651,
    7.3421658641363376E-16, -3.2645222743046805E-14, -7.4756316921436365,
    7.3785772707859251E-16, -3.199357322845382E-14, -7.4719116017523151,
    7.4150014384425152E-16, -3.1342075915078777E-14, -7.4681923802056653,
    7.4514170987611038E-16, -3.0690540551399242E-14, -7.4644729414504143,
    7.4878348859141931E-16, -3.0039043238024193E-14, -207.46075371990256,
    7.5242505462327808E-16, -2.9387507874344658E-14, -7.2279686878035747,
    7.56066833338587E-16, -2.8735983925756472E-14, -6.9951837795055631,
    7.5970854824886086E-16, 1.7350765912975604E-17, 7.6655720108963423E-16,
    -7.5239824411033318, 1.6932312998076405E-17, 7.6572602115984094E-16,
    -7.5202636539824406, 1.6513860083177205E-17, 7.6489484123004774E-16,
    -7.5165448668273047, 1.6095162744846405E-17, 7.64063175797959E-16,
    -7.5128239076032974, 1.5676709829947203E-17, 7.632319958681658E-16,
    -7.5091051204995294, 1.5258256915048003E-17, 7.6240081593837241E-16,
    -7.5053863333786373, 1.48398040001488E-17, 7.6156963600857911E-16,
    -7.5016675462235014, 1.4421106661818004E-17, 7.6073797057649046E-16,
    -7.497946586999495, 1.4002653746918804E-17, 7.5990679064669707E-16,
    -7.4942277998786038, 1.3584200832019601E-17, 7.5907561071690388E-16,
    -7.4905090127234679, 1.3165503493688802E-17, 7.5824394528481522E-16,
    -7.4867880534994615, 1.2747050578789602E-17, 7.5741276535502183E-16,
    -7.4830692663956935, 1.2328597663890403E-17, 7.5658158542522864E-16,
    -7.479350479281651, 1.1910144748991203E-17, 7.5575040549543534E-16,
    -7.4756316921436365, 1.1491545180033043E-17, 7.5491893426426478E-16,
    -7.4719116017523151, 1.1073043380447522E-17, 7.5408765723401239E-16,
    -7.4681923802056653, 1.0654517138518842E-17, 7.5325633165353056E-16,
    -7.4644729414504143, 1.0236015338933321E-17, 7.5242505462327808E-16,
    -207.46075371990256, 9.8174890970046427E-18, 7.5159372904279615E-16,
    -7.2279686878035747, 9.39897018777891E-18, 7.5076241802738315E-16,
    -6.9951837795055631, -7.289195988207168, -4.4833834548748326E-14,
    3.4233007979350409E-17, -7.2855949143456868, -4.3925361320611371E-14,
    3.2797013148784008E-17, -7.2819938404670825, -4.30168880924744E-14,
    3.13610183182176E-17, -7.27839066319173, -4.2107884214087363E-14,
    2.9924184705624005E-17, -7.27478958933025, -4.11994109859504E-14,
    2.8488189875057604E-17, -7.2711885154687668, -4.0290937757813441E-14,
    2.7052195044491203E-17, -7.2675874415901633, -3.938246452967648E-14,
    2.5616200213924802E-17, -7.2639842643148107, -3.8473460651289452E-14,
    2.4179366601331206E-17, -7.26038319045333, -3.7564987423152485E-14,
    2.2743371770764802E-17, -7.2567821165747262, -3.665651419501553E-14,
    2.1307376940198404E-17, -7.2531789392993735, -3.574751031662849E-14,
    1.9870543327604802E-17, -7.2495778654378924, -3.4839037088491522E-14,
    1.84345484970384E-17, -7.24597679157641, -3.3930563860354568E-14,
    1.6998553666472003E-17, -7.2423757177046557, -3.3022090632217606E-14,
    1.5562558835905602E-17, -7.238773381791427, -3.2113299013930596E-14,
    1.4126060736122884E-17, -7.2351718872420347, -3.1204719655743623E-14,
    1.2689898149151041E-17, -7.231570182352967, -3.0296087232531638E-14,
    1.1253651683976483E-17, -7.2279686878035747, -2.9387507874344658E-14,
    9.8174890970046427E-18, -207.2243669829133, -2.8478875451132679E-14,
    8.3812426318300817E-18, -6.9916997477848017, -2.7570258947428195E-14,
    6.945021330116338E-18, -3.9678403132192329E-14, -7.289195988207168,
    6.8908114546172841E-16, -3.9056300013175362E-14, -7.2855949143456868,
    6.927578130476331E-16, -3.8434196894158408E-14, -7.2819938404670825,
    6.9643448063353768E-16, -3.7811730397151363E-14, -7.27839066319173,
    7.00113295805649E-16, -3.7189627278134408E-14, -7.27478958933025,
    7.0378996339155368E-16, -3.6567524159117441E-14, -7.2711885154687668,
    7.0746663097745827E-16, -3.5945421040100487E-14, -7.2675874415901633,
    7.1114329856336306E-16, -3.5322954543093442E-14, -7.2639842643148107,
    7.1482211373547437E-16, -3.4700851424076487E-14, -7.26038319045333,
    7.18498781321379E-16, -3.4078748305059533E-14, -7.2567821165747262,
    7.2217544890728374E-16, -3.3456281808052481E-14, -7.2531789392993735,
    7.2585426407939506E-16, -3.2834178689035521E-14, -7.2495778654378924,
    7.2953093166529954E-16, -3.2212075570018566E-14, -7.24597679157641,
    7.3320759925120423E-16, -3.1589972451001605E-14, -7.2423757177046557,
    7.3688426683710892E-16, -3.09676513051906E-14, -7.238773381791427,
    7.4056222297473752E-16, -3.0345475510575623E-14, -7.2351718872420347,
    7.4423932007788361E-16, -2.972326337816164E-14, -7.231570182352967,
    7.4791663193965026E-16, -2.9101087583546662E-14, -7.2279686878035747,
    7.5159372904279615E-16, -2.8478875451132679E-14, -207.2243669829133,
    7.5527104090456289E-16, -2.7856674220058395E-14, -6.9916997477848017,
    7.5894828833874329E-16, 1.5024325703307205E-17, 7.69594932916526E-16,
    -7.289195988207168, 1.4655320615204805E-17, 7.6879927990767854E-16,
    -7.2855949143456868, 1.4286315527102403E-17, 7.6800362689883138E-16,
    -7.2819938404670825, 1.3917094898644805E-17, 7.6720750913939471E-16,
    -7.27839066319173, 1.3548089810542402E-17, 7.6641185613054735E-16,
    -7.27478958933025, 1.3179084722440001E-17, 7.656162031216999E-16,
    -7.2711885154687668, 1.2810079634337602E-17, 7.6482055011285264E-16,
    -7.2675874415901633, 1.2440859005880003E-17, 7.64024432353416E-16,
    -7.2639842643148107, 1.2071853917777603E-17, 7.6322877934456851E-16,
    -7.26038319045333, 1.17028488296752E-17, 7.6243312633572125E-16,
    -7.2567821165747262, 1.1333628201217601E-17, 7.6163700857628468E-16,
    -7.2531789392993735, 1.0964623113115202E-17, 7.6084135556743722E-16,
    -7.2495778654378924, 1.0595618025012802E-17, 7.6004570255858986E-16,
    -7.24597679157641, 1.0226612936910401E-17, 7.592500495497425E-16,
    -7.2423757177046557, 9.857478524594881E-18, 7.5845411769054164E-16,
    -7.238773381791427, 9.4884303284214412E-18, 7.5765837173157638E-16,
    -7.2351718872420347, 9.1193605782124813E-18, 7.5686257929755226E-16,
    -7.231570182352967, 8.7503123820390416E-18, 7.56066833338587E-16,
    -7.2279686878035747, 8.3812426318300817E-18, 7.5527104090456279E-16,
    -207.2243669829133, 8.0121793478317786E-18, 7.5447526241305643E-16,
    -6.9916997477848017, -7.05441069099154, -4.3683838102227532E-14,
    3.2973766547314689E-17, -7.0509272696833909, -4.2804682454178059E-14,
    3.1587218677204607E-17, -7.0474438483752424, -4.192552680612858E-14,
    3.0200670807094518E-17, -7.0439583923583857, -4.1045857632583755E-14,
    2.8813313037527604E-17, -7.0404749710502372, -4.0166701984534282E-14,
    2.7426765167417521E-17, -7.036991549742087, -3.928754633648481E-14,
    2.6040217297307442E-17, -7.0335081284339385, -3.8408390688435337E-14,
    2.4653669427197363E-17, -7.0300226724170809, -3.7528721514890518E-14,
    2.3266311657630445E-17, -7.0265392511089324, -3.6649565866841046E-14,
    2.1879763787520363E-17, -7.0230558298007839, -3.5770410218791573E-14,
    2.0493215917410281E-17, -7.0195703737839272, -3.4890741045246742E-14,
    1.910585814784336E-17, -7.0160869524757787, -3.4011585397197269E-14,
    1.7719310277733281E-17, -7.0126035311676294, -3.31324297491478E-14,
    1.6332762407623202E-17, -7.00912010985948, -3.2253274101098324E-14,
    1.4946214537513123E-17, -7.0056354677261057, -3.1373810337751641E-14,
    1.3559180727728938E-17, -7.0021516394762164, -3.0494551984603094E-14,
    1.217247087772749E-17, -6.9986676077554542, -2.9615242278905019E-14,
    1.0785680037780358E-17, -6.9951837795055631, -2.8735983925756472E-14,
    9.39897018777891E-18, -6.9916997477848026, -2.7856674220058395E-14,
    8.0121793478317786E-18, -206.9882157771053, -2.6977379920125177E-14,
    6.625412804868351E-18, -3.8241991413041339E-14, -7.05441069099154,
    6.876457981409775E-16, -3.7649205874111854E-14, -7.0509272696833909,
    6.91357992025386E-16, -3.7056420335182388E-14, -7.0474438483752424,
    6.950701859097944E-16, -3.6463288543017558E-14, -7.0439583923583857,
    6.9878454813175211E-16, -3.5870503004088092E-14, -7.0404749710502372,
    7.0249674201616062E-16, -3.5277717465158607E-14, -7.036991549742087,
    7.06208935900569E-16, -3.4684931926229141E-14, -7.0335081284339385,
    7.0992112978497752E-16, -3.4091800134064311E-14, -7.0300226724170809,
    7.1363549200693533E-16, -3.3499014595134845E-14, -7.0265392511089324,
    7.1734768589134373E-16, -3.2906229056205373E-14, -7.0230558298007839,
    7.2105987977575223E-16, -3.2313097264040537E-14, -7.0195703737839272,
    7.2477424199770994E-16, -3.1720311725111065E-14, -7.0160869524757787,
    7.2848643588211825E-16, -3.1127526186181592E-14, -7.0126035311676294,
    7.3219862976652685E-16, -3.0534740647252126E-14, -7.00912010985948,
    7.3591082365093525E-16, -2.9941747356381442E-14, -7.0056354677261057,
    7.3962431853787322E-16, -2.93488925668049E-14, -7.0021516394762164,
    7.4333694608979164E-16, -2.8756003151904819E-14, -6.9986676077554542,
    7.4704979047546483E-16, -2.8163148362328274E-14, -6.9951837795055631,
    7.5076241802738315E-16, -2.7570258947428195E-14, -6.9916997477848026,
    7.5447526241305643E-16, -2.6977379920125177E-14, -206.9882157771053,
    7.5818804174860314E-16, 1.2697926253556044E-17, 7.72632611521462E-16,
    -7.05441069099154, 1.2378368125909965E-17, 7.718724848111184E-16,
    -7.0509272696833909, 1.2058809998263884E-17, 7.71112358100775E-16,
    -7.0474438483752424, 1.1739065212832965E-17, 7.7035178739118468E-16,
    -7.0439583923583857, 1.1419507085186881E-17, 7.6959166068084124E-16,
    -7.0404749710502372, 1.1099948957540802E-17, 7.688315339704976E-16,
    -7.036991549742087, 1.0780390829894721E-17, 7.6807140726015406E-16,
    -7.0335081284339385, 1.0460646044463803E-17, 7.6731083655056378E-16,
    -7.0300226724170809, 1.0141087916817723E-17, 7.6655070984022014E-16,
    -7.0265392511089324, 9.82152978917164E-18, 7.657905831298768E-16,
    -7.0230558298007839, 9.5017850037407208E-18, 7.6503001242028652E-16,
    -7.0195703737839272, 9.1822268760946417E-18, 7.6426988570994288E-16,
    -7.0160869524757787, 8.862668748448561E-18, 7.6350975899959934E-16,
    -7.0126035311676294, 8.5431106208024819E-18, 7.627496322892558E-16,
    -7.00912010985948, 8.223440498485498E-18, 7.6198923917936429E-16,
    -7.0056354677261057, 7.90384503928245E-18, 7.6122902366917126E-16,
    -7.0021516394762164, 7.5842309143009172E-18, 7.6046876375905389E-16,
    -6.9986676077554542, 7.2646354550978688E-18, 7.5970854824886086E-16,
    -6.9951837795055631, 6.945021330116338E-18, 7.5894828833874319E-16,
    -6.9916997477848026, 6.6254128048683518E-18, 7.5818804174860324E-16,
    -206.9882157771053 };

  
   double dv[720] = { 1.0, 0.0, 0.0, 0.016, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.016, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.016, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0,
    0.0, 0.032, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.032, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.032, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.048, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.048, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.048, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0,
    0.064, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.064, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.064, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.08, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.08,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.08, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.096, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.096, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.096, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 1.0, 0.0, 0.0, 0.112, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.112, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.112, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.128, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.128, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.128, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 1.0, 0.0, 0.0, 0.144, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.144, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.144, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.16, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.16, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.16, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0,
    0.0, 0.0, 0.176, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.176, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.176, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.192, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.192, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.192, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0,
    0.0, 0.208, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.208, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.208, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.224, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.224, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.224, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0,
    0.24, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.24, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.24,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.256, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.256,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.256, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.272, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.272, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.272, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 1.0, 0.0, 0.0, 0.288, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.288, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.288, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.304, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.304, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.304, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 1.0, 0.0, 0.0, 0.32, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.32, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.32, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 };

  
   double B[9] = { 0.6, 0.0, 0.0, 0.0, 0.6, 0.0, 0.0, 0.0, 0.6};
  /* double B[9] = { 0.68, 0.0, 0.0, 0.0, 0.332, 0.0, 0.0, 0.0, 0.394
  }; */

   signed char A[400] = { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 };

  double IGA[7200];
  double A_data[3600];
  double a[3600];
  double b_del_lam[240];
  double result[240];
  double varargin_1_data[240];
  double Bineq[120];
  double GAMMA[120];
  double del_lam[120];
  double esig[120];
  double ftol[120];
  double igr2[120];
  double ilam[120];
  double lam[120];
  double mesil[120];
  double H[60];
  double del_z[60];
  double u_max[60];
  double absxk;
  double mu;
  double scale;
  double t;
  int a_tmp;
  int b_i;
  int b_i1;
  int exitflag;
  int i;
  int i1;
  int idx;
  int iter;
  int iy;
  int j2;
  short b_Aineq[7200];
  signed char Aineq[7200];
  signed char At[7200];
  unsigned char ii_data[240];
  boolean_T x[240];

  //  Set Constants - Extract variables from Structure MPCParams
  //  MPCParams=load('MPCParams.mat');
  // MPCParams.X;
  // MPCParams.A;
  // MPCParams.B;
  // MPCParams.P1;
  // MPCParams.PSI;
  // MPCParams.OMEGA;
  // MPCParams.PHI;
  // MPCParams.L1;
  // MPCParams.L2;
  // MPCParams.L3;
  // MPCParams.L4;
  // MPCParams.Fmaxx;
  // MPCParams.Fmaxy;
  // MPCParams.Fmaxz;
  //  MPCParams.ConvThresh;
  // MPCParams.rCollAvoid;
  // MPCParams.maxiter;
  // MPCParams.tol;
  // MPCParams.ObstAvoid;
  // MPCParams.aObst; % Obstacle semi-major axis
  // MPCParams.bObst; % Obstacle semi-minor axis
  // MPCParams.alphaObst; % Orientation of ellipse
  // P_debris = MPCParams.P_debris;
  // MPCParams.pos_of_debris;
  //  MPCParams.AppCone;
  //  MPCParams.phi;
  //  sampT = MPCParams.Ts;
  //  Set Decision Variables - Configure constraints
  //  Compute distance from target
  scale = 3.3121686421112381E-170;
  absxk = std::abs(x0[0]);
  if (absxk > 3.3121686421112381E-170) {
    dr = 1.0;
    scale = absxk;
  } else {
    t = absxk / 3.3121686421112381E-170;
    dr = t * t;
  }

  absxk = std::abs(x0[1]);
  if (absxk > scale) {
    t = scale / absxk;
    dr = dr * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    dr += t * t;
  }

  dr = scale * std::sqrt(dr);

  //  position tolerance to stop control action
  //  Define vectors for Obstacle Avoidance
  //  to current position
  //  to obstacle
  //  to final position
  //  from obstacle to target
  //  from obstacle to FSS
  //  Angle between unit vectors
  //  Approach Cone variables
  //  point where cone is placed
  //  orientation of cone (inertial)
  //  Unit vector of approach axis
  //  Angle between FSS and approach axis
  dock_flag = 0.0;
  CollAvoid_flag = 0.0;

  //  Set obstacle avoidance flag
  //  Define target point for current iteration
  //      xfinal = xfinal;
  //  to final position
  for (i = 0; i < 6; i++) {
    target_state[i] = 0.0;
  }

  //  Define and Solve QP
  //  dock_complete = 0;
  pt_sel = 0.0;
  dock_complete = 0.0;

  //  Compute QP Matrices
  // Dimension of GAMMA is n*horizon x n
  // for i = 1:horizon
  for (b_i = 0; b_i < 120; b_i++) {
    t = 0.0;
    for (i1 = 0; i1 < 6; i1++) {
      t += 2.0 * x0[i1] * dv[i1 + 6 * b_i];
    }

    GAMMA[b_i] = t;
  }

  for (b_i = 0; b_i < 120; b_i++) {
    t = 0.0;
    for (i1 = 0; i1 < 120; i1++) {
      t += GAMMA[i1] * b[i1 + 120 * b_i];
    }

    Bineq[b_i] = t;
  }

  for (b_i = 0; b_i < 60; b_i++) {
    t = 0.0;
    for (i1 = 0; i1 < 120; i1++) {
      t += Bineq[i1] * b_b[i1 + 120 * b_i];
    }

    H[b_i] = t;
  }

  //  Control Consraints
  i = -1;
  for (idx = 0; idx < 20; idx++) {
    for (j2 = 0; j2 < 3; j2++) {
      t = B[3 * j2];
      scale = B[3 * j2 + 1];
      absxk = B[3 * j2 + 2];
      for (b_i1 = 0; b_i1 < 20; b_i1++) {
        i++;
        a_tmp = A[b_i1 + 20 * idx];
        a[i] = static_cast<double>(a_tmp) * t;
        i++;
        a[i] = static_cast<double>(a_tmp) * scale;
        i++;
        a[i] = static_cast<double>(a_tmp) * absxk;
      }
    }
  }

  for (b_i = 0; b_i < 60; b_i++) {
    t = 0.0;
    for (i1 = 0; i1 < 60; i1++) {
      t += a[b_i + 60 * i1];
    }

    u_max[b_i] = t;
  }

  std::memset(&Aineq[0], 0, 7200U * sizeof(signed char));
  std::memset(&a[0], 0, 3600U * sizeof(double));
  for (i = 0; i < 60; i++) {
    a[i + 60 * i] = 1.0;
  }

  for (b_i = 0; b_i < 60; b_i++) {
    for (i1 = 0; i1 < 60; i1++) {
      idx = static_cast<int>(a[i1 + 60 * b_i]);
      i = i1 + 120 * b_i;
      Aineq[i] = static_cast<signed char>(-idx);
      Aineq[i + 60] = static_cast<signed char>(idx);
    }
  }

  std::memset(&Bineq[0], 0, 120U * sizeof(double));

  //  Control Constraints only
  //      Aineq = Aineq;
  //      Bineq = Bineq; % Joh and me comment this on spet 23
  //  Control Constraints only
  //      Aineq = Aineq;
  //      Bineq = Bineq;
  //  Call QP Solver
  for (i = 0; i < 60; i++) {
    t = u_max[i];
    Bineq[i] = t;
    Bineq[i + 60] = t;
  }
   /* if (X_QP[0].empty()){
   for (i = 0; i < 60; i++){
      X_QP[i] = 0.0;

    }
  } */
  //  Solve quadratic programming problem using Wright's (1997) Method
  //  Minimise J(x) = 1/2x'Hx + f'x
  //  Subject to: Ax <= b
  //  Supporting Functions
  //  Reference: S. J. Wright, "Applying New Optimization Algorithms to Model
  //  Predictive Control," in Chemical Process Control-V, CACHE, AIChE
  //  Symposium, 1997, pp. 147-155.
  // Number of decision variables
  //  p = 0;
  // Test for Cold Start
  // Warm Start
  // to tune
  // to tune
  // Default Values
  mu = 10000.0;
  for (i = 0; i < 120; i++) {
    lam[i] = 100.0;
    ftol[i] = 100.0;
    esig[i] = 0.001;
    for (b_i = 0; b_i < 60; b_i++) {
      At[b_i + 60 * i] = Aineq[i + 120 * b_i];
    }
  }

  //  %Linsolve options
  //  opU.UT = true;
  //  opUT.UT = true;
  //  opUT.TRANSA = true;
  // Begin Searching
  //  for iter = 1:maxiter
  iter = 0;
  exitflag = 0;
  while ((iter <= 100) && (exitflag != 1)) {
    boolean_T exitg1;
    boolean_T y;

    // Create common matrices
    for (i = 0; i < 120; i++) {
      t = lam[i];
      scale = 1.0 / t;
      ilam[i] = scale;
      GAMMA[i] = -t / ftol[i];
      mesil[i] = mu * esig[i] * scale;
    }

    // RHS
    for (i = 0; i < 60; i++) {
      for (b_i = 0; b_i < 120; b_i++) {
        idx = b_i + 120 * i;
        IGA[idx] = GAMMA[b_i] * static_cast<double>(Aineq[idx]);
      }

      t = 0.0;
      for (b_i = 0; b_i < 60; b_i++) {
        t += b_a[i + 60 * b_i] * X_QP[b_i];
      }

      scale = 0.0;
      for (b_i = 0; b_i < 120; b_i++) {
        scale += static_cast<double>(At[i + 60 * b_i]) * lam[b_i];
      }

      u_max[i] = (t - scale) - H[i];
    }

    for (b_i = 0; b_i < 7200; b_i++) {
      b_Aineq[b_i] = static_cast<short>(-Aineq[b_i]);
    }

    for (b_i = 0; b_i < 120; b_i++) {
      t = 0.0;
      for (i1 = 0; i1 < 60; i1++) {
        t += static_cast<double>(b_Aineq[b_i + 120 * i1]) * X_QP[i1];
      }

      igr2[b_i] = GAMMA[b_i] * ((t + Bineq[b_i]) - mesil[b_i]);
    }

    // Solve
    for (b_i = 0; b_i < 60; b_i++) {
      for (i1 = 0; i1 < 60; i1++) {
        t = 0.0;
        for (idx = 0; idx < 120; idx++) {
          t += static_cast<double>(At[b_i + 60 * idx]) * IGA[idx + 120 * i1];
        }

        a[b_i + 60 * i1] = t;
      }
    }

    for (b_i = 0; b_i < 3600; b_i++) {
      A_data[b_i] = b_H[b_i] - a[b_i];
    }

    a_tmp = -1;
    idx = 0;
    exitg1 = false;
    while ((!exitg1) && (idx < 60)) {
      int idxA1j;
      int idxAjj;
      int ix;
      idxA1j = idx * 60;
      idxAjj = idxA1j + idx;
      scale = 0.0;
      if (idx >= 1) {
        ix = idxA1j;
        iy = idxA1j;
        for (i = 0; i < idx; i++) {
          scale += A_data[ix] * A_data[iy];
          ix++;
          iy++;
        }
      }

      scale = A_data[idxAjj] - scale;
      if (scale > 0.0) {
        scale = std::sqrt(scale);
        A_data[idxAjj] = scale;
        if (idx + 1 < 60) {
          int idxAjjp1;
          i = idxA1j + 61;
          idxAjjp1 = idxAjj + 61;
          if (idx != 0) {
            iy = idxAjj + 60;
            b_i = (idxA1j + 60 * (58 - idx)) + 61;
            for (j2 = i; j2 <= b_i; j2 += 60) {
              ix = idxA1j;
              absxk = 0.0;
              i1 = (j2 + idx) - 1;
              for (b_i1 = j2; b_i1 <= i1; b_i1++) {
                absxk += A_data[b_i1 - 1] * A_data[ix];
                ix++;
              }

              A_data[iy] += -absxk;
              iy += 60;
            }
          }

          scale = 1.0 / scale;
          b_i = (idxAjj + 60 * (58 - idx)) + 61;
          for (i = idxAjjp1; i <= b_i; i += 60) {
            A_data[i - 1] *= scale;
          }
        }

        idx++;
      } else {
        a_tmp = idx;
        exitg1 = true;
      }
    }

    // [R] = chol(H-At*IGA);
    if (a_tmp + 1 == 0) {
      for (b_i = 0; b_i < 60; b_i++) {
        t = 0.0;
        for (i1 = 0; i1 < 120; i1++) {
          t += static_cast<double>(At[b_i + 60 * i1]) * igr2[i1];
        }

        del_z[b_i] = u_max[b_i] - t;
      }

      for (b_i = 0; b_i < 3600; b_i++) {
        a[b_i] = b_H[b_i] - a[b_i];
      }

      mldivide(a, del_z);

      // old method (LU?)
      //   del_z = linsolve (R, linsolve (R, (r1-At*igr2), opUT), opU); %exploit matrix properties for solving 
    } else {
      // Not Positive Definite (problem? eg infeasible)
      for (b_i = 0; b_i < 60; b_i++) {
        t = 0.0;
        for (i1 = 0; i1 < 120; i1++) {
          t += static_cast<double>(At[b_i + 60 * i1]) * igr2[i1];
        }

        del_z[b_i] = u_max[b_i] - t;
      }

      for (b_i = 0; b_i < 3600; b_i++) {
        a[b_i] = b_H[b_i] - a[b_i];
      }

      mldivide(a, del_z);

      // old method (LU?)
    }

    // Decide on suitable alpha (from Wright's paper)
    // Try Max Increment (alpha = 1)
    // Check lam and ftol > 0
    for (i = 0; i < 120; i++) {
      t = 0.0;
      for (b_i = 0; b_i < 60; b_i++) {
        t += IGA[i + 120 * b_i] * del_z[b_i];
      }

      t = igr2[i] - t;
      del_lam[i] = t;
      scale = ftol[i];
      absxk = (-scale + mesil[i]) - ilam[i] * scale * t;
      mesil[i] = absxk;
      t += lam[i];
      GAMMA[i] = t;
      scale += absxk;
      ilam[i] = scale;
      x[i] = (t < 2.2204460492503131E-16);
      x[i + 120] = (scale < 2.2204460492503131E-16);
    }

    y = false;
    i = 0;
    exitg1 = false;
    while ((!exitg1) && (i < 240)) {
      if (!x[i]) {
        i++;
      } else {
        y = true;
        exitg1 = true;
      }
    }

    if (!y) {
      // KKT met
      std::memcpy(&lam[0], &GAMMA[0], 120U * sizeof(double));
      std::memcpy(&ftol[0], &ilam[0], 120U * sizeof(double));
      for (b_i = 0; b_i < 60; b_i++) {
        X_QP[b_i] += del_z[b_i];
      }
    } else {
      // KKT failed - solve by finding minimum ratio
      for (b_i = 0; b_i < 120; b_i++) {
        result[b_i] = GAMMA[b_i];
        result[b_i + 120] = ilam[b_i];
      }

      idx = 0;
      i = 0;
      exitg1 = false;
      while ((!exitg1) && (i < 240)) {
        if (result[i] < 2.2204460492503131E-16) {
          idx++;
          ii_data[idx - 1] = static_cast<unsigned char>(i + 1);
          if (idx >= 240) {
            exitg1 = true;
          } else {
            i++;
          }
        } else {
          i++;
        }
      }

      if (1 > idx) {
        j2 = 0;
      } else {
        j2 = idx;
      }

      // detects elements breaking KKT condition
      for (b_i = 0; b_i < 120; b_i++) {
        b_del_lam[b_i] = del_lam[b_i];
        b_del_lam[b_i + 120] = mesil[b_i];
      }

      for (b_i = 0; b_i < j2; b_i++) {
        i = ii_data[b_i] - 1;
        varargin_1_data[b_i] = 1.0 - result[i] / b_del_lam[i];
      }

      if (j2 <= 2) {
        if (j2 == 1) {
          scale = varargin_1_data[0];
        } else if ((varargin_1_data[0] > varargin_1_data[1]) || (
                    (varargin_1_data[0]) && (!rtIsNaN(varargin_1_data[1])))) {
          scale = varargin_1_data[1];
        } else {
          scale = varargin_1_data[0];
        }
      } else {
        if (!rtIsNaN(varargin_1_data[0])) {
          idx = 1;
        } else {
          idx = 0;
          i = 2;
          exitg1 = false;
          while ((!exitg1) && (i <= j2)) {
            if (!rtIsNaN(varargin_1_data[i - 1])) {
              idx = i;
              exitg1 = true;
            } else {
              i++;
            }
          }
        }

        if (idx == 0) {
          scale = varargin_1_data[0];
        } else {
          scale = varargin_1_data[idx - 1];
          b_i = idx + 1;
          for (i = b_i; i <= j2; i++) {
            t = varargin_1_data[i - 1];
            if (scale > t) {
              scale = t;
            }
          }
        }
      }

      scale *= 0.995;

      // solves for min ratio (max value of alpha allowed)
      // Increment
      for (b_i = 0; b_i < 120; b_i++) {
        lam[b_i] += scale * del_lam[b_i];
        ftol[b_i] += scale * mesil[b_i];
      }

      for (b_i = 0; b_i < 60; b_i++) {
        X_QP[b_i] += scale * del_z[b_i];
      }
    }

    // Complimentary Gap
    absxk = mu;
    scale = 0.0;
    for (b_i = 0; b_i < 120; b_i++) {
      scale += ftol[b_i] * lam[b_i];
    }

    mu = scale / 120.0;

    //      if(mu < tol)
    //          exitflag = 1;
    //          return
    //      end
    //      %Solve for new Sigma
    //      sigma = mu/mu_old;
    //      if(sigma > 0.1) %to tune
    //          sigma = 0.1;
    //      end
    //      esig = sigma*ones(mc,1);
    if (mu < 0.001) {
      exitflag = 1;
    } else {
      // Solve for new Sigma
      scale = mu / absxk;
      if (scale > 0.1) {
        // to tune
        scale = 0.1;
      }

      for (i = 0; i < 120; i++) {
        esig[i] = scale;
      }
    }

    iter++;
  }

  // Check for failure
  num_iter = iter;

  //  solution to warm-start next iteration
  //  Extract first control
  Fx = X_QP[0];
  Fy = X_QP[1];
  Fz = X_QP[2];
}

template<typename T>
void CoordinatorBase<T>::MPC_Guidance_v3_sand_worst()
{
   double b[14400] = {  100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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
    100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0,
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
    0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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
    0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0,
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
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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
    100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0,
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
    0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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
    0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0,
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
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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
    100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0,
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
    0.0, 0.0, 0.0, 0.0, 200159.068, -1.5235E-9, -5.4232E-8, 78490.927,
    -1.6858E-8, 4.5618E-9, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.5235E-9,
    200159.068, 2.1905E-9, 3.9901E-10, 78490.927, -2.0594E-9, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, -5.4232E-8, 2.1905E-9, 200159.068, -2.562E-8,
    3.0017E-9, 78490.927, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 78490.927,
    3.9901E-10, -2.562E-8, 2.5630794002E+6, 3.1661E-11, 4.0701E-11, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.6858E-8, 78490.927, 3.0017E-9,
    3.1661E-11, 2.5630794002E+6, 2.5296E-7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 4.5618E-9, -2.0594E-9, 78490.927, 4.0701E-11, 2.5296E-7,
    2.5630794002E+6} ;
   double b_b[7200] = {1.0194E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    3.0583E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 5.0972E-5, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 7.136E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 9.1749E-5, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.00011214, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00013253,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00015291, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.0001733, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00019369, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.00021408, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00023447, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.00025486, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00027525,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00029564, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.00031602, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00033641, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.0003568, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00037719, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.00039758, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0,
    1.0194E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 3.0583E-5, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 5.0972E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 7.136E-5, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 9.1749E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00011214,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00013253, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.00015291, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0001733, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.00019369, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00021408, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.00023447, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00025486,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00027525, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.00029564, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00031602, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.00033641, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0003568, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.00037719, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00039758,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0, 1.0194E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    3.0583E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 5.0972E-5, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 7.136E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 9.1749E-5, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.00011214, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00013253,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00015291, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.0001733, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00019369, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.00021408, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00023447, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.00025486, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00027525,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00029564, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.00031602, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00033641, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.0003568, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00037719, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.00039758, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0194E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 3.0583E-5, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 5.0972E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 7.136E-5, 0.0,
    0.0, 0.0012743, 0.0, 0.0, 9.1749E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.00011214, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00013253, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.00015291, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0001733, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.00019369, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00021408,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00023447, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.00025486, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00027525, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.00029564, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00031602, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.00033641, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0003568,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00037719, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0194E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    3.0583E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 5.0972E-5, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 7.136E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 9.1749E-5, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.00011214, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00013253,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00015291, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.0001733, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00019369, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.00021408, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00023447, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.00025486, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00027525,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00029564, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.00031602, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00033641, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.0003568, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00037719, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0194E-5, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 3.0583E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 5.0972E-5,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 7.136E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    9.1749E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00011214, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.00013253, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00015291, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.0001733, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00019369,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00021408, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.00023447, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00025486, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.00027525, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00029564, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.00031602, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00033641,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0003568, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.00037719, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0194E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 3.0583E-5, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 5.0972E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 7.136E-5, 0.0,
    0.0, 0.0012743, 0.0, 0.0, 9.1749E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.00011214, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00013253, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.00015291, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0001733, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.00019369, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00021408,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00023447, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.00025486, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00027525, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.00029564, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00031602, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.00033641, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0003568,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0194E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 3.0583E-5, 0.0,
    0.0, 0.0012743, 0.0, 0.0, 5.0972E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 7.136E-5,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 9.1749E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.00011214, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00013253, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.00015291, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0001733, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.00019369, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00021408,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00023447, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.00025486, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00027525, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.00029564, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00031602, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.00033641, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0003568,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0194E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 3.0583E-5, 0.0,
    0.0, 0.0012743, 0.0, 0.0, 5.0972E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 7.136E-5,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 9.1749E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.00011214, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00013253, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.00015291, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0001733, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.00019369, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00021408,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00023447, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.00025486, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00027525, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.00029564, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00031602, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.00033641, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0003568,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0194E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    3.0583E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 5.0972E-5, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 7.136E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 9.1749E-5, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.00011214, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00013253,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00015291, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.0001733, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00019369, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.00021408, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00023447, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.00025486, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00027525,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00029564, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.00031602, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00033641, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0194E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    3.0583E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 5.0972E-5, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 7.136E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 9.1749E-5, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.00011214, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00013253,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00015291, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.0001733, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00019369, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.00021408, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00023447, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.00025486, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00027525,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00029564, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.00031602, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00033641, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0194E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    3.0583E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 5.0972E-5, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 7.136E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 9.1749E-5, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.00011214, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00013253,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00015291, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.0001733, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00019369, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.00021408, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00023447, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.00025486, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00027525,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00029564, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.00031602, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00033641, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0194E-5, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 3.0583E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 5.0972E-5, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 7.136E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 9.1749E-5, 0.0,
    0.0, 0.0012743, 0.0, 0.0, 0.00011214, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.00013253, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00015291, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.0001733, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00019369, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.00021408, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00023447,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00025486, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.00027525, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00029564, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.00031602, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0194E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 3.0583E-5,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 5.0972E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    7.136E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 9.1749E-5, 0.0, 0.0, 0.0012743, 0.0,
    0.0, 0.00011214, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00013253, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.00015291, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0001733,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00019369, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.00021408, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00023447, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.00025486, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00027525, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.00029564, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00031602,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0194E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 3.0583E-5, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 5.0972E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 7.136E-5, 0.0,
    0.0, 0.0012743, 0.0, 0.0, 9.1749E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.00011214, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00013253, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.00015291, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0001733, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.00019369, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00021408,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00023447, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.00025486, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00027525, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.00029564, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00031602, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0194E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 3.0583E-5, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 5.0972E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 7.136E-5, 0.0,
    0.0, 0.0012743, 0.0, 0.0, 9.1749E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.00011214, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00013253, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.00015291, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0001733, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.00019369, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00021408,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00023447, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.00025486, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00027525, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.00029564, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0194E-5, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 3.0583E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 5.0972E-5,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 7.136E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    9.1749E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00011214, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.00013253, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00015291, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.0001733, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00019369,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00021408, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.00023447, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00025486, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.00027525, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00029564, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0194E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 3.0583E-5,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 5.0972E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    7.136E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 9.1749E-5, 0.0, 0.0, 0.0012743, 0.0,
    0.0, 0.00011214, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00013253, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.00015291, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0001733,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00019369, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.00021408, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00023447, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.00025486, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00027525, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.00029564, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0194E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 3.0583E-5, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 5.0972E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 7.136E-5, 0.0,
    0.0, 0.0012743, 0.0, 0.0, 9.1749E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.00011214, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00013253, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.00015291, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0001733, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.00019369, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00021408,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00023447, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.00025486, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00027525, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0194E-5, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 3.0583E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 5.0972E-5, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 7.136E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 9.1749E-5, 0.0,
    0.0, 0.0012743, 0.0, 0.0, 0.00011214, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.00013253, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00015291, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.0001733, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00019369, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.00021408, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00023447,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00025486, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.00027525, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0194E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 3.0583E-5, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 5.0972E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 7.136E-5, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 9.1749E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00011214,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00013253, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.00015291, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0001733, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.00019369, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00021408, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.00023447, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00025486,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00027525, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0194E-5, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 3.0583E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 5.0972E-5,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 7.136E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    9.1749E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00011214, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.00013253, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00015291, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.0001733, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00019369,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00021408, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.00023447, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00025486, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0194E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 3.0583E-5, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 5.0972E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 7.136E-5, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 9.1749E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00011214,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00013253, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.00015291, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0001733, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.00019369, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00021408, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.00023447, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00025486,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0194E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 3.0583E-5, 0.0,
    0.0, 0.0012743, 0.0, 0.0, 5.0972E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 7.136E-5,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 9.1749E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.00011214, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00013253, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.00015291, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0001733, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.00019369, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00021408,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00023447, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.00025486, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0194E-5, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 3.0583E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 5.0972E-5, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 7.136E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 9.1749E-5, 0.0,
    0.0, 0.0012743, 0.0, 0.0, 0.00011214, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.00013253, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00015291, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.0001733, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00019369, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.00021408, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00023447,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0194E-5, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 3.0583E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 5.0972E-5,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 7.136E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    9.1749E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00011214, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.00013253, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00015291, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.0001733, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00019369,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00021408, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.00023447, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0194E-5, 0.0,
    0.0, 0.0012743, 0.0, 0.0, 3.0583E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    5.0972E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 7.136E-5, 0.0, 0.0, 0.0012743, 0.0,
    0.0, 9.1749E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00011214, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.00013253, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00015291,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0001733, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.00019369, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00021408, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.00023447, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0194E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 3.0583E-5, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 5.0972E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 7.136E-5, 0.0,
    0.0, 0.0012743, 0.0, 0.0, 9.1749E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.00011214, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00013253, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.00015291, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0001733, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.00019369, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00021408,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0194E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 3.0583E-5, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 5.0972E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 7.136E-5, 0.0,
    0.0, 0.0012743, 0.0, 0.0, 9.1749E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.00011214, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00013253, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.00015291, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0001733, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.00019369, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00021408,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0194E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 3.0583E-5, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 5.0972E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 7.136E-5, 0.0,
    0.0, 0.0012743, 0.0, 0.0, 9.1749E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.00011214, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00013253, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.00015291, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0001733, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.00019369, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00021408,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0194E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 3.0583E-5, 0.0,
    0.0, 0.0012743, 0.0, 0.0, 5.0972E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 7.136E-5,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 9.1749E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.00011214, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00013253, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.00015291, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0001733, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.00019369, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0194E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 3.0583E-5, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 5.0972E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 7.136E-5, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 9.1749E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00011214,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00013253, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.00015291, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0001733, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.00019369, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0194E-5, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 3.0583E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 5.0972E-5,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 7.136E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    9.1749E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00011214, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.00013253, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00015291, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.0001733, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00019369,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0194E-5, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 3.0583E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 5.0972E-5,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 7.136E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    9.1749E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00011214, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.00013253, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00015291, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.0001733, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0194E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 3.0583E-5, 0.0,
    0.0, 0.0012743, 0.0, 0.0, 5.0972E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 7.136E-5,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 9.1749E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.00011214, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00013253, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.00015291, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0001733, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0194E-5, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 3.0583E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 5.0972E-5,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 7.136E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    9.1749E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00011214, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.00013253, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00015291, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.0001733, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0194E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    3.0583E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 5.0972E-5, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 7.136E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 9.1749E-5, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.00011214, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00013253,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00015291, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0194E-5,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 3.0583E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    5.0972E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 7.136E-5, 0.0, 0.0, 0.0012743, 0.0,
    0.0, 9.1749E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00011214, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.00013253, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00015291,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0194E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 3.0583E-5, 0.0,
    0.0, 0.0012743, 0.0, 0.0, 5.0972E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 7.136E-5,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 9.1749E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.00011214, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00013253, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.00015291, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0194E-5, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 3.0583E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 5.0972E-5,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 7.136E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    9.1749E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00011214, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.00013253, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0194E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 3.0583E-5, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 5.0972E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 7.136E-5, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 9.1749E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00011214,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00013253, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0194E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 3.0583E-5, 0.0,
    0.0, 0.0012743, 0.0, 0.0, 5.0972E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 7.136E-5,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 9.1749E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    0.00011214, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.00013253, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0194E-5, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 3.0583E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 5.0972E-5, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 7.136E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 9.1749E-5, 0.0,
    0.0, 0.0012743, 0.0, 0.0, 0.00011214, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0194E-5, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 3.0583E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 5.0972E-5, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 7.136E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 9.1749E-5, 0.0,
    0.0, 0.0012743, 0.0, 0.0, 0.00011214, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0194E-5, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 3.0583E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 5.0972E-5, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 7.136E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 9.1749E-5, 0.0,
    0.0, 0.0012743, 0.0, 0.0, 0.00011214, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0194E-5, 0.0,
    0.0, 0.0012743, 0.0, 0.0, 3.0583E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    5.0972E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 7.136E-5, 0.0, 0.0, 0.0012743, 0.0,
    0.0, 9.1749E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0194E-5, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 3.0583E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 5.0972E-5,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 7.136E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    9.1749E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0194E-5, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 3.0583E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 5.0972E-5, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 7.136E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 9.1749E-5, 0.0,
    0.0, 0.0012743, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0194E-5, 0.0, 0.0, 0.0012743,
    0.0, 0.0, 3.0583E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 5.0972E-5, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 7.136E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0194E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 3.0583E-5, 0.0,
    0.0, 0.0012743, 0.0, 0.0, 5.0972E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 7.136E-5,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0194E-5,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 3.0583E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    5.0972E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 7.136E-5, 0.0, 0.0, 0.0012743, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0194E-5, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 3.0583E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 5.0972E-5,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0194E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 3.0583E-5, 0.0,
    0.0, 0.0012743, 0.0, 0.0, 5.0972E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0194E-5, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 3.0583E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 5.0972E-5,
    0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0194E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0,
    3.0583E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0194E-5, 0.0,
    0.0, 0.0012743, 0.0, 0.0, 3.0583E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0194E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 3.0583E-5, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0194E-5, 0.0, 0.0,
    0.0012743, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0194E-5, 0.0, 0.0, 0.0012743, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0194E-5, 0.0, 0.0, 0.0012743 };

  
   double b_H[3600] = { 214.71719164095626, -1.7056254126311141E-14,
    -3.8350397641984218E-14, 14.385084459677893, -1.7052288080175281E-14,
    -3.6139749542427422E-14, 14.052977361645764, -1.7048322034039424E-14,
    -3.3929101442870619E-14, 13.720870430006237, -1.704435598790356E-14,
    -3.1718453343313816E-14, 13.388763747903601, -1.70403899417677E-14,
    -2.9507805243757019E-14, 13.056660989691615, -1.7036425840725581E-14,
    -2.7298241326691417E-14, 12.724555055265856, -1.703245979458972E-14,
    -2.508759322713462E-14, 12.392449619783397, -1.702849374845386E-14,
    -2.287694512757782E-14, 12.060344766374264, -1.7024527702318E-14,
    -2.066629702802102E-14, 11.728240578168471, -1.702056165618214E-14,
    -1.8455648928464221E-14, 11.39613713831646, -1.701659561004628E-14,
    -1.6245000828907424E-14, 11.06403452996865, -1.701262956391042E-14,
    -1.4034352729350621E-14, 10.731932836275458, -1.7008663517774563E-14,
    -1.1823704629793823E-14, 10.399835731596751, -1.700469941673244E-14,
    -9.614140712728219E-15, 10.067736115828115, -1.7000733370596583E-14,
    -7.4034926131714191E-15, 9.7356373050036282, -1.6996767129951346E-14,
    -5.1927360953655005E-15, 9.4035404596997552, -1.6992801278324863E-14,
    -2.9821964140578193E-15, 9.0714453037624718, -1.6988835621207749E-14,
    -7.7176515099926037E-16, 8.7393512020126352, -1.6984869769581263E-14,
    1.4387745303084196E-15, 8.4072585967970479, -1.6980903917954779E-14,
    3.6493142116160993E-15, -1.7056254126311141E-14, 214.71719164095626,
    8.2318069131290087E-13, -1.6155511092337738E-14, 14.385084459677893,
    8.2298918984796686E-13, -1.5254768058364339E-14, 14.052977361645764,
    8.2279768838303274E-13, -1.435402502439094E-14, 13.720870430006237,
    8.2260618691809882E-13, -1.3453281990417541E-14, 13.388763747903601,
    8.224146854531648E-13, -1.2552980713694741E-14, 13.056660989691615,
    8.222232779075368E-13, -1.1652237679721339E-14, 12.724555055265856,
    8.2203177644260278E-13, -1.0751494645747939E-14, 12.392449619783397,
    8.2184027497766887E-13, -9.85075161177454E-15, 12.060344766374264,
    8.2164877351273475E-13, -8.95000857780114E-15, 11.728240578168471,
    8.2145727204780073E-13, -8.04926554382774E-15, 11.39613713831646,
    8.2126577058286682E-13, -7.14852250985434E-15, 11.06403452996865,
    8.210742691179328E-13, -6.2477794758809407E-15, 10.731932836275458,
    8.2088276765299878E-13, -5.3474781991581405E-15, 10.399835731596751,
    8.2069136010737078E-13, -4.4467351651847404E-15, 10.067736115828115,
    8.2049985864243676E-13, -3.5459479554862802E-15, 9.7356373050036282,
    8.2030834778557209E-13, -2.6452490972379395E-15, 9.4035404596997552,
    8.2011685571256873E-13, -1.74459441471466E-15, 9.0714453037624718,
    8.19925373031496E-13, -8.4389555646631984E-16, 8.7393512020126352,
    8.1973388095849255E-13, 5.6803301782020074E-17, 8.4072585967970479,
    8.1954238888548918E-13, -3.8350397641984218E-14, 8.2318069131290067E-13,
    214.71719164095626, -3.7708175572864622E-14, 8.23252194737756E-13,
    14.385084459677893, -3.7065953503745013E-14, 8.2332369816261152E-13,
    14.052977361645764, -3.6423731434625422E-14, 8.2339520158746689E-13,
    13.720870430006237, -3.578150936550582E-14, 8.2346670501232227E-13,
    13.388763747903601, -3.5139602265522621E-14, 8.2353817336928915E-13,
    13.056660989691615, -3.4497380196403018E-14, 8.2360967679414463E-13,
    12.724555055265856, -3.3855158127283415E-14, 8.23681180219E-13,
    12.392449619783397, -3.3212936058163819E-14, 8.2375268364385537E-13,
    12.060344766374264, -3.2570713989044216E-14, 8.2382418706871085E-13,
    11.728240578168471, -3.1928491919924619E-14, 8.2389569049356622E-13,
    11.39613713831646, -3.1286269850805016E-14, 8.2396719391842159E-13,
    11.06403452996865, -3.064404778168542E-14, 8.24038697343277E-13,
    10.731932836275458, -3.0002140681702215E-14, 8.2411016570024385E-13,
    10.399835731596751, -2.9359918612582618E-14, 8.2418166912509923E-13,
    10.067736115828115, -2.8717665046549382E-14, 8.2425317605674349E-13,
    9.7356373050036282, -2.8075474474343419E-14, 8.2432467597481008E-13,
    9.4035404596997552, -2.7433315399051097E-14, 8.2439617238608767E-13,
    9.0714453037624718, -2.679112482684514E-14, 8.2446767230415426E-13,
    8.7393512020126352, -2.6148934254639179E-14, 8.2453917222222064E-13,
    8.4072585967970479, 14.385084459677891, -1.6155511092337741E-14,
    -3.7708175572864622E-14, 214.37791332680061, -1.6152811842850578E-14,
    -3.5542621609082219E-14, 14.045974116808761, -1.6150112593363422E-14,
    -3.3377067645299817E-14, 13.714034990063162, -1.6147413343876259E-14,
    -3.121151368151742E-14, 13.382096029710162, -1.61447140943891E-14,
    -2.9045959717735018E-14, 13.050160828519424, -1.614201616871238E-14,
    -2.6881467820634218E-14, 12.718222449589355, -1.613931691922522E-14,
    -2.4715913856851822E-14, 12.386284486450126, -1.613661766973806E-14,
    -2.2550359893069419E-14, 12.054347022254197, -1.61339184202509E-14,
    -2.038480592928702E-14, 11.722410140131593, -1.613121917076374E-14,
    -1.821925196550462E-14, 11.390473923212333, -1.612851992127658E-14,
    -1.6053698001722224E-14, 11.05853845464685, -1.612582067178942E-14,
    -1.3888144037939821E-14, 10.72660381758557, -1.6123121422302263E-14,
    -1.1722590074157422E-14, 10.394673604804272, -1.6120423496625541E-14,
    -9.5580981770566187E-15, 10.062740879407485, -1.6117724247138381E-14,
    -7.39254421327422E-15, 9.73080888396284, -1.6115024865270176E-14,
    -5.2268840428236609E-15, 9.3988787545878836, -1.6112325748164062E-14,
    -3.0614362857094196E-15, 9.0669502232951515, -1.610962676343899E-14,
    -8.9609473526334018E-16, 8.73502267120601, -1.6106927646332872E-14,
    1.2693530218508995E-15, 8.40309653250478, -1.6104228529226758E-14,
    3.43480077896514E-15, -1.7052288080175281E-14, 14.385084459677891,
    8.2325219473775625E-13, -1.6152811842850578E-14, 214.37791332680061,
    8.2306251468277238E-13, -1.5253335605525879E-14, 14.045974116808761,
    8.2287283462778832E-13, -1.435385936820118E-14, 13.714034990063162,
    8.2268315457280446E-13, -1.345438313087648E-14, 13.382096029710162,
    8.224934745178206E-13, -1.2555348029519081E-14, 13.050160828519424,
    8.2230388748885675E-13, -1.1655871792194379E-14, 12.718222449589355,
    8.2211420743387289E-13, -1.075639555486968E-14, 12.386284486450126,
    8.21924527378889E-13, -9.85691931754498E-15, 12.054347022254197,
    8.2173484732390507E-13, -8.95744308022028E-15, 11.722410140131593,
    8.215451672689212E-13, -8.05796684289558E-15, 11.390473923212333,
    8.2135548721393724E-13, -7.1584906055708821E-15, 11.05853845464685,
    8.2116580715895338E-13, -6.2590143682461807E-15, 10.72660381758557,
    8.2097612710396952E-13, -5.35997926688878E-15, 10.394673604804272,
    8.2078654007500567E-13, -4.46050302956408E-15, 10.062740879407485,
    8.2059686002002171E-13, -3.5609826786426505E-15, 9.73080888396284,
    8.2040717066243584E-13, -2.6615505549146795E-15, 9.3988787545878836,
    8.20217499910054E-13, -1.7621625447834399E-15, 9.0669502232951515,
    8.2002783846027415E-13, -8.627304210554699E-16, 8.73502267120601,
    8.198381677078922E-13, 3.670170267250009E-17, 8.40309653250478,
    8.1964849695551035E-13, -3.6139749542427415E-14, 8.2298918984796665E-13,
    14.385084459677891, -3.5542621609082219E-14, 8.2306251468277218E-13,
    214.37791332680061, -3.4945493675737017E-14, 8.2313583951757771E-13,
    14.045974116808761, -3.4348365742391821E-14, 8.2320916435238314E-13,
    13.714034990063162, -3.3751237809046618E-14, 8.2328248918718867E-13,
    13.382096029710162, -3.3154402729028215E-14, 8.2335577806081982E-13,
    13.050160828519424, -3.2557274795683019E-14, 8.2342910289562534E-13,
    12.718222449589355, -3.1960146862337817E-14, 8.2350242773043077E-13,
    12.386284486450126, -3.1363018928992614E-14, 8.235757525652363E-13,
    12.054347022254197, -3.0765890995647418E-14, 8.2364907740004183E-13,
    11.722410140131593, -3.0168763062302216E-14, 8.2372240223484736E-13,
    11.390473923212333, -2.957163512895702E-14, 8.2379572706965278E-13,
    11.05853845464685, -2.8974507195611817E-14, 8.2386905190445821E-13,
    10.72660381758557, -2.8377672115593414E-14, 8.2394234077808936E-13,
    10.394673604804272, -2.7780544182248218E-14, 8.2401566561289479E-13,
    10.062740879407485, -2.7183386963570337E-14, 8.2408899404381774E-13,
    9.73080888396284, -2.658628831555782E-14, 8.2416231528250585E-13,
    9.3988787545878836, -2.5989218952877975E-14, 8.2423563292507642E-13,
    9.0669502232951515, -2.539212030486546E-14, 8.2430895416376452E-13,
    8.73502267120601, -2.4795021656852936E-14, 8.2438227540245242E-13,
    8.40309653250478, 14.052977361645764, -1.5254768058364339E-14,
    -3.7065953503745019E-14, 14.045974116808761, -1.5253335605525879E-14,
    -3.4945493675737017E-14, 214.03897087206758, -1.525190315268742E-14,
    -3.2825033847729015E-14, 13.70719955021185, -1.525047069984896E-14,
    -3.0704574019721019E-14, 13.375428311602359, -1.5249038247010497E-14,
    -2.8584114191713023E-14, 13.043660667426746, -1.5247606496699179E-14,
    -2.6464694314577019E-14, 12.711889843992372, -1.5246174043860719E-14,
    -2.434423448656902E-14, 12.380119353198417, -1.5244741591022259E-14,
    -2.2223774658561021E-14, 12.048349278195296, -1.5243309138183797E-14,
    -2.0103314830553019E-14, 11.716579702135476, -1.524187668534534E-14,
    -1.7982855002545023E-14, 11.384810708148985, -1.5240444232506877E-14,
    -1.5862395174537024E-14, 11.053042379365833, -1.5239011779668421E-14,
    -1.3741935346529022E-14, 10.721274798936461, -1.5237579326829961E-14,
    -1.1621475518521023E-14, 10.389511478052565, -1.5236147576518639E-14,
    -9.50205564138502E-15, 10.057745643027635, -1.523471512368018E-14,
    -7.38159581337702E-15, 9.72598046296283, -1.5233282600589005E-14,
    -5.2610319902818204E-15, 9.3942170494963975, -1.523185021800326E-14,
    -3.1406761573610196E-15, 9.0624551428278313, -1.523041790567023E-14,
    -1.0204243195274204E-15, 8.7306941403993878, -1.5228985523084482E-14,
    1.0999315133933793E-15, 8.3989344682125129, -1.5227553140498737E-14,
    3.2202873463141794E-15, -1.704832203403942E-14, 14.052977361645764,
    8.2332369816261162E-13, -1.6150112593363418E-14, 14.045974116808761,
    8.2313583951757781E-13, -1.525190315268742E-14, 214.03897087206758,
    8.22947980872544E-13, -1.4353693712011418E-14, 13.70719955021185,
    8.227601222275102E-13, -1.3455484271335421E-14, 13.375428311602359,
    8.2257226358247639E-13, -1.255771534534342E-14, 13.043660667426746,
    8.2238449707017681E-13, -1.1659505904667419E-14, 12.711889843992372,
    8.22196638425143E-13, -1.076129646399142E-14, 12.380119353198417,
    8.2200877978010919E-13, -9.86308702331542E-15, 12.048349278195296,
    8.2182092113507538E-13, -8.96487758263942E-15, 11.716579702135476,
    8.2163306249004158E-13, -8.06666814196342E-15, 11.384810708148985,
    8.2144520384500777E-13, -7.1684587012874218E-15, 11.053042379365833,
    8.21257345199974E-13, -6.2702492606114207E-15, 10.721274798936461,
    8.2106948655494025E-13, -5.37248033461942E-15, 10.389511478052565,
    8.2088172004264057E-13, -4.47427089394342E-15, 10.057745643027635,
    8.2069386139760676E-13, -3.57601740179902E-15, 9.72598046296283,
    8.2050599353929948E-13, -2.6778520125914195E-15, 9.3942170494963975,
    8.2031814410753915E-13, -1.77973067485222E-15, 9.0624551428278313,
    8.2013030388905229E-13, -8.8156528564462E-16, 8.7306941403993878,
    8.1994245445729185E-13, 1.6600103562980056E-17, 8.3989344682125129,
    8.1975460502553141E-13, -3.3929101442870612E-14, 8.2279768838303264E-13,
    14.052977361645764, -3.3377067645299817E-14, 8.2287283462778832E-13,
    14.045974116808761, -3.2825033847729015E-14, 8.229479808725439E-13,
    214.03897087206758, -3.2273000050158219E-14, 8.2302312711729949E-13,
    13.70719955021185, -3.1720966252587417E-14, 8.2309827336205507E-13,
    13.375428311602359, -3.1169203192533816E-14, 8.2317338275235038E-13,
    13.043660667426746, -3.0617169394963014E-14, 8.23248528997106E-13,
    12.711889843992372, -3.0065135597392218E-14, 8.2332367524186164E-13,
    12.380119353198417, -2.9513101799821416E-14, 8.2339882148661723E-13,
    12.048349278195296, -2.8961068002250614E-14, 8.2347396773137281E-13,
    11.716579702135476, -2.8409034204679819E-14, 8.2354911397612839E-13,
    11.384810708148985, -2.7857000407109017E-14, 8.23624260220884E-13,
    11.053042379365833, -2.7304966609538215E-14, 8.2369940646563956E-13,
    10.721274798936461, -2.6753203549484614E-14, 8.2377451585593487E-13,
    10.389511478052565, -2.6201169751913818E-14, 8.2384966210069045E-13,
    10.057745643027635, -2.5649108880591295E-14, 8.23924812030892E-13,
    9.72598046296283, -2.5097102156772217E-14, 8.2399995459020162E-13,
    9.3942170494963975, -2.4545122506704855E-14, 8.2407509346406516E-13,
    9.0624551428278313, -2.3993115782885777E-14, 8.2415023602337478E-13,
    8.7306941403993878, -2.3441109059066697E-14, 8.2422537858268419E-13,
    8.3989344682125129, 13.720870430006235, -1.435402502439094E-14,
    -3.6423731434625416E-14, 13.714034990063162, -1.4353859368201177E-14,
    -3.4348365742391821E-14, 13.707199550211849, -1.4353693712011421E-14,
    -3.2273000050158219E-14, 213.70036411045635, -1.4353528055821658E-14,
    -3.0197634357924617E-14, 13.368760593586316, -1.43533623996319E-14,
    -2.8122268665691019E-14, 13.037160506419706, -1.4353196824685978E-14,
    -2.6047920808519817E-14, 12.705557238474906, -1.4353031168496222E-14,
    -2.3972555116286218E-14, 12.373954220026224, -1.4352865512306459E-14,
    -2.189718942405262E-14, 12.042351534217955, -1.4352699856116696E-14,
    -1.9821823731819021E-14, 11.710749264200526, -1.435253419992694E-14,
    -1.7746458039585419E-14, 11.379147493126398, -1.4352368543737178E-14,
    -1.5671092347351824E-14, 11.047546304125595, -1.4352202887547418E-14,
    -1.3595726655118221E-14, 10.715945780328134, -1.4352037231357662E-14,
    -1.1520360962884622E-14, 10.384349351341639, -1.4351871656411741E-14,
    -9.44601310571342E-15, 10.05275040668856, -1.4351706000221981E-14,
    -7.37064741347982E-15, 9.7211520420035988, -1.4351540335907835E-14,
    -5.2951799377399808E-15, 9.3895553444456912, -1.4351374687842459E-14,
    -3.2199160290126195E-15, 9.0579600623808982, -1.4351209047901468E-14,
    -1.1447539037915006E-15, 8.7263656095927651, -1.4351043399836092E-14,
    9.3051000493585911E-16, 8.3947724039202445, -1.4350877751770716E-14,
    3.0057739136632192E-15, -1.704435598790356E-14, 13.720870430006235,
    8.23395201587467E-13, -1.6147413343876259E-14, 13.714034990063162,
    8.2320916435238334E-13, -1.525047069984896E-14, 13.707199550211849,
    8.2302312711729959E-13, -1.4353528055821658E-14, 213.70036411045635,
    8.2283708988221583E-13, -1.345658541179436E-14, 13.368760593586316,
    8.2265105264713218E-13, -1.256008266116776E-14, 13.037160506419706,
    8.2246510665149676E-13, -1.1663140017140459E-14, 12.705557238474906,
    8.2227906941641311E-13, -1.076619737311316E-14, 12.373954220026224,
    8.2209303218132945E-13, -9.86925472908586E-15, 12.042351534217955,
    8.219069949462456E-13, -8.97231208505856E-15, 11.710749264200526,
    8.2172095771116195E-13, -8.07536944103126E-15, 11.379147493126398,
    8.215349204760783E-13, -7.1784267970039615E-15, 11.047546304125595,
    8.2134888324099454E-13, -6.2814841529766606E-15, 10.715945780328134,
    8.2116284600591089E-13, -5.38498140235006E-15, 10.384349351341639,
    8.2097690001027547E-13, -4.48803875832276E-15, 10.05275040668856,
    8.2079086277519171E-13, -3.59105212495539E-15, 9.7211520420035988,
    8.2060481641616323E-13, -2.6941534702681594E-15, 9.3895553444456912,
    8.2041878830502431E-13, -1.797298804921E-15, 9.0579600623808982,
    8.2023276931783032E-13, -9.0040015023377E-16, 8.7263656095927651,
    8.200467412066914E-13, -3.5014955465399778E-18, 8.3947724039202445,
    8.1986071309555258E-13, -3.1718453343313809E-14, 8.2260618691809872E-13,
    13.720870430006235, -3.1211513681517414E-14, 8.2268315457280436E-13,
    13.714034990063162, -3.0704574019721013E-14, 8.227601222275101E-13,
    13.707199550211849, -3.0197634357924617E-14, 8.2283708988221583E-13,
    213.70036411045635, -2.9690694696128216E-14, 8.2291405753692147E-13,
    13.368760593586316, -2.9184003656039417E-14, 8.22990987443881E-13,
    13.037160506419706, -2.8677063994243015E-14, 8.2306795509858668E-13,
    12.705557238474906, -2.8170124332446614E-14, 8.2314492275329242E-13,
    12.373954220026224, -2.7663184670650212E-14, 8.2322189040799815E-13,
    12.042351534217955, -2.7156245008853817E-14, 8.2329885806270379E-13,
    11.710749264200526, -2.6649305347057415E-14, 8.2337582571740953E-13,
    11.379147493126398, -2.6142365685261014E-14, 8.2345279337211527E-13,
    11.047546304125595, -2.5635426023464616E-14, 8.235297610268208E-13,
    10.715945780328134, -2.5128734983375813E-14, 8.2360669093378037E-13,
    10.384349351341639, -2.4621795321579418E-14, 8.23683658588486E-13,
    10.05275040668856, -2.4114830797612256E-14, 8.2376063001796635E-13,
    9.7211520420035988, -2.3607915997986618E-14, 8.2383759389789749E-13,
    9.3895553444456912, -2.3101026060531735E-14, 8.2391455400305391E-13,
    9.0579600623808982, -2.2594111260906098E-14, 8.2399151788298494E-13,
    8.7263656095927651, -2.2087196461280457E-14, 8.2406848176291587E-13,
    8.3947724039202445, 13.388763747903601, -1.3453281990417542E-14,
    -3.578150936550582E-14, 13.382096029710164, -1.3454383130876481E-14,
    -3.3751237809046618E-14, 13.375428311602359, -1.3455484271335422E-14,
    -3.1720966252587417E-14, 13.368760593586318, -1.345658541179436E-14,
    -2.9690694696128216E-14, 213.36209287566609, -1.3457686552253301E-14,
    -2.7660423139669021E-14, 13.030660345504426, -1.3458787152672781E-14,
    -2.5631147302462618E-14, 12.699224633043077, -1.3459888293131721E-14,
    -2.3600875746003419E-14, 12.367789086933549, -1.346098943359066E-14,
    -2.1570604189544221E-14, 12.036353790320133, -1.3462090574049601E-14,
    -1.954033263308502E-14, 11.704918826347136, -1.346319171450854E-14,
    -1.7510061076625822E-14, 11.373484278164979, -1.346429285496748E-14,
    -1.5479789520166624E-14, 11.042050228926119, -1.3465393995426422E-14,
    -1.3449517963707421E-14, 10.710616761760587, -1.3466495135885361E-14,
    -1.1419246407248223E-14, 10.379187224671494, -1.3467595736304839E-14,
    -9.38997057004182E-15, 10.047755170390262, -1.3468696876763781E-14,
    -7.35969901358262E-15, 9.7163236210851434, -1.3469798071226666E-14,
    -5.32932788519814E-15, 9.3848936394357629, -1.3470899157681661E-14,
    -3.2991559006642194E-15, 9.053464981974745, -1.347200019013271E-14,
    -1.26908348805558E-15, 8.7220370788065313, -1.3473101276587702E-14,
    7.6108849647833972E-16, 8.3906103396279779, -1.3474202363042697E-14,
    2.79126048101226E-15, -1.70403899417677E-14, 13.388763747903601,
    8.2346670501232247E-13, -1.61447140943891E-14, 13.382096029710164,
    8.2328248918718887E-13, -1.52490382470105E-14, 13.375428311602359,
    8.2309827336205517E-13, -1.43533623996319E-14, 13.368760593586318,
    8.2291405753692157E-13, -1.3457686552253301E-14, 213.36209287566609,
    8.22729841711788E-13, -1.25624499769921E-14, 13.030660345504426,
    8.2254571623281681E-13, -1.1666774129613498E-14, 12.699224633043077,
    8.2236150040768321E-13, -1.07710982822349E-14, 12.367789086933549,
    8.2217728458254962E-13, -9.8754224348563E-15, 12.036353790320133,
    8.2199306875741592E-13, -8.9797465874776991E-15, 11.704918826347136,
    8.2180885293228232E-13, -8.0840707400991E-15, 11.373484278164979,
    8.2162463710714882E-13, -7.1883948927205012E-15, 11.042050228926119,
    8.2144042128201522E-13, -6.2927190453419006E-15, 10.710616761760587,
    8.2125620545688163E-13, -5.3974824700807E-15, 10.379187224671494,
    8.2107207997791036E-13, -4.5018066227021E-15, 10.047755170390262,
    8.2088786415277677E-13, -3.6060868481117604E-15, 9.7163236210851434,
    8.2070363929302687E-13, -2.7104549279448994E-15, 9.3848936394357629,
    8.2051943250250957E-13, -1.81486693498978E-15, 9.053464981974745,
    8.2033523474660846E-13, -9.1923501482291989E-16, 8.7220370788065313,
    8.2015102795609105E-13, -2.3603094656059913E-17, 8.3906103396279779,
    8.1996682116557375E-13, -2.9507805243757019E-14, 8.224146854531647E-13,
    13.388763747903601, -2.9045959717735024E-14, 8.224934745178205E-13,
    13.382096029710164, -2.8584114191713017E-14, 8.2257226358247629E-13,
    13.375428311602359, -2.8122268665691022E-14, 8.2265105264713208E-13,
    13.368760593586318, -2.7660423139669021E-14, 8.2272984171178787E-13,
    213.36209287566609, -2.7198804119545018E-14, 8.228085921354116E-13,
    13.030660345504426, -2.6736958593523016E-14, 8.228873812000674E-13,
    12.699224633043077, -2.6275113067501022E-14, 8.2296617026472319E-13,
    12.367789086933549, -2.5813267541479021E-14, 8.23044959329379E-13,
    12.036353790320133, -2.535142201545702E-14, 8.2312374839403477E-13,
    11.704918826347136, -2.4889576489435018E-14, 8.2320253745869067E-13,
    11.373484278164979, -2.4427730963413017E-14, 8.2328132652334646E-13,
    11.042050228926119, -2.396588543739102E-14, 8.2336011558800215E-13,
    10.710616761760587, -2.3504266417267016E-14, 8.2343886601162578E-13,
    10.379187224671494, -2.3042420891245021E-14, 8.2351765507628167E-13,
    10.047755170390262, -2.2580552714633217E-14, 8.2359644800504061E-13,
    9.7163236210851434, -2.2118729839201019E-14, 8.2367523320559326E-13,
    9.3848936394357629, -2.1656929614358619E-14, 8.2375401454204256E-13,
    9.053464981974745, -2.1195106738926421E-14, 8.238327997425952E-13,
    8.7220370788065313, -2.0733283863494218E-14, 8.2391158494314765E-13,
    8.3906103396279779, 13.056660989691615, -1.2552980713694741E-14,
    -3.5139602265522615E-14, 13.050160828519424, -1.2555348029519081E-14,
    -3.3154402729028222E-14, 13.043660667426746, -1.2557715345343419E-14,
    -3.1169203192533816E-14, 13.037160506419706, -1.256008266116776E-14,
    -2.9184003656039417E-14, 13.030660345504424, -1.25624499769921E-14,
    -2.7198804119545018E-14, 213.024163368014, -1.256481613179838E-14,
    -2.5214578197340217E-14, 12.69289512944794, -1.256718344762272E-14,
    -2.322937866084582E-14, 12.361626974087356, -1.2569550763447059E-14,
    -2.1244179124351418E-14, 12.030358985078585, -1.25719180792714E-14,
    -1.9258979587857019E-14, 11.699091245565935, -1.2574285395095739E-14,
    -1.727378005136262E-14, 11.3678238386937, -1.2576652710920079E-14,
    -1.5288580514868223E-14, 11.036556847612303, -1.2579020026744421E-14,
    -1.3303380978373821E-14, 10.705290355474208, -1.2581387342568761E-14,
    -1.1318181441879422E-14, 10.37402762873846, -1.2583753497375038E-14,
    -9.33395551967462E-15, 10.042762383244991, -1.2586120813199381E-14,
    -7.34875598318022E-15, 9.71149756773162, -1.2588488245125525E-14,
    -5.363459085256861E-15, 9.3802342204106779, -1.259085544484806E-14,
    -3.3783569101914196E-15, 9.04897210597735, -1.2593222528468788E-14,
    -1.3933520965549404E-15, 8.71771067084897, -1.2595589728191322E-14,
    5.9175007851049945E-16, 8.3864503165639057, -1.2597956927913856E-14,
    2.5768522535759397E-15, -1.7036425840725581E-14, 13.056660989691615,
    8.2353817336928925E-13, -1.614201616871238E-14, 13.050160828519424,
    8.2335577806081982E-13, -1.5247606496699179E-14, 13.043660667426746,
    8.2317338275235038E-13, -1.4353196824685978E-14, 13.037160506419706,
    8.2299098744388094E-13, -1.345878715267278E-14, 13.030660345504424,
    8.228085921354116E-13, -1.256481613179838E-14, 213.024163368014,
    8.2262628628025683E-13, -1.1670406459785179E-14, 12.69289512944794,
    8.224438909717874E-13, -1.077599678777198E-14, 12.361626974087356,
    8.2226149566331806E-13, -9.88158711575878E-15, 12.030358985078585,
    8.2207910035484852E-13, -8.98717744374558E-15, 11.699091245565935,
    8.2189670504637918E-13, -8.09276777173238E-15, 11.3678238386937,
    8.2171430973790974E-13, -7.19835809971918E-15, 11.036556847612303,
    8.2153191442944041E-13, -6.30394842770598E-15, 10.705290355474208,
    8.21349519120971E-13, -5.4099774068315805E-15, 10.37402762873846,
    8.211672132658162E-13, -4.51556773481838E-15, 10.042762383244991,
    8.2098481795734676E-13, -3.6211141976913E-15, 9.71149756773162,
    8.2080241370354585E-13, -2.7267483907919795E-15, 9.3802342204106779,
    8.20620027340408E-13, -1.83242644900654E-15, 9.04897210597735,
    8.204376499226015E-13, -9.3806064210721983E-16, 8.71771067084897,
    8.2025526355946354E-13, -4.3694835207899946E-17, 8.3864503165639057,
    8.2007287719632558E-13, -2.7298241326691417E-14, 8.222232779075367E-13,
    13.056660989691615, -2.6881467820634218E-14, 8.2230388748885665E-13,
    13.050160828519424, -2.6464694314577013E-14, 8.2238449707017671E-13,
    13.043660667426746, -2.604792080851982E-14, 8.2246510665149666E-13,
    13.037160506419706, -2.5631147302462621E-14, 8.2254571623281671E-13,
    13.030660345504424, -2.5214578197340217E-14, 8.2262628628025683E-13,
    213.024163368014, -2.4797804691283017E-14, 8.2270689586157679E-13,
    12.69289512944794, -2.4381031185225818E-14, 8.2278750544289684E-13,
    12.361626974087356, -2.3964257679168616E-14, 8.2286811502421679E-13,
    12.030358985078585, -2.3547484173111417E-14, 8.2294872460553685E-13,
    11.699091245565935, -2.3130710667054218E-14, 8.230293341868568E-13,
    11.3678238386937, -2.2713937160997019E-14, 8.2310994376817685E-13,
    11.036556847612303, -2.2297163654939817E-14, 8.231905533494967E-13,
    10.705290355474208, -2.1880594549817416E-14, 8.2327112339693682E-13,
    10.37402762873846, -2.146382104376022E-14, 8.2335173297825688E-13,
    10.042762383244991, -2.1047027097609538E-14, 8.2343234651296479E-13,
    9.71149756773162, -2.0630274031645819E-14, 8.2351295214089688E-13,
    9.3802342204106779, -2.0213541405775575E-14, 8.2359355381544081E-13,
    9.04897210597735, -1.9796788339811859E-14, 8.236741594433729E-13,
    8.71771067084897, -1.9380035273848139E-14, 8.2375476507130469E-13,
    8.3864503165639057, 12.724555055265856, -1.165223767972134E-14,
    -3.4497380196403018E-14, 12.718222449589355, -1.165587179219438E-14,
    -3.2557274795683019E-14, 12.711889843992372, -1.165950590466742E-14,
    -3.0617169394963014E-14, 12.705557238474906, -1.166314001714046E-14,
    -2.8677063994243015E-14, 12.699224633043077, -1.16667741296135E-14,
    -2.673695859352302E-14, 12.69289512944794, -1.1670406459785181E-14,
    -2.4797804691283017E-14, 212.68656252361242, -1.1674040572258219E-14,
    -2.2857699290563019E-14, 12.355461840621572, -1.1677674684731259E-14,
    -2.091759388984302E-14, 12.024361240836198, -1.1681308797204299E-14,
    -1.8977488489123021E-14, 11.693260807402639, -1.1684942909677339E-14,
    -1.7037383088403022E-14, 11.362160623465199, -1.1688577022150379E-14,
    -1.5097277687683023E-14, 11.031060772168173, -1.169221113462342E-14,
    -1.3157172286963022E-14, 10.699961336661989, -1.169584524709646E-14,
    -1.1217066886243023E-14, 10.368865501844034, -1.1699477577268138E-14,
    -9.277912984003019E-15, 10.037767146763196, -1.1703111689741181E-14,
    -7.3378075832830187E-15, 9.7066691466704427, -1.1706745980444355E-14,
    -5.3976070327150205E-15, 9.3755725152988045, -1.1710379914687259E-14,
    -3.4575967818430195E-15, 9.04447702551003, -1.1714013670700029E-14,
    -1.5176816808190202E-15, 8.7133821400423468, -1.171764760494293E-14,
    4.2232857005297966E-16, 8.3822882522716373, -1.1721281539185835E-14,
    2.3623388209249795E-15, -1.703245979458972E-14, 12.724555055265856,
    8.2360967679414463E-13, -1.613931691922522E-14, 12.718222449589355,
    8.2342910289562534E-13, -1.5246174043860719E-14, 12.711889843992372,
    8.23248528997106E-13, -1.4353031168496218E-14, 12.705557238474906,
    8.2306795509858668E-13, -1.3459888293131721E-14, 12.699224633043077,
    8.228873812000674E-13, -1.256718344762272E-14, 12.69289512944794,
    8.2270689586157679E-13, -1.1674040572258219E-14, 212.68656252361242,
    8.225263219630575E-13, -1.078089769689372E-14, 12.355461840621572,
    8.2234574806453822E-13, -9.88775482152922E-15, 12.024361240836198,
    8.2216517416601884E-13, -8.99461194616472E-15, 11.693260807402639,
    8.2198460026749955E-13, -8.10146907080022E-15, 11.362160623465199,
    8.2180402636898027E-13, -7.2083261954357218E-15, 11.031060772168173,
    8.21623452470461E-13, -6.31518332007122E-15, 10.699961336661989,
    8.2144287857194171E-13, -5.42247847456222E-15, 10.368865501844034,
    8.212623932334511E-13, -4.52933559919772E-15, 10.037767146763196,
    8.2108181933493171E-13, -3.6361489208476704E-15, 9.7066691466704427,
    8.2090123658040959E-13, -2.7430498484687195E-15, 9.3755725152988045,
    8.2072067153789315E-13, -1.84999457907532E-15, 9.04447702551003,
    8.2054011535137964E-13, -9.5689550669636989E-16, 8.7133821400423468,
    8.203595503088632E-13, -6.3796434317419931E-17, 8.3822882522716373,
    8.2017898526634675E-13, -2.5087593227134614E-14, 8.2203177644260268E-13,
    12.724555055265856, -2.4715913856851815E-14, 8.2211420743387279E-13,
    12.718222449589355, -2.4344234486569017E-14, 8.221966384251429E-13,
    12.711889843992372, -2.3972555116286218E-14, 8.22279069416413E-13,
    12.705557238474906, -2.3600875746003419E-14, 8.2236150040768311E-13,
    12.699224633043077, -2.3229378660845817E-14, 8.224438909717874E-13,
    12.69289512944794, -2.2857699290563019E-14, 8.225263219630575E-13,
    212.68656252361242, -2.2486019920280217E-14, 8.2260875295432761E-13,
    12.355461840621572, -2.2114340549997415E-14, 8.2269118394559772E-13,
    12.024361240836198, -2.1742661179714617E-14, 8.2277361493686783E-13,
    11.693260807402639, -2.1370981809431818E-14, 8.2285604592813793E-13,
    11.362160623465199, -2.0999302439149016E-14, 8.22938476919408E-13,
    11.031060772168173, -2.0627623068866217E-14, 8.23020907910678E-13,
    10.699961336661989, -2.0256125983708615E-14, 8.2310329847478233E-13,
    10.368865501844034, -1.988444661342582E-14, 8.2318572946605244E-13,
    10.037767146763196, -1.9512749014630496E-14, 8.23268164500039E-13,
    9.7066691466704427, -1.914108787286022E-14, 8.2335059144859265E-13,
    9.3755725152988045, -1.8769444959602456E-14, 8.2343301435442956E-13,
    9.04447702551003, -1.8397783817832179E-14, 8.2351544130298306E-13,
    8.7133821400423468, -1.8026122676061897E-14, 8.2359786825153647E-13,
    8.3822882522716373, 12.392449619783395, -1.0751494645747939E-14,
    -3.3855158127283415E-14, 12.386284486450126, -1.075639555486968E-14,
    -3.1960146862337817E-14, 12.380119353198415, -1.0761296463991418E-14,
    -3.0065135597392218E-14, 12.373954220026224, -1.0766197373113158E-14,
    -2.8170124332446617E-14, 12.367789086933547, -1.0771098282234899E-14,
    -2.6275113067501019E-14, 12.361626974087356, -1.077599678777198E-14,
    -2.4381031185225815E-14, 12.355461840621572, -1.0780897696893718E-14,
    -2.248601992028022E-14, 212.34929670721084, -1.0785798606015459E-14,
    -2.0591008655334619E-14, 12.018363496644785, -1.0790699515137199E-14,
    -1.869599739038902E-14, 11.6874303692842, -1.0795600424258937E-14,
    -1.6800986125443422E-14, 11.356497408275434, -1.0800501333380678E-14,
    -1.490597486049782E-14, 11.025564696762782, -1.0805402242502421E-14,
    -1.301096359555222E-14, 10.694632317890548, -1.081030315162416E-14,
    -1.1115952330606622E-14, 10.363703374969997, -1.0815201657161237E-14,
    -9.2218704483314188E-15, 10.0327719102814, -1.0820102566282979E-14,
    -7.3268591833858189E-15, 9.7018407256092676, -1.0825003715763184E-14,
    -5.43175498017318E-15, 9.3709108101869329, -1.082990438452646E-14,
    -3.5368366534946198E-15, 9.03998194504271, -1.0834804812931267E-14,
    -1.6420112650831002E-15, 8.709053609235724, -1.083970548169454E-14,
    2.5290706159545948E-16, 8.37812618797937, -1.0844606150457814E-14,
    2.1478253882740194E-15, -1.702849374845386E-14, 12.392449619783395,
    8.23681180219E-13, -1.613661766973806E-14, 12.386284486450126,
    8.2350242773043087E-13, -1.5244741591022259E-14, 12.380119353198415,
    8.2332367524186154E-13, -1.4352865512306459E-14, 12.373954220026224,
    8.2314492275329242E-13, -1.346098943359066E-14, 12.367789086933547,
    8.2296617026472319E-13, -1.2569550763447061E-14, 12.361626974087356,
    8.2278750544289684E-13, -1.1677674684731259E-14, 12.355461840621572,
    8.2260875295432761E-13, -1.078579860601546E-14, 212.34929670721084,
    8.2243000046575838E-13, -9.89392252729966E-15, 12.018363496644785,
    8.2225124797718915E-13, -9.00204644858386E-15, 11.6874303692842,
    8.2207249548861993E-13, -8.11017036986806E-15, 11.356497408275434,
    8.218937430000508E-13, -7.2182942911522615E-15, 11.025564696762782,
    8.2171499051148157E-13, -6.32641821243646E-15, 10.694632317890548,
    8.2153623802291244E-13, -5.4349795422928604E-15, 10.363703374969997,
    8.21357573201086E-13, -4.54310346357706E-15, 10.0327719102814,
    8.2117882071251677E-13, -3.65118364400404E-15, 9.7018407256092676,
    8.2100005945727324E-13, -2.7593513061454595E-15, 9.3709108101869329,
    8.2082131573537831E-13, -1.8675627091441E-15, 9.03998194504271,
    8.2064258078015778E-13, -9.7573037128552E-16, 8.709053609235724,
    8.2046383705826275E-13, -8.3898033426939964E-17, 8.37812618797937,
    8.2028509333636792E-13, -2.2876945127577817E-14, 8.2184027497766867E-13,
    12.392449619783395, -2.2550359893069419E-14, 8.2192452737888883E-13,
    12.386284486450126, -2.2223774658561015E-14, 8.2200877978010909E-13,
    12.380119353198415, -2.1897189424052616E-14, 8.2209303218132925E-13,
    12.373954220026224, -2.1570604189544218E-14, 8.2217728458254952E-13,
    12.367789086933547, -2.1244179124351418E-14, 8.22261495663318E-13,
    12.361626974087356, -2.091759388984302E-14, 8.2234574806453822E-13,
    12.355461840621572, -2.0591008655334615E-14, 8.2243000046575838E-13,
    212.34929670721084, -2.0264423420826214E-14, 8.2251425286697865E-13,
    12.018363496644785, -1.9937838186317816E-14, 8.2259850526819881E-13,
    11.6874303692842, -1.9611252951809415E-14, 8.22682757669419E-13,
    11.356497408275434, -1.9284667717301016E-14, 8.2276701007063923E-13,
    11.025564696762782, -1.8958082482792615E-14, 8.2285126247185929E-13,
    10.694632317890548, -1.8631657417599815E-14, 8.2293547355262784E-13,
    10.363703374969997, -1.8305072183091417E-14, 8.23019725953848E-13,
    10.0327719102814, -1.7978470931651457E-14, 8.231039824871134E-13,
    9.7018407256092676, -1.7651901714074617E-14, 8.2318823075628842E-13,
    9.3709108101869329, -1.7325348513429336E-14, 8.232724748934183E-13,
    9.03998194504271, -1.6998779295852497E-14, 8.2335672316259332E-13,
    8.709053609235724, -1.6672210078275657E-14, 8.2344097143176814E-13,
    8.37812618797937, 12.060344766374264, -9.85075161177454E-15,
    -3.3212936058163819E-14, 12.054347022254197, -9.85691931754498E-15,
    -3.1363018928992621E-14, 12.048349278195296, -9.86308702331542E-15,
    -2.9513101799821416E-14, 12.042351534217957, -9.86925472908586E-15,
    -2.7663184670650219E-14, 12.036353790320133, -9.8754224348563E-15,
    -2.5813267541479017E-14, 12.030358985078587, -9.88158711575878E-15,
    -2.3964257679168616E-14, 12.024361240836198, -9.8877548215292192E-15,
    -2.2114340549997421E-14, 12.018363496644787, -9.89392252729966E-15,
    -2.026442342082622E-14, 212.01236575250843, -9.9000902330701E-15,
    -1.8414506291655019E-14, 11.681599931216738, -9.9062579388405391E-15,
    -1.6564589162483821E-14, 11.350834193130524, -9.91242564461098E-15,
    -1.4714672033312624E-14, 11.020068621396126, -9.91859335038142E-15,
    -1.2864754904141421E-14, 10.689303299157846, -9.924761056151859E-15,
    -1.1014837774970221E-14, 10.358541248136738, -9.9309257370543388E-15,
    -9.1658279126598186E-15, 10.027776673819995, -9.9370934428247808E-15,
    -7.3159107834886191E-15, 9.69701230454809, -9.9432614510820153E-15,
    -5.4659029276313404E-15, 9.36624910507506, -9.94942885436566E-15,
    -3.6160765251462193E-15, 9.03548686457539, -9.9555959551625079E-15,
    -1.7663408493471802E-15, 8.7047250784291, -9.9617633584461513E-15,
    8.34855531379397E-17, 8.3739641236871023, -9.9679307617297946E-15,
    1.9333119556230596E-15, -1.7024527702318E-14, 12.060344766374264,
    8.2375268364385547E-13, -1.61339184202509E-14, 12.054347022254197,
    8.235757525652364E-13, -1.52433091381838E-14, 12.048349278195296,
    8.2339882148661713E-13, -1.43526998561167E-14, 12.042351534217957,
    8.2322189040799805E-13, -1.3462090574049601E-14, 12.036353790320133,
    8.23044959329379E-13, -1.25719180792714E-14, 12.030358985078587,
    8.2286811502421679E-13, -1.1681308797204299E-14, 12.024361240836198,
    8.2269118394559772E-13, -1.07906995151372E-14, 12.018363496644787,
    8.2251425286697865E-13, -9.9000902330701E-15, 212.01236575250843,
    8.2233732178835947E-13, -9.009480951003E-15, 11.681599931216738,
    8.221603907097404E-13, -8.1188716689359E-15, 11.350834193130524,
    8.2198345963112133E-13, -7.2282623868688012E-15, 11.020068621396126,
    8.2180652855250215E-13, -6.3376531048017E-15, 10.689303299157846,
    8.2162959747388308E-13, -5.4474806100235E-15, 10.358541248136738,
    8.2145275316872089E-13, -4.5568713279564E-15, 10.027776673819995,
    8.2127582209010172E-13, -3.66621836716041E-15, 9.69701230454809,
    8.2109888233413688E-13, -2.7756527638221995E-15, 9.36624910507506,
    8.2092195993286357E-13, -1.88513083921288E-15, 9.03548686457539,
    8.2074504620893591E-13, -9.9456523587466982E-16, 8.7047250784291,
    8.205681238076624E-13, -1.0399963253645995E-16, 8.3739641236871023,
    8.2039120140638908E-13, -2.0666297028021014E-14, 8.2164877351273465E-13,
    12.060344766374264, -2.0384805929287023E-14, 8.21734847323905E-13,
    12.054347022254197, -2.0103314830553016E-14, 8.2182092113507528E-13,
    12.048349278195296, -1.9821823731819018E-14, 8.219069949462456E-13,
    12.042351534217957, -1.9540332633085017E-14, 8.2199306875741592E-13,
    12.036353790320133, -1.9258979587857019E-14, 8.2207910035484862E-13,
    12.030358985078587, -1.8977488489123021E-14, 8.2216517416601894E-13,
    12.024361240836198, -1.8695997390389017E-14, 8.2225124797718915E-13,
    12.018363496644787, -1.8414506291655016E-14, 8.2233732178835947E-13,
    212.01236575250843, -1.8133015192921018E-14, 8.2242339559952979E-13,
    11.681599931216738, -1.7851524094187018E-14, 8.2250946941070011E-13,
    11.350834193130524, -1.7570032995453017E-14, 8.2259554322187042E-13,
    11.020068621396126, -1.7288541896719016E-14, 8.2268161703304064E-13,
    10.689303299157846, -1.7007188851491015E-14, 8.2276764863047334E-13,
    10.358541248136738, -1.672569775275702E-14, 8.2285372244164366E-13,
    10.027776673819995, -1.6444192848672418E-14, 8.2293980047418766E-13,
    9.69701230454809, -1.6162715555289018E-14, 8.2302587006398419E-13,
    9.36624910507506, -1.5881252067256217E-14, 8.2311193543240695E-13,
    9.03548686457539, -1.559977477387282E-14, 8.2319800502220359E-13,
    8.7047250784291, -1.5318297480489417E-14, 8.2328407461199992E-13,
    8.3739641236871023, 11.728240578168471, -8.9500085778011409E-15,
    -3.2570713989044216E-14, 11.722410140131595, -8.9574430802202816E-15,
    -3.0765890995647418E-14, 11.716579702135476, -8.9648775826394208E-15,
    -2.8961068002250621E-14, 11.710749264200528, -8.97231208505856E-15,
    -2.7156245008853817E-14, 11.704918826347136, -8.9797465874777E-15,
    -2.535142201545702E-14, 11.699091245565937, -8.9871774437455809E-15,
    -2.3547484173111417E-14, 11.693260807402639, -8.99461194616472E-15,
    -2.174266117971462E-14, 11.687430369284202, -9.00204644858386E-15,
    -1.9937838186317819E-14, 11.681599931216738, -9.009480951003E-15,
    -1.8133015192921022E-14, 211.67576949320431, -9.0169154534221392E-15,
    -1.6328192199524221E-14, 11.345170978036592, -9.02434995584128E-15,
    -1.4523369206127424E-14, 11.014572546074326, -9.03178445826042E-15,
    -1.271854621273062E-14, 10.683974280463881, -9.03921896067956E-15,
    -1.0913723219333822E-14, 10.353379121342218, -9.0466498169474385E-15,
    -9.10978537698822E-15, 10.022781437399365, -9.05408431936658E-15,
    -7.30496238359142E-15, 9.6921838835073029, -9.0615191864008465E-15,
    -5.5000508750895E-15, 9.3615873999631862, -9.06895332420486E-15,
    -3.69531639679782E-15, 9.03099178410807, -9.0763870973937484E-15,
    -1.89067043361126E-15, 8.7003965476224767, -9.0838212351977611E-15,
    -8.5935955319580086E-17, 8.3698020593948357, -9.0912553730017753E-15,
    1.7187985229720998E-15, -1.702056165618214E-14, 11.728240578168471,
    8.2382418706871085E-13, -1.613121917076374E-14, 11.722410140131595,
    8.2364907740004183E-13, -1.524187668534534E-14, 11.716579702135476,
    8.2347396773137281E-13, -1.435253419992694E-14, 11.710749264200528,
    8.2329885806270379E-13, -1.346319171450854E-14, 11.704918826347136,
    8.2312374839403477E-13, -1.257428539509574E-14, 11.699091245565937,
    8.2294872460553685E-13, -1.1684942909677339E-14, 11.693260807402639,
    8.2277361493686783E-13, -1.0795600424258941E-14, 11.687430369284202,
    8.2259850526819881E-13, -9.90625793884054E-15, 11.681599931216738,
    8.2242339559952979E-13, -9.0169154534221392E-15, 211.67576949320431,
    8.2224828593086077E-13, -8.12757296800374E-15, 11.345170978036592,
    8.2207317626219175E-13, -7.2382304825853409E-15, 11.014572546074326,
    8.2189806659352283E-13, -6.34888799716694E-15, 10.683974280463881,
    8.2172295692485381E-13, -5.45998167775414E-15, 10.353379121342218,
    8.2154793313635579E-13, -4.5706391923357404E-15, 10.022781437399365,
    8.2137282346768677E-13, -3.68125309031678E-15, 9.6921838835073029,
    8.2119770521100063E-13, -2.7919542214989394E-15, 9.3615873999631862,
    8.2102260413034873E-13, -1.90269896928166E-15, 9.03099178410807,
    8.2084751163771405E-13, -1.0134001004638199E-15, 8.7003965476224767,
    8.2067241055706205E-13, -1.2410123164597993E-16, 8.3698020593948357,
    8.2049730947641015E-13, -1.8455648928464217E-14, 8.2145727204780063E-13,
    11.728240578168471, -1.821925196550462E-14, 8.215451672689211E-13,
    11.722410140131595, -1.7982855002545017E-14, 8.2163306249004148E-13,
    11.716579702135476, -1.7746458039585419E-14, 8.2172095771116185E-13,
    11.710749264200528, -1.7510061076625819E-14, 8.2180885293228232E-13,
    11.704918826347136, -1.727378005136262E-14, 8.2189670504637918E-13,
    11.699091245565937, -1.7037383088403022E-14, 8.2198460026749955E-13,
    11.693260807402639, -1.6800986125443419E-14, 8.2207249548862E-13,
    11.687430369284202, -1.6564589162483818E-14, 8.221603907097404E-13,
    11.681599931216738, -1.6328192199524218E-14, 8.2224828593086077E-13,
    211.67576949320431, -1.6091795236564617E-14, 8.2233618115198124E-13,
    11.345170978036592, -1.585539827360502E-14, 8.2242407637310161E-13,
    11.014572546074326, -1.561900131064542E-14, 8.2251197159422188E-13,
    10.683974280463881, -1.5382720285382217E-14, 8.2259982370831885E-13,
    10.353379121342218, -1.514632332242262E-14, 8.2268771892943922E-13,
    10.022781437399365, -1.4909914765693379E-14, 8.2277561846126191E-13,
    9.6921838835073029, -1.4673529396503419E-14, 8.2286350937168007E-13,
    9.3615873999631862, -1.4437155621083097E-14, 8.2295139597139569E-13,
    9.03099178410807, -1.4200770251893141E-14, 8.2303928688181375E-13,
    8.7003965476224767, -1.3964384882703178E-14, 8.231271777922317E-13,
    8.3698020593948357, 11.39613713831646, -8.04926554382774E-15,
    -3.1928491919924619E-14, 11.390473923212333, -8.057966842895581E-15,
    -3.0168763062302216E-14, 11.384810708148985, -8.0666681419634189E-15,
    -2.8409034204679819E-14, 11.379147493126398, -8.07536944103126E-15,
    -2.6649305347057415E-14, 11.373484278164979, -8.0840707400991E-15,
    -2.4889576489435018E-14, 11.367823838693702, -8.09276777173238E-15,
    -2.3130710667054215E-14, 11.362160623465199, -8.10146907080022E-15,
    -2.1370981809431818E-14, 11.356497408275434, -8.1101703698680588E-15,
    -1.9611252951809421E-14, 11.350834193130524, -8.1188716689359E-15,
    -1.7851524094187021E-14, 11.345170978036592, -8.12757296800374E-15,
    -1.6091795236564621E-14, 211.3395077629977, -8.1362742670715787E-15,
    -1.4332066378942224E-14, 11.009076470803505, -8.1449755661394213E-15,
    -1.257233752131982E-14, 10.678645261814772, -8.1536768652072592E-15,
    -1.0812608663697422E-14, 10.348216994586437, -8.1623738968405382E-15,
    -9.05374284131662E-15, 10.017786201017476, -8.17107519590838E-15,
    -7.29401398369422E-15, 9.6873554625072931, -8.1797769217196761E-15,
    -5.53419882254766E-15, 9.3569256948717019, -8.18847779404406E-15,
    -3.7745562684494192E-15, 9.0264967036407491, -8.1971782396249873E-15,
    -2.0150000178753402E-15, 8.696068016815854, -8.2058791119493708E-15,
    -2.5535746377710027E-16, 8.3656399951025673, -8.2145799842737544E-15,
    1.5042850903211397E-15, -1.701659561004628E-14, 11.39613713831646,
    8.2389569049356622E-13, -1.612851992127658E-14, 11.390473923212333,
    8.2372240223484736E-13, -1.524044423250688E-14, 11.384810708148985,
    8.2354911397612839E-13, -1.4352368543737178E-14, 11.379147493126398,
    8.2337582571740953E-13, -1.3464292854967481E-14, 11.373484278164979,
    8.2320253745869056E-13, -1.257665271092008E-14, 11.367823838693702,
    8.230293341868568E-13, -1.1688577022150379E-14, 11.362160623465199,
    8.2285604592813793E-13, -1.0800501333380681E-14, 11.356497408275434,
    8.2268275766941907E-13, -9.91242564461098E-15, 11.350834193130524,
    8.225094694107E-13, -9.02434995584128E-15, 11.345170978036592,
    8.2233618115198114E-13, -8.13627426707158E-15, 211.3395077629977,
    8.2216289289326228E-13, -7.24819857830188E-15, 11.009076470803505,
    8.2198960463454341E-13, -6.36012288953218E-15, 10.678645261814772,
    8.2181631637582455E-13, -5.47248274548478E-15, 10.348216994586437,
    8.2164311310399068E-13, -4.58440705671508E-15, 10.017786201017476,
    8.2146982484527172E-13, -3.6962878134731506E-15, 9.6873554625072931,
    8.2129652808786428E-13, -2.8082556791756794E-15, 9.3569256948717019,
    8.21123248327834E-13, -1.9202670993504398E-15, 9.0264967036407491,
    8.2094997706649209E-13, -1.03223496505297E-15, 8.696068016815854,
    8.207766973064617E-13, -1.4420283075549994E-16, 8.3656399951025673,
    8.2060341754643132E-13, -1.6245000828907414E-14, 8.2126577058286671E-13,
    11.39613713831646, -1.6053698001722218E-14, 8.2135548721393714E-13,
    11.390473923212333, -1.5862395174537014E-14, 8.2144520384500767E-13,
    11.384810708148985, -1.5671092347351818E-14, 8.2153492047607819E-13,
    11.379147493126398, -1.5479789520166618E-14, 8.2162463710714872E-13,
    11.373484278164979, -1.528858051486822E-14, 8.2171430973790974E-13,
    11.367823838693702, -1.5097277687683017E-14, 8.2180402636898027E-13,
    11.362160623465199, -1.4905974860497817E-14, 8.218937430000508E-13,
    11.356497408275434, -1.4714672033312617E-14, 8.2198345963112133E-13,
    11.350834193130524, -1.4523369206127417E-14, 8.2207317626219185E-13,
    11.345170978036592, -1.4332066378942217E-14, 8.2216289289326228E-13,
    211.3395077629977, -1.4140763551757017E-14, 8.2225260952433281E-13,
    11.009076470803505, -1.3949460724571817E-14, 8.2234232615540323E-13,
    10.678645261814772, -1.3758251719273418E-14, 8.2243199878616436E-13,
    10.348216994586437, -1.356694889208822E-14, 8.2252171541723478E-13,
    10.017786201017476, -1.3375636682714337E-14, 8.2261143644833627E-13,
    9.6873554625072931, -1.3184343237717818E-14, 8.2270114867937584E-13,
    9.3569256948717019, -1.2993059174909977E-14, 8.2279085651038444E-13,
    9.0264967036407491, -1.2801765729913459E-14, 8.22880568741424E-13,
    8.696068016815854, -1.2610472284916938E-14, 8.2297028097246337E-13,
    8.3656399951025673, 11.06403452996865, -7.14852250985434E-15,
    -3.1286269850805016E-14, 11.05853845464685, -7.15849060557088E-15,
    -2.957163512895702E-14, 11.053042379365833, -7.16845870128742E-15,
    -2.7857000407109017E-14, 11.047546304125596, -7.17842679700396E-15,
    -2.6142365685261017E-14, 11.042050228926119, -7.1883948927205012E-15,
    -2.4427730963413017E-14, 11.036556847612307, -7.19835809971918E-15,
    -2.2713937160997016E-14, 11.031060772168175, -7.20832619543572E-15,
    -2.0999302439149019E-14, 11.025564696762782, -7.21829429115226E-15,
    -1.928466771730102E-14, 11.020068621396126, -7.2282623868688012E-15,
    -1.757003299545302E-14, 11.014572546074326, -7.23823048258534E-15,
    -1.585539827360502E-14, 11.009076470803505, -7.2481985783018791E-15,
    -1.4140763551757022E-14, 211.00358039558773, -7.2581666740184219E-15,
    -1.2426128829909021E-14, 10.67331624321664, -7.26813476973496E-15,
    -1.0711494108061023E-14, 10.343054867875509, -7.27809797673364E-15,
    -8.99770030564502E-15, 10.012790964674323, -7.28806607245018E-15,
    -7.28306558379702E-15, 9.6825270415460221, -7.2980346570385057E-15,
    -5.56834677000582E-15, 9.3522639898209956, -7.30800226388326E-15,
    -3.8537961401010195E-15, 9.022001623193816, -7.3179693818562278E-15,
    -2.13932960213942E-15, 8.6917394860092312, -7.3279369887009822E-15,
    -4.2477897223462025E-16, 8.3614779308103, -7.3379045955457351E-15,
    1.28977165767018E-15, -1.701262956391042E-14, 11.06403452996865,
    8.2396719391842159E-13, -1.612582067178942E-14, 11.05853845464685,
    8.2379572706965289E-13, -1.5239011779668418E-14, 11.053042379365833,
    8.23624260220884E-13, -1.4352202887547418E-14, 11.047546304125596,
    8.2345279337211517E-13, -1.346539399542642E-14, 11.042050228926119,
    8.2328132652334636E-13, -1.257902002674442E-14, 11.036556847612307,
    8.2310994376817685E-13, -1.1692211134623419E-14, 11.031060772168175,
    8.22938476919408E-13, -1.0805402242502419E-14, 11.025564696762782,
    8.2276701007063923E-13, -9.91859335038142E-15, 11.020068621396126,
    8.2259554322187032E-13, -9.03178445826042E-15, 11.014572546074326,
    8.2242407637310151E-13, -8.14497556613942E-15, 11.009076470803505,
    8.2225260952433281E-13, -7.2581666740184219E-15, 211.00358039558773,
    8.22081142675564E-13, -6.37135778189742E-15, 10.67331624321664,
    8.2190967582679519E-13, -5.48498381321542E-15, 10.343054867875509,
    8.2173829307162558E-13, -4.59817492109442E-15, 10.012790964674323,
    8.2156682622285677E-13, -3.71132253662952E-15, 9.6825270415460221,
    8.21395350964728E-13, -2.8245571368524198E-15, 9.3522639898209956,
    8.2122389252531915E-13, -1.93783522941922E-15, 9.022001623193816,
    8.2105244249527022E-13, -1.0510698296421198E-15, 8.6917394860092312,
    8.2088098405586125E-13, -1.6430442986501993E-16, 8.3614779308103,
    8.2070952561645249E-13, -1.4034352729350619E-14, 8.210742691179327E-13,
    11.06403452996865, -1.3888144037939823E-14, 8.2116580715895328E-13,
    11.05853845464685, -1.3741935346529019E-14, 8.2125734519997386E-13,
    11.053042379365833, -1.3595726655118219E-14, 8.2134888324099444E-13,
    11.047546304125596, -1.344951796370742E-14, 8.2144042128201512E-13,
    11.042050228926119, -1.330338097837382E-14, 8.2153191442944041E-13,
    11.036556847612307, -1.315717228696302E-14, 8.21623452470461E-13,
    11.031060772168175, -1.301096359555222E-14, 8.2171499051148157E-13,
    11.025564696762782, -1.2864754904141419E-14, 8.2180652855250225E-13,
    11.020068621396126, -1.271854621273062E-14, 8.2189806659352283E-13,
    11.014572546074326, -1.257233752131982E-14, 8.2198960463454341E-13,
    11.009076470803505, -1.2426128829909021E-14, 8.22081142675564E-13,
    211.00358039558773, -1.227992013849822E-14, 8.2217268071658448E-13,
    10.67331624321664, -1.213378315316462E-14, 8.2226417386400986E-13,
    10.343054867875509, -1.198757446175382E-14, 8.2235571190503044E-13,
    10.012790964674323, -1.1841358599735299E-14, 8.2244725443541052E-13,
    9.6825270415460221, -1.1695157078932219E-14, 8.2253878798707161E-13,
    9.3522639898209956, -1.1548962728736858E-14, 8.2263031704937319E-13,
    9.022001623193816, -1.140276120793378E-14, 8.2272185060103427E-13,
    8.6917394860092312, -1.12565596871307E-14, 8.2281338415269515E-13,
    8.3614779308103, 10.731932836275458, -6.24777947588094E-15,
    -3.064404778168542E-14, 10.726603817585572, -6.25901436824618E-15,
    -2.8974507195611817E-14, 10.721274798936461, -6.27024926061142E-15,
    -2.7304966609538218E-14, 10.715945780328134, -6.28148415297666E-15,
    -2.5635426023464619E-14, 10.710616761760587, -6.2927190453419E-15,
    -2.396588543739102E-14, 10.705290355474208, -6.30394842770598E-15,
    -2.2297163654939817E-14, 10.699961336661989, -6.3151833200712195E-15,
    -2.0627623068866221E-14, 10.694632317890548, -6.3264182124364595E-15,
    -1.8958082482792618E-14, 10.689303299157846, -6.3376531048017E-15,
    -1.7288541896719019E-14, 10.683974280463879, -6.3488879971669394E-15,
    -1.561900131064542E-14, 10.678645261814772, -6.3601228895321794E-15,
    -1.3949460724571822E-14, 10.673316243216638, -6.37135778189742E-15,
    -1.2279920138498221E-14, 210.66798722467357, -6.38259267426266E-15,
    -1.0610379552424622E-14, 10.337892741215558, -6.39382205662674E-15,
    -8.9416577699734192E-15, 10.007795728376026, -6.40505694899198E-15,
    -7.27211718389982E-15, 9.6776986206234881, -6.4162923923573352E-15,
    -5.60249471746398E-15, 9.3476022848090281, -6.4275267337224606E-15,
    -3.93303601175262E-15, 9.0175065427876628, -6.4387605240874682E-15,
    -2.2636591864035E-15, 8.6874109552229974, -6.4499948654525912E-15,
    -5.9420048069214023E-16, 8.3573158665180323, -6.4612292068177158E-15,
    1.0752582250192198E-15, -1.7008663517774559E-14, 10.731932836275458,
    8.2403869734327707E-13, -1.612312142230226E-14, 10.726603817585572,
    8.2386905190445831E-13, -1.5237579326829958E-14, 10.721274798936461,
    8.2369940646563956E-13, -1.4352037231357659E-14, 10.715945780328134,
    8.235297610268209E-13, -1.346649513588536E-14, 10.710616761760587,
    8.2336011558800215E-13, -1.2581387342568761E-14, 10.705290355474208,
    8.231905533494968E-13, -1.1695845247096459E-14, 10.699961336661989,
    8.2302090791067815E-13, -1.081030315162416E-14, 10.694632317890548,
    8.2285126247185939E-13, -9.92476105615186E-15, 10.689303299157846,
    8.2268161703304064E-13, -9.03921896067956E-15, 10.683974280463879,
    8.22511971594222E-13, -8.1536768652072592E-15, 10.678645261814772,
    8.2234232615540323E-13, -7.2681347697349616E-15, 10.673316243216638,
    8.2217268071658458E-13, -6.38259267426266E-15, 210.66798722467357,
    8.2200303527776592E-13, -5.4974848809460606E-15, 10.337892741215558,
    8.2183347303926048E-13, -4.61194278547376E-15, 10.007795728376026,
    8.2166382760044172E-13, -3.72635725978589E-15, 9.6776986206234881,
    8.2149417384159167E-13, -2.8408585945291597E-15, 9.3476022848090281,
    8.2132453672280431E-13, -1.955403359488E-15, 9.0175065427876628,
    8.2115490792404836E-13, -1.0699046942312699E-15, 8.6874109552229974,
    8.2098527080526091E-13, -1.8440602897453994E-16, 8.3573158665180323,
    8.2081563368647355E-13, -1.1823704629793819E-14, 8.2088276765299868E-13,
    10.731932836275458, -1.172259007415742E-14, 8.2097612710396942E-13,
    10.726603817585572, -1.1621475518521016E-14, 8.2106948655494005E-13,
    10.721274798936461, -1.1520360962884617E-14, 8.2116284600591079E-13,
    10.715945780328134, -1.1419246407248218E-14, 8.2125620545688153E-13,
    10.710616761760587, -1.131818144187942E-14, 8.21349519120971E-13,
    10.705290355474208, -1.1217066886243018E-14, 8.2144287857194171E-13,
    10.699961336661989, -1.1115952330606619E-14, 8.2153623802291244E-13,
    10.694632317890548, -1.1014837774970218E-14, 8.2162959747388308E-13,
    10.689303299157846, -1.0913723219333818E-14, 8.2172295692485381E-13,
    10.683974280463879, -1.0812608663697419E-14, 8.2181631637582455E-13,
    10.678645261814772, -1.0711494108061018E-14, 8.2190967582679519E-13,
    10.673316243216638, -1.0610379552424619E-14, 8.2200303527776582E-13,
    210.66798722467357, -1.0509314587055818E-14, 8.2209634894185537E-13,
    10.337892741215558, -1.0408200031419419E-14, 8.22189708392826E-13,
    10.007795728376026, -1.0307080516756259E-14, 8.2228307242248478E-13,
    9.6776986206234881, -1.0205970920146619E-14, 8.2237642729476748E-13,
    9.3476022848090281, -1.0104866282563738E-14, 8.2246977758836183E-13,
    9.0175065427876628, -1.0003756685954099E-14, 8.2256313246064443E-13,
    8.6874109552229974, -9.9026470893444588E-15, 8.2265648733292693E-13,
    8.3573158665180323, 10.399835731596752, -5.3474781991581405E-15,
    -3.0002140681702221E-14, 10.394673604804272, -5.3599792668887808E-15,
    -2.8377672115593417E-14, 10.389511478052567, -5.37248033461942E-15,
    -2.6753203549484617E-14, 10.38434935134164, -5.38498140235006E-15,
    -2.5128734983375817E-14, 10.379187224671494, -5.3974824700807E-15,
    -2.3504266417267019E-14, 10.374027628738462, -5.4099774068315805E-15,
    -2.1880594549817416E-14, 10.368865501844036, -5.42247847456222E-15,
    -2.0256125983708619E-14, 10.36370337497, -5.43497954229286E-15,
    -1.8631657417599818E-14, 10.35854124813674, -5.4474806100235008E-15,
    -1.7007188851491021E-14, 10.353379121342218, -5.45998167775414E-15,
    -1.538272028538222E-14, 10.348216994586437, -5.47248274548478E-15,
    -1.3758251719273422E-14, 10.343054867875509, -5.48498381321542E-15,
    -1.213378315316462E-14, 10.33789274121556, -5.4974848809460606E-15,
    -1.0509314587055822E-14, 210.332733145307, -5.50997981769694E-15,
    -8.8856427196062188E-15, 10.002802941240949, -5.52248088542758E-15,
    -7.2611741534974184E-15, 9.6728725672699643, -5.5349825662561957E-15,
    -5.6366259175227E-15, 9.3429428657798645, -5.5474830208888611E-15,
    -4.0122370212798192E-15, 9.013013666788229, -5.5599828624235483E-15,
    -2.3879277949028598E-15, 8.6830845472654339, -5.5724833170562113E-15,
    -7.635388986599801E-16, 8.3531558434539619, -5.5849837716888766E-15,
    8.6084999758289989E-16, -1.700469941673244E-14, 10.399835731596752,
    8.2411016570024385E-13, -1.6120423496625538E-14, 10.394673604804272,
    8.2394234077808936E-13, -1.5236147576518639E-14, 10.389511478052567,
    8.2377451585593477E-13, -1.4351871656411738E-14, 10.38434935134164,
    8.2360669093378027E-13, -1.346759573630484E-14, 10.379187224671494,
    8.2343886601162578E-13, -1.258375349737504E-14, 10.374027628738462,
    8.2327112339693682E-13, -1.1699477577268139E-14, 10.368865501844036,
    8.2310329847478233E-13, -1.0815201657161241E-14, 10.36370337497,
    8.2293547355262784E-13, -9.93092573705434E-15, 10.35854124813674,
    8.2276764863047324E-13, -9.04664981694744E-15, 10.353379121342218,
    8.2259982370831875E-13, -8.16237389684054E-15, 10.348216994586437,
    8.2243199878616425E-13, -7.278097976733641E-15, 10.343054867875509,
    8.2226417386400976E-13, -6.3938220566267407E-15, 10.33789274121556,
    8.2209634894185527E-13, -5.50997981769694E-15, 210.332733145307,
    8.2192860632716631E-13, -4.6257038975900405E-15, 10.002802941240949,
    8.2176078140501172E-13, -3.74138460936543E-15, 9.6728725672699643,
    8.2159294825211064E-13, -2.8571520573762395E-15, 9.3429428657798645,
    8.2142513156070273E-13, -1.9729628735047598E-15, 9.013013666788229,
    8.2125732310004141E-13, -1.0887303215155698E-15, 8.6830845472654339,
    8.210895064086334E-13, -2.0449776952637994E-16, 8.3531558434539619,
    8.2092168971722549E-13, -9.61414071272822E-15, 8.2069136010737068E-13,
    10.399835731596752, -9.55809817705662E-15, 8.2078654007500557E-13,
    10.394673604804272, -9.5020556413850185E-15, 8.2088172004264047E-13,
    10.389511478052567, -9.4460131057134214E-15, 8.2097690001027537E-13,
    10.38434935134164, -9.38997057004182E-15, 8.2107207997791026E-13,
    10.379187224671494, -9.33395551967462E-15, 8.211672132658162E-13,
    10.374027628738462, -9.27791298400302E-15, 8.212623932334511E-13,
    10.368865501844036, -9.22187044833142E-15, 8.21357573201086E-13,
    10.36370337497, -9.16582791265982E-15, 8.2145275316872089E-13,
    10.35854124813674, -9.1097853769882183E-15, 8.2154793313635579E-13,
    10.353379121342218, -9.05374284131662E-15, 8.2164311310399068E-13,
    10.348216994586437, -8.99770030564502E-15, 8.2173829307162558E-13,
    10.343054867875509, -8.9416577699734192E-15, 8.2183347303926048E-13,
    10.33789274121556, -8.8856427196062188E-15, 8.2192860632716631E-13,
    210.332733145307, -8.82960018393462E-15, 8.2202378629480121E-13,
    10.002802941240949, -8.7735548997325784E-15, 8.22118970930409E-13,
    9.6728725672699643, -8.71751511259142E-15, 8.22214146230071E-13,
    9.3429428657798645, -8.6614780739807E-15, 8.2230931686176009E-13,
    9.013013666788229, -8.60543828683954E-15, 8.2240449216142213E-13,
    8.6830845472654339, -8.549398499698379E-15, 8.22499667461084E-13,
    8.3531558434539619, 10.067736115828113, -4.44673516518474E-15,
    -2.9359918612582618E-14, 10.062740879407485, -4.46050302956408E-15,
    -2.7780544182248218E-14, 10.057745643027634, -4.47427089394342E-15,
    -2.6201169751913818E-14, 10.05275040668856, -4.4880387583227591E-15,
    -2.4621795321579418E-14, 10.047755170390262, -4.5018066227021E-15,
    -2.3042420891245018E-14, 10.042762383244991, -4.51556773481838E-15,
    -2.1463821043760217E-14, 10.037767146763196, -4.5293355991977193E-15,
    -1.988444661342582E-14, 10.032771910281403, -4.5431034635770592E-15,
    -1.830507218309142E-14, 10.027776673819995, -4.5568713279564E-15,
    -1.672569775275702E-14, 10.022781437399365, -4.57063919233574E-15,
    -1.514632332242262E-14, 10.017786201017476, -4.5844070567150795E-15,
    -1.3566948892088222E-14, 10.012790964674322, -4.5981749210944193E-15,
    -1.198757446175382E-14, 10.007795728376026, -4.61194278547376E-15,
    -1.0408200031419422E-14, 10.002802941240949, -4.62570389759004E-15,
    -8.8296001839346186E-15, 209.99780770478361, -4.63947176196938E-15,
    -7.2502257536002186E-15, 9.6680441462189837, -4.6532403015750253E-15,
    -5.6707738649808604E-15, 9.338281160667993, -4.66700749072806E-15,
    -4.0914768929314195E-15, 9.008518586316832, -4.680774004654788E-15,
    -2.51225737916694E-15, 8.6787560164567736, -4.6945411938078211E-15,
    -9.3296040711750028E-16, 8.3489937791616935, -4.7083083829608558E-15,
    6.4633656493193974E-16, -1.700073337059658E-14, 10.067736115828113,
    8.2418166912509923E-13, -1.6117724247138378E-14, 10.062740879407485,
    8.2401566561289489E-13, -1.523471512368018E-14, 10.057745643027634,
    8.2384966210069035E-13, -1.4351706000221978E-14, 10.05275040668856,
    8.23683658588486E-13, -1.346869687676378E-14, 10.047755170390262,
    8.2351765507628157E-13, -1.2586120813199381E-14, 10.042762383244991,
    8.2335173297825678E-13, -1.1703111689741179E-14, 10.037767146763196,
    8.2318572946605244E-13, -1.0820102566282981E-14, 10.032771910281403,
    8.23019725953848E-13, -9.93709344282478E-15, 10.027776673819995,
    8.2285372244164356E-13, -9.0540843193665792E-15, 10.022781437399365,
    8.2268771892943912E-13, -8.1710751959083792E-15, 10.017786201017476,
    8.2252171541723478E-13, -7.28806607245018E-15, 10.012790964674322,
    8.2235571190503034E-13, -6.4050569489919807E-15, 10.007795728376026,
    8.22189708392826E-13, -5.52248088542758E-15, 10.002802941240949,
    8.2202378629480121E-13, -4.63947176196938E-15, 209.99780770478361,
    8.2185778278259677E-13, -3.7564193325218E-15, 9.6680441462189837,
    8.2169177112897439E-13, -2.8734535150529794E-15, 9.338281160667993,
    8.21525775758188E-13, -1.9905310035735397E-15, 9.008518586316832,
    8.2135978852881954E-13, -1.1075651861047199E-15, 8.6787560164567736,
    8.21193793158033E-13, -2.2459936863589995E-16, 8.3489937791616935,
    8.2102779778724665E-13, -7.4034926131714175E-15, 8.2049985864243666E-13,
    10.067736115828113, -7.39254421327422E-15, 8.2059686002002171E-13,
    10.062740879407485, -7.3815958133770164E-15, 8.2069386139760666E-13,
    10.057745643027634, -7.37064741347982E-15, 8.2079086277519171E-13,
    10.05275040668856, -7.3596990135826183E-15, 8.2088786415277666E-13,
    10.047755170390262, -7.3487559831802184E-15, 8.2098481795734676E-13,
    10.042762383244991, -7.3378075832830187E-15, 8.2108181933493181E-13,
    10.037767146763196, -7.3268591833858189E-15, 8.2117882071251677E-13,
    10.032771910281403, -7.3159107834886191E-15, 8.2127582209010182E-13,
    10.027776673819995, -7.3049623835914177E-15, 8.2137282346768677E-13,
    10.022781437399365, -7.2940139836942179E-15, 8.2146982484527182E-13,
    10.017786201017476, -7.2830655837970181E-15, 8.2156682622285687E-13,
    10.012790964674322, -7.2721171838998183E-15, 8.2166382760044172E-13,
    10.007795728376026, -7.2611741534974184E-15, 8.2176078140501182E-13,
    10.002802941240949, -7.2502257536002186E-15, 8.2185778278259687E-13,
    209.99780770478361, -7.2392768167535378E-15, 8.2195478891748332E-13,
    9.6680441462189837, -7.2283289538058191E-15, 8.2205178553776687E-13,
    9.338281160667993, -7.2173816278075783E-15, 8.2214877740074883E-13,
    9.008518586316832, -7.20643376485986E-15, 8.2224577402103239E-13,
    8.6787560164567736, -7.19548590191214E-15, 8.2234277064131575E-13,
    8.3489937791616935, 9.7356373050036282, -3.5459479554862806E-15,
    -2.8717665046549376E-14, 9.7308088839628422, -3.5609826786426505E-15,
    -2.7183386963570337E-14, 9.72598046296283, -3.57601740179902E-15,
    -2.5649108880591298E-14, 9.7211520420035988, -3.59105212495539E-15,
    -2.4114830797612256E-14, 9.7163236210851434, -3.6060868481117604E-15,
    -2.2580552714633217E-14, 9.71149756773162, -3.6211141976913E-15,
    -2.1047027097609535E-14, 9.7066691466704427, -3.6361489208476704E-15,
    -1.95127490146305E-14, 9.7018407256092676, -3.65118364400404E-15,
    -1.797847093165146E-14, 9.69701230454809, -3.66621836716041E-15,
    -1.6444192848672421E-14, 9.6921838835073011, -3.68125309031678E-15,
    -1.4909914765693379E-14, 9.6873554625072931, -3.69628781347315E-15,
    -1.3375636682714341E-14, 9.6825270415460221, -3.71132253662952E-15,
    -1.18413585997353E-14, 9.6776986206234881, -3.72635725978589E-15,
    -1.0307080516756262E-14, 9.6728725672699643, -3.74138460936543E-15,
    -8.7735548997325784E-15, 9.6680441462189837, -3.7564193325218E-15,
    -7.23927681675354E-15, 209.66321548842987, -3.7714547930358532E-15,
    -5.7049234871789638E-15, 9.3336192269719067, -3.7864887788345408E-15,
    -4.17072065079546E-15, 9.0040232854127122, -3.8015220272755441E-15,
    -2.6365930610074917E-15, 8.6744272733672823, -3.8165560130742309E-15,
    -1.102390224623988E-15, 8.3448315107466069, -3.8315899988729185E-15,
    4.31812611759516E-16, -1.6996767129951346E-14, 9.7356373050036282,
    8.2425317605674349E-13, -1.6115024865270176E-14, 9.7308088839628422,
    8.2408899404381785E-13, -1.5233282600589005E-14, 9.72598046296283,
    8.23924812030892E-13, -1.4351540335907835E-14, 9.7211520420035988,
    8.2376063001796635E-13, -1.3469798071226666E-14, 9.7163236210851434,
    8.2359644800504061E-13, -1.2588488245125527E-14, 9.71149756773162,
    8.2343234651296479E-13, -1.1706745980444355E-14, 9.7066691466704427,
    8.2326816450003915E-13, -1.0825003715763186E-14, 9.7018407256092676,
    8.231039824871134E-13, -9.9432614510820153E-15, 9.69701230454809,
    8.2293980047418766E-13, -9.0615191864008465E-15, 9.6921838835073011,
    8.2277561846126191E-13, -8.1797769217196761E-15, 9.6873554625072931,
    8.2261143644833627E-13, -7.2980346570385072E-15, 9.6825270415460221,
    8.2244725443541052E-13, -6.416292392357336E-15, 9.6776986206234881,
    8.2228307242248488E-13, -5.5349825662561965E-15, 9.6728725672699643,
    8.2211897093040906E-13, -4.6532403015750261E-15, 9.6680441462189837,
    8.2195478891748322E-13, -3.7714547930358532E-15, 209.66321548842987,
    8.2179059885247251E-13, -2.8897557722126856E-15, 9.3336192269719067,
    8.2162642489163183E-13, -2.0080999952475221E-15, 9.0040232854127122,
    8.2146225898287611E-13, -1.1264009744243549E-15, 8.6744272733672823,
    8.2129808502203533E-13, -2.4470195360118795E-16, 8.3448315107466069,
    8.2113391106119464E-13, -5.1927360953655E-15, 8.2030834778557209E-13,
    9.7356373050036282, -5.2268840428236609E-15, 8.2040717066243574E-13,
    9.7308088839628422, -5.26103199028182E-15, 8.2050599353929948E-13,
    9.72598046296283, -5.2951799377399808E-15, 8.2060481641616313E-13,
    9.7211520420035988, -5.3293278851981395E-15, 8.2070363929302687E-13,
    9.7163236210851434, -5.3634590852568594E-15, 8.2080241370354595E-13,
    9.71149756773162, -5.3976070327150205E-15, 8.2090123658040959E-13,
    9.7066691466704427, -5.43175498017318E-15, 8.2100005945727334E-13,
    9.7018407256092676, -5.4659029276313404E-15, 8.21098882334137E-13,
    9.69701230454809, -5.5000508750894992E-15, 8.2119770521100073E-13,
    9.6921838835073011, -5.5341988225476587E-15, 8.2129652808786438E-13,
    9.6873554625072931, -5.56834677000582E-15, 8.2139535096472812E-13,
    9.6825270415460221, -5.6024947174639794E-15, 8.2149417384159167E-13,
    9.6776986206234881, -5.6366259175226985E-15, 8.2159294825211074E-13,
    9.6728725672699643, -5.6707738649808612E-15, 8.2169177112897449E-13,
    9.6680441462189837, -5.704923487178963E-15, 8.2179059885247261E-13,
    209.66321548842987, -5.73906975989718E-15, 8.2188941688270188E-13,
    9.3336192269719067, -5.7732143578754514E-15, 8.2198823006629657E-13,
    9.0040232854127122, -5.807360630593668E-15, 8.2208704809652584E-13,
    8.6744272733672823, -5.8415069033118829E-15, 8.2218586612675491E-13,
    8.3448315107466069, 9.4035404596997552, -2.64524909723794E-15,
    -2.8075474474343419E-14, 9.3988787545878836, -2.66155055491468E-15,
    -2.6586288315557816E-14, 9.3942170494963975, -2.67785201259142E-15,
    -2.5097102156772217E-14, 9.389555344445693, -2.69415347026816E-15,
    -2.3607915997986618E-14, 9.3848936394357629, -2.7104549279448994E-15,
    -2.2118729839201016E-14, 9.3802342204106779, -2.7267483907919795E-15,
    -2.0630274031645815E-14, 9.3755725152988045, -2.74304984846872E-15,
    -1.914108787286022E-14, 9.3709108101869329, -2.75935130614546E-15,
    -1.765190171407462E-14, 9.36624910507506, -2.7756527638221995E-15,
    -1.6162715555289018E-14, 9.3615873999631862, -2.79195422149894E-15,
    -1.4673529396503419E-14, 9.3569256948717019, -2.8082556791756794E-15,
    -1.3184343237717822E-14, 9.3522639898209938, -2.82455713685242E-15,
    -1.1695157078932219E-14, 9.3476022848090281, -2.84085859452916E-15,
    -1.0205970920146622E-14, 9.3429428657798645, -2.85715205737624E-15,
    -8.71751511259142E-15, 9.338281160667993, -2.87345351505298E-15,
    -7.2283289538058191E-15, 9.3336192269719067, -2.8897557722126852E-15,
    -5.7390697598971795E-15, 209.32895775049727, -2.9060564304064598E-15,
    -4.2499566362346193E-15, 8.9995284254066572, -2.9223562891172677E-15,
    -2.7609165476950996E-15, 8.6700989548496441, -2.9386569473110415E-15,
    -1.2718034240325402E-15, 8.3406696505771585, -2.954957605504816E-15,
    2.1730969963001983E-16, -1.6992801278324859E-14, 9.4035404596997552,
    8.2432467597481008E-13, -1.6112325748164058E-14, 9.3988787545878836,
    8.2416231528250585E-13, -1.523185021800326E-14, 9.3942170494963975,
    8.2399995459020162E-13, -1.4351374687842459E-14, 9.389555344445693,
    8.2383759389789739E-13, -1.347089915768166E-14, 9.3848936394357629,
    8.2367523320559316E-13, -1.259085544484806E-14, 9.3802342204106779,
    8.2351295214089678E-13, -1.1710379914687259E-14, 9.3755725152988045,
    8.2335059144859265E-13, -1.082990438452646E-14, 9.3709108101869329,
    8.2318823075628842E-13, -9.94942885436566E-15, 9.36624910507506,
    8.2302587006398419E-13, -9.0689533242048591E-15, 9.3615873999631862,
    8.2286350937168E-13, -8.18847779404406E-15, 9.3569256948717019,
    8.2270114867937573E-13, -7.3080022638832617E-15, 9.3522639898209938,
    8.2253878798707161E-13, -6.4275267337224606E-15, 9.3476022848090281,
    8.2237642729476738E-13, -5.54748302088886E-15, 9.3429428657798645,
    8.22214146230071E-13, -4.66700749072806E-15, 9.338281160667993,
    8.2205178553776677E-13, -3.78648877883454E-15, 9.3336192269719067,
    8.2188941688270178E-13, -2.9060564304064594E-15, 209.32895775049727,
    8.2172706415315831E-13, -2.0256672637111E-15, 8.9995284254066572,
    8.2156471938637582E-13, -1.1452349152830198E-15, 8.6700989548496441,
    8.2140236665683225E-13, -2.6480256685493997E-16, 8.3406696505771585,
    8.2124001392728889E-13, -2.9821964140578193E-15, 8.2011685571256862E-13,
    9.4035404596997552, -3.0614362857094204E-15, 8.2021749991005389E-13,
    9.3988787545878836, -3.1406761573610176E-15, 8.20318144107539E-13,
    9.3942170494963975, -3.2199160290126203E-15, 8.2041878830502431E-13,
    9.389555344445693, -3.299155900664219E-15, 8.2051943250250947E-13,
    9.3848936394357629, -3.3783569101914192E-15, 8.20620027340408E-13,
    9.3802342204106779, -3.4575967818430195E-15, 8.2072067153789315E-13,
    9.3755725152988045, -3.5368366534946194E-15, 8.2082131573537841E-13,
    9.3709108101869329, -3.6160765251462193E-15, 8.2092195993286357E-13,
    9.36624910507506, -3.6953163967978189E-15, 8.2102260413034883E-13,
    9.3615873999631862, -3.7745562684494192E-15, 8.21123248327834E-13,
    9.3569256948717019, -3.8537961401010187E-15, 8.2122389252531925E-13,
    9.3522639898209938, -3.9330360117526182E-15, 8.2132453672280431E-13,
    9.3476022848090281, -4.0122370212798184E-15, 8.2142513156070283E-13,
    9.3429428657798645, -4.0914768929314195E-15, 8.21525775758188E-13,
    9.338281160667993, -4.1707206507954589E-15, 8.2162642489163183E-13,
    9.3336192269719067, -4.2499566362346193E-15, 8.2172706415315841E-13,
    209.32895775049727, -4.3291887354613391E-15, 8.2182769847872623E-13,
    8.9995284254066572, -4.4084247209004996E-15, 8.2192833774025281E-13,
    8.6700989548496441, -4.4876607063396585E-15, 8.220289770017792E-13,
    8.3406696505771585, 9.0714453037624718, -1.7445944147146601E-15,
    -2.7433315399051097E-14, 9.0669502232951515, -1.7621625447834399E-15,
    -2.5989218952877978E-14, 9.0624551428278313, -1.77973067485222E-15,
    -2.4545122506704858E-14, 9.0579600623809, -1.797298804921E-15,
    -2.3101026060531739E-14, 9.053464981974745, -1.81486693498978E-15,
    -2.1656929614358619E-14, 9.0489721059773522, -1.8324264490065397E-15,
    -2.0213541405775575E-14, 9.04447702551003, -1.84999457907532E-15,
    -1.8769444959602459E-14, 9.0399819450427117, -1.8675627091441E-15,
    -1.7325348513429339E-14, 9.03548686457539, -1.88513083921288E-15,
    -1.588125206725622E-14, 9.03099178410807, -1.90269896928166E-15,
    -1.44371556210831E-14, 9.0264967036407491, -1.9202670993504398E-15,
    -1.2993059174909982E-14, 9.022001623193816, -1.93783522941922E-15,
    -1.1548962728736859E-14, 9.0175065427876628, -1.9554033594880004E-15,
    -1.0104866282563741E-14, 9.013013666788229, -1.9729628735047598E-15,
    -8.6614780739807E-15, 9.008518586316832, -1.99053100357354E-15,
    -7.2173816278075783E-15, 9.004023285412714, -2.0080999952475221E-15,
    -5.7732143578754514E-15, 8.9995284254066572, -2.0256672637111E-15,
    -4.3291887354613391E-15, 208.99503378584149, -2.0432336705694759E-15,
    -2.8852339368062357E-15, 8.665770848614871, -2.0608009390330539E-15,
    -1.441208314392124E-15, 8.33650799453053, -2.0783682074966322E-15,
    2.8173080219879851E-18, -1.6988835621207746E-14, 9.0714453037624718,
    8.2439617238608767E-13, -1.6109626763438987E-14, 9.0669502232951515,
    8.2423563292507642E-13, -1.5230417905670227E-14, 9.0624551428278313,
    8.2407509346406506E-13, -1.4351209047901468E-14, 9.0579600623809,
    8.2391455400305381E-13, -1.3472000190132708E-14, 9.053464981974745,
    8.2375401454204256E-13, -1.2593222528468788E-14, 9.0489721059773522,
    8.2359355381544081E-13, -1.1714013670700027E-14, 9.04447702551003,
    8.2343301435442956E-13, -1.0834804812931267E-14, 9.0399819450427117,
    8.232724748934183E-13, -9.9555959551625079E-15, 9.03548686457539,
    8.2311193543240695E-13, -9.0763870973937484E-15, 9.03099178410807,
    8.2295139597139559E-13, -8.1971782396249873E-15, 9.0264967036407491,
    8.2279085651038434E-13, -7.31796938185623E-15, 9.022001623193816,
    8.2263031704937309E-13, -6.4387605240874682E-15, 9.0175065427876628,
    8.2246977758836183E-13, -5.5599828624235483E-15, 9.013013666788229,
    8.2230931686176009E-13, -4.680774004654788E-15, 9.008518586316832,
    8.2214877740074873E-13, -3.8015220272755441E-15, 9.004023285412714,
    8.2198823006629647E-13, -2.9223562891172677E-15, 8.9995284254066572,
    8.2182769847872623E-13, -2.0432336705694759E-15, 208.99503378584149,
    8.2166717476459689E-13, -1.1640679324111999E-15, 8.665770848614871,
    8.2150664317702645E-13, -2.8490219425292393E-16, 8.33650799453053,
    8.2134611158945621E-13, -7.7176515099925918E-16, 8.1992537303149591E-13,
    9.0714453037624718, -8.9609473526334057E-16, 8.2002783846027395E-13,
    9.0669502232951515, -1.0204243195274188E-15, 8.2013030388905208E-13,
    9.0624551428278313, -1.1447539037915002E-15, 8.2023276931783022E-13,
    9.0579600623809, -1.2690834880555792E-15, 8.2033523474660836E-13,
    9.053464981974745, -1.3933520965549396E-15, 8.204376499226015E-13,
    9.0489721059773522, -1.5176816808190202E-15, 8.2054011535137964E-13,
    9.04447702551003, -1.6420112650830996E-15, 8.2064258078015778E-13,
    9.0399819450427117, -1.7663408493471798E-15, 8.2074504620893591E-13,
    9.03548686457539, -1.8906704336112588E-15, 8.2084751163771405E-13,
    9.03099178410807, -2.01500001787534E-15, 8.2094997706649209E-13,
    9.0264967036407491, -2.13932960213942E-15, 8.2105244249527022E-13,
    9.022001623193816, -2.2636591864034994E-15, 8.2115490792404826E-13,
    9.0175065427876628, -2.3879277949028594E-15, 8.2125732310004141E-13,
    9.013013666788229, -2.51225737916694E-15, 8.2135978852881954E-13,
    9.008518586316832, -2.6365930610074913E-15, 8.2146225898287611E-13,
    9.004023285412714, -2.7609165476950996E-15, 8.2156471938637582E-13,
    8.9995284254066572, -2.8852339368062357E-15, 8.2166717476459689E-13,
    208.99503378584149, -3.0095574234938437E-15, 8.217696351680965E-13,
    8.665770848614871, -3.1338809101814516E-15, 8.21872095571596E-13,
    8.33650799453053, 8.7393512020126352, -8.4389555646632023E-16,
    -2.6791124826845137E-14, 8.7350226712060124, -8.627304210554701E-16,
    -2.5392120304865457E-14, 8.7306941403993878, -8.8156528564462016E-16,
    -2.3993115782885777E-14, 8.7263656095927669, -9.0040015023377E-16,
    -2.2594111260906098E-14, 8.7220370788065313, -9.1923501482291989E-16,
    -2.1195106738926418E-14, 8.71771067084897, -9.3806064210722022E-16,
    -1.9796788339811856E-14, 8.7133821400423468, -9.5689550669637E-16,
    -1.8397783817832179E-14, 8.709053609235724, -9.7573037128552E-16,
    -1.69987792958525E-14, 8.7047250784291, -9.9456523587467E-16,
    -1.559977477387282E-14, 8.7003965476224767, -1.0134001004638199E-15,
    -1.4200770251893141E-14, 8.696068016815854, -1.0322349650529701E-15,
    -1.2801765729913461E-14, 8.6917394860092312, -1.0510698296421202E-15,
    -1.140276120793378E-14, 8.6874109552229974, -1.0699046942312703E-15,
    -1.0003756685954102E-14, 8.6830845472654339, -1.08873032151557E-15,
    -8.6054382868395392E-15, 8.6787560164567736, -1.10756518610472E-15,
    -7.20643376485986E-15, 8.6744272733672823, -1.1264009744243549E-15,
    -5.807360630593668E-15, 8.6700989548496441, -1.1452349152830202E-15,
    -4.4084247209004988E-15, 8.665770848614871, -1.1640679324112003E-15,
    -3.0095574234938437E-15, 208.6614425300952, -1.1829018732698649E-15,
    -1.6106215138006761E-15, 8.3323461343610816, -1.20173581412853E-15,
    -2.1168560410750796E-16, -1.6984869769581263E-14, 8.7393512020126352,
    8.2446767230415426E-13, -1.6106927646332872E-14, 8.7350226712060124,
    8.2430895416376452E-13, -1.5228985523084482E-14, 8.7306941403993878,
    8.2415023602337468E-13, -1.4351043399836092E-14, 8.7263656095927669,
    8.2399151788298494E-13, -1.3473101276587702E-14, 8.7220370788065313,
    8.238327997425951E-13, -1.2595589728191322E-14, 8.71771067084897,
    8.236741594433728E-13, -1.171764760494293E-14, 8.7133821400423468,
    8.2351544130298306E-13, -1.0839705481694541E-14, 8.709053609235724,
    8.2335672316259332E-13, -9.9617633584461513E-15, 8.7047250784291,
    8.2319800502220348E-13, -9.0838212351977626E-15, 8.7003965476224767,
    8.2303928688181364E-13, -8.2058791119493724E-15, 8.696068016815854,
    8.2288056874142391E-13, -7.3279369887009838E-15, 8.6917394860092312,
    8.2272185060103417E-13, -6.449994865452592E-15, 8.6874109552229974,
    8.2256313246064443E-13, -5.5724833170562121E-15, 8.6830845472654339,
    8.2240449216142213E-13, -4.6945411938078219E-15, 8.6787560164567736,
    8.2224577402103219E-13, -3.8165560130742309E-15, 8.6744272733672823,
    8.2208704809652574E-13, -2.9386569473110415E-15, 8.6700989548496441,
    8.2192833774025271E-13, -2.0608009390330539E-15, 8.665770848614871,
    8.217696351680965E-13, -1.1829018732698649E-15, 208.6614425300952,
    8.2161092481182338E-13, -3.0500280750667595E-16, 8.3323461343610816,
    8.2145221445555035E-13, 1.43877453030842E-15, 8.1973388095849244E-13,
    8.7393512020126352, 1.2693530218508995E-15, 8.198381677078921E-13,
    8.7350226712060124, 1.09993151339338E-15, 8.1994245445729175E-13,
    8.7306941403993878, 9.305100049358595E-16, 8.200467412066914E-13,
    8.7263656095927669, 7.6108849647834011E-16, 8.2015102795609105E-13,
    8.7220370788065313, 5.9175007851049984E-16, 8.2025526355946354E-13,
    8.71771067084897, 4.2232857005297966E-16, 8.203595503088632E-13,
    8.7133821400423468, 2.5290706159545988E-16, 8.2046383705826285E-13,
    8.709053609235724, 8.34855531379397E-17, 8.205681238076625E-13,
    8.7047250784291, -8.5935955319579692E-17, 8.2067241055706215E-13,
    8.7003965476224767, -2.5535746377710027E-16, 8.207766973064617E-13,
    8.696068016815854, -4.2477897223462005E-16, 8.2088098405586136E-13,
    8.6917394860092312, -5.9420048069213984E-16, 8.2098527080526091E-13,
    8.6874109552229974, -7.6353889865997991E-16, 8.210895064086335E-13,
    8.6830845472654339, -9.329604071175E-16, 8.2119379315803315E-13,
    8.6787560164567736, -1.102390224623988E-15, 8.2129808502203543E-13,
    8.6744272733672823, -1.27180342403254E-15, 8.2140236665683235E-13,
    8.6700989548496441, -1.4412083143921236E-15, 8.2150664317702655E-13,
    8.665770848614871, -1.6106215138006761E-15, 8.2161092481182348E-13,
    208.6614425300952, -1.780034713209228E-15, 8.217152064466203E-13,
    8.3323461343610816, 8.4072585967970479, 5.6803301782019976E-17,
    -2.6148934254639176E-14, 8.4030965325047813, 3.6701702672499991E-17,
    -2.4795021656852936E-14, 8.3989344682125129, 1.6600103562980007E-17,
    -2.3441109059066697E-14, 8.3947724039202445, -3.5014955465399778E-18,
    -2.2087196461280457E-14, 8.3906103396279761, -2.3603094656059913E-17,
    -2.0733283863494218E-14, 8.3864503165639075, -4.3694835207899946E-17,
    -1.9380035273848136E-14, 8.3822882522716373, -6.3796434317419931E-17,
    -1.80261226760619E-14, 8.37812618797937, -8.3898033426940014E-17,
    -1.667221007827566E-14, 8.3739641236871023, -1.0399963253646005E-16,
    -1.531829748048942E-14, 8.3698020593948339, -1.2410123164597988E-16,
    -1.3964384882703179E-14, 8.3656399951025673, -1.4420283075549997E-16,
    -1.2610472284916941E-14, 8.3614779308102989, -1.6430442986502E-16,
    -1.12565596871307E-14, 8.3573158665180323, -1.8440602897454003E-16,
    -9.9026470893444619E-15, 8.35315584345396, -2.0449776952637997E-16,
    -8.549398499698379E-15, 8.3489937791616935, -2.245993686359E-16,
    -7.19548590191214E-15, 8.3448315107466069, -2.4470195360118795E-16,
    -5.8415069033118837E-15, 8.3406696505771585, -2.6480256685493997E-16,
    -4.4876607063396593E-15, 8.33650799453053, -2.8490219425292403E-16,
    -3.1338809101814516E-15, 8.3323461343610816, -3.0500280750667595E-16,
    -1.780034713209228E-15, 208.32818427419164, -3.2510342076042796E-16,
    -4.2618851623700396E-16, -1.6980903917954776E-14, 8.4072585967970479,
    8.2453917222222084E-13, -1.6104228529226755E-14, 8.4030965325047813,
    8.2438227540245252E-13, -1.5227553140498734E-14, 8.3989344682125129,
    8.2422537858268419E-13, -1.4350877751770716E-14, 8.3947724039202445,
    8.24068481762916E-13, -1.3474202363042697E-14, 8.3906103396279761,
    8.2391158494314775E-13, -1.2597956927913856E-14, 8.3864503165639075,
    8.2375476507130479E-13, -1.1721281539185835E-14, 8.3822882522716373,
    8.2359786825153657E-13, -1.0844606150457815E-14, 8.37812618797937,
    8.2344097143176834E-13, -9.9679307617297962E-15, 8.3739641236871023,
    8.2328407461199992E-13, -9.0912553730017753E-15, 8.3698020593948339,
    8.231271777922317E-13, -8.214579984273756E-15, 8.3656399951025673,
    8.2297028097246347E-13, -7.3379045955457367E-15, 8.3614779308102989,
    8.2281338415269525E-13, -6.4612292068177166E-15, 8.3573158665180323,
    8.22656487332927E-13, -5.5849837716888766E-15, 8.35315584345396,
    8.2249966746108407E-13, -4.7083083829608558E-15, 8.3489937791616935,
    8.2234277064131575E-13, -3.8315899988729185E-15, 8.3448315107466069,
    8.2218586612675491E-13, -2.9549576055048156E-15, 8.3406696505771585,
    8.220289770017792E-13, -2.0783682074966318E-15, 8.33650799453053,
    8.2187209557159611E-13, -1.2017358141285298E-15, 8.3323461343610816,
    8.217152064466203E-13, -3.2510342076042796E-16, 208.32818427419164,
    8.215583173216446E-13, 3.6493142116161E-15, 8.1954238888548908E-13,
    8.4072585967970479, 3.4348007789651395E-15, 8.1964849695551025E-13,
    8.4030965325047813, 3.22028734631418E-15, 8.1975460502553131E-13,
    8.3989344682125129, 3.00577391366322E-15, 8.1986071309555248E-13,
    8.3947724039202445, 2.7912604810122602E-15, 8.1996682116557365E-13,
    8.3906103396279761, 2.57685225357594E-15, 8.2007287719632558E-13,
    8.3864503165639075, 2.3623388209249795E-15, 8.2017898526634675E-13,
    8.3822882522716373, 2.1478253882740197E-15, 8.2028509333636792E-13,
    8.37812618797937, 1.93331195562306E-15, 8.2039120140638908E-13,
    8.3739641236871023, 1.7187985229721004E-15, 8.2049730947641025E-13,
    8.3698020593948339, 1.5042850903211397E-15, 8.2060341754643132E-13,
    8.3656399951025673, 1.28977165767018E-15, 8.2070952561645249E-13,
    8.3614779308102989, 1.07525822501922E-15, 8.2081563368647355E-13,
    8.3573158665180323, 8.6084999758289989E-16, 8.2092168971722549E-13,
    8.35315584345396, 6.4633656493193984E-16, 8.2102779778724665E-13,
    8.3489937791616935, 4.3181261175951612E-16, 8.2113391106119474E-13,
    8.3448315107466069, 2.1730969963002012E-16, 8.21240013927289E-13,
    8.3406696505771585, 2.8173080219881823E-18, 8.2134611158945621E-13,
    8.33650799453053, -2.1168560410750791E-16, 8.2145221445555045E-13,
    8.3323461343610816, -4.2618851623700386E-16, 8.215583173216446E-13,
    208.32818427419164 };
   double b_a[3600] = { -214.71719164095626, 1.7056254126311141E-14,
    3.8350397641984218E-14, -14.385084459677893, 1.7052288080175281E-14,
    3.6139749542427422E-14, -14.052977361645764, 1.7048322034039424E-14,
    3.3929101442870619E-14, -13.720870430006237, 1.704435598790356E-14,
    3.1718453343313816E-14, -13.388763747903601, 1.70403899417677E-14,
    2.9507805243757019E-14, -13.056660989691615, 1.7036425840725581E-14,
    2.7298241326691417E-14, -12.724555055265856, 1.703245979458972E-14,
    2.508759322713462E-14, -12.392449619783397, 1.702849374845386E-14,
    2.287694512757782E-14, -12.060344766374264, 1.7024527702318E-14,
    2.066629702802102E-14, -11.728240578168471, 1.702056165618214E-14,
    1.8455648928464221E-14, -11.39613713831646, 1.701659561004628E-14,
    1.6245000828907424E-14, -11.06403452996865, 1.701262956391042E-14,
    1.4034352729350621E-14, -10.731932836275458, 1.7008663517774563E-14,
    1.1823704629793823E-14, -10.399835731596751, 1.700469941673244E-14,
    9.614140712728219E-15, -10.067736115828115, 1.7000733370596583E-14,
    7.4034926131714191E-15, -9.7356373050036282, 1.6996767129951346E-14,
    5.1927360953655005E-15, -9.4035404596997552, 1.6992801278324863E-14,
    2.9821964140578193E-15, -9.0714453037624718, 1.6988835621207749E-14,
    7.7176515099926037E-16, -8.7393512020126352, 1.6984869769581263E-14,
    -1.4387745303084196E-15, -8.4072585967970479, 1.6980903917954779E-14,
    -3.6493142116160993E-15, 1.7056254126311141E-14, -214.71719164095626,
    -8.2318069131290087E-13, 1.6155511092337738E-14, -14.385084459677893,
    -8.2298918984796686E-13, 1.5254768058364339E-14, -14.052977361645764,
    -8.2279768838303274E-13, 1.435402502439094E-14, -13.720870430006237,
    -8.2260618691809882E-13, 1.3453281990417541E-14, -13.388763747903601,
    -8.224146854531648E-13, 1.2552980713694741E-14, -13.056660989691615,
    -8.222232779075368E-13, 1.1652237679721339E-14, -12.724555055265856,
    -8.2203177644260278E-13, 1.0751494645747939E-14, -12.392449619783397,
    -8.2184027497766887E-13, 9.85075161177454E-15, -12.060344766374264,
    -8.2164877351273475E-13, 8.95000857780114E-15, -11.728240578168471,
    -8.2145727204780073E-13, 8.04926554382774E-15, -11.39613713831646,
    -8.2126577058286682E-13, 7.14852250985434E-15, -11.06403452996865,
    -8.210742691179328E-13, 6.2477794758809407E-15, -10.731932836275458,
    -8.2088276765299878E-13, 5.3474781991581405E-15, -10.399835731596751,
    -8.2069136010737078E-13, 4.4467351651847404E-15, -10.067736115828115,
    -8.2049985864243676E-13, 3.5459479554862802E-15, -9.7356373050036282,
    -8.2030834778557209E-13, 2.6452490972379395E-15, -9.4035404596997552,
    -8.2011685571256873E-13, 1.74459441471466E-15, -9.0714453037624718,
    -8.19925373031496E-13, 8.4389555646631984E-16, -8.7393512020126352,
    -8.1973388095849255E-13, -5.6803301782020074E-17, -8.4072585967970479,
    -8.1954238888548918E-13, 3.8350397641984218E-14, -8.2318069131290067E-13,
    -214.71719164095626, 3.7708175572864622E-14, -8.23252194737756E-13,
    -14.385084459677893, 3.7065953503745013E-14, -8.2332369816261152E-13,
    -14.052977361645764, 3.6423731434625422E-14, -8.2339520158746689E-13,
    -13.720870430006237, 3.578150936550582E-14, -8.2346670501232227E-13,
    -13.388763747903601, 3.5139602265522621E-14, -8.2353817336928915E-13,
    -13.056660989691615, 3.4497380196403018E-14, -8.2360967679414463E-13,
    -12.724555055265856, 3.3855158127283415E-14, -8.23681180219E-13,
    -12.392449619783397, 3.3212936058163819E-14, -8.2375268364385537E-13,
    -12.060344766374264, 3.2570713989044216E-14, -8.2382418706871085E-13,
    -11.728240578168471, 3.1928491919924619E-14, -8.2389569049356622E-13,
    -11.39613713831646, 3.1286269850805016E-14, -8.2396719391842159E-13,
    -11.06403452996865, 3.064404778168542E-14, -8.24038697343277E-13,
    -10.731932836275458, 3.0002140681702215E-14, -8.2411016570024385E-13,
    -10.399835731596751, 2.9359918612582618E-14, -8.2418166912509923E-13,
    -10.067736115828115, 2.8717665046549382E-14, -8.2425317605674349E-13,
    -9.7356373050036282, 2.8075474474343419E-14, -8.2432467597481008E-13,
    -9.4035404596997552, 2.7433315399051097E-14, -8.2439617238608767E-13,
    -9.0714453037624718, 2.679112482684514E-14, -8.2446767230415426E-13,
    -8.7393512020126352, 2.6148934254639179E-14, -8.2453917222222064E-13,
    -8.4072585967970479, -14.385084459677891, 1.6155511092337741E-14,
    3.7708175572864622E-14, -214.37791332680061, 1.6152811842850578E-14,
    3.5542621609082219E-14, -14.045974116808761, 1.6150112593363422E-14,
    3.3377067645299817E-14, -13.714034990063162, 1.6147413343876259E-14,
    3.121151368151742E-14, -13.382096029710162, 1.61447140943891E-14,
    2.9045959717735018E-14, -13.050160828519424, 1.614201616871238E-14,
    2.6881467820634218E-14, -12.718222449589355, 1.613931691922522E-14,
    2.4715913856851822E-14, -12.386284486450126, 1.613661766973806E-14,
    2.2550359893069419E-14, -12.054347022254197, 1.61339184202509E-14,
    2.038480592928702E-14, -11.722410140131593, 1.613121917076374E-14,
    1.821925196550462E-14, -11.390473923212333, 1.612851992127658E-14,
    1.6053698001722224E-14, -11.05853845464685, 1.612582067178942E-14,
    1.3888144037939821E-14, -10.72660381758557, 1.6123121422302263E-14,
    1.1722590074157422E-14, -10.394673604804272, 1.6120423496625541E-14,
    9.5580981770566187E-15, -10.062740879407485, 1.6117724247138381E-14,
    7.39254421327422E-15, -9.73080888396284, 1.6115024865270176E-14,
    5.2268840428236609E-15, -9.3988787545878836, 1.6112325748164062E-14,
    3.0614362857094196E-15, -9.0669502232951515, 1.610962676343899E-14,
    8.9609473526334018E-16, -8.73502267120601, 1.6106927646332872E-14,
    -1.2693530218508995E-15, -8.40309653250478, 1.6104228529226758E-14,
    -3.43480077896514E-15, 1.7052288080175281E-14, -14.385084459677891,
    -8.2325219473775625E-13, 1.6152811842850578E-14, -214.37791332680061,
    -8.2306251468277238E-13, 1.5253335605525879E-14, -14.045974116808761,
    -8.2287283462778832E-13, 1.435385936820118E-14, -13.714034990063162,
    -8.2268315457280446E-13, 1.345438313087648E-14, -13.382096029710162,
    -8.224934745178206E-13, 1.2555348029519081E-14, -13.050160828519424,
    -8.2230388748885675E-13, 1.1655871792194379E-14, -12.718222449589355,
    -8.2211420743387289E-13, 1.075639555486968E-14, -12.386284486450126,
    -8.21924527378889E-13, 9.85691931754498E-15, -12.054347022254197,
    -8.2173484732390507E-13, 8.95744308022028E-15, -11.722410140131593,
    -8.215451672689212E-13, 8.05796684289558E-15, -11.390473923212333,
    -8.2135548721393724E-13, 7.1584906055708821E-15, -11.05853845464685,
    -8.2116580715895338E-13, 6.2590143682461807E-15, -10.72660381758557,
    -8.2097612710396952E-13, 5.35997926688878E-15, -10.394673604804272,
    -8.2078654007500567E-13, 4.46050302956408E-15, -10.062740879407485,
    -8.2059686002002171E-13, 3.5609826786426505E-15, -9.73080888396284,
    -8.2040717066243584E-13, 2.6615505549146795E-15, -9.3988787545878836,
    -8.20217499910054E-13, 1.7621625447834399E-15, -9.0669502232951515,
    -8.2002783846027415E-13, 8.627304210554699E-16, -8.73502267120601,
    -8.198381677078922E-13, -3.670170267250009E-17, -8.40309653250478,
    -8.1964849695551035E-13, 3.6139749542427415E-14, -8.2298918984796665E-13,
    -14.385084459677891, 3.5542621609082219E-14, -8.2306251468277218E-13,
    -214.37791332680061, 3.4945493675737017E-14, -8.2313583951757771E-13,
    -14.045974116808761, 3.4348365742391821E-14, -8.2320916435238314E-13,
    -13.714034990063162, 3.3751237809046618E-14, -8.2328248918718867E-13,
    -13.382096029710162, 3.3154402729028215E-14, -8.2335577806081982E-13,
    -13.050160828519424, 3.2557274795683019E-14, -8.2342910289562534E-13,
    -12.718222449589355, 3.1960146862337817E-14, -8.2350242773043077E-13,
    -12.386284486450126, 3.1363018928992614E-14, -8.235757525652363E-13,
    -12.054347022254197, 3.0765890995647418E-14, -8.2364907740004183E-13,
    -11.722410140131593, 3.0168763062302216E-14, -8.2372240223484736E-13,
    -11.390473923212333, 2.957163512895702E-14, -8.2379572706965278E-13,
    -11.05853845464685, 2.8974507195611817E-14, -8.2386905190445821E-13,
    -10.72660381758557, 2.8377672115593414E-14, -8.2394234077808936E-13,
    -10.394673604804272, 2.7780544182248218E-14, -8.2401566561289479E-13,
    -10.062740879407485, 2.7183386963570337E-14, -8.2408899404381774E-13,
    -9.73080888396284, 2.658628831555782E-14, -8.2416231528250585E-13,
    -9.3988787545878836, 2.5989218952877975E-14, -8.2423563292507642E-13,
    -9.0669502232951515, 2.539212030486546E-14, -8.2430895416376452E-13,
    -8.73502267120601, 2.4795021656852936E-14, -8.2438227540245242E-13,
    -8.40309653250478, -14.052977361645764, 1.5254768058364339E-14,
    3.7065953503745019E-14, -14.045974116808761, 1.5253335605525879E-14,
    3.4945493675737017E-14, -214.03897087206758, 1.525190315268742E-14,
    3.2825033847729015E-14, -13.70719955021185, 1.525047069984896E-14,
    3.0704574019721019E-14, -13.375428311602359, 1.5249038247010497E-14,
    2.8584114191713023E-14, -13.043660667426746, 1.5247606496699179E-14,
    2.6464694314577019E-14, -12.711889843992372, 1.5246174043860719E-14,
    2.434423448656902E-14, -12.380119353198417, 1.5244741591022259E-14,
    2.2223774658561021E-14, -12.048349278195296, 1.5243309138183797E-14,
    2.0103314830553019E-14, -11.716579702135476, 1.524187668534534E-14,
    1.7982855002545023E-14, -11.384810708148985, 1.5240444232506877E-14,
    1.5862395174537024E-14, -11.053042379365833, 1.5239011779668421E-14,
    1.3741935346529022E-14, -10.721274798936461, 1.5237579326829961E-14,
    1.1621475518521023E-14, -10.389511478052565, 1.5236147576518639E-14,
    9.50205564138502E-15, -10.057745643027635, 1.523471512368018E-14,
    7.38159581337702E-15, -9.72598046296283, 1.5233282600589005E-14,
    5.2610319902818204E-15, -9.3942170494963975, 1.523185021800326E-14,
    3.1406761573610196E-15, -9.0624551428278313, 1.523041790567023E-14,
    1.0204243195274204E-15, -8.7306941403993878, 1.5228985523084482E-14,
    -1.0999315133933793E-15, -8.3989344682125129, 1.5227553140498737E-14,
    -3.2202873463141794E-15, 1.704832203403942E-14, -14.052977361645764,
    -8.2332369816261162E-13, 1.6150112593363418E-14, -14.045974116808761,
    -8.2313583951757781E-13, 1.525190315268742E-14, -214.03897087206758,
    -8.22947980872544E-13, 1.4353693712011418E-14, -13.70719955021185,
    -8.227601222275102E-13, 1.3455484271335421E-14, -13.375428311602359,
    -8.2257226358247639E-13, 1.255771534534342E-14, -13.043660667426746,
    -8.2238449707017681E-13, 1.1659505904667419E-14, -12.711889843992372,
    -8.22196638425143E-13, 1.076129646399142E-14, -12.380119353198417,
    -8.2200877978010919E-13, 9.86308702331542E-15, -12.048349278195296,
    -8.2182092113507538E-13, 8.96487758263942E-15, -11.716579702135476,
    -8.2163306249004158E-13, 8.06666814196342E-15, -11.384810708148985,
    -8.2144520384500777E-13, 7.1684587012874218E-15, -11.053042379365833,
    -8.21257345199974E-13, 6.2702492606114207E-15, -10.721274798936461,
    -8.2106948655494025E-13, 5.37248033461942E-15, -10.389511478052565,
    -8.2088172004264057E-13, 4.47427089394342E-15, -10.057745643027635,
    -8.2069386139760676E-13, 3.57601740179902E-15, -9.72598046296283,
    -8.2050599353929948E-13, 2.6778520125914195E-15, -9.3942170494963975,
    -8.2031814410753915E-13, 1.77973067485222E-15, -9.0624551428278313,
    -8.2013030388905229E-13, 8.8156528564462E-16, -8.7306941403993878,
    -8.1994245445729185E-13, -1.6600103562980056E-17, -8.3989344682125129,
    -8.1975460502553141E-13, 3.3929101442870612E-14, -8.2279768838303264E-13,
    -14.052977361645764, 3.3377067645299817E-14, -8.2287283462778832E-13,
    -14.045974116808761, 3.2825033847729015E-14, -8.229479808725439E-13,
    -214.03897087206758, 3.2273000050158219E-14, -8.2302312711729949E-13,
    -13.70719955021185, 3.1720966252587417E-14, -8.2309827336205507E-13,
    -13.375428311602359, 3.1169203192533816E-14, -8.2317338275235038E-13,
    -13.043660667426746, 3.0617169394963014E-14, -8.23248528997106E-13,
    -12.711889843992372, 3.0065135597392218E-14, -8.2332367524186164E-13,
    -12.380119353198417, 2.9513101799821416E-14, -8.2339882148661723E-13,
    -12.048349278195296, 2.8961068002250614E-14, -8.2347396773137281E-13,
    -11.716579702135476, 2.8409034204679819E-14, -8.2354911397612839E-13,
    -11.384810708148985, 2.7857000407109017E-14, -8.23624260220884E-13,
    -11.053042379365833, 2.7304966609538215E-14, -8.2369940646563956E-13,
    -10.721274798936461, 2.6753203549484614E-14, -8.2377451585593487E-13,
    -10.389511478052565, 2.6201169751913818E-14, -8.2384966210069045E-13,
    -10.057745643027635, 2.5649108880591295E-14, -8.23924812030892E-13,
    -9.72598046296283, 2.5097102156772217E-14, -8.2399995459020162E-13,
    -9.3942170494963975, 2.4545122506704855E-14, -8.2407509346406516E-13,
    -9.0624551428278313, 2.3993115782885777E-14, -8.2415023602337478E-13,
    -8.7306941403993878, 2.3441109059066697E-14, -8.2422537858268419E-13,
    -8.3989344682125129, -13.720870430006235, 1.435402502439094E-14,
    3.6423731434625416E-14, -13.714034990063162, 1.4353859368201177E-14,
    3.4348365742391821E-14, -13.707199550211849, 1.4353693712011421E-14,
    3.2273000050158219E-14, -213.70036411045635, 1.4353528055821658E-14,
    3.0197634357924617E-14, -13.368760593586316, 1.43533623996319E-14,
    2.8122268665691019E-14, -13.037160506419706, 1.4353196824685978E-14,
    2.6047920808519817E-14, -12.705557238474906, 1.4353031168496222E-14,
    2.3972555116286218E-14, -12.373954220026224, 1.4352865512306459E-14,
    2.189718942405262E-14, -12.042351534217955, 1.4352699856116696E-14,
    1.9821823731819021E-14, -11.710749264200526, 1.435253419992694E-14,
    1.7746458039585419E-14, -11.379147493126398, 1.4352368543737178E-14,
    1.5671092347351824E-14, -11.047546304125595, 1.4352202887547418E-14,
    1.3595726655118221E-14, -10.715945780328134, 1.4352037231357662E-14,
    1.1520360962884622E-14, -10.384349351341639, 1.4351871656411741E-14,
    9.44601310571342E-15, -10.05275040668856, 1.4351706000221981E-14,
    7.37064741347982E-15, -9.7211520420035988, 1.4351540335907835E-14,
    5.2951799377399808E-15, -9.3895553444456912, 1.4351374687842459E-14,
    3.2199160290126195E-15, -9.0579600623808982, 1.4351209047901468E-14,
    1.1447539037915006E-15, -8.7263656095927651, 1.4351043399836092E-14,
    -9.3051000493585911E-16, -8.3947724039202445, 1.4350877751770716E-14,
    -3.0057739136632192E-15, 1.704435598790356E-14, -13.720870430006235,
    -8.23395201587467E-13, 1.6147413343876259E-14, -13.714034990063162,
    -8.2320916435238334E-13, 1.525047069984896E-14, -13.707199550211849,
    -8.2302312711729959E-13, 1.4353528055821658E-14, -213.70036411045635,
    -8.2283708988221583E-13, 1.345658541179436E-14, -13.368760593586316,
    -8.2265105264713218E-13, 1.256008266116776E-14, -13.037160506419706,
    -8.2246510665149676E-13, 1.1663140017140459E-14, -12.705557238474906,
    -8.2227906941641311E-13, 1.076619737311316E-14, -12.373954220026224,
    -8.2209303218132945E-13, 9.86925472908586E-15, -12.042351534217955,
    -8.219069949462456E-13, 8.97231208505856E-15, -11.710749264200526,
    -8.2172095771116195E-13, 8.07536944103126E-15, -11.379147493126398,
    -8.215349204760783E-13, 7.1784267970039615E-15, -11.047546304125595,
    -8.2134888324099454E-13, 6.2814841529766606E-15, -10.715945780328134,
    -8.2116284600591089E-13, 5.38498140235006E-15, -10.384349351341639,
    -8.2097690001027547E-13, 4.48803875832276E-15, -10.05275040668856,
    -8.2079086277519171E-13, 3.59105212495539E-15, -9.7211520420035988,
    -8.2060481641616323E-13, 2.6941534702681594E-15, -9.3895553444456912,
    -8.2041878830502431E-13, 1.797298804921E-15, -9.0579600623808982,
    -8.2023276931783032E-13, 9.0040015023377E-16, -8.7263656095927651,
    -8.200467412066914E-13, 3.5014955465399778E-18, -8.3947724039202445,
    -8.1986071309555258E-13, 3.1718453343313809E-14, -8.2260618691809872E-13,
    -13.720870430006235, 3.1211513681517414E-14, -8.2268315457280436E-13,
    -13.714034990063162, 3.0704574019721013E-14, -8.227601222275101E-13,
    -13.707199550211849, 3.0197634357924617E-14, -8.2283708988221583E-13,
    -213.70036411045635, 2.9690694696128216E-14, -8.2291405753692147E-13,
    -13.368760593586316, 2.9184003656039417E-14, -8.22990987443881E-13,
    -13.037160506419706, 2.8677063994243015E-14, -8.2306795509858668E-13,
    -12.705557238474906, 2.8170124332446614E-14, -8.2314492275329242E-13,
    -12.373954220026224, 2.7663184670650212E-14, -8.2322189040799815E-13,
    -12.042351534217955, 2.7156245008853817E-14, -8.2329885806270379E-13,
    -11.710749264200526, 2.6649305347057415E-14, -8.2337582571740953E-13,
    -11.379147493126398, 2.6142365685261014E-14, -8.2345279337211527E-13,
    -11.047546304125595, 2.5635426023464616E-14, -8.235297610268208E-13,
    -10.715945780328134, 2.5128734983375813E-14, -8.2360669093378037E-13,
    -10.384349351341639, 2.4621795321579418E-14, -8.23683658588486E-13,
    -10.05275040668856, 2.4114830797612256E-14, -8.2376063001796635E-13,
    -9.7211520420035988, 2.3607915997986618E-14, -8.2383759389789749E-13,
    -9.3895553444456912, 2.3101026060531735E-14, -8.2391455400305391E-13,
    -9.0579600623808982, 2.2594111260906098E-14, -8.2399151788298494E-13,
    -8.7263656095927651, 2.2087196461280457E-14, -8.2406848176291587E-13,
    -8.3947724039202445, -13.388763747903601, 1.3453281990417542E-14,
    3.578150936550582E-14, -13.382096029710164, 1.3454383130876481E-14,
    3.3751237809046618E-14, -13.375428311602359, 1.3455484271335422E-14,
    3.1720966252587417E-14, -13.368760593586318, 1.345658541179436E-14,
    2.9690694696128216E-14, -213.36209287566609, 1.3457686552253301E-14,
    2.7660423139669021E-14, -13.030660345504426, 1.3458787152672781E-14,
    2.5631147302462618E-14, -12.699224633043077, 1.3459888293131721E-14,
    2.3600875746003419E-14, -12.367789086933549, 1.346098943359066E-14,
    2.1570604189544221E-14, -12.036353790320133, 1.3462090574049601E-14,
    1.954033263308502E-14, -11.704918826347136, 1.346319171450854E-14,
    1.7510061076625822E-14, -11.373484278164979, 1.346429285496748E-14,
    1.5479789520166624E-14, -11.042050228926119, 1.3465393995426422E-14,
    1.3449517963707421E-14, -10.710616761760587, 1.3466495135885361E-14,
    1.1419246407248223E-14, -10.379187224671494, 1.3467595736304839E-14,
    9.38997057004182E-15, -10.047755170390262, 1.3468696876763781E-14,
    7.35969901358262E-15, -9.7163236210851434, 1.3469798071226666E-14,
    5.32932788519814E-15, -9.3848936394357629, 1.3470899157681661E-14,
    3.2991559006642194E-15, -9.053464981974745, 1.347200019013271E-14,
    1.26908348805558E-15, -8.7220370788065313, 1.3473101276587702E-14,
    -7.6108849647833972E-16, -8.3906103396279779, 1.3474202363042697E-14,
    -2.79126048101226E-15, 1.70403899417677E-14, -13.388763747903601,
    -8.2346670501232247E-13, 1.61447140943891E-14, -13.382096029710164,
    -8.2328248918718887E-13, 1.52490382470105E-14, -13.375428311602359,
    -8.2309827336205517E-13, 1.43533623996319E-14, -13.368760593586318,
    -8.2291405753692157E-13, 1.3457686552253301E-14, -213.36209287566609,
    -8.22729841711788E-13, 1.25624499769921E-14, -13.030660345504426,
    -8.2254571623281681E-13, 1.1666774129613498E-14, -12.699224633043077,
    -8.2236150040768321E-13, 1.07710982822349E-14, -12.367789086933549,
    -8.2217728458254962E-13, 9.8754224348563E-15, -12.036353790320133,
    -8.2199306875741592E-13, 8.9797465874776991E-15, -11.704918826347136,
    -8.2180885293228232E-13, 8.0840707400991E-15, -11.373484278164979,
    -8.2162463710714882E-13, 7.1883948927205012E-15, -11.042050228926119,
    -8.2144042128201522E-13, 6.2927190453419006E-15, -10.710616761760587,
    -8.2125620545688163E-13, 5.3974824700807E-15, -10.379187224671494,
    -8.2107207997791036E-13, 4.5018066227021E-15, -10.047755170390262,
    -8.2088786415277677E-13, 3.6060868481117604E-15, -9.7163236210851434,
    -8.2070363929302687E-13, 2.7104549279448994E-15, -9.3848936394357629,
    -8.2051943250250957E-13, 1.81486693498978E-15, -9.053464981974745,
    -8.2033523474660846E-13, 9.1923501482291989E-16, -8.7220370788065313,
    -8.2015102795609105E-13, 2.3603094656059913E-17, -8.3906103396279779,
    -8.1996682116557375E-13, 2.9507805243757019E-14, -8.224146854531647E-13,
    -13.388763747903601, 2.9045959717735024E-14, -8.224934745178205E-13,
    -13.382096029710164, 2.8584114191713017E-14, -8.2257226358247629E-13,
    -13.375428311602359, 2.8122268665691022E-14, -8.2265105264713208E-13,
    -13.368760593586318, 2.7660423139669021E-14, -8.2272984171178787E-13,
    -213.36209287566609, 2.7198804119545018E-14, -8.228085921354116E-13,
    -13.030660345504426, 2.6736958593523016E-14, -8.228873812000674E-13,
    -12.699224633043077, 2.6275113067501022E-14, -8.2296617026472319E-13,
    -12.367789086933549, 2.5813267541479021E-14, -8.23044959329379E-13,
    -12.036353790320133, 2.535142201545702E-14, -8.2312374839403477E-13,
    -11.704918826347136, 2.4889576489435018E-14, -8.2320253745869067E-13,
    -11.373484278164979, 2.4427730963413017E-14, -8.2328132652334646E-13,
    -11.042050228926119, 2.396588543739102E-14, -8.2336011558800215E-13,
    -10.710616761760587, 2.3504266417267016E-14, -8.2343886601162578E-13,
    -10.379187224671494, 2.3042420891245021E-14, -8.2351765507628167E-13,
    -10.047755170390262, 2.2580552714633217E-14, -8.2359644800504061E-13,
    -9.7163236210851434, 2.2118729839201019E-14, -8.2367523320559326E-13,
    -9.3848936394357629, 2.1656929614358619E-14, -8.2375401454204256E-13,
    -9.053464981974745, 2.1195106738926421E-14, -8.238327997425952E-13,
    -8.7220370788065313, 2.0733283863494218E-14, -8.2391158494314765E-13,
    -8.3906103396279779, -13.056660989691615, 1.2552980713694741E-14,
    3.5139602265522615E-14, -13.050160828519424, 1.2555348029519081E-14,
    3.3154402729028222E-14, -13.043660667426746, 1.2557715345343419E-14,
    3.1169203192533816E-14, -13.037160506419706, 1.256008266116776E-14,
    2.9184003656039417E-14, -13.030660345504424, 1.25624499769921E-14,
    2.7198804119545018E-14, -213.024163368014, 1.256481613179838E-14,
    2.5214578197340217E-14, -12.69289512944794, 1.256718344762272E-14,
    2.322937866084582E-14, -12.361626974087356, 1.2569550763447059E-14,
    2.1244179124351418E-14, -12.030358985078585, 1.25719180792714E-14,
    1.9258979587857019E-14, -11.699091245565935, 1.2574285395095739E-14,
    1.727378005136262E-14, -11.3678238386937, 1.2576652710920079E-14,
    1.5288580514868223E-14, -11.036556847612303, 1.2579020026744421E-14,
    1.3303380978373821E-14, -10.705290355474208, 1.2581387342568761E-14,
    1.1318181441879422E-14, -10.37402762873846, 1.2583753497375038E-14,
    9.33395551967462E-15, -10.042762383244991, 1.2586120813199381E-14,
    7.34875598318022E-15, -9.71149756773162, 1.2588488245125525E-14,
    5.363459085256861E-15, -9.3802342204106779, 1.259085544484806E-14,
    3.3783569101914196E-15, -9.04897210597735, 1.2593222528468788E-14,
    1.3933520965549404E-15, -8.71771067084897, 1.2595589728191322E-14,
    -5.9175007851049945E-16, -8.3864503165639057, 1.2597956927913856E-14,
    -2.5768522535759397E-15, 1.7036425840725581E-14, -13.056660989691615,
    -8.2353817336928925E-13, 1.614201616871238E-14, -13.050160828519424,
    -8.2335577806081982E-13, 1.5247606496699179E-14, -13.043660667426746,
    -8.2317338275235038E-13, 1.4353196824685978E-14, -13.037160506419706,
    -8.2299098744388094E-13, 1.345878715267278E-14, -13.030660345504424,
    -8.228085921354116E-13, 1.256481613179838E-14, -213.024163368014,
    -8.2262628628025683E-13, 1.1670406459785179E-14, -12.69289512944794,
    -8.224438909717874E-13, 1.077599678777198E-14, -12.361626974087356,
    -8.2226149566331806E-13, 9.88158711575878E-15, -12.030358985078585,
    -8.2207910035484852E-13, 8.98717744374558E-15, -11.699091245565935,
    -8.2189670504637918E-13, 8.09276777173238E-15, -11.3678238386937,
    -8.2171430973790974E-13, 7.19835809971918E-15, -11.036556847612303,
    -8.2153191442944041E-13, 6.30394842770598E-15, -10.705290355474208,
    -8.21349519120971E-13, 5.4099774068315805E-15, -10.37402762873846,
    -8.211672132658162E-13, 4.51556773481838E-15, -10.042762383244991,
    -8.2098481795734676E-13, 3.6211141976913E-15, -9.71149756773162,
    -8.2080241370354585E-13, 2.7267483907919795E-15, -9.3802342204106779,
    -8.20620027340408E-13, 1.83242644900654E-15, -9.04897210597735,
    -8.204376499226015E-13, 9.3806064210721983E-16, -8.71771067084897,
    -8.2025526355946354E-13, 4.3694835207899946E-17, -8.3864503165639057,
    -8.2007287719632558E-13, 2.7298241326691417E-14, -8.222232779075367E-13,
    -13.056660989691615, 2.6881467820634218E-14, -8.2230388748885665E-13,
    -13.050160828519424, 2.6464694314577013E-14, -8.2238449707017671E-13,
    -13.043660667426746, 2.604792080851982E-14, -8.2246510665149666E-13,
    -13.037160506419706, 2.5631147302462621E-14, -8.2254571623281671E-13,
    -13.030660345504424, 2.5214578197340217E-14, -8.2262628628025683E-13,
    -213.024163368014, 2.4797804691283017E-14, -8.2270689586157679E-13,
    -12.69289512944794, 2.4381031185225818E-14, -8.2278750544289684E-13,
    -12.361626974087356, 2.3964257679168616E-14, -8.2286811502421679E-13,
    -12.030358985078585, 2.3547484173111417E-14, -8.2294872460553685E-13,
    -11.699091245565935, 2.3130710667054218E-14, -8.230293341868568E-13,
    -11.3678238386937, 2.2713937160997019E-14, -8.2310994376817685E-13,
    -11.036556847612303, 2.2297163654939817E-14, -8.231905533494967E-13,
    -10.705290355474208, 2.1880594549817416E-14, -8.2327112339693682E-13,
    -10.37402762873846, 2.146382104376022E-14, -8.2335173297825688E-13,
    -10.042762383244991, 2.1047027097609538E-14, -8.2343234651296479E-13,
    -9.71149756773162, 2.0630274031645819E-14, -8.2351295214089688E-13,
    -9.3802342204106779, 2.0213541405775575E-14, -8.2359355381544081E-13,
    -9.04897210597735, 1.9796788339811859E-14, -8.236741594433729E-13,
    -8.71771067084897, 1.9380035273848139E-14, -8.2375476507130469E-13,
    -8.3864503165639057, -12.724555055265856, 1.165223767972134E-14,
    3.4497380196403018E-14, -12.718222449589355, 1.165587179219438E-14,
    3.2557274795683019E-14, -12.711889843992372, 1.165950590466742E-14,
    3.0617169394963014E-14, -12.705557238474906, 1.166314001714046E-14,
    2.8677063994243015E-14, -12.699224633043077, 1.16667741296135E-14,
    2.673695859352302E-14, -12.69289512944794, 1.1670406459785181E-14,
    2.4797804691283017E-14, -212.68656252361242, 1.1674040572258219E-14,
    2.2857699290563019E-14, -12.355461840621572, 1.1677674684731259E-14,
    2.091759388984302E-14, -12.024361240836198, 1.1681308797204299E-14,
    1.8977488489123021E-14, -11.693260807402639, 1.1684942909677339E-14,
    1.7037383088403022E-14, -11.362160623465199, 1.1688577022150379E-14,
    1.5097277687683023E-14, -11.031060772168173, 1.169221113462342E-14,
    1.3157172286963022E-14, -10.699961336661989, 1.169584524709646E-14,
    1.1217066886243023E-14, -10.368865501844034, 1.1699477577268138E-14,
    9.277912984003019E-15, -10.037767146763196, 1.1703111689741181E-14,
    7.3378075832830187E-15, -9.7066691466704427, 1.1706745980444355E-14,
    5.3976070327150205E-15, -9.3755725152988045, 1.1710379914687259E-14,
    3.4575967818430195E-15, -9.04447702551003, 1.1714013670700029E-14,
    1.5176816808190202E-15, -8.7133821400423468, 1.171764760494293E-14,
    -4.2232857005297966E-16, -8.3822882522716373, 1.1721281539185835E-14,
    -2.3623388209249795E-15, 1.703245979458972E-14, -12.724555055265856,
    -8.2360967679414463E-13, 1.613931691922522E-14, -12.718222449589355,
    -8.2342910289562534E-13, 1.5246174043860719E-14, -12.711889843992372,
    -8.23248528997106E-13, 1.4353031168496218E-14, -12.705557238474906,
    -8.2306795509858668E-13, 1.3459888293131721E-14, -12.699224633043077,
    -8.228873812000674E-13, 1.256718344762272E-14, -12.69289512944794,
    -8.2270689586157679E-13, 1.1674040572258219E-14, -212.68656252361242,
    -8.225263219630575E-13, 1.078089769689372E-14, -12.355461840621572,
    -8.2234574806453822E-13, 9.88775482152922E-15, -12.024361240836198,
    -8.2216517416601884E-13, 8.99461194616472E-15, -11.693260807402639,
    -8.2198460026749955E-13, 8.10146907080022E-15, -11.362160623465199,
    -8.2180402636898027E-13, 7.2083261954357218E-15, -11.031060772168173,
    -8.21623452470461E-13, 6.31518332007122E-15, -10.699961336661989,
    -8.2144287857194171E-13, 5.42247847456222E-15, -10.368865501844034,
    -8.212623932334511E-13, 4.52933559919772E-15, -10.037767146763196,
    -8.2108181933493171E-13, 3.6361489208476704E-15, -9.7066691466704427,
    -8.2090123658040959E-13, 2.7430498484687195E-15, -9.3755725152988045,
    -8.2072067153789315E-13, 1.84999457907532E-15, -9.04447702551003,
    -8.2054011535137964E-13, 9.5689550669636989E-16, -8.7133821400423468,
    -8.203595503088632E-13, 6.3796434317419931E-17, -8.3822882522716373,
    -8.2017898526634675E-13, 2.5087593227134614E-14, -8.2203177644260268E-13,
    -12.724555055265856, 2.4715913856851815E-14, -8.2211420743387279E-13,
    -12.718222449589355, 2.4344234486569017E-14, -8.221966384251429E-13,
    -12.711889843992372, 2.3972555116286218E-14, -8.22279069416413E-13,
    -12.705557238474906, 2.3600875746003419E-14, -8.2236150040768311E-13,
    -12.699224633043077, 2.3229378660845817E-14, -8.224438909717874E-13,
    -12.69289512944794, 2.2857699290563019E-14, -8.225263219630575E-13,
    -212.68656252361242, 2.2486019920280217E-14, -8.2260875295432761E-13,
    -12.355461840621572, 2.2114340549997415E-14, -8.2269118394559772E-13,
    -12.024361240836198, 2.1742661179714617E-14, -8.2277361493686783E-13,
    -11.693260807402639, 2.1370981809431818E-14, -8.2285604592813793E-13,
    -11.362160623465199, 2.0999302439149016E-14, -8.22938476919408E-13,
    -11.031060772168173, 2.0627623068866217E-14, -8.23020907910678E-13,
    -10.699961336661989, 2.0256125983708615E-14, -8.2310329847478233E-13,
    -10.368865501844034, 1.988444661342582E-14, -8.2318572946605244E-13,
    -10.037767146763196, 1.9512749014630496E-14, -8.23268164500039E-13,
    -9.7066691466704427, 1.914108787286022E-14, -8.2335059144859265E-13,
    -9.3755725152988045, 1.8769444959602456E-14, -8.2343301435442956E-13,
    -9.04447702551003, 1.8397783817832179E-14, -8.2351544130298306E-13,
    -8.7133821400423468, 1.8026122676061897E-14, -8.2359786825153647E-13,
    -8.3822882522716373, -12.392449619783395, 1.0751494645747939E-14,
    3.3855158127283415E-14, -12.386284486450126, 1.075639555486968E-14,
    3.1960146862337817E-14, -12.380119353198415, 1.0761296463991418E-14,
    3.0065135597392218E-14, -12.373954220026224, 1.0766197373113158E-14,
    2.8170124332446617E-14, -12.367789086933547, 1.0771098282234899E-14,
    2.6275113067501019E-14, -12.361626974087356, 1.077599678777198E-14,
    2.4381031185225815E-14, -12.355461840621572, 1.0780897696893718E-14,
    2.248601992028022E-14, -212.34929670721084, 1.0785798606015459E-14,
    2.0591008655334619E-14, -12.018363496644785, 1.0790699515137199E-14,
    1.869599739038902E-14, -11.6874303692842, 1.0795600424258937E-14,
    1.6800986125443422E-14, -11.356497408275434, 1.0800501333380678E-14,
    1.490597486049782E-14, -11.025564696762782, 1.0805402242502421E-14,
    1.301096359555222E-14, -10.694632317890548, 1.081030315162416E-14,
    1.1115952330606622E-14, -10.363703374969997, 1.0815201657161237E-14,
    9.2218704483314188E-15, -10.0327719102814, 1.0820102566282979E-14,
    7.3268591833858189E-15, -9.7018407256092676, 1.0825003715763184E-14,
    5.43175498017318E-15, -9.3709108101869329, 1.082990438452646E-14,
    3.5368366534946198E-15, -9.03998194504271, 1.0834804812931267E-14,
    1.6420112650831002E-15, -8.709053609235724, 1.083970548169454E-14,
    -2.5290706159545948E-16, -8.37812618797937, 1.0844606150457814E-14,
    -2.1478253882740194E-15, 1.702849374845386E-14, -12.392449619783395,
    -8.23681180219E-13, 1.613661766973806E-14, -12.386284486450126,
    -8.2350242773043087E-13, 1.5244741591022259E-14, -12.380119353198415,
    -8.2332367524186154E-13, 1.4352865512306459E-14, -12.373954220026224,
    -8.2314492275329242E-13, 1.346098943359066E-14, -12.367789086933547,
    -8.2296617026472319E-13, 1.2569550763447061E-14, -12.361626974087356,
    -8.2278750544289684E-13, 1.1677674684731259E-14, -12.355461840621572,
    -8.2260875295432761E-13, 1.078579860601546E-14, -212.34929670721084,
    -8.2243000046575838E-13, 9.89392252729966E-15, -12.018363496644785,
    -8.2225124797718915E-13, 9.00204644858386E-15, -11.6874303692842,
    -8.2207249548861993E-13, 8.11017036986806E-15, -11.356497408275434,
    -8.218937430000508E-13, 7.2182942911522615E-15, -11.025564696762782,
    -8.2171499051148157E-13, 6.32641821243646E-15, -10.694632317890548,
    -8.2153623802291244E-13, 5.4349795422928604E-15, -10.363703374969997,
    -8.21357573201086E-13, 4.54310346357706E-15, -10.0327719102814,
    -8.2117882071251677E-13, 3.65118364400404E-15, -9.7018407256092676,
    -8.2100005945727324E-13, 2.7593513061454595E-15, -9.3709108101869329,
    -8.2082131573537831E-13, 1.8675627091441E-15, -9.03998194504271,
    -8.2064258078015778E-13, 9.7573037128552E-16, -8.709053609235724,
    -8.2046383705826275E-13, 8.3898033426939964E-17, -8.37812618797937,
    -8.2028509333636792E-13, 2.2876945127577817E-14, -8.2184027497766867E-13,
    -12.392449619783395, 2.2550359893069419E-14, -8.2192452737888883E-13,
    -12.386284486450126, 2.2223774658561015E-14, -8.2200877978010909E-13,
    -12.380119353198415, 2.1897189424052616E-14, -8.2209303218132925E-13,
    -12.373954220026224, 2.1570604189544218E-14, -8.2217728458254952E-13,
    -12.367789086933547, 2.1244179124351418E-14, -8.22261495663318E-13,
    -12.361626974087356, 2.091759388984302E-14, -8.2234574806453822E-13,
    -12.355461840621572, 2.0591008655334615E-14, -8.2243000046575838E-13,
    -212.34929670721084, 2.0264423420826214E-14, -8.2251425286697865E-13,
    -12.018363496644785, 1.9937838186317816E-14, -8.2259850526819881E-13,
    -11.6874303692842, 1.9611252951809415E-14, -8.22682757669419E-13,
    -11.356497408275434, 1.9284667717301016E-14, -8.2276701007063923E-13,
    -11.025564696762782, 1.8958082482792615E-14, -8.2285126247185929E-13,
    -10.694632317890548, 1.8631657417599815E-14, -8.2293547355262784E-13,
    -10.363703374969997, 1.8305072183091417E-14, -8.23019725953848E-13,
    -10.0327719102814, 1.7978470931651457E-14, -8.231039824871134E-13,
    -9.7018407256092676, 1.7651901714074617E-14, -8.2318823075628842E-13,
    -9.3709108101869329, 1.7325348513429336E-14, -8.232724748934183E-13,
    -9.03998194504271, 1.6998779295852497E-14, -8.2335672316259332E-13,
    -8.709053609235724, 1.6672210078275657E-14, -8.2344097143176814E-13,
    -8.37812618797937, -12.060344766374264, 9.85075161177454E-15,
    3.3212936058163819E-14, -12.054347022254197, 9.85691931754498E-15,
    3.1363018928992621E-14, -12.048349278195296, 9.86308702331542E-15,
    2.9513101799821416E-14, -12.042351534217957, 9.86925472908586E-15,
    2.7663184670650219E-14, -12.036353790320133, 9.8754224348563E-15,
    2.5813267541479017E-14, -12.030358985078587, 9.88158711575878E-15,
    2.3964257679168616E-14, -12.024361240836198, 9.8877548215292192E-15,
    2.2114340549997421E-14, -12.018363496644787, 9.89392252729966E-15,
    2.026442342082622E-14, -212.01236575250843, 9.9000902330701E-15,
    1.8414506291655019E-14, -11.681599931216738, 9.9062579388405391E-15,
    1.6564589162483821E-14, -11.350834193130524, 9.91242564461098E-15,
    1.4714672033312624E-14, -11.020068621396126, 9.91859335038142E-15,
    1.2864754904141421E-14, -10.689303299157846, 9.924761056151859E-15,
    1.1014837774970221E-14, -10.358541248136738, 9.9309257370543388E-15,
    9.1658279126598186E-15, -10.027776673819995, 9.9370934428247808E-15,
    7.3159107834886191E-15, -9.69701230454809, 9.9432614510820153E-15,
    5.4659029276313404E-15, -9.36624910507506, 9.94942885436566E-15,
    3.6160765251462193E-15, -9.03548686457539, 9.9555959551625079E-15,
    1.7663408493471802E-15, -8.7047250784291, 9.9617633584461513E-15,
    -8.34855531379397E-17, -8.3739641236871023, 9.9679307617297946E-15,
    -1.9333119556230596E-15, 1.7024527702318E-14, -12.060344766374264,
    -8.2375268364385547E-13, 1.61339184202509E-14, -12.054347022254197,
    -8.235757525652364E-13, 1.52433091381838E-14, -12.048349278195296,
    -8.2339882148661713E-13, 1.43526998561167E-14, -12.042351534217957,
    -8.2322189040799805E-13, 1.3462090574049601E-14, -12.036353790320133,
    -8.23044959329379E-13, 1.25719180792714E-14, -12.030358985078587,
    -8.2286811502421679E-13, 1.1681308797204299E-14, -12.024361240836198,
    -8.2269118394559772E-13, 1.07906995151372E-14, -12.018363496644787,
    -8.2251425286697865E-13, 9.9000902330701E-15, -212.01236575250843,
    -8.2233732178835947E-13, 9.009480951003E-15, -11.681599931216738,
    -8.221603907097404E-13, 8.1188716689359E-15, -11.350834193130524,
    -8.2198345963112133E-13, 7.2282623868688012E-15, -11.020068621396126,
    -8.2180652855250215E-13, 6.3376531048017E-15, -10.689303299157846,
    -8.2162959747388308E-13, 5.4474806100235E-15, -10.358541248136738,
    -8.2145275316872089E-13, 4.5568713279564E-15, -10.027776673819995,
    -8.2127582209010172E-13, 3.66621836716041E-15, -9.69701230454809,
    -8.2109888233413688E-13, 2.7756527638221995E-15, -9.36624910507506,
    -8.2092195993286357E-13, 1.88513083921288E-15, -9.03548686457539,
    -8.2074504620893591E-13, 9.9456523587466982E-16, -8.7047250784291,
    -8.205681238076624E-13, 1.0399963253645995E-16, -8.3739641236871023,
    -8.2039120140638908E-13, 2.0666297028021014E-14, -8.2164877351273465E-13,
    -12.060344766374264, 2.0384805929287023E-14, -8.21734847323905E-13,
    -12.054347022254197, 2.0103314830553016E-14, -8.2182092113507528E-13,
    -12.048349278195296, 1.9821823731819018E-14, -8.219069949462456E-13,
    -12.042351534217957, 1.9540332633085017E-14, -8.2199306875741592E-13,
    -12.036353790320133, 1.9258979587857019E-14, -8.2207910035484862E-13,
    -12.030358985078587, 1.8977488489123021E-14, -8.2216517416601894E-13,
    -12.024361240836198, 1.8695997390389017E-14, -8.2225124797718915E-13,
    -12.018363496644787, 1.8414506291655016E-14, -8.2233732178835947E-13,
    -212.01236575250843, 1.8133015192921018E-14, -8.2242339559952979E-13,
    -11.681599931216738, 1.7851524094187018E-14, -8.2250946941070011E-13,
    -11.350834193130524, 1.7570032995453017E-14, -8.2259554322187042E-13,
    -11.020068621396126, 1.7288541896719016E-14, -8.2268161703304064E-13,
    -10.689303299157846, 1.7007188851491015E-14, -8.2276764863047334E-13,
    -10.358541248136738, 1.672569775275702E-14, -8.2285372244164366E-13,
    -10.027776673819995, 1.6444192848672418E-14, -8.2293980047418766E-13,
    -9.69701230454809, 1.6162715555289018E-14, -8.2302587006398419E-13,
    -9.36624910507506, 1.5881252067256217E-14, -8.2311193543240695E-13,
    -9.03548686457539, 1.559977477387282E-14, -8.2319800502220359E-13,
    -8.7047250784291, 1.5318297480489417E-14, -8.2328407461199992E-13,
    -8.3739641236871023, -11.728240578168471, 8.9500085778011409E-15,
    3.2570713989044216E-14, -11.722410140131595, 8.9574430802202816E-15,
    3.0765890995647418E-14, -11.716579702135476, 8.9648775826394208E-15,
    2.8961068002250621E-14, -11.710749264200528, 8.97231208505856E-15,
    2.7156245008853817E-14, -11.704918826347136, 8.9797465874777E-15,
    2.535142201545702E-14, -11.699091245565937, 8.9871774437455809E-15,
    2.3547484173111417E-14, -11.693260807402639, 8.99461194616472E-15,
    2.174266117971462E-14, -11.687430369284202, 9.00204644858386E-15,
    1.9937838186317819E-14, -11.681599931216738, 9.009480951003E-15,
    1.8133015192921022E-14, -211.67576949320431, 9.0169154534221392E-15,
    1.6328192199524221E-14, -11.345170978036592, 9.02434995584128E-15,
    1.4523369206127424E-14, -11.014572546074326, 9.03178445826042E-15,
    1.271854621273062E-14, -10.683974280463881, 9.03921896067956E-15,
    1.0913723219333822E-14, -10.353379121342218, 9.0466498169474385E-15,
    9.10978537698822E-15, -10.022781437399365, 9.05408431936658E-15,
    7.30496238359142E-15, -9.6921838835073029, 9.0615191864008465E-15,
    5.5000508750895E-15, -9.3615873999631862, 9.06895332420486E-15,
    3.69531639679782E-15, -9.03099178410807, 9.0763870973937484E-15,
    1.89067043361126E-15, -8.7003965476224767, 9.0838212351977611E-15,
    8.5935955319580086E-17, -8.3698020593948357, 9.0912553730017753E-15,
    -1.7187985229720998E-15, 1.702056165618214E-14, -11.728240578168471,
    -8.2382418706871085E-13, 1.613121917076374E-14, -11.722410140131595,
    -8.2364907740004183E-13, 1.524187668534534E-14, -11.716579702135476,
    -8.2347396773137281E-13, 1.435253419992694E-14, -11.710749264200528,
    -8.2329885806270379E-13, 1.346319171450854E-14, -11.704918826347136,
    -8.2312374839403477E-13, 1.257428539509574E-14, -11.699091245565937,
    -8.2294872460553685E-13, 1.1684942909677339E-14, -11.693260807402639,
    -8.2277361493686783E-13, 1.0795600424258941E-14, -11.687430369284202,
    -8.2259850526819881E-13, 9.90625793884054E-15, -11.681599931216738,
    -8.2242339559952979E-13, 9.0169154534221392E-15, -211.67576949320431,
    -8.2224828593086077E-13, 8.12757296800374E-15, -11.345170978036592,
    -8.2207317626219175E-13, 7.2382304825853409E-15, -11.014572546074326,
    -8.2189806659352283E-13, 6.34888799716694E-15, -10.683974280463881,
    -8.2172295692485381E-13, 5.45998167775414E-15, -10.353379121342218,
    -8.2154793313635579E-13, 4.5706391923357404E-15, -10.022781437399365,
    -8.2137282346768677E-13, 3.68125309031678E-15, -9.6921838835073029,
    -8.2119770521100063E-13, 2.7919542214989394E-15, -9.3615873999631862,
    -8.2102260413034873E-13, 1.90269896928166E-15, -9.03099178410807,
    -8.2084751163771405E-13, 1.0134001004638199E-15, -8.7003965476224767,
    -8.2067241055706205E-13, 1.2410123164597993E-16, -8.3698020593948357,
    -8.2049730947641015E-13, 1.8455648928464217E-14, -8.2145727204780063E-13,
    -11.728240578168471, 1.821925196550462E-14, -8.215451672689211E-13,
    -11.722410140131595, 1.7982855002545017E-14, -8.2163306249004148E-13,
    -11.716579702135476, 1.7746458039585419E-14, -8.2172095771116185E-13,
    -11.710749264200528, 1.7510061076625819E-14, -8.2180885293228232E-13,
    -11.704918826347136, 1.727378005136262E-14, -8.2189670504637918E-13,
    -11.699091245565937, 1.7037383088403022E-14, -8.2198460026749955E-13,
    -11.693260807402639, 1.6800986125443419E-14, -8.2207249548862E-13,
    -11.687430369284202, 1.6564589162483818E-14, -8.221603907097404E-13,
    -11.681599931216738, 1.6328192199524218E-14, -8.2224828593086077E-13,
    -211.67576949320431, 1.6091795236564617E-14, -8.2233618115198124E-13,
    -11.345170978036592, 1.585539827360502E-14, -8.2242407637310161E-13,
    -11.014572546074326, 1.561900131064542E-14, -8.2251197159422188E-13,
    -10.683974280463881, 1.5382720285382217E-14, -8.2259982370831885E-13,
    -10.353379121342218, 1.514632332242262E-14, -8.2268771892943922E-13,
    -10.022781437399365, 1.4909914765693379E-14, -8.2277561846126191E-13,
    -9.6921838835073029, 1.4673529396503419E-14, -8.2286350937168007E-13,
    -9.3615873999631862, 1.4437155621083097E-14, -8.2295139597139569E-13,
    -9.03099178410807, 1.4200770251893141E-14, -8.2303928688181375E-13,
    -8.7003965476224767, 1.3964384882703178E-14, -8.231271777922317E-13,
    -8.3698020593948357, -11.39613713831646, 8.04926554382774E-15,
    3.1928491919924619E-14, -11.390473923212333, 8.057966842895581E-15,
    3.0168763062302216E-14, -11.384810708148985, 8.0666681419634189E-15,
    2.8409034204679819E-14, -11.379147493126398, 8.07536944103126E-15,
    2.6649305347057415E-14, -11.373484278164979, 8.0840707400991E-15,
    2.4889576489435018E-14, -11.367823838693702, 8.09276777173238E-15,
    2.3130710667054215E-14, -11.362160623465199, 8.10146907080022E-15,
    2.1370981809431818E-14, -11.356497408275434, 8.1101703698680588E-15,
    1.9611252951809421E-14, -11.350834193130524, 8.1188716689359E-15,
    1.7851524094187021E-14, -11.345170978036592, 8.12757296800374E-15,
    1.6091795236564621E-14, -211.3395077629977, 8.1362742670715787E-15,
    1.4332066378942224E-14, -11.009076470803505, 8.1449755661394213E-15,
    1.257233752131982E-14, -10.678645261814772, 8.1536768652072592E-15,
    1.0812608663697422E-14, -10.348216994586437, 8.1623738968405382E-15,
    9.05374284131662E-15, -10.017786201017476, 8.17107519590838E-15,
    7.29401398369422E-15, -9.6873554625072931, 8.1797769217196761E-15,
    5.53419882254766E-15, -9.3569256948717019, 8.18847779404406E-15,
    3.7745562684494192E-15, -9.0264967036407491, 8.1971782396249873E-15,
    2.0150000178753402E-15, -8.696068016815854, 8.2058791119493708E-15,
    2.5535746377710027E-16, -8.3656399951025673, 8.2145799842737544E-15,
    -1.5042850903211397E-15, 1.701659561004628E-14, -11.39613713831646,
    -8.2389569049356622E-13, 1.612851992127658E-14, -11.390473923212333,
    -8.2372240223484736E-13, 1.524044423250688E-14, -11.384810708148985,
    -8.2354911397612839E-13, 1.4352368543737178E-14, -11.379147493126398,
    -8.2337582571740953E-13, 1.3464292854967481E-14, -11.373484278164979,
    -8.2320253745869056E-13, 1.257665271092008E-14, -11.367823838693702,
    -8.230293341868568E-13, 1.1688577022150379E-14, -11.362160623465199,
    -8.2285604592813793E-13, 1.0800501333380681E-14, -11.356497408275434,
    -8.2268275766941907E-13, 9.91242564461098E-15, -11.350834193130524,
    -8.225094694107E-13, 9.02434995584128E-15, -11.345170978036592,
    -8.2233618115198114E-13, 8.13627426707158E-15, -211.3395077629977,
    -8.2216289289326228E-13, 7.24819857830188E-15, -11.009076470803505,
    -8.2198960463454341E-13, 6.36012288953218E-15, -10.678645261814772,
    -8.2181631637582455E-13, 5.47248274548478E-15, -10.348216994586437,
    -8.2164311310399068E-13, 4.58440705671508E-15, -10.017786201017476,
    -8.2146982484527172E-13, 3.6962878134731506E-15, -9.6873554625072931,
    -8.2129652808786428E-13, 2.8082556791756794E-15, -9.3569256948717019,
    -8.21123248327834E-13, 1.9202670993504398E-15, -9.0264967036407491,
    -8.2094997706649209E-13, 1.03223496505297E-15, -8.696068016815854,
    -8.207766973064617E-13, 1.4420283075549994E-16, -8.3656399951025673,
    -8.2060341754643132E-13, 1.6245000828907414E-14, -8.2126577058286671E-13,
    -11.39613713831646, 1.6053698001722218E-14, -8.2135548721393714E-13,
    -11.390473923212333, 1.5862395174537014E-14, -8.2144520384500767E-13,
    -11.384810708148985, 1.5671092347351818E-14, -8.2153492047607819E-13,
    -11.379147493126398, 1.5479789520166618E-14, -8.2162463710714872E-13,
    -11.373484278164979, 1.528858051486822E-14, -8.2171430973790974E-13,
    -11.367823838693702, 1.5097277687683017E-14, -8.2180402636898027E-13,
    -11.362160623465199, 1.4905974860497817E-14, -8.218937430000508E-13,
    -11.356497408275434, 1.4714672033312617E-14, -8.2198345963112133E-13,
    -11.350834193130524, 1.4523369206127417E-14, -8.2207317626219185E-13,
    -11.345170978036592, 1.4332066378942217E-14, -8.2216289289326228E-13,
    -211.3395077629977, 1.4140763551757017E-14, -8.2225260952433281E-13,
    -11.009076470803505, 1.3949460724571817E-14, -8.2234232615540323E-13,
    -10.678645261814772, 1.3758251719273418E-14, -8.2243199878616436E-13,
    -10.348216994586437, 1.356694889208822E-14, -8.2252171541723478E-13,
    -10.017786201017476, 1.3375636682714337E-14, -8.2261143644833627E-13,
    -9.6873554625072931, 1.3184343237717818E-14, -8.2270114867937584E-13,
    -9.3569256948717019, 1.2993059174909977E-14, -8.2279085651038444E-13,
    -9.0264967036407491, 1.2801765729913459E-14, -8.22880568741424E-13,
    -8.696068016815854, 1.2610472284916938E-14, -8.2297028097246337E-13,
    -8.3656399951025673, -11.06403452996865, 7.14852250985434E-15,
    3.1286269850805016E-14, -11.05853845464685, 7.15849060557088E-15,
    2.957163512895702E-14, -11.053042379365833, 7.16845870128742E-15,
    2.7857000407109017E-14, -11.047546304125596, 7.17842679700396E-15,
    2.6142365685261017E-14, -11.042050228926119, 7.1883948927205012E-15,
    2.4427730963413017E-14, -11.036556847612307, 7.19835809971918E-15,
    2.2713937160997016E-14, -11.031060772168175, 7.20832619543572E-15,
    2.0999302439149019E-14, -11.025564696762782, 7.21829429115226E-15,
    1.928466771730102E-14, -11.020068621396126, 7.2282623868688012E-15,
    1.757003299545302E-14, -11.014572546074326, 7.23823048258534E-15,
    1.585539827360502E-14, -11.009076470803505, 7.2481985783018791E-15,
    1.4140763551757022E-14, -211.00358039558773, 7.2581666740184219E-15,
    1.2426128829909021E-14, -10.67331624321664, 7.26813476973496E-15,
    1.0711494108061023E-14, -10.343054867875509, 7.27809797673364E-15,
    8.99770030564502E-15, -10.012790964674323, 7.28806607245018E-15,
    7.28306558379702E-15, -9.6825270415460221, 7.2980346570385057E-15,
    5.56834677000582E-15, -9.3522639898209956, 7.30800226388326E-15,
    3.8537961401010195E-15, -9.022001623193816, 7.3179693818562278E-15,
    2.13932960213942E-15, -8.6917394860092312, 7.3279369887009822E-15,
    4.2477897223462025E-16, -8.3614779308103, 7.3379045955457351E-15,
    -1.28977165767018E-15, 1.701262956391042E-14, -11.06403452996865,
    -8.2396719391842159E-13, 1.612582067178942E-14, -11.05853845464685,
    -8.2379572706965289E-13, 1.5239011779668418E-14, -11.053042379365833,
    -8.23624260220884E-13, 1.4352202887547418E-14, -11.047546304125596,
    -8.2345279337211517E-13, 1.346539399542642E-14, -11.042050228926119,
    -8.2328132652334636E-13, 1.257902002674442E-14, -11.036556847612307,
    -8.2310994376817685E-13, 1.1692211134623419E-14, -11.031060772168175,
    -8.22938476919408E-13, 1.0805402242502419E-14, -11.025564696762782,
    -8.2276701007063923E-13, 9.91859335038142E-15, -11.020068621396126,
    -8.2259554322187032E-13, 9.03178445826042E-15, -11.014572546074326,
    -8.2242407637310151E-13, 8.14497556613942E-15, -11.009076470803505,
    -8.2225260952433281E-13, 7.2581666740184219E-15, -211.00358039558773,
    -8.22081142675564E-13, 6.37135778189742E-15, -10.67331624321664,
    -8.2190967582679519E-13, 5.48498381321542E-15, -10.343054867875509,
    -8.2173829307162558E-13, 4.59817492109442E-15, -10.012790964674323,
    -8.2156682622285677E-13, 3.71132253662952E-15, -9.6825270415460221,
    -8.21395350964728E-13, 2.8245571368524198E-15, -9.3522639898209956,
    -8.2122389252531915E-13, 1.93783522941922E-15, -9.022001623193816,
    -8.2105244249527022E-13, 1.0510698296421198E-15, -8.6917394860092312,
    -8.2088098405586125E-13, 1.6430442986501993E-16, -8.3614779308103,
    -8.2070952561645249E-13, 1.4034352729350619E-14, -8.210742691179327E-13,
    -11.06403452996865, 1.3888144037939823E-14, -8.2116580715895328E-13,
    -11.05853845464685, 1.3741935346529019E-14, -8.2125734519997386E-13,
    -11.053042379365833, 1.3595726655118219E-14, -8.2134888324099444E-13,
    -11.047546304125596, 1.344951796370742E-14, -8.2144042128201512E-13,
    -11.042050228926119, 1.330338097837382E-14, -8.2153191442944041E-13,
    -11.036556847612307, 1.315717228696302E-14, -8.21623452470461E-13,
    -11.031060772168175, 1.301096359555222E-14, -8.2171499051148157E-13,
    -11.025564696762782, 1.2864754904141419E-14, -8.2180652855250225E-13,
    -11.020068621396126, 1.271854621273062E-14, -8.2189806659352283E-13,
    -11.014572546074326, 1.257233752131982E-14, -8.2198960463454341E-13,
    -11.009076470803505, 1.2426128829909021E-14, -8.22081142675564E-13,
    -211.00358039558773, 1.227992013849822E-14, -8.2217268071658448E-13,
    -10.67331624321664, 1.213378315316462E-14, -8.2226417386400986E-13,
    -10.343054867875509, 1.198757446175382E-14, -8.2235571190503044E-13,
    -10.012790964674323, 1.1841358599735299E-14, -8.2244725443541052E-13,
    -9.6825270415460221, 1.1695157078932219E-14, -8.2253878798707161E-13,
    -9.3522639898209956, 1.1548962728736858E-14, -8.2263031704937319E-13,
    -9.022001623193816, 1.140276120793378E-14, -8.2272185060103427E-13,
    -8.6917394860092312, 1.12565596871307E-14, -8.2281338415269515E-13,
    -8.3614779308103, -10.731932836275458, 6.24777947588094E-15,
    3.064404778168542E-14, -10.726603817585572, 6.25901436824618E-15,
    2.8974507195611817E-14, -10.721274798936461, 6.27024926061142E-15,
    2.7304966609538218E-14, -10.715945780328134, 6.28148415297666E-15,
    2.5635426023464619E-14, -10.710616761760587, 6.2927190453419E-15,
    2.396588543739102E-14, -10.705290355474208, 6.30394842770598E-15,
    2.2297163654939817E-14, -10.699961336661989, 6.3151833200712195E-15,
    2.0627623068866221E-14, -10.694632317890548, 6.3264182124364595E-15,
    1.8958082482792618E-14, -10.689303299157846, 6.3376531048017E-15,
    1.7288541896719019E-14, -10.683974280463879, 6.3488879971669394E-15,
    1.561900131064542E-14, -10.678645261814772, 6.3601228895321794E-15,
    1.3949460724571822E-14, -10.673316243216638, 6.37135778189742E-15,
    1.2279920138498221E-14, -210.66798722467357, 6.38259267426266E-15,
    1.0610379552424622E-14, -10.337892741215558, 6.39382205662674E-15,
    8.9416577699734192E-15, -10.007795728376026, 6.40505694899198E-15,
    7.27211718389982E-15, -9.6776986206234881, 6.4162923923573352E-15,
    5.60249471746398E-15, -9.3476022848090281, 6.4275267337224606E-15,
    3.93303601175262E-15, -9.0175065427876628, 6.4387605240874682E-15,
    2.2636591864035E-15, -8.6874109552229974, 6.4499948654525912E-15,
    5.9420048069214023E-16, -8.3573158665180323, 6.4612292068177158E-15,
    -1.0752582250192198E-15, 1.7008663517774559E-14, -10.731932836275458,
    -8.2403869734327707E-13, 1.612312142230226E-14, -10.726603817585572,
    -8.2386905190445831E-13, 1.5237579326829958E-14, -10.721274798936461,
    -8.2369940646563956E-13, 1.4352037231357659E-14, -10.715945780328134,
    -8.235297610268209E-13, 1.346649513588536E-14, -10.710616761760587,
    -8.2336011558800215E-13, 1.2581387342568761E-14, -10.705290355474208,
    -8.231905533494968E-13, 1.1695845247096459E-14, -10.699961336661989,
    -8.2302090791067815E-13, 1.081030315162416E-14, -10.694632317890548,
    -8.2285126247185939E-13, 9.92476105615186E-15, -10.689303299157846,
    -8.2268161703304064E-13, 9.03921896067956E-15, -10.683974280463879,
    -8.22511971594222E-13, 8.1536768652072592E-15, -10.678645261814772,
    -8.2234232615540323E-13, 7.2681347697349616E-15, -10.673316243216638,
    -8.2217268071658458E-13, 6.38259267426266E-15, -210.66798722467357,
    -8.2200303527776592E-13, 5.4974848809460606E-15, -10.337892741215558,
    -8.2183347303926048E-13, 4.61194278547376E-15, -10.007795728376026,
    -8.2166382760044172E-13, 3.72635725978589E-15, -9.6776986206234881,
    -8.2149417384159167E-13, 2.8408585945291597E-15, -9.3476022848090281,
    -8.2132453672280431E-13, 1.955403359488E-15, -9.0175065427876628,
    -8.2115490792404836E-13, 1.0699046942312699E-15, -8.6874109552229974,
    -8.2098527080526091E-13, 1.8440602897453994E-16, -8.3573158665180323,
    -8.2081563368647355E-13, 1.1823704629793819E-14, -8.2088276765299868E-13,
    -10.731932836275458, 1.172259007415742E-14, -8.2097612710396942E-13,
    -10.726603817585572, 1.1621475518521016E-14, -8.2106948655494005E-13,
    -10.721274798936461, 1.1520360962884617E-14, -8.2116284600591079E-13,
    -10.715945780328134, 1.1419246407248218E-14, -8.2125620545688153E-13,
    -10.710616761760587, 1.131818144187942E-14, -8.21349519120971E-13,
    -10.705290355474208, 1.1217066886243018E-14, -8.2144287857194171E-13,
    -10.699961336661989, 1.1115952330606619E-14, -8.2153623802291244E-13,
    -10.694632317890548, 1.1014837774970218E-14, -8.2162959747388308E-13,
    -10.689303299157846, 1.0913723219333818E-14, -8.2172295692485381E-13,
    -10.683974280463879, 1.0812608663697419E-14, -8.2181631637582455E-13,
    -10.678645261814772, 1.0711494108061018E-14, -8.2190967582679519E-13,
    -10.673316243216638, 1.0610379552424619E-14, -8.2200303527776582E-13,
    -210.66798722467357, 1.0509314587055818E-14, -8.2209634894185537E-13,
    -10.337892741215558, 1.0408200031419419E-14, -8.22189708392826E-13,
    -10.007795728376026, 1.0307080516756259E-14, -8.2228307242248478E-13,
    -9.6776986206234881, 1.0205970920146619E-14, -8.2237642729476748E-13,
    -9.3476022848090281, 1.0104866282563738E-14, -8.2246977758836183E-13,
    -9.0175065427876628, 1.0003756685954099E-14, -8.2256313246064443E-13,
    -8.6874109552229974, 9.9026470893444588E-15, -8.2265648733292693E-13,
    -8.3573158665180323, -10.399835731596752, 5.3474781991581405E-15,
    3.0002140681702221E-14, -10.394673604804272, 5.3599792668887808E-15,
    2.8377672115593417E-14, -10.389511478052567, 5.37248033461942E-15,
    2.6753203549484617E-14, -10.38434935134164, 5.38498140235006E-15,
    2.5128734983375817E-14, -10.379187224671494, 5.3974824700807E-15,
    2.3504266417267019E-14, -10.374027628738462, 5.4099774068315805E-15,
    2.1880594549817416E-14, -10.368865501844036, 5.42247847456222E-15,
    2.0256125983708619E-14, -10.36370337497, 5.43497954229286E-15,
    1.8631657417599818E-14, -10.35854124813674, 5.4474806100235008E-15,
    1.7007188851491021E-14, -10.353379121342218, 5.45998167775414E-15,
    1.538272028538222E-14, -10.348216994586437, 5.47248274548478E-15,
    1.3758251719273422E-14, -10.343054867875509, 5.48498381321542E-15,
    1.213378315316462E-14, -10.33789274121556, 5.4974848809460606E-15,
    1.0509314587055822E-14, -210.332733145307, 5.50997981769694E-15,
    8.8856427196062188E-15, -10.002802941240949, 5.52248088542758E-15,
    7.2611741534974184E-15, -9.6728725672699643, 5.5349825662561957E-15,
    5.6366259175227E-15, -9.3429428657798645, 5.5474830208888611E-15,
    4.0122370212798192E-15, -9.013013666788229, 5.5599828624235483E-15,
    2.3879277949028598E-15, -8.6830845472654339, 5.5724833170562113E-15,
    7.635388986599801E-16, -8.3531558434539619, 5.5849837716888766E-15,
    -8.6084999758289989E-16, 1.700469941673244E-14, -10.399835731596752,
    -8.2411016570024385E-13, 1.6120423496625538E-14, -10.394673604804272,
    -8.2394234077808936E-13, 1.5236147576518639E-14, -10.389511478052567,
    -8.2377451585593477E-13, 1.4351871656411738E-14, -10.38434935134164,
    -8.2360669093378027E-13, 1.346759573630484E-14, -10.379187224671494,
    -8.2343886601162578E-13, 1.258375349737504E-14, -10.374027628738462,
    -8.2327112339693682E-13, 1.1699477577268139E-14, -10.368865501844036,
    -8.2310329847478233E-13, 1.0815201657161241E-14, -10.36370337497,
    -8.2293547355262784E-13, 9.93092573705434E-15, -10.35854124813674,
    -8.2276764863047324E-13, 9.04664981694744E-15, -10.353379121342218,
    -8.2259982370831875E-13, 8.16237389684054E-15, -10.348216994586437,
    -8.2243199878616425E-13, 7.278097976733641E-15, -10.343054867875509,
    -8.2226417386400976E-13, 6.3938220566267407E-15, -10.33789274121556,
    -8.2209634894185527E-13, 5.50997981769694E-15, -210.332733145307,
    -8.2192860632716631E-13, 4.6257038975900405E-15, -10.002802941240949,
    -8.2176078140501172E-13, 3.74138460936543E-15, -9.6728725672699643,
    -8.2159294825211064E-13, 2.8571520573762395E-15, -9.3429428657798645,
    -8.2142513156070273E-13, 1.9729628735047598E-15, -9.013013666788229,
    -8.2125732310004141E-13, 1.0887303215155698E-15, -8.6830845472654339,
    -8.210895064086334E-13, 2.0449776952637994E-16, -8.3531558434539619,
    -8.2092168971722549E-13, 9.61414071272822E-15, -8.2069136010737068E-13,
    -10.399835731596752, 9.55809817705662E-15, -8.2078654007500557E-13,
    -10.394673604804272, 9.5020556413850185E-15, -8.2088172004264047E-13,
    -10.389511478052567, 9.4460131057134214E-15, -8.2097690001027537E-13,
    -10.38434935134164, 9.38997057004182E-15, -8.2107207997791026E-13,
    -10.379187224671494, 9.33395551967462E-15, -8.211672132658162E-13,
    -10.374027628738462, 9.27791298400302E-15, -8.212623932334511E-13,
    -10.368865501844036, 9.22187044833142E-15, -8.21357573201086E-13,
    -10.36370337497, 9.16582791265982E-15, -8.2145275316872089E-13,
    -10.35854124813674, 9.1097853769882183E-15, -8.2154793313635579E-13,
    -10.353379121342218, 9.05374284131662E-15, -8.2164311310399068E-13,
    -10.348216994586437, 8.99770030564502E-15, -8.2173829307162558E-13,
    -10.343054867875509, 8.9416577699734192E-15, -8.2183347303926048E-13,
    -10.33789274121556, 8.8856427196062188E-15, -8.2192860632716631E-13,
    -210.332733145307, 8.82960018393462E-15, -8.2202378629480121E-13,
    -10.002802941240949, 8.7735548997325784E-15, -8.22118970930409E-13,
    -9.6728725672699643, 8.71751511259142E-15, -8.22214146230071E-13,
    -9.3429428657798645, 8.6614780739807E-15, -8.2230931686176009E-13,
    -9.013013666788229, 8.60543828683954E-15, -8.2240449216142213E-13,
    -8.6830845472654339, 8.549398499698379E-15, -8.22499667461084E-13,
    -8.3531558434539619, -10.067736115828113, 4.44673516518474E-15,
    2.9359918612582618E-14, -10.062740879407485, 4.46050302956408E-15,
    2.7780544182248218E-14, -10.057745643027634, 4.47427089394342E-15,
    2.6201169751913818E-14, -10.05275040668856, 4.4880387583227591E-15,
    2.4621795321579418E-14, -10.047755170390262, 4.5018066227021E-15,
    2.3042420891245018E-14, -10.042762383244991, 4.51556773481838E-15,
    2.1463821043760217E-14, -10.037767146763196, 4.5293355991977193E-15,
    1.988444661342582E-14, -10.032771910281403, 4.5431034635770592E-15,
    1.830507218309142E-14, -10.027776673819995, 4.5568713279564E-15,
    1.672569775275702E-14, -10.022781437399365, 4.57063919233574E-15,
    1.514632332242262E-14, -10.017786201017476, 4.5844070567150795E-15,
    1.3566948892088222E-14, -10.012790964674322, 4.5981749210944193E-15,
    1.198757446175382E-14, -10.007795728376026, 4.61194278547376E-15,
    1.0408200031419422E-14, -10.002802941240949, 4.62570389759004E-15,
    8.8296001839346186E-15, -209.99780770478361, 4.63947176196938E-15,
    7.2502257536002186E-15, -9.6680441462189837, 4.6532403015750253E-15,
    5.6707738649808604E-15, -9.338281160667993, 4.66700749072806E-15,
    4.0914768929314195E-15, -9.008518586316832, 4.680774004654788E-15,
    2.51225737916694E-15, -8.6787560164567736, 4.6945411938078211E-15,
    9.3296040711750028E-16, -8.3489937791616935, 4.7083083829608558E-15,
    -6.4633656493193974E-16, 1.700073337059658E-14, -10.067736115828113,
    -8.2418166912509923E-13, 1.6117724247138378E-14, -10.062740879407485,
    -8.2401566561289489E-13, 1.523471512368018E-14, -10.057745643027634,
    -8.2384966210069035E-13, 1.4351706000221978E-14, -10.05275040668856,
    -8.23683658588486E-13, 1.346869687676378E-14, -10.047755170390262,
    -8.2351765507628157E-13, 1.2586120813199381E-14, -10.042762383244991,
    -8.2335173297825678E-13, 1.1703111689741179E-14, -10.037767146763196,
    -8.2318572946605244E-13, 1.0820102566282981E-14, -10.032771910281403,
    -8.23019725953848E-13, 9.93709344282478E-15, -10.027776673819995,
    -8.2285372244164356E-13, 9.0540843193665792E-15, -10.022781437399365,
    -8.2268771892943912E-13, 8.1710751959083792E-15, -10.017786201017476,
    -8.2252171541723478E-13, 7.28806607245018E-15, -10.012790964674322,
    -8.2235571190503034E-13, 6.4050569489919807E-15, -10.007795728376026,
    -8.22189708392826E-13, 5.52248088542758E-15, -10.002802941240949,
    -8.2202378629480121E-13, 4.63947176196938E-15, -209.99780770478361,
    -8.2185778278259677E-13, 3.7564193325218E-15, -9.6680441462189837,
    -8.2169177112897439E-13, 2.8734535150529794E-15, -9.338281160667993,
    -8.21525775758188E-13, 1.9905310035735397E-15, -9.008518586316832,
    -8.2135978852881954E-13, 1.1075651861047199E-15, -8.6787560164567736,
    -8.21193793158033E-13, 2.2459936863589995E-16, -8.3489937791616935,
    -8.2102779778724665E-13, 7.4034926131714175E-15, -8.2049985864243666E-13,
    -10.067736115828113, 7.39254421327422E-15, -8.2059686002002171E-13,
    -10.062740879407485, 7.3815958133770164E-15, -8.2069386139760666E-13,
    -10.057745643027634, 7.37064741347982E-15, -8.2079086277519171E-13,
    -10.05275040668856, 7.3596990135826183E-15, -8.2088786415277666E-13,
    -10.047755170390262, 7.3487559831802184E-15, -8.2098481795734676E-13,
    -10.042762383244991, 7.3378075832830187E-15, -8.2108181933493181E-13,
    -10.037767146763196, 7.3268591833858189E-15, -8.2117882071251677E-13,
    -10.032771910281403, 7.3159107834886191E-15, -8.2127582209010182E-13,
    -10.027776673819995, 7.3049623835914177E-15, -8.2137282346768677E-13,
    -10.022781437399365, 7.2940139836942179E-15, -8.2146982484527182E-13,
    -10.017786201017476, 7.2830655837970181E-15, -8.2156682622285687E-13,
    -10.012790964674322, 7.2721171838998183E-15, -8.2166382760044172E-13,
    -10.007795728376026, 7.2611741534974184E-15, -8.2176078140501182E-13,
    -10.002802941240949, 7.2502257536002186E-15, -8.2185778278259687E-13,
    -209.99780770478361, 7.2392768167535378E-15, -8.2195478891748332E-13,
    -9.6680441462189837, 7.2283289538058191E-15, -8.2205178553776687E-13,
    -9.338281160667993, 7.2173816278075783E-15, -8.2214877740074883E-13,
    -9.008518586316832, 7.20643376485986E-15, -8.2224577402103239E-13,
    -8.6787560164567736, 7.19548590191214E-15, -8.2234277064131575E-13,
    -8.3489937791616935, -9.7356373050036282, 3.5459479554862806E-15,
    2.8717665046549376E-14, -9.7308088839628422, 3.5609826786426505E-15,
    2.7183386963570337E-14, -9.72598046296283, 3.57601740179902E-15,
    2.5649108880591298E-14, -9.7211520420035988, 3.59105212495539E-15,
    2.4114830797612256E-14, -9.7163236210851434, 3.6060868481117604E-15,
    2.2580552714633217E-14, -9.71149756773162, 3.6211141976913E-15,
    2.1047027097609535E-14, -9.7066691466704427, 3.6361489208476704E-15,
    1.95127490146305E-14, -9.7018407256092676, 3.65118364400404E-15,
    1.797847093165146E-14, -9.69701230454809, 3.66621836716041E-15,
    1.6444192848672421E-14, -9.6921838835073011, 3.68125309031678E-15,
    1.4909914765693379E-14, -9.6873554625072931, 3.69628781347315E-15,
    1.3375636682714341E-14, -9.6825270415460221, 3.71132253662952E-15,
    1.18413585997353E-14, -9.6776986206234881, 3.72635725978589E-15,
    1.0307080516756262E-14, -9.6728725672699643, 3.74138460936543E-15,
    8.7735548997325784E-15, -9.6680441462189837, 3.7564193325218E-15,
    7.23927681675354E-15, -209.66321548842987, 3.7714547930358532E-15,
    5.7049234871789638E-15, -9.3336192269719067, 3.7864887788345408E-15,
    4.17072065079546E-15, -9.0040232854127122, 3.8015220272755441E-15,
    2.6365930610074917E-15, -8.6744272733672823, 3.8165560130742309E-15,
    1.102390224623988E-15, -8.3448315107466069, 3.8315899988729185E-15,
    -4.31812611759516E-16, 1.6996767129951346E-14, -9.7356373050036282,
    -8.2425317605674349E-13, 1.6115024865270176E-14, -9.7308088839628422,
    -8.2408899404381785E-13, 1.5233282600589005E-14, -9.72598046296283,
    -8.23924812030892E-13, 1.4351540335907835E-14, -9.7211520420035988,
    -8.2376063001796635E-13, 1.3469798071226666E-14, -9.7163236210851434,
    -8.2359644800504061E-13, 1.2588488245125527E-14, -9.71149756773162,
    -8.2343234651296479E-13, 1.1706745980444355E-14, -9.7066691466704427,
    -8.2326816450003915E-13, 1.0825003715763186E-14, -9.7018407256092676,
    -8.231039824871134E-13, 9.9432614510820153E-15, -9.69701230454809,
    -8.2293980047418766E-13, 9.0615191864008465E-15, -9.6921838835073011,
    -8.2277561846126191E-13, 8.1797769217196761E-15, -9.6873554625072931,
    -8.2261143644833627E-13, 7.2980346570385072E-15, -9.6825270415460221,
    -8.2244725443541052E-13, 6.416292392357336E-15, -9.6776986206234881,
    -8.2228307242248488E-13, 5.5349825662561965E-15, -9.6728725672699643,
    -8.2211897093040906E-13, 4.6532403015750261E-15, -9.6680441462189837,
    -8.2195478891748322E-13, 3.7714547930358532E-15, -209.66321548842987,
    -8.2179059885247251E-13, 2.8897557722126856E-15, -9.3336192269719067,
    -8.2162642489163183E-13, 2.0080999952475221E-15, -9.0040232854127122,
    -8.2146225898287611E-13, 1.1264009744243549E-15, -8.6744272733672823,
    -8.2129808502203533E-13, 2.4470195360118795E-16, -8.3448315107466069,
    -8.2113391106119464E-13, 5.1927360953655E-15, -8.2030834778557209E-13,
    -9.7356373050036282, 5.2268840428236609E-15, -8.2040717066243574E-13,
    -9.7308088839628422, 5.26103199028182E-15, -8.2050599353929948E-13,
    -9.72598046296283, 5.2951799377399808E-15, -8.2060481641616313E-13,
    -9.7211520420035988, 5.3293278851981395E-15, -8.2070363929302687E-13,
    -9.7163236210851434, 5.3634590852568594E-15, -8.2080241370354595E-13,
    -9.71149756773162, 5.3976070327150205E-15, -8.2090123658040959E-13,
    -9.7066691466704427, 5.43175498017318E-15, -8.2100005945727334E-13,
    -9.7018407256092676, 5.4659029276313404E-15, -8.21098882334137E-13,
    -9.69701230454809, 5.5000508750894992E-15, -8.2119770521100073E-13,
    -9.6921838835073011, 5.5341988225476587E-15, -8.2129652808786438E-13,
    -9.6873554625072931, 5.56834677000582E-15, -8.2139535096472812E-13,
    -9.6825270415460221, 5.6024947174639794E-15, -8.2149417384159167E-13,
    -9.6776986206234881, 5.6366259175226985E-15, -8.2159294825211074E-13,
    -9.6728725672699643, 5.6707738649808612E-15, -8.2169177112897449E-13,
    -9.6680441462189837, 5.704923487178963E-15, -8.2179059885247261E-13,
    -209.66321548842987, 5.73906975989718E-15, -8.2188941688270188E-13,
    -9.3336192269719067, 5.7732143578754514E-15, -8.2198823006629657E-13,
    -9.0040232854127122, 5.807360630593668E-15, -8.2208704809652584E-13,
    -8.6744272733672823, 5.8415069033118829E-15, -8.2218586612675491E-13,
    -8.3448315107466069, -9.4035404596997552, 2.64524909723794E-15,
    2.8075474474343419E-14, -9.3988787545878836, 2.66155055491468E-15,
    2.6586288315557816E-14, -9.3942170494963975, 2.67785201259142E-15,
    2.5097102156772217E-14, -9.389555344445693, 2.69415347026816E-15,
    2.3607915997986618E-14, -9.3848936394357629, 2.7104549279448994E-15,
    2.2118729839201016E-14, -9.3802342204106779, 2.7267483907919795E-15,
    2.0630274031645815E-14, -9.3755725152988045, 2.74304984846872E-15,
    1.914108787286022E-14, -9.3709108101869329, 2.75935130614546E-15,
    1.765190171407462E-14, -9.36624910507506, 2.7756527638221995E-15,
    1.6162715555289018E-14, -9.3615873999631862, 2.79195422149894E-15,
    1.4673529396503419E-14, -9.3569256948717019, 2.8082556791756794E-15,
    1.3184343237717822E-14, -9.3522639898209938, 2.82455713685242E-15,
    1.1695157078932219E-14, -9.3476022848090281, 2.84085859452916E-15,
    1.0205970920146622E-14, -9.3429428657798645, 2.85715205737624E-15,
    8.71751511259142E-15, -9.338281160667993, 2.87345351505298E-15,
    7.2283289538058191E-15, -9.3336192269719067, 2.8897557722126852E-15,
    5.7390697598971795E-15, -209.32895775049727, 2.9060564304064598E-15,
    4.2499566362346193E-15, -8.9995284254066572, 2.9223562891172677E-15,
    2.7609165476950996E-15, -8.6700989548496441, 2.9386569473110415E-15,
    1.2718034240325402E-15, -8.3406696505771585, 2.954957605504816E-15,
    -2.1730969963001983E-16, 1.6992801278324859E-14, -9.4035404596997552,
    -8.2432467597481008E-13, 1.6112325748164058E-14, -9.3988787545878836,
    -8.2416231528250585E-13, 1.523185021800326E-14, -9.3942170494963975,
    -8.2399995459020162E-13, 1.4351374687842459E-14, -9.389555344445693,
    -8.2383759389789739E-13, 1.347089915768166E-14, -9.3848936394357629,
    -8.2367523320559316E-13, 1.259085544484806E-14, -9.3802342204106779,
    -8.2351295214089678E-13, 1.1710379914687259E-14, -9.3755725152988045,
    -8.2335059144859265E-13, 1.082990438452646E-14, -9.3709108101869329,
    -8.2318823075628842E-13, 9.94942885436566E-15, -9.36624910507506,
    -8.2302587006398419E-13, 9.0689533242048591E-15, -9.3615873999631862,
    -8.2286350937168E-13, 8.18847779404406E-15, -9.3569256948717019,
    -8.2270114867937573E-13, 7.3080022638832617E-15, -9.3522639898209938,
    -8.2253878798707161E-13, 6.4275267337224606E-15, -9.3476022848090281,
    -8.2237642729476738E-13, 5.54748302088886E-15, -9.3429428657798645,
    -8.22214146230071E-13, 4.66700749072806E-15, -9.338281160667993,
    -8.2205178553776677E-13, 3.78648877883454E-15, -9.3336192269719067,
    -8.2188941688270178E-13, 2.9060564304064594E-15, -209.32895775049727,
    -8.2172706415315831E-13, 2.0256672637111E-15, -8.9995284254066572,
    -8.2156471938637582E-13, 1.1452349152830198E-15, -8.6700989548496441,
    -8.2140236665683225E-13, 2.6480256685493997E-16, -8.3406696505771585,
    -8.2124001392728889E-13, 2.9821964140578193E-15, -8.2011685571256862E-13,
    -9.4035404596997552, 3.0614362857094204E-15, -8.2021749991005389E-13,
    -9.3988787545878836, 3.1406761573610176E-15, -8.20318144107539E-13,
    -9.3942170494963975, 3.2199160290126203E-15, -8.2041878830502431E-13,
    -9.389555344445693, 3.299155900664219E-15, -8.2051943250250947E-13,
    -9.3848936394357629, 3.3783569101914192E-15, -8.20620027340408E-13,
    -9.3802342204106779, 3.4575967818430195E-15, -8.2072067153789315E-13,
    -9.3755725152988045, 3.5368366534946194E-15, -8.2082131573537841E-13,
    -9.3709108101869329, 3.6160765251462193E-15, -8.2092195993286357E-13,
    -9.36624910507506, 3.6953163967978189E-15, -8.2102260413034883E-13,
    -9.3615873999631862, 3.7745562684494192E-15, -8.21123248327834E-13,
    -9.3569256948717019, 3.8537961401010187E-15, -8.2122389252531925E-13,
    -9.3522639898209938, 3.9330360117526182E-15, -8.2132453672280431E-13,
    -9.3476022848090281, 4.0122370212798184E-15, -8.2142513156070283E-13,
    -9.3429428657798645, 4.0914768929314195E-15, -8.21525775758188E-13,
    -9.338281160667993, 4.1707206507954589E-15, -8.2162642489163183E-13,
    -9.3336192269719067, 4.2499566362346193E-15, -8.2172706415315841E-13,
    -209.32895775049727, 4.3291887354613391E-15, -8.2182769847872623E-13,
    -8.9995284254066572, 4.4084247209004996E-15, -8.2192833774025281E-13,
    -8.6700989548496441, 4.4876607063396585E-15, -8.220289770017792E-13,
    -8.3406696505771585, -9.0714453037624718, 1.7445944147146601E-15,
    2.7433315399051097E-14, -9.0669502232951515, 1.7621625447834399E-15,
    2.5989218952877978E-14, -9.0624551428278313, 1.77973067485222E-15,
    2.4545122506704858E-14, -9.0579600623809, 1.797298804921E-15,
    2.3101026060531739E-14, -9.053464981974745, 1.81486693498978E-15,
    2.1656929614358619E-14, -9.0489721059773522, 1.8324264490065397E-15,
    2.0213541405775575E-14, -9.04447702551003, 1.84999457907532E-15,
    1.8769444959602459E-14, -9.0399819450427117, 1.8675627091441E-15,
    1.7325348513429339E-14, -9.03548686457539, 1.88513083921288E-15,
    1.588125206725622E-14, -9.03099178410807, 1.90269896928166E-15,
    1.44371556210831E-14, -9.0264967036407491, 1.9202670993504398E-15,
    1.2993059174909982E-14, -9.022001623193816, 1.93783522941922E-15,
    1.1548962728736859E-14, -9.0175065427876628, 1.9554033594880004E-15,
    1.0104866282563741E-14, -9.013013666788229, 1.9729628735047598E-15,
    8.6614780739807E-15, -9.008518586316832, 1.99053100357354E-15,
    7.2173816278075783E-15, -9.004023285412714, 2.0080999952475221E-15,
    5.7732143578754514E-15, -8.9995284254066572, 2.0256672637111E-15,
    4.3291887354613391E-15, -208.99503378584149, 2.0432336705694759E-15,
    2.8852339368062357E-15, -8.665770848614871, 2.0608009390330539E-15,
    1.441208314392124E-15, -8.33650799453053, 2.0783682074966322E-15,
    -2.8173080219879851E-18, 1.6988835621207746E-14, -9.0714453037624718,
    -8.2439617238608767E-13, 1.6109626763438987E-14, -9.0669502232951515,
    -8.2423563292507642E-13, 1.5230417905670227E-14, -9.0624551428278313,
    -8.2407509346406506E-13, 1.4351209047901468E-14, -9.0579600623809,
    -8.2391455400305381E-13, 1.3472000190132708E-14, -9.053464981974745,
    -8.2375401454204256E-13, 1.2593222528468788E-14, -9.0489721059773522,
    -8.2359355381544081E-13, 1.1714013670700027E-14, -9.04447702551003,
    -8.2343301435442956E-13, 1.0834804812931267E-14, -9.0399819450427117,
    -8.232724748934183E-13, 9.9555959551625079E-15, -9.03548686457539,
    -8.2311193543240695E-13, 9.0763870973937484E-15, -9.03099178410807,
    -8.2295139597139559E-13, 8.1971782396249873E-15, -9.0264967036407491,
    -8.2279085651038434E-13, 7.31796938185623E-15, -9.022001623193816,
    -8.2263031704937309E-13, 6.4387605240874682E-15, -9.0175065427876628,
    -8.2246977758836183E-13, 5.5599828624235483E-15, -9.013013666788229,
    -8.2230931686176009E-13, 4.680774004654788E-15, -9.008518586316832,
    -8.2214877740074873E-13, 3.8015220272755441E-15, -9.004023285412714,
    -8.2198823006629647E-13, 2.9223562891172677E-15, -8.9995284254066572,
    -8.2182769847872623E-13, 2.0432336705694759E-15, -208.99503378584149,
    -8.2166717476459689E-13, 1.1640679324111999E-15, -8.665770848614871,
    -8.2150664317702645E-13, 2.8490219425292393E-16, -8.33650799453053,
    -8.2134611158945621E-13, 7.7176515099925918E-16, -8.1992537303149591E-13,
    -9.0714453037624718, 8.9609473526334057E-16, -8.2002783846027395E-13,
    -9.0669502232951515, 1.0204243195274188E-15, -8.2013030388905208E-13,
    -9.0624551428278313, 1.1447539037915002E-15, -8.2023276931783022E-13,
    -9.0579600623809, 1.2690834880555792E-15, -8.2033523474660836E-13,
    -9.053464981974745, 1.3933520965549396E-15, -8.204376499226015E-13,
    -9.0489721059773522, 1.5176816808190202E-15, -8.2054011535137964E-13,
    -9.04447702551003, 1.6420112650830996E-15, -8.2064258078015778E-13,
    -9.0399819450427117, 1.7663408493471798E-15, -8.2074504620893591E-13,
    -9.03548686457539, 1.8906704336112588E-15, -8.2084751163771405E-13,
    -9.03099178410807, 2.01500001787534E-15, -8.2094997706649209E-13,
    -9.0264967036407491, 2.13932960213942E-15, -8.2105244249527022E-13,
    -9.022001623193816, 2.2636591864034994E-15, -8.2115490792404826E-13,
    -9.0175065427876628, 2.3879277949028594E-15, -8.2125732310004141E-13,
    -9.013013666788229, 2.51225737916694E-15, -8.2135978852881954E-13,
    -9.008518586316832, 2.6365930610074913E-15, -8.2146225898287611E-13,
    -9.004023285412714, 2.7609165476950996E-15, -8.2156471938637582E-13,
    -8.9995284254066572, 2.8852339368062357E-15, -8.2166717476459689E-13,
    -208.99503378584149, 3.0095574234938437E-15, -8.217696351680965E-13,
    -8.665770848614871, 3.1338809101814516E-15, -8.21872095571596E-13,
    -8.33650799453053, -8.7393512020126352, 8.4389555646632023E-16,
    2.6791124826845137E-14, -8.7350226712060124, 8.627304210554701E-16,
    2.5392120304865457E-14, -8.7306941403993878, 8.8156528564462016E-16,
    2.3993115782885777E-14, -8.7263656095927669, 9.0040015023377E-16,
    2.2594111260906098E-14, -8.7220370788065313, 9.1923501482291989E-16,
    2.1195106738926418E-14, -8.71771067084897, 9.3806064210722022E-16,
    1.9796788339811856E-14, -8.7133821400423468, 9.5689550669637E-16,
    1.8397783817832179E-14, -8.709053609235724, 9.7573037128552E-16,
    1.69987792958525E-14, -8.7047250784291, 9.9456523587467E-16,
    1.559977477387282E-14, -8.7003965476224767, 1.0134001004638199E-15,
    1.4200770251893141E-14, -8.696068016815854, 1.0322349650529701E-15,
    1.2801765729913461E-14, -8.6917394860092312, 1.0510698296421202E-15,
    1.140276120793378E-14, -8.6874109552229974, 1.0699046942312703E-15,
    1.0003756685954102E-14, -8.6830845472654339, 1.08873032151557E-15,
    8.6054382868395392E-15, -8.6787560164567736, 1.10756518610472E-15,
    7.20643376485986E-15, -8.6744272733672823, 1.1264009744243549E-15,
    5.807360630593668E-15, -8.6700989548496441, 1.1452349152830202E-15,
    4.4084247209004988E-15, -8.665770848614871, 1.1640679324112003E-15,
    3.0095574234938437E-15, -208.6614425300952, 1.1829018732698649E-15,
    1.6106215138006761E-15, -8.3323461343610816, 1.20173581412853E-15,
    2.1168560410750796E-16, 1.6984869769581263E-14, -8.7393512020126352,
    -8.2446767230415426E-13, 1.6106927646332872E-14, -8.7350226712060124,
    -8.2430895416376452E-13, 1.5228985523084482E-14, -8.7306941403993878,
    -8.2415023602337468E-13, 1.4351043399836092E-14, -8.7263656095927669,
    -8.2399151788298494E-13, 1.3473101276587702E-14, -8.7220370788065313,
    -8.238327997425951E-13, 1.2595589728191322E-14, -8.71771067084897,
    -8.236741594433728E-13, 1.171764760494293E-14, -8.7133821400423468,
    -8.2351544130298306E-13, 1.0839705481694541E-14, -8.709053609235724,
    -8.2335672316259332E-13, 9.9617633584461513E-15, -8.7047250784291,
    -8.2319800502220348E-13, 9.0838212351977626E-15, -8.7003965476224767,
    -8.2303928688181364E-13, 8.2058791119493724E-15, -8.696068016815854,
    -8.2288056874142391E-13, 7.3279369887009838E-15, -8.6917394860092312,
    -8.2272185060103417E-13, 6.449994865452592E-15, -8.6874109552229974,
    -8.2256313246064443E-13, 5.5724833170562121E-15, -8.6830845472654339,
    -8.2240449216142213E-13, 4.6945411938078219E-15, -8.6787560164567736,
    -8.2224577402103219E-13, 3.8165560130742309E-15, -8.6744272733672823,
    -8.2208704809652574E-13, 2.9386569473110415E-15, -8.6700989548496441,
    -8.2192833774025271E-13, 2.0608009390330539E-15, -8.665770848614871,
    -8.217696351680965E-13, 1.1829018732698649E-15, -208.6614425300952,
    -8.2161092481182338E-13, 3.0500280750667595E-16, -8.3323461343610816,
    -8.2145221445555035E-13, -1.43877453030842E-15, -8.1973388095849244E-13,
    -8.7393512020126352, -1.2693530218508995E-15, -8.198381677078921E-13,
    -8.7350226712060124, -1.09993151339338E-15, -8.1994245445729175E-13,
    -8.7306941403993878, -9.305100049358595E-16, -8.200467412066914E-13,
    -8.7263656095927669, -7.6108849647834011E-16, -8.2015102795609105E-13,
    -8.7220370788065313, -5.9175007851049984E-16, -8.2025526355946354E-13,
    -8.71771067084897, -4.2232857005297966E-16, -8.203595503088632E-13,
    -8.7133821400423468, -2.5290706159545988E-16, -8.2046383705826285E-13,
    -8.709053609235724, -8.34855531379397E-17, -8.205681238076625E-13,
    -8.7047250784291, 8.5935955319579692E-17, -8.2067241055706215E-13,
    -8.7003965476224767, 2.5535746377710027E-16, -8.207766973064617E-13,
    -8.696068016815854, 4.2477897223462005E-16, -8.2088098405586136E-13,
    -8.6917394860092312, 5.9420048069213984E-16, -8.2098527080526091E-13,
    -8.6874109552229974, 7.6353889865997991E-16, -8.210895064086335E-13,
    -8.6830845472654339, 9.329604071175E-16, -8.2119379315803315E-13,
    -8.6787560164567736, 1.102390224623988E-15, -8.2129808502203543E-13,
    -8.6744272733672823, 1.27180342403254E-15, -8.2140236665683235E-13,
    -8.6700989548496441, 1.4412083143921236E-15, -8.2150664317702655E-13,
    -8.665770848614871, 1.6106215138006761E-15, -8.2161092481182348E-13,
    -208.6614425300952, 1.780034713209228E-15, -8.217152064466203E-13,
    -8.3323461343610816, -8.4072585967970479, -5.6803301782019976E-17,
    2.6148934254639176E-14, -8.4030965325047813, -3.6701702672499991E-17,
    2.4795021656852936E-14, -8.3989344682125129, -1.6600103562980007E-17,
    2.3441109059066697E-14, -8.3947724039202445, 3.5014955465399778E-18,
    2.2087196461280457E-14, -8.3906103396279761, 2.3603094656059913E-17,
    2.0733283863494218E-14, -8.3864503165639075, 4.3694835207899946E-17,
    1.9380035273848136E-14, -8.3822882522716373, 6.3796434317419931E-17,
    1.80261226760619E-14, -8.37812618797937, 8.3898033426940014E-17,
    1.667221007827566E-14, -8.3739641236871023, 1.0399963253646005E-16,
    1.531829748048942E-14, -8.3698020593948339, 1.2410123164597988E-16,
    1.3964384882703179E-14, -8.3656399951025673, 1.4420283075549997E-16,
    1.2610472284916941E-14, -8.3614779308102989, 1.6430442986502E-16,
    1.12565596871307E-14, -8.3573158665180323, 1.8440602897454003E-16,
    9.9026470893444619E-15, -8.35315584345396, 2.0449776952637997E-16,
    8.549398499698379E-15, -8.3489937791616935, 2.245993686359E-16,
    7.19548590191214E-15, -8.3448315107466069, 2.4470195360118795E-16,
    5.8415069033118837E-15, -8.3406696505771585, 2.6480256685493997E-16,
    4.4876607063396593E-15, -8.33650799453053, 2.8490219425292403E-16,
    3.1338809101814516E-15, -8.3323461343610816, 3.0500280750667595E-16,
    1.780034713209228E-15, -208.32818427419164, 3.2510342076042796E-16,
    4.2618851623700396E-16, 1.6980903917954776E-14, -8.4072585967970479,
    -8.2453917222222084E-13, 1.6104228529226755E-14, -8.4030965325047813,
    -8.2438227540245252E-13, 1.5227553140498734E-14, -8.3989344682125129,
    -8.2422537858268419E-13, 1.4350877751770716E-14, -8.3947724039202445,
    -8.24068481762916E-13, 1.3474202363042697E-14, -8.3906103396279761,
    -8.2391158494314775E-13, 1.2597956927913856E-14, -8.3864503165639075,
    -8.2375476507130479E-13, 1.1721281539185835E-14, -8.3822882522716373,
    -8.2359786825153657E-13, 1.0844606150457815E-14, -8.37812618797937,
    -8.2344097143176834E-13, 9.9679307617297962E-15, -8.3739641236871023,
    -8.2328407461199992E-13, 9.0912553730017753E-15, -8.3698020593948339,
    -8.231271777922317E-13, 8.214579984273756E-15, -8.3656399951025673,
    -8.2297028097246347E-13, 7.3379045955457367E-15, -8.3614779308102989,
    -8.2281338415269525E-13, 6.4612292068177166E-15, -8.3573158665180323,
    -8.22656487332927E-13, 5.5849837716888766E-15, -8.35315584345396,
    -8.2249966746108407E-13, 4.7083083829608558E-15, -8.3489937791616935,
    -8.2234277064131575E-13, 3.8315899988729185E-15, -8.3448315107466069,
    -8.2218586612675491E-13, 2.9549576055048156E-15, -8.3406696505771585,
    -8.220289770017792E-13, 2.0783682074966318E-15, -8.33650799453053,
    -8.2187209557159611E-13, 1.2017358141285298E-15, -8.3323461343610816,
    -8.217152064466203E-13, 3.2510342076042796E-16, -208.32818427419164,
    -8.215583173216446E-13, -3.6493142116161E-15, -8.1954238888548908E-13,
    -8.4072585967970479, -3.4348007789651395E-15, -8.1964849695551025E-13,
    -8.4030965325047813, -3.22028734631418E-15, -8.1975460502553131E-13,
    -8.3989344682125129, -3.00577391366322E-15, -8.1986071309555248E-13,
    -8.3947724039202445, -2.7912604810122602E-15, -8.1996682116557365E-13,
    -8.3906103396279761, -2.57685225357594E-15, -8.2007287719632558E-13,
    -8.3864503165639075, -2.3623388209249795E-15, -8.2017898526634675E-13,
    -8.3822882522716373, -2.1478253882740197E-15, -8.2028509333636792E-13,
    -8.37812618797937, -1.93331195562306E-15, -8.2039120140638908E-13,
    -8.3739641236871023, -1.7187985229721004E-15, -8.2049730947641025E-13,
    -8.3698020593948339, -1.5042850903211397E-15, -8.2060341754643132E-13,
    -8.3656399951025673, -1.28977165767018E-15, -8.2070952561645249E-13,
    -8.3614779308102989, -1.07525822501922E-15, -8.2081563368647355E-13,
    -8.3573158665180323, -8.6084999758289989E-16, -8.2092168971722549E-13,
    -8.35315584345396, -6.4633656493193984E-16, -8.2102779778724665E-13,
    -8.3489937791616935, -4.3181261175951612E-16, -8.2113391106119474E-13,
    -8.3448315107466069, -2.1730969963002012E-16, -8.21240013927289E-13,
    -8.3406696505771585, -2.8173080219881823E-18, -8.2134611158945621E-13,
    -8.33650799453053, 2.1168560410750791E-16, -8.2145221445555045E-13,
    -8.3323461343610816, 4.2618851623700386E-16, -8.215583173216446E-13,
    -208.32818427419164 };

  
   double dv[720] = { 1.0, 0.0, 0.0, 0.016, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.016, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.016, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0,
    0.0, 0.032, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.032, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.032, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.048, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.048, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.048, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0,
    0.064, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.064, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.064, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.08, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.08,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.08, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.096, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.096, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.096, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 1.0, 0.0, 0.0, 0.112, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.112, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.112, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.128, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.128, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.128, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 1.0, 0.0, 0.0, 0.144, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.144, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.144, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.16, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.16, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.16, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0,
    0.0, 0.0, 0.176, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.176, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.176, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.192, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.192, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.192, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0,
    0.0, 0.208, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.208, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.208, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.224, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.224, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.224, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0,
    0.24, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.24, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.24,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.256, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.256,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.256, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.272, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.272, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.272, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 1.0, 0.0, 0.0, 0.288, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.288, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.288, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.304, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.304, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.304, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 1.0, 0.0, 0.0, 0.32, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.32, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.32, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 };

  
   double B[9] = {0.6, 0.0, 0.0, 0.0, 0.6, 0.0, 0.0, 0.0, 0.6};
  /* double B[9] = { 0.68, 0.0, 0.0, 0.0, 0.332, 0.0, 0.0, 0.0, 0.394
  }; */

   signed char A[400] = { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 };

  
  double IGA[7200];
  double A_data[3600];
  double a[3600];
  double b_del_lam[240];
  double result[240];
  double varargin_1_data[240];
  double Bineq[120];
  double GAMMA[120];
  double del_lam[120];
  double esig[120];
  double ftol[120];
  double igr2[120];
  double ilam[120];
  double lam[120];
  double mesil[120];
  double H[60];
  double del_z[60];
  double u_max[60];
  double absxk;
  double mu;
  double scale;
  double t;
  int a_tmp;
  int b_i;
  int b_i1;
  int exitflag;
  int i;
  int i1;
  int idx;
  int iter;
  int iy;
  int j2;
  short b_Aineq[7200];
  signed char Aineq[7200];
  signed char At[7200];
  unsigned char ii_data[240];
  boolean_T x[240];

  //  Set Constants - Extract variables from Structure MPCParams
  //  MPCParams=load('MPCParams.mat');
  // MPCParams.X;
  // MPCParams.A;
  // MPCParams.B;
  // MPCParams.P1;
  // MPCParams.PSI;
  // MPCParams.OMEGA;
  // MPCParams.PHI;
  // MPCParams.L1;
  // MPCParams.L2;
  // MPCParams.L3;
  // MPCParams.L4;
  // MPCParams.Fmaxx;
  // MPCParams.Fmaxy;
  // MPCParams.Fmaxz;
  //  MPCParams.ConvThresh;
  // MPCParams.rCollAvoid;
  // MPCParams.maxiter;
  // MPCParams.tol;
  // MPCParams.ObstAvoid;
  // MPCParams.aObst; % Obstacle semi-major axis
  // MPCParams.bObst; % Obstacle semi-minor axis
  // MPCParams.alphaObst; % Orientation of ellipse
  // P_debris = MPCParams.P_debris;
  // MPCParams.pos_of_debris;
  //  MPCParams.AppCone;
  //  MPCParams.phi;
  //  sampT = MPCParams.Ts;
  //  Set Decision Variables - Configure constraints
  //  Compute distance from target
  scale = 3.3121686421112381E-170;
  absxk = std::abs(x0[0]);
  if (absxk > 3.3121686421112381E-170) {
    dr = 1.0;
    scale = absxk;
  } else {
    t = absxk / 3.3121686421112381E-170;
    dr = t * t;
  }

  absxk = std::abs(x0[1]);
  if (absxk > scale) {
    t = scale / absxk;
    dr = dr * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    dr += t * t;
  }

  dr = scale * std::sqrt(dr);

  //  position tolerance to stop control action
  //  Define vectors for Obstacle Avoidance
  //  to current position
  //  to obstacle
  //  to final position
  //  from obstacle to target
  //  from obstacle to FSS
  //  Angle between unit vectors
  //  Approach Cone variables
  //  point where cone is placed
  //  orientation of cone (inertial)
  //  Unit vector of approach axis
  //  Angle between FSS and approach axis
  dock_flag = 0.0;
  CollAvoid_flag = 0.0;

  //  Set obstacle avoidance flag
  //  Define target point for current iteration
  //      xfinal = xfinal;
  //  to final position
  for (i = 0; i < 6; i++) {
    target_state[i] = 0.0;
  }

  //  Define and Solve QP
  //  dock_complete = 0;
  pt_sel = 0.0;
  dock_complete = 0.0;

  //  Compute QP Matrices
  // Dimension of GAMMA is n*horizon x n
  // for i = 1:horizon
  for (b_i = 0; b_i < 120; b_i++) {
    t = 0.0;
    for (i1 = 0; i1 < 6; i1++) {
      t += 2.0 * x0[i1] * dv[i1 + 6 * b_i];
    }

    GAMMA[b_i] = t;
  }

  for (b_i = 0; b_i < 120; b_i++) {
    t = 0.0;
    for (i1 = 0; i1 < 120; i1++) {
      t += GAMMA[i1] * b[i1 + 120 * b_i];
    }

    Bineq[b_i] = t;
  }

  for (b_i = 0; b_i < 60; b_i++) {
    t = 0.0;
    for (i1 = 0; i1 < 120; i1++) {
      t += Bineq[i1] * b_b[i1 + 120 * b_i];
    }

    H[b_i] = t;
  }

  //  Control Consraints
  i = -1;
  for (idx = 0; idx < 20; idx++) {
    for (j2 = 0; j2 < 3; j2++) {
      t = B[3 * j2];
      scale = B[3 * j2 + 1];
      absxk = B[3 * j2 + 2];
      for (b_i1 = 0; b_i1 < 20; b_i1++) {
        i++;
        a_tmp = A[b_i1 + 20 * idx];
        a[i] = static_cast<double>(a_tmp) * t;
        i++;
        a[i] = static_cast<double>(a_tmp) * scale;
        i++;
        a[i] = static_cast<double>(a_tmp) * absxk;
      }
    }
  }

  for (b_i = 0; b_i < 60; b_i++) {
    t = 0.0;
    for (i1 = 0; i1 < 60; i1++) {
      t += a[b_i + 60 * i1];
    }

    u_max[b_i] = t;
  }

  std::memset(&Aineq[0], 0, 7200U * sizeof(signed char));
  std::memset(&a[0], 0, 3600U * sizeof(double));
  for (i = 0; i < 60; i++) {
    a[i + 60 * i] = 1.0;
  }

  for (b_i = 0; b_i < 60; b_i++) {
    for (i1 = 0; i1 < 60; i1++) {
      idx = static_cast<int>(a[i1 + 60 * b_i]);
      i = i1 + 120 * b_i;
      Aineq[i] = static_cast<signed char>(-idx);
      Aineq[i + 60] = static_cast<signed char>(idx);
    }
  }

  std::memset(&Bineq[0], 0, 120U * sizeof(double));

  //  Control Constraints only
  //      Aineq = Aineq;
  //      Bineq = Bineq; % Joh and me comment this on spet 23
  //  Control Constraints only
  //      Aineq = Aineq;
  //      Bineq = Bineq;
  //  Call QP Solver
  for (i = 0; i < 60; i++) {
    t = u_max[i];
    Bineq[i] = t;
    Bineq[i + 60] = t;
  }
   /* if (X_QP[0].empty()){
   for (i = 0; i < 60; i++){
      X_QP[i] = 0.0;

    }
  } */
  //  Solve quadratic programming problem using Wright's (1997) Method
  //  Minimise J(x) = 1/2x'Hx + f'x
  //  Subject to: Ax <= b
  //  Supporting Functions
  //  Reference: S. J. Wright, "Applying New Optimization Algorithms to Model
  //  Predictive Control," in Chemical Process Control-V, CACHE, AIChE
  //  Symposium, 1997, pp. 147-155.
  // Number of decision variables
  //  p = 0;
  // Test for Cold Start
  // Warm Start
  // to tune
  // to tune
  // Default Values
  mu = 10000.0;
  for (i = 0; i < 120; i++) {
    lam[i] = 100.0;
    ftol[i] = 100.0;
    esig[i] = 0.001;
    for (b_i = 0; b_i < 60; b_i++) {
      At[b_i + 60 * i] = Aineq[i + 120 * b_i];
    }
  }

  //  %Linsolve options
  //  opU.UT = true;
  //  opUT.UT = true;
  //  opUT.TRANSA = true;
  // Begin Searching
  //  for iter = 1:maxiter
  iter = 0;
  exitflag = 0;
  while ((iter <= 100) && (exitflag != 1)) {
    boolean_T exitg1;
    boolean_T y;

    // Create common matrices
    for (i = 0; i < 120; i++) {
      t = lam[i];
      scale = 1.0 / t;
      ilam[i] = scale;
      GAMMA[i] = -t / ftol[i];
      mesil[i] = mu * esig[i] * scale;
    }

    // RHS
    for (i = 0; i < 60; i++) {
      for (b_i = 0; b_i < 120; b_i++) {
        idx = b_i + 120 * i;
        IGA[idx] = GAMMA[b_i] * static_cast<double>(Aineq[idx]);
      }

      t = 0.0;
      for (b_i = 0; b_i < 60; b_i++) {
        t += b_a[i + 60 * b_i] * X_QP[b_i];
      }

      scale = 0.0;
      for (b_i = 0; b_i < 120; b_i++) {
        scale += static_cast<double>(At[i + 60 * b_i]) * lam[b_i];
      }

      u_max[i] = (t - scale) - H[i];
    }

    for (b_i = 0; b_i < 7200; b_i++) {
      b_Aineq[b_i] = static_cast<short>(-Aineq[b_i]);
    }

    for (b_i = 0; b_i < 120; b_i++) {
      t = 0.0;
      for (i1 = 0; i1 < 60; i1++) {
        t += static_cast<double>(b_Aineq[b_i + 120 * i1]) * X_QP[i1];
      }

      igr2[b_i] = GAMMA[b_i] * ((t + Bineq[b_i]) - mesil[b_i]);
    }

    // Solve
    for (b_i = 0; b_i < 60; b_i++) {
      for (i1 = 0; i1 < 60; i1++) {
        t = 0.0;
        for (idx = 0; idx < 120; idx++) {
          t += static_cast<double>(At[b_i + 60 * idx]) * IGA[idx + 120 * i1];
        }

        a[b_i + 60 * i1] = t;
      }
    }

    for (b_i = 0; b_i < 3600; b_i++) {
      A_data[b_i] = b_H[b_i] - a[b_i];
    }

    a_tmp = -1;
    idx = 0;
    exitg1 = false;
    while ((!exitg1) && (idx < 60)) {
      int idxA1j;
      int idxAjj;
      int ix;
      idxA1j = idx * 60;
      idxAjj = idxA1j + idx;
      scale = 0.0;
      if (idx >= 1) {
        ix = idxA1j;
        iy = idxA1j;
        for (i = 0; i < idx; i++) {
          scale += A_data[ix] * A_data[iy];
          ix++;
          iy++;
        }
      }

      scale = A_data[idxAjj] - scale;
      if (scale > 0.0) {
        scale = std::sqrt(scale);
        A_data[idxAjj] = scale;
        if (idx + 1 < 60) {
          int idxAjjp1;
          i = idxA1j + 61;
          idxAjjp1 = idxAjj + 61;
          if (idx != 0) {
            iy = idxAjj + 60;
            b_i = (idxA1j + 60 * (58 - idx)) + 61;
            for (j2 = i; j2 <= b_i; j2 += 60) {
              ix = idxA1j;
              absxk = 0.0;
              i1 = (j2 + idx) - 1;
              for (b_i1 = j2; b_i1 <= i1; b_i1++) {
                absxk += A_data[b_i1 - 1] * A_data[ix];
                ix++;
              }

              A_data[iy] += -absxk;
              iy += 60;
            }
          }

          scale = 1.0 / scale;
          b_i = (idxAjj + 60 * (58 - idx)) + 61;
          for (i = idxAjjp1; i <= b_i; i += 60) {
            A_data[i - 1] *= scale;
          }
        }

        idx++;
      } else {
        a_tmp = idx;
        exitg1 = true;
      }
    }

    // [R] = chol(H-At*IGA);
    if (a_tmp + 1 == 0) {
      for (b_i = 0; b_i < 60; b_i++) {
        t = 0.0;
        for (i1 = 0; i1 < 120; i1++) {
          t += static_cast<double>(At[b_i + 60 * i1]) * igr2[i1];
        }

        del_z[b_i] = u_max[b_i] - t;
      }

      for (b_i = 0; b_i < 3600; b_i++) {
        a[b_i] = b_H[b_i] - a[b_i];
      }

      mldivide(a, del_z);

      // old method (LU?)
      //   del_z = linsolve (R, linsolve (R, (r1-At*igr2), opUT), opU); %exploit matrix properties for solving 
    } else {
      // Not Positive Definite (problem? eg infeasible)
      for (b_i = 0; b_i < 60; b_i++) {
        t = 0.0;
        for (i1 = 0; i1 < 120; i1++) {
          t += static_cast<double>(At[b_i + 60 * i1]) * igr2[i1];
        }

        del_z[b_i] = u_max[b_i] - t;
      }

      for (b_i = 0; b_i < 3600; b_i++) {
        a[b_i] = b_H[b_i] - a[b_i];
      }

      mldivide(a, del_z);

      // old method (LU?)
    }

    // Decide on suitable alpha (from Wright's paper)
    // Try Max Increment (alpha = 1)
    // Check lam and ftol > 0
    for (i = 0; i < 120; i++) {
      t = 0.0;
      for (b_i = 0; b_i < 60; b_i++) {
        t += IGA[i + 120 * b_i] * del_z[b_i];
      }

      t = igr2[i] - t;
      del_lam[i] = t;
      scale = ftol[i];
      absxk = (-scale + mesil[i]) - ilam[i] * scale * t;
      mesil[i] = absxk;
      t += lam[i];
      GAMMA[i] = t;
      scale += absxk;
      ilam[i] = scale;
      x[i] = (t < 2.2204460492503131E-16);
      x[i + 120] = (scale < 2.2204460492503131E-16);
    }

    y = false;
    i = 0;
    exitg1 = false;
    while ((!exitg1) && (i < 240)) {
      if (!x[i]) {
        i++;
      } else {
        y = true;
        exitg1 = true;
      }
    }

    if (!y) {
      // KKT met
      std::memcpy(&lam[0], &GAMMA[0], 120U * sizeof(double));
      std::memcpy(&ftol[0], &ilam[0], 120U * sizeof(double));
      for (b_i = 0; b_i < 60; b_i++) {
        X_QP[b_i] += del_z[b_i];
      }
    } else {
      // KKT failed - solve by finding minimum ratio
      for (b_i = 0; b_i < 120; b_i++) {
        result[b_i] = GAMMA[b_i];
        result[b_i + 120] = ilam[b_i];
      }

      idx = 0;
      i = 0;
      exitg1 = false;
      while ((!exitg1) && (i < 240)) {
        if (result[i] < 2.2204460492503131E-16) {
          idx++;
          ii_data[idx - 1] = static_cast<unsigned char>(i + 1);
          if (idx >= 240) {
            exitg1 = true;
          } else {
            i++;
          }
        } else {
          i++;
        }
      }

      if (1 > idx) {
        j2 = 0;
      } else {
        j2 = idx;
      }

      // detects elements breaking KKT condition
      for (b_i = 0; b_i < 120; b_i++) {
        b_del_lam[b_i] = del_lam[b_i];
        b_del_lam[b_i + 120] = mesil[b_i];
      }

      for (b_i = 0; b_i < j2; b_i++) {
        i = ii_data[b_i] - 1;
        varargin_1_data[b_i] = 1.0 - result[i] / b_del_lam[i];
      }

      if (j2 <= 2) {
        if (j2 == 1) {
          scale = varargin_1_data[0];
        } else if ((varargin_1_data[0] > varargin_1_data[1]) || (
                    (varargin_1_data[0]) && (!rtIsNaN(varargin_1_data[1])))) {
          scale = varargin_1_data[1];
        } else {
          scale = varargin_1_data[0];
        }
      } else {
        if (!rtIsNaN(varargin_1_data[0])) {
          idx = 1;
        } else {
          idx = 0;
          i = 2;
          exitg1 = false;
          while ((!exitg1) && (i <= j2)) {
            if (!rtIsNaN(varargin_1_data[i - 1])) {
              idx = i;
              exitg1 = true;
            } else {
              i++;
            }
          }
        }

        if (idx == 0) {
          scale = varargin_1_data[0];
        } else {
          scale = varargin_1_data[idx - 1];
          b_i = idx + 1;
          for (i = b_i; i <= j2; i++) {
            t = varargin_1_data[i - 1];
            if (scale > t) {
              scale = t;
            }
          }
        }
      }

      scale *= 0.995;

      // solves for min ratio (max value of alpha allowed)
      // Increment
      for (b_i = 0; b_i < 120; b_i++) {
        lam[b_i] += scale * del_lam[b_i];
        ftol[b_i] += scale * mesil[b_i];
      }

      for (b_i = 0; b_i < 60; b_i++) {
        X_QP[b_i] += scale * del_z[b_i];
      }
    }

    // Complimentary Gap
    absxk = mu;
    scale = 0.0;
    for (b_i = 0; b_i < 120; b_i++) {
      scale += ftol[b_i] * lam[b_i];
    }

    mu = scale / 120.0;

    //      if(mu < tol)
    //          exitflag = 1;
    //          return
    //      end
    //      %Solve for new Sigma
    //      sigma = mu/mu_old;
    //      if(sigma > 0.1) %to tune
    //          sigma = 0.1;
    //      end
    //      esig = sigma*ones(mc,1);
    if (mu < 0.001) {
      exitflag = 1;
    } else {
      // Solve for new Sigma
      scale = mu / absxk;
      if (scale > 0.1) {
        // to tune
        scale = 0.1;
      }

      for (i = 0; i < 120; i++) {
        esig[i] = scale;
      }
    }

    iter++;
  }

  // Check for failure
  num_iter = iter;

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

