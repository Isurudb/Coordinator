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
#include <ff_common/ff_names.h>
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
#include <coordinator/Prediction.h>
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





static std::string TOPIC_ASAP_STATUS = "/asap/status";
static std::string TOPIC_ASAP_TEST_NUMBER = "/asap/test_number";
static std::string TOPIC_GNC_CTL_CMD = "/gnc/ctl/command";
static std::string TOPIC_GNC_EKF_ = "/gnc/ekf";
static std::string SERVICE_GNC_CTL_ENABLE_ = "/gnc/ctl/enable" ;

static std::string TOPIC_ASAP_STATUS_s = "/bumble/asap/status";
static std::string TOPIC_ASAP_TEST_NUMBER_s = "/bumble/asap/test_number";
static std::string TOPIC_GNC_CTL_CMD_s = "/bumble/gnc/ctl/command";



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

  ros::Publisher VL_status;

  ros::Subscriber sub_flight_mode_;
  ros::Subscriber sub_ekf_;
  ros::Subscriber sub_ekf_VF;  // virtual follower >> dds topic
  ros::Subscriber sub_ekf_VL;  // virtual leader >> dds topic
  ros::Subscriber sub_test_number_;
  ros::Subscriber sub_VL_status;

  ros::ServiceClient serv_ctl_enable_;

  ros::Timer status_timer_;
  ros::Timer ctl_disable_timer_;

  ff_msgs::FlightMode flight_mode_;
  ff_msgs::EkfState ekf_state_;
  ff_msgs::FamCommand gnc_setpoint;

  coordinator::Prediction mpc_pred;

  geometry_msgs::Wrench ctl_input;
  geometry_msgs::Quaternion attitude;
  geometry_msgs::Vector3 omega,velocity_, position_, position_error, position_ref,velocity;
  geometry_msgs::Vector3 pos_ref2, vel_ref_2,position_error_2, vl_pose;

  tf2::Quaternion attitude_,q_ref,q_e,q_ref_inv;

  // Parameters
  bool ground_ = false;  // whether or not this is a ground test
  bool sim_ = false;
  std::string flight_mode_name_;

  // Stored status parameters
  std::string stored_control_mode_ = "track";  // stored control_mode, set by parameter inputs
  std::string robot ;
  std::string Estimate_status = "Best";
  bool publishflag = false;

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
  void VF_callback(const ff_msgs::EkfState::ConstPtr msg);
  void VL_callback(const ff_msgs::EkfState::ConstPtr msg);
  //void VL_callback(const coordinator::Prediction::ConstPtr  msg);

  void debug();

  // Astrobee GNC interface
  void disable_default_ctl_callback(const ros::TimerEvent&);
  void disable_default_ctl();
  void enable_default_ctl();
  

  // Virtual test list: to be replaced on each derived coordinator
  virtual void RunTest0(ros::NodeHandle *nh) {};
  virtual void RunTest1(ros::NodeHandle *nh) {};
  virtual void RunTest2(ros::NodeHandle *nh) {}; // <--X
  virtual void RunTest3(ros::NodeHandle *nh) {}; // <--X
  virtual void RunTest4(ros::NodeHandle *nh) {}; // <--X
  virtual void RunTest5(ros::NodeHandle *nh) {}; // <--X
  virtual void RunTest6(ros::NodeHandle *nh) {}; // <--X

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

 Eigen::Vector3d x0_;
 Eigen::Vector4d a0_;

double x0[6];
double x0_vl[6];
double x_pred[120]={0};
double Fx;
double Fy;
double Fz;
double target_state[6];
double dock_flag;
double CollAvoid_flag;
double dock_complete;
double num_iter;
double X_QP[60]={0};
double X_QP_vl[60]={0};
double pt_sel;
double dr;
double inti_e_x=0;
double inti_e_y=0;
double inti_e_z=0;
float q0_x = 0;
float q0_y = 0;
float q0_z = 0;

double L=0.5,L0=0.5;

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
void mldivide(double A[3600], double B[60]);
bool rtIsNaN(double value);
void nominal_dynamics();
void nominal_dynamics_good();
void nominal_dynamics_worst();
void tubing_mpc();
void tubing_mpc_good();
void tubing_mpc_worst();
bool coupled=false;
void MPC();
void MPC_rev_1();

double F_r[3];
double F_l[3];
void LSS_primary(double u_primary[3], double Tau[3]);
void LSS_secondary(double u_secondary[3], double Tau[3]);

void MPC_vl();
void vl_pred();


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
     /*  else if(base_status_.test_number  == 2) {
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
      } */
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

  if (qx != 0 || qy != 0 || qz != 0 || qw != 0)
   {
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
    if ((prod_val==0))
    {
      deno=1;
    }
    if(prod_val/deno <-0.9 )
    {
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
    double r=3.14159265, p=0, y=0  ;// for wanna bee upside down
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

  if(initialzation)
  {
        // For the primary
        position_error.x = position_.x - position_ref.x;
        position_error.y = position_.y - position_ref.y;
        position_error.z = position_.z - position_ref.z;

        /* velocity_.x=vx - velocity.x;
        velocity_.y=vy - velocity.y;
        velocity_.z=vz - velocity.z; */

        // for the secondary
        position_error_2.x = position_.x - pos_ref2.x + L;// neg x direction
        position_error_2.y = position_.y - pos_ref2.y ;//0.5; // position off set
        position_error_2.z = position_.z - pos_ref2.z;

        velocity_.x=vx - vel_ref_2.x;
        velocity_.y=vy - vel_ref_2.y;
        velocity_.z=vz - vel_ref_2.z;

      

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
        
        step_PID_worst();


        // MPC Controller inbound <<<<<<<<<<<<<<<<<<<<<<<ID
        
          if(robot=="Primary")
         {  
            x0[0]=position_error.x;
            x0[1]=position_error.y;
            x0[2]=position_error.z;
            x0[3]=vx;
            x0[4]=vy;
            x0[5]=vz;

            // not using
            x0_vl[0]=position_.x;
            x0_vl[1]=position_.y  -0.8;
            x0_vl[2]=position_.z;
            x0_vl[3]=velocity_.x;
            x0_vl[4]=velocity_.y;
            x0_vl[5]=velocity_.z;
            if(publishflag)
            {
                vl_pred();
                
                mpc_pred.x0.x = x0_vl[0];         mpc_pred.x0.y = x0_vl[1];       mpc_pred.x0.z = x0_vl[2];
                mpc_pred.v0.x = x0_vl[3];         mpc_pred.v0.y = x0_vl[4];       mpc_pred.v0.z = x0_vl[5];


                mpc_pred.x1.x = x_pred[0];         mpc_pred.x1.y = x_pred[1];       mpc_pred.x1.z = x_pred[2];
                mpc_pred.v1.x = x_pred[3];         mpc_pred.v1.y = x_pred[4];       mpc_pred.v1.z = x_pred[5];

                mpc_pred.x2.x = x_pred[6];         mpc_pred.x2.y = x_pred[7];       mpc_pred.x2.z = x_pred[8];
                mpc_pred.v2.x = x_pred[9];         mpc_pred.v2.y = x_pred[10];       mpc_pred.v2.z = x_pred[11];

                mpc_pred.x3.x = x_pred[12];         mpc_pred.x3.y = x_pred[13];       mpc_pred.x3.z = x_pred[14];        
                mpc_pred.v3.x = x_pred[15];         mpc_pred.v3.y = x_pred[16];       mpc_pred.v3.z = x_pred[17];

                mpc_pred.x4.x = x_pred[18];         mpc_pred.x4.y = x_pred[19];       mpc_pred.x4.z = x_pred[20];
                mpc_pred.v4.x = x_pred[21];         mpc_pred.v4.y = x_pred[22];       mpc_pred.v4.z = x_pred[23];

                mpc_pred.x5.x = x_pred[24];         mpc_pred.x5.y = x_pred[25];       mpc_pred.x5.z = x_pred[26];
                mpc_pred.v5.x = x_pred[27];         mpc_pred.v5.y = x_pred[28];       mpc_pred.v5.z = x_pred[29];

                mpc_pred.x6.x = x_pred[30];         mpc_pred.x6.y = x_pred[31];       mpc_pred.x6.z = x_pred[32];
                mpc_pred.v6.x = x_pred[33];         mpc_pred.v6.y = x_pred[34];       mpc_pred.v6.z = x_pred[35];

            }

            MPC_rev_1();

          }
          

         else
          {
            /* x0[0]=position_.x - pos_ref2.x;
            x0[1]=position_.y - pos_ref2.y;
            x0[2]=position_.z - pos_ref2.z;
            x0[3]=vx - vel_ref_2.x;
            x0[4]=vy - vel_ref_2.y;
            x0[5]=vz - vel_ref_2.z; */
            x0[0]=position_error_2.x;
            x0[1]=position_error_2.y;
            x0[2]=position_error_2.z;
            x0[3]=velocity_.x;
            x0[4]=velocity_.y;
            x0[5]=velocity_.z;

            MPC_rev_1();
          
          }
          
       /*  count+=1;
        
        if (count==5)
        {

          count =0;
        }  */
        if (coupled && robot=="Primary")
        {
          double u_primary[3];
          u_primary[0]=X_QP[0];
          u_primary[1]=X_QP[1];
          u_primary[2]=X_QP[2];

          double Tau[3];
          Tau[0]=arg_tau_x;
          Tau[1]=arg_tau_y;
          Tau[2]=arg_tau_z;
          LSS_primary(u_primary,Tau);
          X_QP[0]=F_r[0];
          X_QP[1]=F_r[1];
          X_QP[2]=F_r[2];

        }
        else if (coupled && robot=="Secondary")
        {
          double u_primary[3];
          u_primary[0]=X_QP[0];
          u_primary[1]=X_QP[1];
          u_primary[2]=X_QP[2];

          double Tau[3];
          Tau[0]=arg_tau_x;
          Tau[1]=arg_tau_y;
          Tau[2]=arg_tau_z;
          LSS_secondary(u_primary,Tau);
          X_QP[0]=F_l[0];
          X_QP[1]=F_l[1];
          X_QP[2]=F_l[2];
        }
        


  }
}


/* ************************************************************************** */
template<typename T>
void CoordinatorBase<T>::VL_callback(const ff_msgs::EkfState::ConstPtr msg) {
  /**
   * @brief The `gnc/ekf` subscriber callback. from the leader Called at 62.5 Hz.
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

  if (qx != 0 || qy != 0 || qz != 0 || qw != 0)
   {
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
  

  //pos_ref2, vel_ref_2;
  pos_ref2.x = px;
  pos_ref2.y = py;
  pos_ref2.z = pz;
  vel_ref_2.x = vx;
  vel_ref_2.y = vy;
  vel_ref_2.z = vz;

}

/* ************************************************************************** */
template<typename T>
void CoordinatorBase<T>::VF_callback(const ff_msgs::EkfState::ConstPtr msg) {
  /**
   * @brief The `gnc/ekf` subscriber callback from the follower. Called at 62.5 Hz.
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

  if (qx != 0 || qy != 0 || qz != 0 || qw != 0)
   {
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

  pos_ref2.x = px;
  pos_ref2.y = py;
  pos_ref2.z = pz;
  vel_ref_2.x = vx;
  vel_ref_2.y = vy;
  vel_ref_2.z = vz;

}


/* template<typename T>
void CoordinatorBase<T>::VL_callback(const coordinator::Prediction::ConstPtr msg){
  //pos_ref2, vel_ref_2;
  pos_ref2.x = msg->x0.x;
  pos_ref2.y = msg->x0.y;
  pos_ref2.z = msg->x0.z;
  vel_ref_2.x = msg->v0.x;
  vel_ref_2.y = msg->v0.y;
  vel_ref_2.z = msg->v0.z;

  //ROS_INFO("fx: [%f]  fy: [%f] fz: [%f] tau_x: [%f] tau_y: [%f] tau_y: [%f]", pos_ref2.x,pos_ref2.y,pos_ref2.z, vel_ref_2.x,vel_ref_2.y,vel_ref_2.z);

} */
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
    //ROS_INFO(" Quaternion sign change detected >>>>>>>>>>>>>>>>>>>>> ");
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
void CoordinatorBase<T>::vl_pred()
{
  static const double b_a[7200] = { 0.125, 0.0, 0.0, 0.5, 0.0, 0.0, 0.375, 0.0,
    0.0, 0.5, 0.0, 0.0, 0.625, 0.0, 0.0, 0.5, 0.0, 0.0, 0.875, 0.0, 0.0, 0.5,
    0.0, 0.0, 1.125, 0.0, 0.0, 0.5, 0.0, 0.0, 1.375, 0.0, 0.0, 0.5, 0.0, 0.0,
    1.625, 0.0, 0.0, 0.5, 0.0, 0.0, 1.875, 0.0, 0.0, 0.5, 0.0, 0.0, 2.125, 0.0,
    0.0, 0.5, 0.0, 0.0, 2.375, 0.0, 0.0, 0.5, 0.0, 0.0, 2.625, 0.0, 0.0, 0.5,
    0.0, 0.0, 2.875, 0.0, 0.0, 0.5, 0.0, 0.0, 3.125, 0.0, 0.0, 0.5, 0.0, 0.0,
    3.375, 0.0, 0.0, 0.5, 0.0, 0.0, 3.625, 0.0, 0.0, 0.5, 0.0, 0.0, 3.875, 0.0,
    0.0, 0.5, 0.0, 0.0, 4.125, 0.0, 0.0, 0.5, 0.0, 0.0, 4.375, 0.0, 0.0, 0.5,
    0.0, 0.0, 4.625, 0.0, 0.0, 0.5, 0.0, 0.0, 4.875, 0.0, 0.0, 0.5, 0.0, 0.0,
    0.0, 0.125, 0.0, 0.0, 0.5, 0.0, 0.0, 0.375, 0.0, 0.0, 0.5, 0.0, 0.0, 0.625,
    0.0, 0.0, 0.5, 0.0, 0.0, 0.875, 0.0, 0.0, 0.5, 0.0, 0.0, 1.125, 0.0, 0.0,
    0.5, 0.0, 0.0, 1.375, 0.0, 0.0, 0.5, 0.0, 0.0, 1.625, 0.0, 0.0, 0.5, 0.0,
    0.0, 1.875, 0.0, 0.0, 0.5, 0.0, 0.0, 2.125, 0.0, 0.0, 0.5, 0.0, 0.0, 2.375,
    0.0, 0.0, 0.5, 0.0, 0.0, 2.625, 0.0, 0.0, 0.5, 0.0, 0.0, 2.875, 0.0, 0.0,
    0.5, 0.0, 0.0, 3.125, 0.0, 0.0, 0.5, 0.0, 0.0, 3.375, 0.0, 0.0, 0.5, 0.0,
    0.0, 3.625, 0.0, 0.0, 0.5, 0.0, 0.0, 3.875, 0.0, 0.0, 0.5, 0.0, 0.0, 4.125,
    0.0, 0.0, 0.5, 0.0, 0.0, 4.375, 0.0, 0.0, 0.5, 0.0, 0.0, 4.625, 0.0, 0.0,
    0.5, 0.0, 0.0, 4.875, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.125, 0.0, 0.0, 0.5,
    0.0, 0.0, 0.375, 0.0, 0.0, 0.5, 0.0, 0.0, 0.625, 0.0, 0.0, 0.5, 0.0, 0.0,
    0.875, 0.0, 0.0, 0.5, 0.0, 0.0, 1.125, 0.0, 0.0, 0.5, 0.0, 0.0, 1.375, 0.0,
    0.0, 0.5, 0.0, 0.0, 1.625, 0.0, 0.0, 0.5, 0.0, 0.0, 1.875, 0.0, 0.0, 0.5,
    0.0, 0.0, 2.125, 0.0, 0.0, 0.5, 0.0, 0.0, 2.375, 0.0, 0.0, 0.5, 0.0, 0.0,
    2.625, 0.0, 0.0, 0.5, 0.0, 0.0, 2.875, 0.0, 0.0, 0.5, 0.0, 0.0, 3.125, 0.0,
    0.0, 0.5, 0.0, 0.0, 3.375, 0.0, 0.0, 0.5, 0.0, 0.0, 3.625, 0.0, 0.0, 0.5,
    0.0, 0.0, 3.875, 0.0, 0.0, 0.5, 0.0, 0.0, 4.125, 0.0, 0.0, 0.5, 0.0, 0.0,
    4.375, 0.0, 0.0, 0.5, 0.0, 0.0, 4.625, 0.0, 0.0, 0.5, 0.0, 0.0, 4.875, 0.0,
    0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.125, 0.0, 0.0, 0.5, 0.0, 0.0,
    0.375, 0.0, 0.0, 0.5, 0.0, 0.0, 0.625, 0.0, 0.0, 0.5, 0.0, 0.0, 0.875, 0.0,
    0.0, 0.5, 0.0, 0.0, 1.125, 0.0, 0.0, 0.5, 0.0, 0.0, 1.375, 0.0, 0.0, 0.5,
    0.0, 0.0, 1.625, 0.0, 0.0, 0.5, 0.0, 0.0, 1.875, 0.0, 0.0, 0.5, 0.0, 0.0,
    2.125, 0.0, 0.0, 0.5, 0.0, 0.0, 2.375, 0.0, 0.0, 0.5, 0.0, 0.0, 2.625, 0.0,
    0.0, 0.5, 0.0, 0.0, 2.875, 0.0, 0.0, 0.5, 0.0, 0.0, 3.125, 0.0, 0.0, 0.5,
    0.0, 0.0, 3.375, 0.0, 0.0, 0.5, 0.0, 0.0, 3.625, 0.0, 0.0, 0.5, 0.0, 0.0,
    3.875, 0.0, 0.0, 0.5, 0.0, 0.0, 4.125, 0.0, 0.0, 0.5, 0.0, 0.0, 4.375, 0.0,
    0.0, 0.5, 0.0, 0.0, 4.625, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.125, 0.0, 0.0, 0.5, 0.0, 0.0, 0.375, 0.0, 0.0, 0.5, 0.0, 0.0,
    0.625, 0.0, 0.0, 0.5, 0.0, 0.0, 0.875, 0.0, 0.0, 0.5, 0.0, 0.0, 1.125, 0.0,
    0.0, 0.5, 0.0, 0.0, 1.375, 0.0, 0.0, 0.5, 0.0, 0.0, 1.625, 0.0, 0.0, 0.5,
    0.0, 0.0, 1.875, 0.0, 0.0, 0.5, 0.0, 0.0, 2.125, 0.0, 0.0, 0.5, 0.0, 0.0,
    2.375, 0.0, 0.0, 0.5, 0.0, 0.0, 2.625, 0.0, 0.0, 0.5, 0.0, 0.0, 2.875, 0.0,
    0.0, 0.5, 0.0, 0.0, 3.125, 0.0, 0.0, 0.5, 0.0, 0.0, 3.375, 0.0, 0.0, 0.5,
    0.0, 0.0, 3.625, 0.0, 0.0, 0.5, 0.0, 0.0, 3.875, 0.0, 0.0, 0.5, 0.0, 0.0,
    4.125, 0.0, 0.0, 0.5, 0.0, 0.0, 4.375, 0.0, 0.0, 0.5, 0.0, 0.0, 4.625, 0.0,
    0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.125, 0.0, 0.0, 0.5,
    0.0, 0.0, 0.375, 0.0, 0.0, 0.5, 0.0, 0.0, 0.625, 0.0, 0.0, 0.5, 0.0, 0.0,
    0.875, 0.0, 0.0, 0.5, 0.0, 0.0, 1.125, 0.0, 0.0, 0.5, 0.0, 0.0, 1.375, 0.0,
    0.0, 0.5, 0.0, 0.0, 1.625, 0.0, 0.0, 0.5, 0.0, 0.0, 1.875, 0.0, 0.0, 0.5,
    0.0, 0.0, 2.125, 0.0, 0.0, 0.5, 0.0, 0.0, 2.375, 0.0, 0.0, 0.5, 0.0, 0.0,
    2.625, 0.0, 0.0, 0.5, 0.0, 0.0, 2.875, 0.0, 0.0, 0.5, 0.0, 0.0, 3.125, 0.0,
    0.0, 0.5, 0.0, 0.0, 3.375, 0.0, 0.0, 0.5, 0.0, 0.0, 3.625, 0.0, 0.0, 0.5,
    0.0, 0.0, 3.875, 0.0, 0.0, 0.5, 0.0, 0.0, 4.125, 0.0, 0.0, 0.5, 0.0, 0.0,
    4.375, 0.0, 0.0, 0.5, 0.0, 0.0, 4.625, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.125, 0.0, 0.0, 0.5, 0.0, 0.0,
    0.375, 0.0, 0.0, 0.5, 0.0, 0.0, 0.625, 0.0, 0.0, 0.5, 0.0, 0.0, 0.875, 0.0,
    0.0, 0.5, 0.0, 0.0, 1.125, 0.0, 0.0, 0.5, 0.0, 0.0, 1.375, 0.0, 0.0, 0.5,
    0.0, 0.0, 1.625, 0.0, 0.0, 0.5, 0.0, 0.0, 1.875, 0.0, 0.0, 0.5, 0.0, 0.0,
    2.125, 0.0, 0.0, 0.5, 0.0, 0.0, 2.375, 0.0, 0.0, 0.5, 0.0, 0.0, 2.625, 0.0,
    0.0, 0.5, 0.0, 0.0, 2.875, 0.0, 0.0, 0.5, 0.0, 0.0, 3.125, 0.0, 0.0, 0.5,
    0.0, 0.0, 3.375, 0.0, 0.0, 0.5, 0.0, 0.0, 3.625, 0.0, 0.0, 0.5, 0.0, 0.0,
    3.875, 0.0, 0.0, 0.5, 0.0, 0.0, 4.125, 0.0, 0.0, 0.5, 0.0, 0.0, 4.375, 0.0,
    0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.125, 0.0, 0.0, 0.5, 0.0, 0.0, 0.375, 0.0, 0.0, 0.5, 0.0, 0.0,
    0.625, 0.0, 0.0, 0.5, 0.0, 0.0, 0.875, 0.0, 0.0, 0.5, 0.0, 0.0, 1.125, 0.0,
    0.0, 0.5, 0.0, 0.0, 1.375, 0.0, 0.0, 0.5, 0.0, 0.0, 1.625, 0.0, 0.0, 0.5,
    0.0, 0.0, 1.875, 0.0, 0.0, 0.5, 0.0, 0.0, 2.125, 0.0, 0.0, 0.5, 0.0, 0.0,
    2.375, 0.0, 0.0, 0.5, 0.0, 0.0, 2.625, 0.0, 0.0, 0.5, 0.0, 0.0, 2.875, 0.0,
    0.0, 0.5, 0.0, 0.0, 3.125, 0.0, 0.0, 0.5, 0.0, 0.0, 3.375, 0.0, 0.0, 0.5,
    0.0, 0.0, 3.625, 0.0, 0.0, 0.5, 0.0, 0.0, 3.875, 0.0, 0.0, 0.5, 0.0, 0.0,
    4.125, 0.0, 0.0, 0.5, 0.0, 0.0, 4.375, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.125, 0.0, 0.0, 0.5,
    0.0, 0.0, 0.375, 0.0, 0.0, 0.5, 0.0, 0.0, 0.625, 0.0, 0.0, 0.5, 0.0, 0.0,
    0.875, 0.0, 0.0, 0.5, 0.0, 0.0, 1.125, 0.0, 0.0, 0.5, 0.0, 0.0, 1.375, 0.0,
    0.0, 0.5, 0.0, 0.0, 1.625, 0.0, 0.0, 0.5, 0.0, 0.0, 1.875, 0.0, 0.0, 0.5,
    0.0, 0.0, 2.125, 0.0, 0.0, 0.5, 0.0, 0.0, 2.375, 0.0, 0.0, 0.5, 0.0, 0.0,
    2.625, 0.0, 0.0, 0.5, 0.0, 0.0, 2.875, 0.0, 0.0, 0.5, 0.0, 0.0, 3.125, 0.0,
    0.0, 0.5, 0.0, 0.0, 3.375, 0.0, 0.0, 0.5, 0.0, 0.0, 3.625, 0.0, 0.0, 0.5,
    0.0, 0.0, 3.875, 0.0, 0.0, 0.5, 0.0, 0.0, 4.125, 0.0, 0.0, 0.5, 0.0, 0.0,
    4.375, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.125, 0.0, 0.0, 0.5, 0.0, 0.0, 0.375,
    0.0, 0.0, 0.5, 0.0, 0.0, 0.625, 0.0, 0.0, 0.5, 0.0, 0.0, 0.875, 0.0, 0.0,
    0.5, 0.0, 0.0, 1.125, 0.0, 0.0, 0.5, 0.0, 0.0, 1.375, 0.0, 0.0, 0.5, 0.0,
    0.0, 1.625, 0.0, 0.0, 0.5, 0.0, 0.0, 1.875, 0.0, 0.0, 0.5, 0.0, 0.0, 2.125,
    0.0, 0.0, 0.5, 0.0, 0.0, 2.375, 0.0, 0.0, 0.5, 0.0, 0.0, 2.625, 0.0, 0.0,
    0.5, 0.0, 0.0, 2.875, 0.0, 0.0, 0.5, 0.0, 0.0, 3.125, 0.0, 0.0, 0.5, 0.0,
    0.0, 3.375, 0.0, 0.0, 0.5, 0.0, 0.0, 3.625, 0.0, 0.0, 0.5, 0.0, 0.0, 3.875,
    0.0, 0.0, 0.5, 0.0, 0.0, 4.125, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.125, 0.0, 0.0, 0.5, 0.0, 0.0, 0.375, 0.0, 0.0, 0.5, 0.0, 0.0, 0.625, 0.0,
    0.0, 0.5, 0.0, 0.0, 0.875, 0.0, 0.0, 0.5, 0.0, 0.0, 1.125, 0.0, 0.0, 0.5,
    0.0, 0.0, 1.375, 0.0, 0.0, 0.5, 0.0, 0.0, 1.625, 0.0, 0.0, 0.5, 0.0, 0.0,
    1.875, 0.0, 0.0, 0.5, 0.0, 0.0, 2.125, 0.0, 0.0, 0.5, 0.0, 0.0, 2.375, 0.0,
    0.0, 0.5, 0.0, 0.0, 2.625, 0.0, 0.0, 0.5, 0.0, 0.0, 2.875, 0.0, 0.0, 0.5,
    0.0, 0.0, 3.125, 0.0, 0.0, 0.5, 0.0, 0.0, 3.375, 0.0, 0.0, 0.5, 0.0, 0.0,
    3.625, 0.0, 0.0, 0.5, 0.0, 0.0, 3.875, 0.0, 0.0, 0.5, 0.0, 0.0, 4.125, 0.0,
    0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.125, 0.0, 0.0, 0.5, 0.0, 0.0,
    0.375, 0.0, 0.0, 0.5, 0.0, 0.0, 0.625, 0.0, 0.0, 0.5, 0.0, 0.0, 0.875, 0.0,
    0.0, 0.5, 0.0, 0.0, 1.125, 0.0, 0.0, 0.5, 0.0, 0.0, 1.375, 0.0, 0.0, 0.5,
    0.0, 0.0, 1.625, 0.0, 0.0, 0.5, 0.0, 0.0, 1.875, 0.0, 0.0, 0.5, 0.0, 0.0,
    2.125, 0.0, 0.0, 0.5, 0.0, 0.0, 2.375, 0.0, 0.0, 0.5, 0.0, 0.0, 2.625, 0.0,
    0.0, 0.5, 0.0, 0.0, 2.875, 0.0, 0.0, 0.5, 0.0, 0.0, 3.125, 0.0, 0.0, 0.5,
    0.0, 0.0, 3.375, 0.0, 0.0, 0.5, 0.0, 0.0, 3.625, 0.0, 0.0, 0.5, 0.0, 0.0,
    3.875, 0.0, 0.0, 0.5, 0.0, 0.0, 4.125, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.125, 0.0, 0.0, 0.5, 0.0, 0.0, 0.375, 0.0, 0.0,
    0.5, 0.0, 0.0, 0.625, 0.0, 0.0, 0.5, 0.0, 0.0, 0.875, 0.0, 0.0, 0.5, 0.0,
    0.0, 1.125, 0.0, 0.0, 0.5, 0.0, 0.0, 1.375, 0.0, 0.0, 0.5, 0.0, 0.0, 1.625,
    0.0, 0.0, 0.5, 0.0, 0.0, 1.875, 0.0, 0.0, 0.5, 0.0, 0.0, 2.125, 0.0, 0.0,
    0.5, 0.0, 0.0, 2.375, 0.0, 0.0, 0.5, 0.0, 0.0, 2.625, 0.0, 0.0, 0.5, 0.0,
    0.0, 2.875, 0.0, 0.0, 0.5, 0.0, 0.0, 3.125, 0.0, 0.0, 0.5, 0.0, 0.0, 3.375,
    0.0, 0.0, 0.5, 0.0, 0.0, 3.625, 0.0, 0.0, 0.5, 0.0, 0.0, 3.875, 0.0, 0.0,
    0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.125, 0.0,
    0.0, 0.5, 0.0, 0.0, 0.375, 0.0, 0.0, 0.5, 0.0, 0.0, 0.625, 0.0, 0.0, 0.5,
    0.0, 0.0, 0.875, 0.0, 0.0, 0.5, 0.0, 0.0, 1.125, 0.0, 0.0, 0.5, 0.0, 0.0,
    1.375, 0.0, 0.0, 0.5, 0.0, 0.0, 1.625, 0.0, 0.0, 0.5, 0.0, 0.0, 1.875, 0.0,
    0.0, 0.5, 0.0, 0.0, 2.125, 0.0, 0.0, 0.5, 0.0, 0.0, 2.375, 0.0, 0.0, 0.5,
    0.0, 0.0, 2.625, 0.0, 0.0, 0.5, 0.0, 0.0, 2.875, 0.0, 0.0, 0.5, 0.0, 0.0,
    3.125, 0.0, 0.0, 0.5, 0.0, 0.0, 3.375, 0.0, 0.0, 0.5, 0.0, 0.0, 3.625, 0.0,
    0.0, 0.5, 0.0, 0.0, 3.875, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.125, 0.0, 0.0, 0.5, 0.0, 0.0, 0.375, 0.0, 0.0,
    0.5, 0.0, 0.0, 0.625, 0.0, 0.0, 0.5, 0.0, 0.0, 0.875, 0.0, 0.0, 0.5, 0.0,
    0.0, 1.125, 0.0, 0.0, 0.5, 0.0, 0.0, 1.375, 0.0, 0.0, 0.5, 0.0, 0.0, 1.625,
    0.0, 0.0, 0.5, 0.0, 0.0, 1.875, 0.0, 0.0, 0.5, 0.0, 0.0, 2.125, 0.0, 0.0,
    0.5, 0.0, 0.0, 2.375, 0.0, 0.0, 0.5, 0.0, 0.0, 2.625, 0.0, 0.0, 0.5, 0.0,
    0.0, 2.875, 0.0, 0.0, 0.5, 0.0, 0.0, 3.125, 0.0, 0.0, 0.5, 0.0, 0.0, 3.375,
    0.0, 0.0, 0.5, 0.0, 0.0, 3.625, 0.0, 0.0, 0.5, 0.0, 0.0, 3.875, 0.0, 0.0,
    0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.125, 0.0, 0.0, 0.5, 0.0, 0.0, 0.375, 0.0, 0.0, 0.5, 0.0, 0.0, 0.625,
    0.0, 0.0, 0.5, 0.0, 0.0, 0.875, 0.0, 0.0, 0.5, 0.0, 0.0, 1.125, 0.0, 0.0,
    0.5, 0.0, 0.0, 1.375, 0.0, 0.0, 0.5, 0.0, 0.0, 1.625, 0.0, 0.0, 0.5, 0.0,
    0.0, 1.875, 0.0, 0.0, 0.5, 0.0, 0.0, 2.125, 0.0, 0.0, 0.5, 0.0, 0.0, 2.375,
    0.0, 0.0, 0.5, 0.0, 0.0, 2.625, 0.0, 0.0, 0.5, 0.0, 0.0, 2.875, 0.0, 0.0,
    0.5, 0.0, 0.0, 3.125, 0.0, 0.0, 0.5, 0.0, 0.0, 3.375, 0.0, 0.0, 0.5, 0.0,
    0.0, 3.625, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.125, 0.0, 0.0, 0.5, 0.0, 0.0,
    0.375, 0.0, 0.0, 0.5, 0.0, 0.0, 0.625, 0.0, 0.0, 0.5, 0.0, 0.0, 0.875, 0.0,
    0.0, 0.5, 0.0, 0.0, 1.125, 0.0, 0.0, 0.5, 0.0, 0.0, 1.375, 0.0, 0.0, 0.5,
    0.0, 0.0, 1.625, 0.0, 0.0, 0.5, 0.0, 0.0, 1.875, 0.0, 0.0, 0.5, 0.0, 0.0,
    2.125, 0.0, 0.0, 0.5, 0.0, 0.0, 2.375, 0.0, 0.0, 0.5, 0.0, 0.0, 2.625, 0.0,
    0.0, 0.5, 0.0, 0.0, 2.875, 0.0, 0.0, 0.5, 0.0, 0.0, 3.125, 0.0, 0.0, 0.5,
    0.0, 0.0, 3.375, 0.0, 0.0, 0.5, 0.0, 0.0, 3.625, 0.0, 0.0, 0.5, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.125, 0.0, 0.0, 0.5, 0.0, 0.0, 0.375, 0.0, 0.0, 0.5, 0.0, 0.0, 0.625,
    0.0, 0.0, 0.5, 0.0, 0.0, 0.875, 0.0, 0.0, 0.5, 0.0, 0.0, 1.125, 0.0, 0.0,
    0.5, 0.0, 0.0, 1.375, 0.0, 0.0, 0.5, 0.0, 0.0, 1.625, 0.0, 0.0, 0.5, 0.0,
    0.0, 1.875, 0.0, 0.0, 0.5, 0.0, 0.0, 2.125, 0.0, 0.0, 0.5, 0.0, 0.0, 2.375,
    0.0, 0.0, 0.5, 0.0, 0.0, 2.625, 0.0, 0.0, 0.5, 0.0, 0.0, 2.875, 0.0, 0.0,
    0.5, 0.0, 0.0, 3.125, 0.0, 0.0, 0.5, 0.0, 0.0, 3.375, 0.0, 0.0, 0.5, 0.0,
    0.0, 3.625, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.125, 0.0, 0.0, 0.5,
    0.0, 0.0, 0.375, 0.0, 0.0, 0.5, 0.0, 0.0, 0.625, 0.0, 0.0, 0.5, 0.0, 0.0,
    0.875, 0.0, 0.0, 0.5, 0.0, 0.0, 1.125, 0.0, 0.0, 0.5, 0.0, 0.0, 1.375, 0.0,
    0.0, 0.5, 0.0, 0.0, 1.625, 0.0, 0.0, 0.5, 0.0, 0.0, 1.875, 0.0, 0.0, 0.5,
    0.0, 0.0, 2.125, 0.0, 0.0, 0.5, 0.0, 0.0, 2.375, 0.0, 0.0, 0.5, 0.0, 0.0,
    2.625, 0.0, 0.0, 0.5, 0.0, 0.0, 2.875, 0.0, 0.0, 0.5, 0.0, 0.0, 3.125, 0.0,
    0.0, 0.5, 0.0, 0.0, 3.375, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.125, 0.0, 0.0, 0.5, 0.0, 0.0, 0.375, 0.0, 0.0, 0.5, 0.0, 0.0,
    0.625, 0.0, 0.0, 0.5, 0.0, 0.0, 0.875, 0.0, 0.0, 0.5, 0.0, 0.0, 1.125, 0.0,
    0.0, 0.5, 0.0, 0.0, 1.375, 0.0, 0.0, 0.5, 0.0, 0.0, 1.625, 0.0, 0.0, 0.5,
    0.0, 0.0, 1.875, 0.0, 0.0, 0.5, 0.0, 0.0, 2.125, 0.0, 0.0, 0.5, 0.0, 0.0,
    2.375, 0.0, 0.0, 0.5, 0.0, 0.0, 2.625, 0.0, 0.0, 0.5, 0.0, 0.0, 2.875, 0.0,
    0.0, 0.5, 0.0, 0.0, 3.125, 0.0, 0.0, 0.5, 0.0, 0.0, 3.375, 0.0, 0.0, 0.5,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.125, 0.0, 0.0, 0.5, 0.0, 0.0,
    0.375, 0.0, 0.0, 0.5, 0.0, 0.0, 0.625, 0.0, 0.0, 0.5, 0.0, 0.0, 0.875, 0.0,
    0.0, 0.5, 0.0, 0.0, 1.125, 0.0, 0.0, 0.5, 0.0, 0.0, 1.375, 0.0, 0.0, 0.5,
    0.0, 0.0, 1.625, 0.0, 0.0, 0.5, 0.0, 0.0, 1.875, 0.0, 0.0, 0.5, 0.0, 0.0,
    2.125, 0.0, 0.0, 0.5, 0.0, 0.0, 2.375, 0.0, 0.0, 0.5, 0.0, 0.0, 2.625, 0.0,
    0.0, 0.5, 0.0, 0.0, 2.875, 0.0, 0.0, 0.5, 0.0, 0.0, 3.125, 0.0, 0.0, 0.5,
    0.0, 0.0, 3.375, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.125, 0.0, 0.0, 0.5, 0.0, 0.0, 0.375, 0.0, 0.0, 0.5, 0.0,
    0.0, 0.625, 0.0, 0.0, 0.5, 0.0, 0.0, 0.875, 0.0, 0.0, 0.5, 0.0, 0.0, 1.125,
    0.0, 0.0, 0.5, 0.0, 0.0, 1.375, 0.0, 0.0, 0.5, 0.0, 0.0, 1.625, 0.0, 0.0,
    0.5, 0.0, 0.0, 1.875, 0.0, 0.0, 0.5, 0.0, 0.0, 2.125, 0.0, 0.0, 0.5, 0.0,
    0.0, 2.375, 0.0, 0.0, 0.5, 0.0, 0.0, 2.625, 0.0, 0.0, 0.5, 0.0, 0.0, 2.875,
    0.0, 0.0, 0.5, 0.0, 0.0, 3.125, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.125, 0.0, 0.0, 0.5, 0.0, 0.0,
    0.375, 0.0, 0.0, 0.5, 0.0, 0.0, 0.625, 0.0, 0.0, 0.5, 0.0, 0.0, 0.875, 0.0,
    0.0, 0.5, 0.0, 0.0, 1.125, 0.0, 0.0, 0.5, 0.0, 0.0, 1.375, 0.0, 0.0, 0.5,
    0.0, 0.0, 1.625, 0.0, 0.0, 0.5, 0.0, 0.0, 1.875, 0.0, 0.0, 0.5, 0.0, 0.0,
    2.125, 0.0, 0.0, 0.5, 0.0, 0.0, 2.375, 0.0, 0.0, 0.5, 0.0, 0.0, 2.625, 0.0,
    0.0, 0.5, 0.0, 0.0, 2.875, 0.0, 0.0, 0.5, 0.0, 0.0, 3.125, 0.0, 0.0, 0.5,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.125, 0.0, 0.0, 0.5, 0.0, 0.0, 0.375, 0.0, 0.0, 0.5, 0.0, 0.0, 0.625, 0.0,
    0.0, 0.5, 0.0, 0.0, 0.875, 0.0, 0.0, 0.5, 0.0, 0.0, 1.125, 0.0, 0.0, 0.5,
    0.0, 0.0, 1.375, 0.0, 0.0, 0.5, 0.0, 0.0, 1.625, 0.0, 0.0, 0.5, 0.0, 0.0,
    1.875, 0.0, 0.0, 0.5, 0.0, 0.0, 2.125, 0.0, 0.0, 0.5, 0.0, 0.0, 2.375, 0.0,
    0.0, 0.5, 0.0, 0.0, 2.625, 0.0, 0.0, 0.5, 0.0, 0.0, 2.875, 0.0, 0.0, 0.5,
    0.0, 0.0, 3.125, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.125, 0.0, 0.0, 0.5, 0.0, 0.0,
    0.375, 0.0, 0.0, 0.5, 0.0, 0.0, 0.625, 0.0, 0.0, 0.5, 0.0, 0.0, 0.875, 0.0,
    0.0, 0.5, 0.0, 0.0, 1.125, 0.0, 0.0, 0.5, 0.0, 0.0, 1.375, 0.0, 0.0, 0.5,
    0.0, 0.0, 1.625, 0.0, 0.0, 0.5, 0.0, 0.0, 1.875, 0.0, 0.0, 0.5, 0.0, 0.0,
    2.125, 0.0, 0.0, 0.5, 0.0, 0.0, 2.375, 0.0, 0.0, 0.5, 0.0, 0.0, 2.625, 0.0,
    0.0, 0.5, 0.0, 0.0, 2.875, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.125,
    0.0, 0.0, 0.5, 0.0, 0.0, 0.375, 0.0, 0.0, 0.5, 0.0, 0.0, 0.625, 0.0, 0.0,
    0.5, 0.0, 0.0, 0.875, 0.0, 0.0, 0.5, 0.0, 0.0, 1.125, 0.0, 0.0, 0.5, 0.0,
    0.0, 1.375, 0.0, 0.0, 0.5, 0.0, 0.0, 1.625, 0.0, 0.0, 0.5, 0.0, 0.0, 1.875,
    0.0, 0.0, 0.5, 0.0, 0.0, 2.125, 0.0, 0.0, 0.5, 0.0, 0.0, 2.375, 0.0, 0.0,
    0.5, 0.0, 0.0, 2.625, 0.0, 0.0, 0.5, 0.0, 0.0, 2.875, 0.0, 0.0, 0.5, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.125, 0.0, 0.0, 0.5, 0.0, 0.0, 0.375, 0.0, 0.0,
    0.5, 0.0, 0.0, 0.625, 0.0, 0.0, 0.5, 0.0, 0.0, 0.875, 0.0, 0.0, 0.5, 0.0,
    0.0, 1.125, 0.0, 0.0, 0.5, 0.0, 0.0, 1.375, 0.0, 0.0, 0.5, 0.0, 0.0, 1.625,
    0.0, 0.0, 0.5, 0.0, 0.0, 1.875, 0.0, 0.0, 0.5, 0.0, 0.0, 2.125, 0.0, 0.0,
    0.5, 0.0, 0.0, 2.375, 0.0, 0.0, 0.5, 0.0, 0.0, 2.625, 0.0, 0.0, 0.5, 0.0,
    0.0, 2.875, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.125,
    0.0, 0.0, 0.5, 0.0, 0.0, 0.375, 0.0, 0.0, 0.5, 0.0, 0.0, 0.625, 0.0, 0.0,
    0.5, 0.0, 0.0, 0.875, 0.0, 0.0, 0.5, 0.0, 0.0, 1.125, 0.0, 0.0, 0.5, 0.0,
    0.0, 1.375, 0.0, 0.0, 0.5, 0.0, 0.0, 1.625, 0.0, 0.0, 0.5, 0.0, 0.0, 1.875,
    0.0, 0.0, 0.5, 0.0, 0.0, 2.125, 0.0, 0.0, 0.5, 0.0, 0.0, 2.375, 0.0, 0.0,
    0.5, 0.0, 0.0, 2.625, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.125, 0.0, 0.0, 0.5, 0.0, 0.0, 0.375, 0.0, 0.0, 0.5,
    0.0, 0.0, 0.625, 0.0, 0.0, 0.5, 0.0, 0.0, 0.875, 0.0, 0.0, 0.5, 0.0, 0.0,
    1.125, 0.0, 0.0, 0.5, 0.0, 0.0, 1.375, 0.0, 0.0, 0.5, 0.0, 0.0, 1.625, 0.0,
    0.0, 0.5, 0.0, 0.0, 1.875, 0.0, 0.0, 0.5, 0.0, 0.0, 2.125, 0.0, 0.0, 0.5,
    0.0, 0.0, 2.375, 0.0, 0.0, 0.5, 0.0, 0.0, 2.625, 0.0, 0.0, 0.5, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.125, 0.0, 0.0, 0.5, 0.0,
    0.0, 0.375, 0.0, 0.0, 0.5, 0.0, 0.0, 0.625, 0.0, 0.0, 0.5, 0.0, 0.0, 0.875,
    0.0, 0.0, 0.5, 0.0, 0.0, 1.125, 0.0, 0.0, 0.5, 0.0, 0.0, 1.375, 0.0, 0.0,
    0.5, 0.0, 0.0, 1.625, 0.0, 0.0, 0.5, 0.0, 0.0, 1.875, 0.0, 0.0, 0.5, 0.0,
    0.0, 2.125, 0.0, 0.0, 0.5, 0.0, 0.0, 2.375, 0.0, 0.0, 0.5, 0.0, 0.0, 2.625,
    0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.125, 0.0, 0.0, 0.5, 0.0, 0.0, 0.375, 0.0, 0.0, 0.5, 0.0,
    0.0, 0.625, 0.0, 0.0, 0.5, 0.0, 0.0, 0.875, 0.0, 0.0, 0.5, 0.0, 0.0, 1.125,
    0.0, 0.0, 0.5, 0.0, 0.0, 1.375, 0.0, 0.0, 0.5, 0.0, 0.0, 1.625, 0.0, 0.0,
    0.5, 0.0, 0.0, 1.875, 0.0, 0.0, 0.5, 0.0, 0.0, 2.125, 0.0, 0.0, 0.5, 0.0,
    0.0, 2.375, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.125, 0.0, 0.0, 0.5, 0.0, 0.0,
    0.375, 0.0, 0.0, 0.5, 0.0, 0.0, 0.625, 0.0, 0.0, 0.5, 0.0, 0.0, 0.875, 0.0,
    0.0, 0.5, 0.0, 0.0, 1.125, 0.0, 0.0, 0.5, 0.0, 0.0, 1.375, 0.0, 0.0, 0.5,
    0.0, 0.0, 1.625, 0.0, 0.0, 0.5, 0.0, 0.0, 1.875, 0.0, 0.0, 0.5, 0.0, 0.0,
    2.125, 0.0, 0.0, 0.5, 0.0, 0.0, 2.375, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.125,
    0.0, 0.0, 0.5, 0.0, 0.0, 0.375, 0.0, 0.0, 0.5, 0.0, 0.0, 0.625, 0.0, 0.0,
    0.5, 0.0, 0.0, 0.875, 0.0, 0.0, 0.5, 0.0, 0.0, 1.125, 0.0, 0.0, 0.5, 0.0,
    0.0, 1.375, 0.0, 0.0, 0.5, 0.0, 0.0, 1.625, 0.0, 0.0, 0.5, 0.0, 0.0, 1.875,
    0.0, 0.0, 0.5, 0.0, 0.0, 2.125, 0.0, 0.0, 0.5, 0.0, 0.0, 2.375, 0.0, 0.0,
    0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.125, 0.0, 0.0, 0.5, 0.0, 0.0, 0.375,
    0.0, 0.0, 0.5, 0.0, 0.0, 0.625, 0.0, 0.0, 0.5, 0.0, 0.0, 0.875, 0.0, 0.0,
    0.5, 0.0, 0.0, 1.125, 0.0, 0.0, 0.5, 0.0, 0.0, 1.375, 0.0, 0.0, 0.5, 0.0,
    0.0, 1.625, 0.0, 0.0, 0.5, 0.0, 0.0, 1.875, 0.0, 0.0, 0.5, 0.0, 0.0, 2.125,
    0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.125, 0.0, 0.0,
    0.5, 0.0, 0.0, 0.375, 0.0, 0.0, 0.5, 0.0, 0.0, 0.625, 0.0, 0.0, 0.5, 0.0,
    0.0, 0.875, 0.0, 0.0, 0.5, 0.0, 0.0, 1.125, 0.0, 0.0, 0.5, 0.0, 0.0, 1.375,
    0.0, 0.0, 0.5, 0.0, 0.0, 1.625, 0.0, 0.0, 0.5, 0.0, 0.0, 1.875, 0.0, 0.0,
    0.5, 0.0, 0.0, 2.125, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.125, 0.0, 0.0, 0.5, 0.0, 0.0, 0.375, 0.0, 0.0, 0.5, 0.0, 0.0, 0.625,
    0.0, 0.0, 0.5, 0.0, 0.0, 0.875, 0.0, 0.0, 0.5, 0.0, 0.0, 1.125, 0.0, 0.0,
    0.5, 0.0, 0.0, 1.375, 0.0, 0.0, 0.5, 0.0, 0.0, 1.625, 0.0, 0.0, 0.5, 0.0,
    0.0, 1.875, 0.0, 0.0, 0.5, 0.0, 0.0, 2.125, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.125, 0.0, 0.0, 0.5, 0.0, 0.0,
    0.375, 0.0, 0.0, 0.5, 0.0, 0.0, 0.625, 0.0, 0.0, 0.5, 0.0, 0.0, 0.875, 0.0,
    0.0, 0.5, 0.0, 0.0, 1.125, 0.0, 0.0, 0.5, 0.0, 0.0, 1.375, 0.0, 0.0, 0.5,
    0.0, 0.0, 1.625, 0.0, 0.0, 0.5, 0.0, 0.0, 1.875, 0.0, 0.0, 0.5, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.125, 0.0,
    0.0, 0.5, 0.0, 0.0, 0.375, 0.0, 0.0, 0.5, 0.0, 0.0, 0.625, 0.0, 0.0, 0.5,
    0.0, 0.0, 0.875, 0.0, 0.0, 0.5, 0.0, 0.0, 1.125, 0.0, 0.0, 0.5, 0.0, 0.0,
    1.375, 0.0, 0.0, 0.5, 0.0, 0.0, 1.625, 0.0, 0.0, 0.5, 0.0, 0.0, 1.875, 0.0,
    0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.125, 0.0, 0.0, 0.5, 0.0, 0.0, 0.375, 0.0, 0.0, 0.5, 0.0, 0.0,
    0.625, 0.0, 0.0, 0.5, 0.0, 0.0, 0.875, 0.0, 0.0, 0.5, 0.0, 0.0, 1.125, 0.0,
    0.0, 0.5, 0.0, 0.0, 1.375, 0.0, 0.0, 0.5, 0.0, 0.0, 1.625, 0.0, 0.0, 0.5,
    0.0, 0.0, 1.875, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.125, 0.0, 0.0, 0.5, 0.0, 0.0,
    0.375, 0.0, 0.0, 0.5, 0.0, 0.0, 0.625, 0.0, 0.0, 0.5, 0.0, 0.0, 0.875, 0.0,
    0.0, 0.5, 0.0, 0.0, 1.125, 0.0, 0.0, 0.5, 0.0, 0.0, 1.375, 0.0, 0.0, 0.5,
    0.0, 0.0, 1.625, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.125, 0.0, 0.0,
    0.5, 0.0, 0.0, 0.375, 0.0, 0.0, 0.5, 0.0, 0.0, 0.625, 0.0, 0.0, 0.5, 0.0,
    0.0, 0.875, 0.0, 0.0, 0.5, 0.0, 0.0, 1.125, 0.0, 0.0, 0.5, 0.0, 0.0, 1.375,
    0.0, 0.0, 0.5, 0.0, 0.0, 1.625, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.125, 0.0, 0.0, 0.5, 0.0, 0.0, 0.375, 0.0, 0.0, 0.5, 0.0, 0.0, 0.625, 0.0,
    0.0, 0.5, 0.0, 0.0, 0.875, 0.0, 0.0, 0.5, 0.0, 0.0, 1.125, 0.0, 0.0, 0.5,
    0.0, 0.0, 1.375, 0.0, 0.0, 0.5, 0.0, 0.0, 1.625, 0.0, 0.0, 0.5, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.125, 0.0, 0.0, 0.5, 0.0, 0.0, 0.375,
    0.0, 0.0, 0.5, 0.0, 0.0, 0.625, 0.0, 0.0, 0.5, 0.0, 0.0, 0.875, 0.0, 0.0,
    0.5, 0.0, 0.0, 1.125, 0.0, 0.0, 0.5, 0.0, 0.0, 1.375, 0.0, 0.0, 0.5, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.125, 0.0, 0.0, 0.5,
    0.0, 0.0, 0.375, 0.0, 0.0, 0.5, 0.0, 0.0, 0.625, 0.0, 0.0, 0.5, 0.0, 0.0,
    0.875, 0.0, 0.0, 0.5, 0.0, 0.0, 1.125, 0.0, 0.0, 0.5, 0.0, 0.0, 1.375, 0.0,
    0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.125,
    0.0, 0.0, 0.5, 0.0, 0.0, 0.375, 0.0, 0.0, 0.5, 0.0, 0.0, 0.625, 0.0, 0.0,
    0.5, 0.0, 0.0, 0.875, 0.0, 0.0, 0.5, 0.0, 0.0, 1.125, 0.0, 0.0, 0.5, 0.0,
    0.0, 1.375, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.125, 0.0, 0.0, 0.5, 0.0, 0.0, 0.375, 0.0, 0.0,
    0.5, 0.0, 0.0, 0.625, 0.0, 0.0, 0.5, 0.0, 0.0, 0.875, 0.0, 0.0, 0.5, 0.0,
    0.0, 1.125, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.125, 0.0, 0.0, 0.5, 0.0, 0.0,
    0.375, 0.0, 0.0, 0.5, 0.0, 0.0, 0.625, 0.0, 0.0, 0.5, 0.0, 0.0, 0.875, 0.0,
    0.0, 0.5, 0.0, 0.0, 1.125, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.125, 0.0, 0.0, 0.5,
    0.0, 0.0, 0.375, 0.0, 0.0, 0.5, 0.0, 0.0, 0.625, 0.0, 0.0, 0.5, 0.0, 0.0,
    0.875, 0.0, 0.0, 0.5, 0.0, 0.0, 1.125, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.125, 0.0, 0.0, 0.5, 0.0, 0.0, 0.375, 0.0, 0.0, 0.5, 0.0, 0.0,
    0.625, 0.0, 0.0, 0.5, 0.0, 0.0, 0.875, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.125, 0.0, 0.0, 0.5, 0.0, 0.0, 0.375, 0.0, 0.0,
    0.5, 0.0, 0.0, 0.625, 0.0, 0.0, 0.5, 0.0, 0.0, 0.875, 0.0, 0.0, 0.5, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.125, 0.0, 0.0, 0.5, 0.0, 0.0,
    0.375, 0.0, 0.0, 0.5, 0.0, 0.0, 0.625, 0.0, 0.0, 0.5, 0.0, 0.0, 0.875, 0.0,
    0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.125,
    0.0, 0.0, 0.5, 0.0, 0.0, 0.375, 0.0, 0.0, 0.5, 0.0, 0.0, 0.625, 0.0, 0.0,
    0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.125, 0.0, 0.0, 0.5, 0.0, 0.0, 0.375, 0.0, 0.0, 0.5, 0.0, 0.0, 0.625,
    0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.125, 0.0, 0.0, 0.5, 0.0, 0.0, 0.375, 0.0, 0.0, 0.5, 0.0,
    0.0, 0.625, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.125, 0.0, 0.0, 0.5, 0.0, 0.0,
    0.375, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.125, 0.0, 0.0, 0.5, 0.0,
    0.0, 0.375, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.125, 0.0, 0.0, 0.5,
    0.0, 0.0, 0.375, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.125, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.125, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.125, 0.0, 0.0, 0.5 };

  static const double a[720] = { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 1.5, 0.0, 0.0, 1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    2.5, 0.0, 0.0, 1.0, 0.0, 0.0, 3.0, 0.0, 0.0, 1.0, 0.0, 0.0, 3.5, 0.0, 0.0,
    1.0, 0.0, 0.0, 4.0, 0.0, 0.0, 1.0, 0.0, 0.0, 4.5, 0.0, 0.0, 1.0, 0.0, 0.0,
    5.0, 0.0, 0.0, 1.0, 0.0, 0.0, 5.5, 0.0, 0.0, 1.0, 0.0, 0.0, 6.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 6.5, 0.0, 0.0, 1.0, 0.0, 0.0, 7.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    7.5, 0.0, 0.0, 1.0, 0.0, 0.0, 8.0, 0.0, 0.0, 1.0, 0.0, 0.0, 8.5, 0.0, 0.0,
    1.0, 0.0, 0.0, 9.0, 0.0, 0.0, 1.0, 0.0, 0.0, 9.5, 0.0, 0.0, 1.0, 0.0, 0.0,
    10.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 1.5, 0.0, 0.0, 1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 2.5, 0.0, 0.0, 1.0, 0.0, 0.0, 3.0, 0.0, 0.0, 1.0, 0.0, 0.0, 3.5, 0.0,
    0.0, 1.0, 0.0, 0.0, 4.0, 0.0, 0.0, 1.0, 0.0, 0.0, 4.5, 0.0, 0.0, 1.0, 0.0,
    0.0, 5.0, 0.0, 0.0, 1.0, 0.0, 0.0, 5.5, 0.0, 0.0, 1.0, 0.0, 0.0, 6.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 6.5, 0.0, 0.0, 1.0, 0.0, 0.0, 7.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 7.5, 0.0, 0.0, 1.0, 0.0, 0.0, 8.0, 0.0, 0.0, 1.0, 0.0, 0.0, 8.5, 0.0,
    0.0, 1.0, 0.0, 0.0, 9.0, 0.0, 0.0, 1.0, 0.0, 0.0, 9.5, 0.0, 0.0, 1.0, 0.0,
    0.0, 10.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 1.5, 0.0, 0.0, 1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 2.5, 0.0, 0.0, 1.0, 0.0, 0.0, 3.0, 0.0, 0.0, 1.0, 0.0, 0.0, 3.5,
    0.0, 0.0, 1.0, 0.0, 0.0, 4.0, 0.0, 0.0, 1.0, 0.0, 0.0, 4.5, 0.0, 0.0, 1.0,
    0.0, 0.0, 5.0, 0.0, 0.0, 1.0, 0.0, 0.0, 5.5, 0.0, 0.0, 1.0, 0.0, 0.0, 6.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 6.5, 0.0, 0.0, 1.0, 0.0, 0.0, 7.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 7.5, 0.0, 0.0, 1.0, 0.0, 0.0, 8.0, 0.0, 0.0, 1.0, 0.0, 0.0, 8.5,
    0.0, 0.0, 1.0, 0.0, 0.0, 9.0, 0.0, 0.0, 1.0, 0.0, 0.0, 9.5, 0.0, 0.0, 1.0,
    0.0, 0.0, 10.0, 0.0, 0.0, 1.0 };

  double d;
  double d1;
  for (int i = 0; i < 120; i++) {
    int i1;
    d = 0.0;
    for (i1 = 0; i1 < 6; i1++) {
      d += a[i + (120 * i1)] * x0_vl[i1];
    }

    d1 = 0.0;
    for (i1 = 0; i1 < 60; i1++) {
      d1 += b_a[i + (120 * i1)] * X_QP[i1];
    }

    x_pred[i] = d + d1;
  }
}

template<typename T>
void CoordinatorBase<T>::MPC_vl()
{
 static const double b_H[3600] = { 14.5253, 0.0, 0.0, 11.8592, 0.0, 0.0,
    11.1948, 0.0, 0.0, 10.5338, 0.0, 0.0, 9.8779, 0.0, 0.0, 9.2289, 0.0, 0.0,
    8.5883, 0.0, 0.0, 7.958, 0.0, 0.0, 7.3395, 0.0, 0.0, 6.7347, 0.0, 0.0,
    6.1452, 0.0, 0.0, 5.5728, 0.0, 0.0, 5.0191, 0.0, 0.0, 4.4858, 0.0, 0.0,
    3.9746, 0.0, 0.0, 3.4873, 0.0, 0.0, 3.0255, 0.0, 0.0, 2.5909, 0.0, 0.0,
    2.1853, 0.0, 0.0, 1.8103, 0.0, 0.0, 0.0, 14.5253, 0.0, 0.0, 11.8592, 0.0,
    0.0, 11.1948, 0.0, 0.0, 10.5338, 0.0, 0.0, 9.8779, 0.0, 0.0, 9.2289, 0.0,
    0.0, 8.5883, 0.0, 0.0, 7.958, 0.0, 0.0, 7.3395, 0.0, 0.0, 6.7347, 0.0, 0.0,
    6.1452, 0.0, 0.0, 5.5728, 0.0, 0.0, 5.0191, 0.0, 0.0, 4.4858, 0.0, 0.0,
    3.9746, 0.0, 0.0, 3.4873, 0.0, 0.0, 3.0255, 0.0, 0.0, 2.5909, 0.0, 0.0,
    2.1853, 0.0, 0.0, 1.8103, 0.0, 0.0, 0.0, 14.5253, 0.0, 0.0, 11.8592, 0.0,
    0.0, 11.1948, 0.0, 0.0, 10.5338, 0.0, 0.0, 9.8779, 0.0, 0.0, 9.2289, 0.0,
    0.0, 8.5883, 0.0, 0.0, 7.958, 0.0, 0.0, 7.3395, 0.0, 0.0, 6.7347, 0.0, 0.0,
    6.1452, 0.0, 0.0, 5.5728, 0.0, 0.0, 5.0191, 0.0, 0.0, 4.4858, 0.0, 0.0,
    3.9746, 0.0, 0.0, 3.4873, 0.0, 0.0, 3.0255, 0.0, 0.0, 2.5909, 0.0, 0.0,
    2.1853, 0.0, 0.0, 1.8103, 11.8592, 0.0, 0.0, 13.2435, 0.0, 0.0, 10.6241, 0.0,
    0.0, 10.0064, 0.0, 0.0, 9.3921, 0.0, 0.0, 8.783, 0.0, 0.0, 8.1806, 0.0, 0.0,
    7.5868, 0.0, 0.0, 7.0031, 0.0, 0.0, 6.4314, 0.0, 0.0, 5.8733, 0.0, 0.0,
    5.3306, 0.0, 0.0, 4.8048, 0.0, 0.0, 4.2978, 0.0, 0.0, 3.8112, 0.0, 0.0,
    3.3468, 0.0, 0.0, 2.9062, 0.0, 0.0, 2.4911, 0.0, 0.0, 2.1033, 0.0, 0.0,
    1.7444, 0.0, 0.0, 0.0, 11.8592, 0.0, 0.0, 13.2435, 0.0, 0.0, 10.6241, 0.0,
    0.0, 10.0064, 0.0, 0.0, 9.3921, 0.0, 0.0, 8.783, 0.0, 0.0, 8.1806, 0.0, 0.0,
    7.5868, 0.0, 0.0, 7.0031, 0.0, 0.0, 6.4314, 0.0, 0.0, 5.8733, 0.0, 0.0,
    5.3306, 0.0, 0.0, 4.8048, 0.0, 0.0, 4.2978, 0.0, 0.0, 3.8112, 0.0, 0.0,
    3.3468, 0.0, 0.0, 2.9062, 0.0, 0.0, 2.4911, 0.0, 0.0, 2.1033, 0.0, 0.0,
    1.7444, 0.0, 0.0, 0.0, 11.8592, 0.0, 0.0, 13.2435, 0.0, 0.0, 10.6241, 0.0,
    0.0, 10.0064, 0.0, 0.0, 9.3921, 0.0, 0.0, 8.783, 0.0, 0.0, 8.1806, 0.0, 0.0,
    7.5868, 0.0, 0.0, 7.0031, 0.0, 0.0, 6.4314, 0.0, 0.0, 5.8733, 0.0, 0.0,
    5.3306, 0.0, 0.0, 4.8048, 0.0, 0.0, 4.2978, 0.0, 0.0, 3.8112, 0.0, 0.0,
    3.3468, 0.0, 0.0, 2.9062, 0.0, 0.0, 2.4911, 0.0, 0.0, 2.1033, 0.0, 0.0,
    1.7444, 11.1948, 0.0, 0.0, 10.6241, 0.0, 0.0, 12.0534, 0.0, 0.0, 9.479, 0.0,
    0.0, 8.9063, 0.0, 0.0, 8.3371, 0.0, 0.0, 7.7729, 0.0, 0.0, 7.2156, 0.0, 0.0,
    6.6668, 0.0, 0.0, 6.1281, 0.0, 0.0, 5.6014, 0.0, 0.0, 5.0884, 0.0, 0.0,
    4.5906, 0.0, 0.0, 4.1099, 0.0, 0.0, 3.6479, 0.0, 0.0, 3.2063, 0.0, 0.0,
    2.7869, 0.0, 0.0, 2.3913, 0.0, 0.0, 2.0212, 0.0, 0.0, 1.6784, 0.0, 0.0, 0.0,
    11.1948, 0.0, 0.0, 10.6241, 0.0, 0.0, 12.0534, 0.0, 0.0, 9.479, 0.0, 0.0,
    8.9063, 0.0, 0.0, 8.3371, 0.0, 0.0, 7.7729, 0.0, 0.0, 7.2156, 0.0, 0.0,
    6.6668, 0.0, 0.0, 6.1281, 0.0, 0.0, 5.6014, 0.0, 0.0, 5.0884, 0.0, 0.0,
    4.5906, 0.0, 0.0, 4.1099, 0.0, 0.0, 3.6479, 0.0, 0.0, 3.2063, 0.0, 0.0,
    2.7869, 0.0, 0.0, 2.3913, 0.0, 0.0, 2.0212, 0.0, 0.0, 1.6784, 0.0, 0.0, 0.0,
    11.1948, 0.0, 0.0, 10.6241, 0.0, 0.0, 12.0534, 0.0, 0.0, 9.479, 0.0, 0.0,
    8.9063, 0.0, 0.0, 8.3371, 0.0, 0.0, 7.7729, 0.0, 0.0, 7.2156, 0.0, 0.0,
    6.6668, 0.0, 0.0, 6.1281, 0.0, 0.0, 5.6014, 0.0, 0.0, 5.0884, 0.0, 0.0,
    4.5906, 0.0, 0.0, 4.1099, 0.0, 0.0, 3.6479, 0.0, 0.0, 3.2063, 0.0, 0.0,
    2.7869, 0.0, 0.0, 2.3913, 0.0, 0.0, 2.0212, 0.0, 0.0, 1.6784, 10.5338, 0.0,
    0.0, 10.0064, 0.0, 0.0, 9.479, 0.0, 0.0, 10.9516, 0.0, 0.0, 8.4206, 0.0, 0.0,
    7.8912, 0.0, 0.0, 7.3653, 0.0, 0.0, 6.8444, 0.0, 0.0, 6.3304, 0.0, 0.0,
    5.8249, 0.0, 0.0, 5.3296, 0.0, 0.0, 4.8462, 0.0, 0.0, 4.3764, 0.0, 0.0,
    3.922, 0.0, 0.0, 3.4845, 0.0, 0.0, 3.0658, 0.0, 0.0, 2.6676, 0.0, 0.0,
    2.2915, 0.0, 0.0, 1.9392, 0.0, 0.0, 1.6124, 0.0, 0.0, 0.0, 10.5338, 0.0, 0.0,
    10.0064, 0.0, 0.0, 9.479, 0.0, 0.0, 10.9516, 0.0, 0.0, 8.4206, 0.0, 0.0,
    7.8912, 0.0, 0.0, 7.3653, 0.0, 0.0, 6.8444, 0.0, 0.0, 6.3304, 0.0, 0.0,
    5.8249, 0.0, 0.0, 5.3296, 0.0, 0.0, 4.8462, 0.0, 0.0, 4.3764, 0.0, 0.0,
    3.922, 0.0, 0.0, 3.4845, 0.0, 0.0, 3.0658, 0.0, 0.0, 2.6676, 0.0, 0.0,
    2.2915, 0.0, 0.0, 1.9392, 0.0, 0.0, 1.6124, 0.0, 0.0, 0.0, 10.5338, 0.0, 0.0,
    10.0064, 0.0, 0.0, 9.479, 0.0, 0.0, 10.9516, 0.0, 0.0, 8.4206, 0.0, 0.0,
    7.8912, 0.0, 0.0, 7.3653, 0.0, 0.0, 6.8444, 0.0, 0.0, 6.3304, 0.0, 0.0,
    5.8249, 0.0, 0.0, 5.3296, 0.0, 0.0, 4.8462, 0.0, 0.0, 4.3764, 0.0, 0.0,
    3.922, 0.0, 0.0, 3.4845, 0.0, 0.0, 3.0658, 0.0, 0.0, 2.6676, 0.0, 0.0,
    2.2915, 0.0, 0.0, 1.9392, 0.0, 0.0, 1.6124, 9.8779, 0.0, 0.0, 9.3921, 0.0,
    0.0, 8.9063, 0.0, 0.0, 8.4206, 0.0, 0.0, 9.9348, 0.0, 0.0, 7.4453, 0.0, 0.0,
    6.9576, 0.0, 0.0, 6.4732, 0.0, 0.0, 5.994, 0.0, 0.0, 5.5216, 0.0, 0.0,
    5.0577, 0.0, 0.0, 4.604, 0.0, 0.0, 4.1622, 0.0, 0.0, 3.734, 0.0, 0.0, 3.3212,
    0.0, 0.0, 2.9254, 0.0, 0.0, 2.5483, 0.0, 0.0, 2.1916, 0.0, 0.0, 1.8571, 0.0,
    0.0, 1.5464, 0.0, 0.0, 0.0, 9.8779, 0.0, 0.0, 9.3921, 0.0, 0.0, 8.9063, 0.0,
    0.0, 8.4206, 0.0, 0.0, 9.9348, 0.0, 0.0, 7.4453, 0.0, 0.0, 6.9576, 0.0, 0.0,
    6.4732, 0.0, 0.0, 5.994, 0.0, 0.0, 5.5216, 0.0, 0.0, 5.0577, 0.0, 0.0, 4.604,
    0.0, 0.0, 4.1622, 0.0, 0.0, 3.734, 0.0, 0.0, 3.3212, 0.0, 0.0, 2.9254, 0.0,
    0.0, 2.5483, 0.0, 0.0, 2.1916, 0.0, 0.0, 1.8571, 0.0, 0.0, 1.5464, 0.0, 0.0,
    0.0, 9.8779, 0.0, 0.0, 9.3921, 0.0, 0.0, 8.9063, 0.0, 0.0, 8.4206, 0.0, 0.0,
    9.9348, 0.0, 0.0, 7.4453, 0.0, 0.0, 6.9576, 0.0, 0.0, 6.4732, 0.0, 0.0,
    5.994, 0.0, 0.0, 5.5216, 0.0, 0.0, 5.0577, 0.0, 0.0, 4.604, 0.0, 0.0, 4.1622,
    0.0, 0.0, 3.734, 0.0, 0.0, 3.3212, 0.0, 0.0, 2.9254, 0.0, 0.0, 2.5483, 0.0,
    0.0, 2.1916, 0.0, 0.0, 1.8571, 0.0, 0.0, 1.5464, 9.2289, 0.0, 0.0, 8.783,
    0.0, 0.0, 8.3371, 0.0, 0.0, 7.8912, 0.0, 0.0, 7.4453, 0.0, 0.0, 8.9994, 0.0,
    0.0, 6.5499, 0.0, 0.0, 6.102, 0.0, 0.0, 5.6576, 0.0, 0.0, 5.2183, 0.0, 0.0,
    4.7858, 0.0, 0.0, 4.3618, 0.0, 0.0, 3.948, 0.0, 0.0, 3.5461, 0.0, 0.0,
    3.1578, 0.0, 0.0, 2.7849, 0.0, 0.0, 2.429, 0.0, 0.0, 2.0918, 0.0, 0.0,
    1.7751, 0.0, 0.0, 1.4805, 0.0, 0.0, 0.0, 9.2289, 0.0, 0.0, 8.783, 0.0, 0.0,
    8.3371, 0.0, 0.0, 7.8912, 0.0, 0.0, 7.4453, 0.0, 0.0, 8.9994, 0.0, 0.0,
    6.5499, 0.0, 0.0, 6.102, 0.0, 0.0, 5.6576, 0.0, 0.0, 5.2183, 0.0, 0.0,
    4.7858, 0.0, 0.0, 4.3618, 0.0, 0.0, 3.948, 0.0, 0.0, 3.5461, 0.0, 0.0,
    3.1578, 0.0, 0.0, 2.7849, 0.0, 0.0, 2.429, 0.0, 0.0, 2.0918, 0.0, 0.0,
    1.7751, 0.0, 0.0, 1.4805, 0.0, 0.0, 0.0, 9.2289, 0.0, 0.0, 8.783, 0.0, 0.0,
    8.3371, 0.0, 0.0, 7.8912, 0.0, 0.0, 7.4453, 0.0, 0.0, 8.9994, 0.0, 0.0,
    6.5499, 0.0, 0.0, 6.102, 0.0, 0.0, 5.6576, 0.0, 0.0, 5.2183, 0.0, 0.0,
    4.7858, 0.0, 0.0, 4.3618, 0.0, 0.0, 3.948, 0.0, 0.0, 3.5461, 0.0, 0.0,
    3.1578, 0.0, 0.0, 2.7849, 0.0, 0.0, 2.429, 0.0, 0.0, 2.0918, 0.0, 0.0,
    1.7751, 0.0, 0.0, 1.4805, 8.5883, 0.0, 0.0, 8.1806, 0.0, 0.0, 7.7729, 0.0,
    0.0, 7.3653, 0.0, 0.0, 6.9576, 0.0, 0.0, 6.5499, 0.0, 0.0, 8.1422, 0.0, 0.0,
    5.7309, 0.0, 0.0, 5.3212, 0.0, 0.0, 4.915, 0.0, 0.0, 4.5139, 0.0, 0.0,
    4.1195, 0.0, 0.0, 3.7337, 0.0, 0.0, 3.3582, 0.0, 0.0, 2.9945, 0.0, 0.0,
    2.6444, 0.0, 0.0, 2.3097, 0.0, 0.0, 1.992, 0.0, 0.0, 1.693, 0.0, 0.0, 1.4145,
    0.0, 0.0, 0.0, 8.5883, 0.0, 0.0, 8.1806, 0.0, 0.0, 7.7729, 0.0, 0.0, 7.3653,
    0.0, 0.0, 6.9576, 0.0, 0.0, 6.5499, 0.0, 0.0, 8.1422, 0.0, 0.0, 5.7309, 0.0,
    0.0, 5.3212, 0.0, 0.0, 4.915, 0.0, 0.0, 4.5139, 0.0, 0.0, 4.1195, 0.0, 0.0,
    3.7337, 0.0, 0.0, 3.3582, 0.0, 0.0, 2.9945, 0.0, 0.0, 2.6444, 0.0, 0.0,
    2.3097, 0.0, 0.0, 1.992, 0.0, 0.0, 1.693, 0.0, 0.0, 1.4145, 0.0, 0.0, 0.0,
    8.5883, 0.0, 0.0, 8.1806, 0.0, 0.0, 7.7729, 0.0, 0.0, 7.3653, 0.0, 0.0,
    6.9576, 0.0, 0.0, 6.5499, 0.0, 0.0, 8.1422, 0.0, 0.0, 5.7309, 0.0, 0.0,
    5.3212, 0.0, 0.0, 4.915, 0.0, 0.0, 4.5139, 0.0, 0.0, 4.1195, 0.0, 0.0,
    3.7337, 0.0, 0.0, 3.3582, 0.0, 0.0, 2.9945, 0.0, 0.0, 2.6444, 0.0, 0.0,
    2.3097, 0.0, 0.0, 1.992, 0.0, 0.0, 1.693, 0.0, 0.0, 1.4145, 7.958, 0.0, 0.0,
    7.5868, 0.0, 0.0, 7.2156, 0.0, 0.0, 6.8444, 0.0, 0.0, 6.4732, 0.0, 0.0,
    6.102, 0.0, 0.0, 5.7309, 0.0, 0.0, 7.3597, 0.0, 0.0, 4.9848, 0.0, 0.0,
    4.6117, 0.0, 0.0, 4.242, 0.0, 0.0, 3.8773, 0.0, 0.0, 3.5195, 0.0, 0.0,
    3.1702, 0.0, 0.0, 2.8311, 0.0, 0.0, 2.504, 0.0, 0.0, 2.1904, 0.0, 0.0,
    1.8922, 0.0, 0.0, 1.611, 0.0, 0.0, 1.3485, 0.0, 0.0, 0.0, 7.958, 0.0, 0.0,
    7.5868, 0.0, 0.0, 7.2156, 0.0, 0.0, 6.8444, 0.0, 0.0, 6.4732, 0.0, 0.0,
    6.102, 0.0, 0.0, 5.7309, 0.0, 0.0, 7.3597, 0.0, 0.0, 4.9848, 0.0, 0.0,
    4.6117, 0.0, 0.0, 4.242, 0.0, 0.0, 3.8773, 0.0, 0.0, 3.5195, 0.0, 0.0,
    3.1702, 0.0, 0.0, 2.8311, 0.0, 0.0, 2.504, 0.0, 0.0, 2.1904, 0.0, 0.0,
    1.8922, 0.0, 0.0, 1.611, 0.0, 0.0, 1.3485, 0.0, 0.0, 0.0, 7.958, 0.0, 0.0,
    7.5868, 0.0, 0.0, 7.2156, 0.0, 0.0, 6.8444, 0.0, 0.0, 6.4732, 0.0, 0.0,
    6.102, 0.0, 0.0, 5.7309, 0.0, 0.0, 7.3597, 0.0, 0.0, 4.9848, 0.0, 0.0,
    4.6117, 0.0, 0.0, 4.242, 0.0, 0.0, 3.8773, 0.0, 0.0, 3.5195, 0.0, 0.0,
    3.1702, 0.0, 0.0, 2.8311, 0.0, 0.0, 2.504, 0.0, 0.0, 2.1904, 0.0, 0.0,
    1.8922, 0.0, 0.0, 1.611, 0.0, 0.0, 1.3485, 7.3395, 0.0, 0.0, 7.0031, 0.0,
    0.0, 6.6668, 0.0, 0.0, 6.3304, 0.0, 0.0, 5.994, 0.0, 0.0, 5.6576, 0.0, 0.0,
    5.3212, 0.0, 0.0, 4.9848, 0.0, 0.0, 6.6485, 0.0, 0.0, 4.3084, 0.0, 0.0,
    3.9701, 0.0, 0.0, 3.6351, 0.0, 0.0, 3.3053, 0.0, 0.0, 2.9823, 0.0, 0.0,
    2.6678, 0.0, 0.0, 2.3635, 0.0, 0.0, 2.0711, 0.0, 0.0, 1.7924, 0.0, 0.0,
    1.5289, 0.0, 0.0, 1.2825, 0.0, 0.0, 0.0, 7.3395, 0.0, 0.0, 7.0031, 0.0, 0.0,
    6.6668, 0.0, 0.0, 6.3304, 0.0, 0.0, 5.994, 0.0, 0.0, 5.6576, 0.0, 0.0,
    5.3212, 0.0, 0.0, 4.9848, 0.0, 0.0, 6.6485, 0.0, 0.0, 4.3084, 0.0, 0.0,
    3.9701, 0.0, 0.0, 3.6351, 0.0, 0.0, 3.3053, 0.0, 0.0, 2.9823, 0.0, 0.0,
    2.6678, 0.0, 0.0, 2.3635, 0.0, 0.0, 2.0711, 0.0, 0.0, 1.7924, 0.0, 0.0,
    1.5289, 0.0, 0.0, 1.2825, 0.0, 0.0, 0.0, 7.3395, 0.0, 0.0, 7.0031, 0.0, 0.0,
    6.6668, 0.0, 0.0, 6.3304, 0.0, 0.0, 5.994, 0.0, 0.0, 5.6576, 0.0, 0.0,
    5.3212, 0.0, 0.0, 4.9848, 0.0, 0.0, 6.6485, 0.0, 0.0, 4.3084, 0.0, 0.0,
    3.9701, 0.0, 0.0, 3.6351, 0.0, 0.0, 3.3053, 0.0, 0.0, 2.9823, 0.0, 0.0,
    2.6678, 0.0, 0.0, 2.3635, 0.0, 0.0, 2.0711, 0.0, 0.0, 1.7924, 0.0, 0.0,
    1.5289, 0.0, 0.0, 1.2825, 6.7347, 0.0, 0.0, 6.4314, 0.0, 0.0, 6.1281, 0.0,
    0.0, 5.8249, 0.0, 0.0, 5.5216, 0.0, 0.0, 5.2183, 0.0, 0.0, 4.915, 0.0, 0.0,
    4.6117, 0.0, 0.0, 4.3084, 0.0, 0.0, 6.0051, 0.0, 0.0, 3.6982, 0.0, 0.0,
    3.3929, 0.0, 0.0, 3.0911, 0.0, 0.0, 2.7944, 0.0, 0.0, 2.5044, 0.0, 0.0,
    2.223, 0.0, 0.0, 1.9518, 0.0, 0.0, 1.6925, 0.0, 0.0, 1.4469, 0.0, 0.0,
    1.2165, 0.0, 0.0, 0.0, 6.7347, 0.0, 0.0, 6.4314, 0.0, 0.0, 6.1281, 0.0, 0.0,
    5.8249, 0.0, 0.0, 5.5216, 0.0, 0.0, 5.2183, 0.0, 0.0, 4.915, 0.0, 0.0,
    4.6117, 0.0, 0.0, 4.3084, 0.0, 0.0, 6.0051, 0.0, 0.0, 3.6982, 0.0, 0.0,
    3.3929, 0.0, 0.0, 3.0911, 0.0, 0.0, 2.7944, 0.0, 0.0, 2.5044, 0.0, 0.0,
    2.223, 0.0, 0.0, 1.9518, 0.0, 0.0, 1.6925, 0.0, 0.0, 1.4469, 0.0, 0.0,
    1.2165, 0.0, 0.0, 0.0, 6.7347, 0.0, 0.0, 6.4314, 0.0, 0.0, 6.1281, 0.0, 0.0,
    5.8249, 0.0, 0.0, 5.5216, 0.0, 0.0, 5.2183, 0.0, 0.0, 4.915, 0.0, 0.0,
    4.6117, 0.0, 0.0, 4.3084, 0.0, 0.0, 6.0051, 0.0, 0.0, 3.6982, 0.0, 0.0,
    3.3929, 0.0, 0.0, 3.0911, 0.0, 0.0, 2.7944, 0.0, 0.0, 2.5044, 0.0, 0.0,
    2.223, 0.0, 0.0, 1.9518, 0.0, 0.0, 1.6925, 0.0, 0.0, 1.4469, 0.0, 0.0,
    1.2165, 6.1452, 0.0, 0.0, 5.8733, 0.0, 0.0, 5.6014, 0.0, 0.0, 5.3296, 0.0,
    0.0, 5.0577, 0.0, 0.0, 4.7858, 0.0, 0.0, 4.5139, 0.0, 0.0, 4.242, 0.0, 0.0,
    3.9701, 0.0, 0.0, 3.6982, 0.0, 0.0, 5.4263, 0.0, 0.0, 3.1507, 0.0, 0.0,
    2.8769, 0.0, 0.0, 2.6064, 0.0, 0.0, 2.3411, 0.0, 0.0, 2.0825, 0.0, 0.0,
    1.8325, 0.0, 0.0, 1.5927, 0.0, 0.0, 1.3648, 0.0, 0.0, 1.1506, 0.0, 0.0, 0.0,
    6.1452, 0.0, 0.0, 5.8733, 0.0, 0.0, 5.6014, 0.0, 0.0, 5.3296, 0.0, 0.0,
    5.0577, 0.0, 0.0, 4.7858, 0.0, 0.0, 4.5139, 0.0, 0.0, 4.242, 0.0, 0.0,
    3.9701, 0.0, 0.0, 3.6982, 0.0, 0.0, 5.4263, 0.0, 0.0, 3.1507, 0.0, 0.0,
    2.8769, 0.0, 0.0, 2.6064, 0.0, 0.0, 2.3411, 0.0, 0.0, 2.0825, 0.0, 0.0,
    1.8325, 0.0, 0.0, 1.5927, 0.0, 0.0, 1.3648, 0.0, 0.0, 1.1506, 0.0, 0.0, 0.0,
    6.1452, 0.0, 0.0, 5.8733, 0.0, 0.0, 5.6014, 0.0, 0.0, 5.3296, 0.0, 0.0,
    5.0577, 0.0, 0.0, 4.7858, 0.0, 0.0, 4.5139, 0.0, 0.0, 4.242, 0.0, 0.0,
    3.9701, 0.0, 0.0, 3.6982, 0.0, 0.0, 5.4263, 0.0, 0.0, 3.1507, 0.0, 0.0,
    2.8769, 0.0, 0.0, 2.6064, 0.0, 0.0, 2.3411, 0.0, 0.0, 2.0825, 0.0, 0.0,
    1.8325, 0.0, 0.0, 1.5927, 0.0, 0.0, 1.3648, 0.0, 0.0, 1.1506, 5.5728, 0.0,
    0.0, 5.3306, 0.0, 0.0, 5.0884, 0.0, 0.0, 4.8462, 0.0, 0.0, 4.604, 0.0, 0.0,
    4.3618, 0.0, 0.0, 4.1195, 0.0, 0.0, 3.8773, 0.0, 0.0, 3.6351, 0.0, 0.0,
    3.3929, 0.0, 0.0, 3.1507, 0.0, 0.0, 4.9085, 0.0, 0.0, 2.6626, 0.0, 0.0,
    2.4185, 0.0, 0.0, 2.1777, 0.0, 0.0, 1.9421, 0.0, 0.0, 1.7132, 0.0, 0.0,
    1.4929, 0.0, 0.0, 1.2828, 0.0, 0.0, 1.0846, 0.0, 0.0, 0.0, 5.5728, 0.0, 0.0,
    5.3306, 0.0, 0.0, 5.0884, 0.0, 0.0, 4.8462, 0.0, 0.0, 4.604, 0.0, 0.0,
    4.3618, 0.0, 0.0, 4.1195, 0.0, 0.0, 3.8773, 0.0, 0.0, 3.6351, 0.0, 0.0,
    3.3929, 0.0, 0.0, 3.1507, 0.0, 0.0, 4.9085, 0.0, 0.0, 2.6626, 0.0, 0.0,
    2.4185, 0.0, 0.0, 2.1777, 0.0, 0.0, 1.9421, 0.0, 0.0, 1.7132, 0.0, 0.0,
    1.4929, 0.0, 0.0, 1.2828, 0.0, 0.0, 1.0846, 0.0, 0.0, 0.0, 5.5728, 0.0, 0.0,
    5.3306, 0.0, 0.0, 5.0884, 0.0, 0.0, 4.8462, 0.0, 0.0, 4.604, 0.0, 0.0,
    4.3618, 0.0, 0.0, 4.1195, 0.0, 0.0, 3.8773, 0.0, 0.0, 3.6351, 0.0, 0.0,
    3.3929, 0.0, 0.0, 3.1507, 0.0, 0.0, 4.9085, 0.0, 0.0, 2.6626, 0.0, 0.0,
    2.4185, 0.0, 0.0, 2.1777, 0.0, 0.0, 1.9421, 0.0, 0.0, 1.7132, 0.0, 0.0,
    1.4929, 0.0, 0.0, 1.2828, 0.0, 0.0, 1.0846, 5.0191, 0.0, 0.0, 4.8048, 0.0,
    0.0, 4.5906, 0.0, 0.0, 4.3764, 0.0, 0.0, 4.1622, 0.0, 0.0, 3.948, 0.0, 0.0,
    3.7337, 0.0, 0.0, 3.5195, 0.0, 0.0, 3.3053, 0.0, 0.0, 3.0911, 0.0, 0.0,
    2.8769, 0.0, 0.0, 2.6626, 0.0, 0.0, 4.4484, 0.0, 0.0, 2.2305, 0.0, 0.0,
    2.0144, 0.0, 0.0, 1.8016, 0.0, 0.0, 1.5939, 0.0, 0.0, 1.3931, 0.0, 0.0,
    1.2007, 0.0, 0.0, 1.0186, 0.0, 0.0, 0.0, 5.0191, 0.0, 0.0, 4.8048, 0.0, 0.0,
    4.5906, 0.0, 0.0, 4.3764, 0.0, 0.0, 4.1622, 0.0, 0.0, 3.948, 0.0, 0.0,
    3.7337, 0.0, 0.0, 3.5195, 0.0, 0.0, 3.3053, 0.0, 0.0, 3.0911, 0.0, 0.0,
    2.8769, 0.0, 0.0, 2.6626, 0.0, 0.0, 4.4484, 0.0, 0.0, 2.2305, 0.0, 0.0,
    2.0144, 0.0, 0.0, 1.8016, 0.0, 0.0, 1.5939, 0.0, 0.0, 1.3931, 0.0, 0.0,
    1.2007, 0.0, 0.0, 1.0186, 0.0, 0.0, 0.0, 5.0191, 0.0, 0.0, 4.8048, 0.0, 0.0,
    4.5906, 0.0, 0.0, 4.3764, 0.0, 0.0, 4.1622, 0.0, 0.0, 3.948, 0.0, 0.0,
    3.7337, 0.0, 0.0, 3.5195, 0.0, 0.0, 3.3053, 0.0, 0.0, 3.0911, 0.0, 0.0,
    2.8769, 0.0, 0.0, 2.6626, 0.0, 0.0, 4.4484, 0.0, 0.0, 2.2305, 0.0, 0.0,
    2.0144, 0.0, 0.0, 1.8016, 0.0, 0.0, 1.5939, 0.0, 0.0, 1.3931, 0.0, 0.0,
    1.2007, 0.0, 0.0, 1.0186, 4.4858, 0.0, 0.0, 4.2978, 0.0, 0.0, 4.1099, 0.0,
    0.0, 3.922, 0.0, 0.0, 3.734, 0.0, 0.0, 3.5461, 0.0, 0.0, 3.3582, 0.0, 0.0,
    3.1702, 0.0, 0.0, 2.9823, 0.0, 0.0, 2.7944, 0.0, 0.0, 2.6064, 0.0, 0.0,
    2.4185, 0.0, 0.0, 2.2305, 0.0, 0.0, 4.0426, 0.0, 0.0, 1.851, 0.0, 0.0,
    1.6611, 0.0, 0.0, 1.4746, 0.0, 0.0, 1.2933, 0.0, 0.0, 1.1187, 0.0, 0.0,
    0.95264, 0.0, 0.0, 0.0, 4.4858, 0.0, 0.0, 4.2978, 0.0, 0.0, 4.1099, 0.0, 0.0,
    3.922, 0.0, 0.0, 3.734, 0.0, 0.0, 3.5461, 0.0, 0.0, 3.3582, 0.0, 0.0, 3.1702,
    0.0, 0.0, 2.9823, 0.0, 0.0, 2.7944, 0.0, 0.0, 2.6064, 0.0, 0.0, 2.4185, 0.0,
    0.0, 2.2305, 0.0, 0.0, 4.0426, 0.0, 0.0, 1.851, 0.0, 0.0, 1.6611, 0.0, 0.0,
    1.4746, 0.0, 0.0, 1.2933, 0.0, 0.0, 1.1187, 0.0, 0.0, 0.95264, 0.0, 0.0, 0.0,
    4.4858, 0.0, 0.0, 4.2978, 0.0, 0.0, 4.1099, 0.0, 0.0, 3.922, 0.0, 0.0, 3.734,
    0.0, 0.0, 3.5461, 0.0, 0.0, 3.3582, 0.0, 0.0, 3.1702, 0.0, 0.0, 2.9823, 0.0,
    0.0, 2.7944, 0.0, 0.0, 2.6064, 0.0, 0.0, 2.4185, 0.0, 0.0, 2.2305, 0.0, 0.0,
    4.0426, 0.0, 0.0, 1.851, 0.0, 0.0, 1.6611, 0.0, 0.0, 1.4746, 0.0, 0.0,
    1.2933, 0.0, 0.0, 1.1187, 0.0, 0.0, 0.95264, 3.9746, 0.0, 0.0, 3.8112, 0.0,
    0.0, 3.6479, 0.0, 0.0, 3.4845, 0.0, 0.0, 3.3212, 0.0, 0.0, 3.1578, 0.0, 0.0,
    2.9945, 0.0, 0.0, 2.8311, 0.0, 0.0, 2.6678, 0.0, 0.0, 2.5044, 0.0, 0.0,
    2.3411, 0.0, 0.0, 2.1777, 0.0, 0.0, 2.0144, 0.0, 0.0, 1.851, 0.0, 0.0,
    3.6877, 0.0, 0.0, 1.5207, 0.0, 0.0, 1.3553, 0.0, 0.0, 1.1934, 0.0, 0.0,
    1.0366, 0.0, 0.0, 0.88666, 0.0, 0.0, 0.0, 3.9746, 0.0, 0.0, 3.8112, 0.0, 0.0,
    3.6479, 0.0, 0.0, 3.4845, 0.0, 0.0, 3.3212, 0.0, 0.0, 3.1578, 0.0, 0.0,
    2.9945, 0.0, 0.0, 2.8311, 0.0, 0.0, 2.6678, 0.0, 0.0, 2.5044, 0.0, 0.0,
    2.3411, 0.0, 0.0, 2.1777, 0.0, 0.0, 2.0144, 0.0, 0.0, 1.851, 0.0, 0.0,
    3.6877, 0.0, 0.0, 1.5207, 0.0, 0.0, 1.3553, 0.0, 0.0, 1.1934, 0.0, 0.0,
    1.0366, 0.0, 0.0, 0.88666, 0.0, 0.0, 0.0, 3.9746, 0.0, 0.0, 3.8112, 0.0, 0.0,
    3.6479, 0.0, 0.0, 3.4845, 0.0, 0.0, 3.3212, 0.0, 0.0, 3.1578, 0.0, 0.0,
    2.9945, 0.0, 0.0, 2.8311, 0.0, 0.0, 2.6678, 0.0, 0.0, 2.5044, 0.0, 0.0,
    2.3411, 0.0, 0.0, 2.1777, 0.0, 0.0, 2.0144, 0.0, 0.0, 1.851, 0.0, 0.0,
    3.6877, 0.0, 0.0, 1.5207, 0.0, 0.0, 1.3553, 0.0, 0.0, 1.1934, 0.0, 0.0,
    1.0366, 0.0, 0.0, 0.88666, 3.4873, 0.0, 0.0, 3.3468, 0.0, 0.0, 3.2063, 0.0,
    0.0, 3.0658, 0.0, 0.0, 2.9254, 0.0, 0.0, 2.7849, 0.0, 0.0, 2.6444, 0.0, 0.0,
    2.504, 0.0, 0.0, 2.3635, 0.0, 0.0, 2.223, 0.0, 0.0, 2.0825, 0.0, 0.0, 1.9421,
    0.0, 0.0, 1.8016, 0.0, 0.0, 1.6611, 0.0, 0.0, 1.5207, 0.0, 0.0, 3.3802, 0.0,
    0.0, 1.2361, 0.0, 0.0, 1.0936, 0.0, 0.0, 0.9546, 0.0, 0.0, 0.82068, 0.0, 0.0,
    0.0, 3.4873, 0.0, 0.0, 3.3468, 0.0, 0.0, 3.2063, 0.0, 0.0, 3.0658, 0.0, 0.0,
    2.9254, 0.0, 0.0, 2.7849, 0.0, 0.0, 2.6444, 0.0, 0.0, 2.504, 0.0, 0.0,
    2.3635, 0.0, 0.0, 2.223, 0.0, 0.0, 2.0825, 0.0, 0.0, 1.9421, 0.0, 0.0,
    1.8016, 0.0, 0.0, 1.6611, 0.0, 0.0, 1.5207, 0.0, 0.0, 3.3802, 0.0, 0.0,
    1.2361, 0.0, 0.0, 1.0936, 0.0, 0.0, 0.9546, 0.0, 0.0, 0.82068, 0.0, 0.0, 0.0,
    3.4873, 0.0, 0.0, 3.3468, 0.0, 0.0, 3.2063, 0.0, 0.0, 3.0658, 0.0, 0.0,
    2.9254, 0.0, 0.0, 2.7849, 0.0, 0.0, 2.6444, 0.0, 0.0, 2.504, 0.0, 0.0,
    2.3635, 0.0, 0.0, 2.223, 0.0, 0.0, 2.0825, 0.0, 0.0, 1.9421, 0.0, 0.0,
    1.8016, 0.0, 0.0, 1.6611, 0.0, 0.0, 1.5207, 0.0, 0.0, 3.3802, 0.0, 0.0,
    1.2361, 0.0, 0.0, 1.0936, 0.0, 0.0, 0.9546, 0.0, 0.0, 0.82068, 3.0255, 0.0,
    0.0, 2.9062, 0.0, 0.0, 2.7869, 0.0, 0.0, 2.6676, 0.0, 0.0, 2.5483, 0.0, 0.0,
    2.429, 0.0, 0.0, 2.3097, 0.0, 0.0, 2.1904, 0.0, 0.0, 2.0711, 0.0, 0.0,
    1.9518, 0.0, 0.0, 1.8325, 0.0, 0.0, 1.7132, 0.0, 0.0, 1.5939, 0.0, 0.0,
    1.4746, 0.0, 0.0, 1.3553, 0.0, 0.0, 1.2361, 0.0, 0.0, 3.1168, 0.0, 0.0,
    0.99381, 0.0, 0.0, 0.87255, 0.0, 0.0, 0.75471, 0.0, 0.0, 0.0, 3.0255, 0.0,
    0.0, 2.9062, 0.0, 0.0, 2.7869, 0.0, 0.0, 2.6676, 0.0, 0.0, 2.5483, 0.0, 0.0,
    2.429, 0.0, 0.0, 2.3097, 0.0, 0.0, 2.1904, 0.0, 0.0, 2.0711, 0.0, 0.0,
    1.9518, 0.0, 0.0, 1.8325, 0.0, 0.0, 1.7132, 0.0, 0.0, 1.5939, 0.0, 0.0,
    1.4746, 0.0, 0.0, 1.3553, 0.0, 0.0, 1.2361, 0.0, 0.0, 3.1168, 0.0, 0.0,
    0.99381, 0.0, 0.0, 0.87255, 0.0, 0.0, 0.75471, 0.0, 0.0, 0.0, 3.0255, 0.0,
    0.0, 2.9062, 0.0, 0.0, 2.7869, 0.0, 0.0, 2.6676, 0.0, 0.0, 2.5483, 0.0, 0.0,
    2.429, 0.0, 0.0, 2.3097, 0.0, 0.0, 2.1904, 0.0, 0.0, 2.0711, 0.0, 0.0,
    1.9518, 0.0, 0.0, 1.8325, 0.0, 0.0, 1.7132, 0.0, 0.0, 1.5939, 0.0, 0.0,
    1.4746, 0.0, 0.0, 1.3553, 0.0, 0.0, 1.2361, 0.0, 0.0, 3.1168, 0.0, 0.0,
    0.99381, 0.0, 0.0, 0.87255, 0.0, 0.0, 0.75471, 2.5909, 0.0, 0.0, 2.4911, 0.0,
    0.0, 2.3913, 0.0, 0.0, 2.2915, 0.0, 0.0, 2.1916, 0.0, 0.0, 2.0918, 0.0, 0.0,
    1.992, 0.0, 0.0, 1.8922, 0.0, 0.0, 1.7924, 0.0, 0.0, 1.6925, 0.0, 0.0,
    1.5927, 0.0, 0.0, 1.4929, 0.0, 0.0, 1.3931, 0.0, 0.0, 1.2933, 0.0, 0.0,
    1.1934, 0.0, 0.0, 1.0936, 0.0, 0.0, 0.99381, 0.0, 0.0, 2.894, 0.0, 0.0,
    0.79051, 0.0, 0.0, 0.68873, 0.0, 0.0, 0.0, 2.5909, 0.0, 0.0, 2.4911, 0.0,
    0.0, 2.3913, 0.0, 0.0, 2.2915, 0.0, 0.0, 2.1916, 0.0, 0.0, 2.0918, 0.0, 0.0,
    1.992, 0.0, 0.0, 1.8922, 0.0, 0.0, 1.7924, 0.0, 0.0, 1.6925, 0.0, 0.0,
    1.5927, 0.0, 0.0, 1.4929, 0.0, 0.0, 1.3931, 0.0, 0.0, 1.2933, 0.0, 0.0,
    1.1934, 0.0, 0.0, 1.0936, 0.0, 0.0, 0.99381, 0.0, 0.0, 2.894, 0.0, 0.0,
    0.79051, 0.0, 0.0, 0.68873, 0.0, 0.0, 0.0, 2.5909, 0.0, 0.0, 2.4911, 0.0,
    0.0, 2.3913, 0.0, 0.0, 2.2915, 0.0, 0.0, 2.1916, 0.0, 0.0, 2.0918, 0.0, 0.0,
    1.992, 0.0, 0.0, 1.8922, 0.0, 0.0, 1.7924, 0.0, 0.0, 1.6925, 0.0, 0.0,
    1.5927, 0.0, 0.0, 1.4929, 0.0, 0.0, 1.3931, 0.0, 0.0, 1.2933, 0.0, 0.0,
    1.1934, 0.0, 0.0, 1.0936, 0.0, 0.0, 0.99381, 0.0, 0.0, 2.894, 0.0, 0.0,
    0.79051, 0.0, 0.0, 0.68873, 2.1853, 0.0, 0.0, 2.1033, 0.0, 0.0, 2.0212, 0.0,
    0.0, 1.9392, 0.0, 0.0, 1.8571, 0.0, 0.0, 1.7751, 0.0, 0.0, 1.693, 0.0, 0.0,
    1.611, 0.0, 0.0, 1.5289, 0.0, 0.0, 1.4469, 0.0, 0.0, 1.3648, 0.0, 0.0,
    1.2828, 0.0, 0.0, 1.2007, 0.0, 0.0, 1.1187, 0.0, 0.0, 1.0366, 0.0, 0.0,
    0.9546, 0.0, 0.0, 0.87255, 0.0, 0.0, 0.79051, 0.0, 0.0, 2.7085, 0.0, 0.0,
    0.62275, 0.0, 0.0, 0.0, 2.1853, 0.0, 0.0, 2.1033, 0.0, 0.0, 2.0212, 0.0, 0.0,
    1.9392, 0.0, 0.0, 1.8571, 0.0, 0.0, 1.7751, 0.0, 0.0, 1.693, 0.0, 0.0, 1.611,
    0.0, 0.0, 1.5289, 0.0, 0.0, 1.4469, 0.0, 0.0, 1.3648, 0.0, 0.0, 1.2828, 0.0,
    0.0, 1.2007, 0.0, 0.0, 1.1187, 0.0, 0.0, 1.0366, 0.0, 0.0, 0.9546, 0.0, 0.0,
    0.87255, 0.0, 0.0, 0.79051, 0.0, 0.0, 2.7085, 0.0, 0.0, 0.62275, 0.0, 0.0,
    0.0, 2.1853, 0.0, 0.0, 2.1033, 0.0, 0.0, 2.0212, 0.0, 0.0, 1.9392, 0.0, 0.0,
    1.8571, 0.0, 0.0, 1.7751, 0.0, 0.0, 1.693, 0.0, 0.0, 1.611, 0.0, 0.0, 1.5289,
    0.0, 0.0, 1.4469, 0.0, 0.0, 1.3648, 0.0, 0.0, 1.2828, 0.0, 0.0, 1.2007, 0.0,
    0.0, 1.1187, 0.0, 0.0, 1.0366, 0.0, 0.0, 0.9546, 0.0, 0.0, 0.87255, 0.0, 0.0,
    0.79051, 0.0, 0.0, 2.7085, 0.0, 0.0, 0.62275, 1.8103, 0.0, 0.0, 1.7444, 0.0,
    0.0, 1.6784, 0.0, 0.0, 1.6124, 0.0, 0.0, 1.5464, 0.0, 0.0, 1.4805, 0.0, 0.0,
    1.4145, 0.0, 0.0, 1.3485, 0.0, 0.0, 1.2825, 0.0, 0.0, 1.2165, 0.0, 0.0,
    1.1506, 0.0, 0.0, 1.0846, 0.0, 0.0, 1.0186, 0.0, 0.0, 0.95264, 0.0, 0.0,
    0.88666, 0.0, 0.0, 0.82068, 0.0, 0.0, 0.75471, 0.0, 0.0, 0.68873, 0.0, 0.0,
    0.62275, 0.0, 0.0, 2.5568, 0.0, 0.0, 0.0, 1.8103, 0.0, 0.0, 1.7444, 0.0, 0.0,
    1.6784, 0.0, 0.0, 1.6124, 0.0, 0.0, 1.5464, 0.0, 0.0, 1.4805, 0.0, 0.0,
    1.4145, 0.0, 0.0, 1.3485, 0.0, 0.0, 1.2825, 0.0, 0.0, 1.2165, 0.0, 0.0,
    1.1506, 0.0, 0.0, 1.0846, 0.0, 0.0, 1.0186, 0.0, 0.0, 0.95264, 0.0, 0.0,
    0.88666, 0.0, 0.0, 0.82068, 0.0, 0.0, 0.75471, 0.0, 0.0, 0.68873, 0.0, 0.0,
    0.62275, 0.0, 0.0, 2.5568, 0.0, 0.0, 0.0, 1.8103, 0.0, 0.0, 1.7444, 0.0, 0.0,
    1.6784, 0.0, 0.0, 1.6124, 0.0, 0.0, 1.5464, 0.0, 0.0, 1.4805, 0.0, 0.0,
    1.4145, 0.0, 0.0, 1.3485, 0.0, 0.0, 1.2825, 0.0, 0.0, 1.2165, 0.0, 0.0,
    1.1506, 0.0, 0.0, 1.0846, 0.0, 0.0, 1.0186, 0.0, 0.0, 0.95264, 0.0, 0.0,
    0.88666, 0.0, 0.0, 0.82068, 0.0, 0.0, 0.75471, 0.0, 0.0, 0.68873, 0.0, 0.0,
    0.62275, 0.0, 0.0, 2.5568 };

  static const double b_a[3600] = { -14.5253, -0.0, -0.0, -11.8592, -0.0, -0.0,
    -11.1948, -0.0, -0.0, -10.5338, -0.0, -0.0, -9.8779, -0.0, -0.0, -9.2289,
    -0.0, -0.0, -8.5883, -0.0, -0.0, -7.958, -0.0, -0.0, -7.3395, -0.0, -0.0,
    -6.7347, -0.0, -0.0, -6.1452, -0.0, -0.0, -5.5728, -0.0, -0.0, -5.0191, -0.0,
    -0.0, -4.4858, -0.0, -0.0, -3.9746, -0.0, -0.0, -3.4873, -0.0, -0.0, -3.0255,
    -0.0, -0.0, -2.5909, -0.0, -0.0, -2.1853, -0.0, -0.0, -1.8103, -0.0, -0.0,
    -0.0, -14.5253, -0.0, -0.0, -11.8592, -0.0, -0.0, -11.1948, -0.0, -0.0,
    -10.5338, -0.0, -0.0, -9.8779, -0.0, -0.0, -9.2289, -0.0, -0.0, -8.5883,
    -0.0, -0.0, -7.958, -0.0, -0.0, -7.3395, -0.0, -0.0, -6.7347, -0.0, -0.0,
    -6.1452, -0.0, -0.0, -5.5728, -0.0, -0.0, -5.0191, -0.0, -0.0, -4.4858, -0.0,
    -0.0, -3.9746, -0.0, -0.0, -3.4873, -0.0, -0.0, -3.0255, -0.0, -0.0, -2.5909,
    -0.0, -0.0, -2.1853, -0.0, -0.0, -1.8103, -0.0, -0.0, -0.0, -14.5253, -0.0,
    -0.0, -11.8592, -0.0, -0.0, -11.1948, -0.0, -0.0, -10.5338, -0.0, -0.0,
    -9.8779, -0.0, -0.0, -9.2289, -0.0, -0.0, -8.5883, -0.0, -0.0, -7.958, -0.0,
    -0.0, -7.3395, -0.0, -0.0, -6.7347, -0.0, -0.0, -6.1452, -0.0, -0.0, -5.5728,
    -0.0, -0.0, -5.0191, -0.0, -0.0, -4.4858, -0.0, -0.0, -3.9746, -0.0, -0.0,
    -3.4873, -0.0, -0.0, -3.0255, -0.0, -0.0, -2.5909, -0.0, -0.0, -2.1853, -0.0,
    -0.0, -1.8103, -11.8592, -0.0, -0.0, -13.2435, -0.0, -0.0, -10.6241, -0.0,
    -0.0, -10.0064, -0.0, -0.0, -9.3921, -0.0, -0.0, -8.783, -0.0, -0.0, -8.1806,
    -0.0, -0.0, -7.5868, -0.0, -0.0, -7.0031, -0.0, -0.0, -6.4314, -0.0, -0.0,
    -5.8733, -0.0, -0.0, -5.3306, -0.0, -0.0, -4.8048, -0.0, -0.0, -4.2978, -0.0,
    -0.0, -3.8112, -0.0, -0.0, -3.3468, -0.0, -0.0, -2.9062, -0.0, -0.0, -2.4911,
    -0.0, -0.0, -2.1033, -0.0, -0.0, -1.7444, -0.0, -0.0, -0.0, -11.8592, -0.0,
    -0.0, -13.2435, -0.0, -0.0, -10.6241, -0.0, -0.0, -10.0064, -0.0, -0.0,
    -9.3921, -0.0, -0.0, -8.783, -0.0, -0.0, -8.1806, -0.0, -0.0, -7.5868, -0.0,
    -0.0, -7.0031, -0.0, -0.0, -6.4314, -0.0, -0.0, -5.8733, -0.0, -0.0, -5.3306,
    -0.0, -0.0, -4.8048, -0.0, -0.0, -4.2978, -0.0, -0.0, -3.8112, -0.0, -0.0,
    -3.3468, -0.0, -0.0, -2.9062, -0.0, -0.0, -2.4911, -0.0, -0.0, -2.1033, -0.0,
    -0.0, -1.7444, -0.0, -0.0, -0.0, -11.8592, -0.0, -0.0, -13.2435, -0.0, -0.0,
    -10.6241, -0.0, -0.0, -10.0064, -0.0, -0.0, -9.3921, -0.0, -0.0, -8.783,
    -0.0, -0.0, -8.1806, -0.0, -0.0, -7.5868, -0.0, -0.0, -7.0031, -0.0, -0.0,
    -6.4314, -0.0, -0.0, -5.8733, -0.0, -0.0, -5.3306, -0.0, -0.0, -4.8048, -0.0,
    -0.0, -4.2978, -0.0, -0.0, -3.8112, -0.0, -0.0, -3.3468, -0.0, -0.0, -2.9062,
    -0.0, -0.0, -2.4911, -0.0, -0.0, -2.1033, -0.0, -0.0, -1.7444, -11.1948,
    -0.0, -0.0, -10.6241, -0.0, -0.0, -12.0534, -0.0, -0.0, -9.479, -0.0, -0.0,
    -8.9063, -0.0, -0.0, -8.3371, -0.0, -0.0, -7.7729, -0.0, -0.0, -7.2156, -0.0,
    -0.0, -6.6668, -0.0, -0.0, -6.1281, -0.0, -0.0, -5.6014, -0.0, -0.0, -5.0884,
    -0.0, -0.0, -4.5906, -0.0, -0.0, -4.1099, -0.0, -0.0, -3.6479, -0.0, -0.0,
    -3.2063, -0.0, -0.0, -2.7869, -0.0, -0.0, -2.3913, -0.0, -0.0, -2.0212, -0.0,
    -0.0, -1.6784, -0.0, -0.0, -0.0, -11.1948, -0.0, -0.0, -10.6241, -0.0, -0.0,
    -12.0534, -0.0, -0.0, -9.479, -0.0, -0.0, -8.9063, -0.0, -0.0, -8.3371, -0.0,
    -0.0, -7.7729, -0.0, -0.0, -7.2156, -0.0, -0.0, -6.6668, -0.0, -0.0, -6.1281,
    -0.0, -0.0, -5.6014, -0.0, -0.0, -5.0884, -0.0, -0.0, -4.5906, -0.0, -0.0,
    -4.1099, -0.0, -0.0, -3.6479, -0.0, -0.0, -3.2063, -0.0, -0.0, -2.7869, -0.0,
    -0.0, -2.3913, -0.0, -0.0, -2.0212, -0.0, -0.0, -1.6784, -0.0, -0.0, -0.0,
    -11.1948, -0.0, -0.0, -10.6241, -0.0, -0.0, -12.0534, -0.0, -0.0, -9.479,
    -0.0, -0.0, -8.9063, -0.0, -0.0, -8.3371, -0.0, -0.0, -7.7729, -0.0, -0.0,
    -7.2156, -0.0, -0.0, -6.6668, -0.0, -0.0, -6.1281, -0.0, -0.0, -5.6014, -0.0,
    -0.0, -5.0884, -0.0, -0.0, -4.5906, -0.0, -0.0, -4.1099, -0.0, -0.0, -3.6479,
    -0.0, -0.0, -3.2063, -0.0, -0.0, -2.7869, -0.0, -0.0, -2.3913, -0.0, -0.0,
    -2.0212, -0.0, -0.0, -1.6784, -10.5338, -0.0, -0.0, -10.0064, -0.0, -0.0,
    -9.479, -0.0, -0.0, -10.9516, -0.0, -0.0, -8.4206, -0.0, -0.0, -7.8912, -0.0,
    -0.0, -7.3653, -0.0, -0.0, -6.8444, -0.0, -0.0, -6.3304, -0.0, -0.0, -5.8249,
    -0.0, -0.0, -5.3296, -0.0, -0.0, -4.8462, -0.0, -0.0, -4.3764, -0.0, -0.0,
    -3.922, -0.0, -0.0, -3.4845, -0.0, -0.0, -3.0658, -0.0, -0.0, -2.6676, -0.0,
    -0.0, -2.2915, -0.0, -0.0, -1.9392, -0.0, -0.0, -1.6124, -0.0, -0.0, -0.0,
    -10.5338, -0.0, -0.0, -10.0064, -0.0, -0.0, -9.479, -0.0, -0.0, -10.9516,
    -0.0, -0.0, -8.4206, -0.0, -0.0, -7.8912, -0.0, -0.0, -7.3653, -0.0, -0.0,
    -6.8444, -0.0, -0.0, -6.3304, -0.0, -0.0, -5.8249, -0.0, -0.0, -5.3296, -0.0,
    -0.0, -4.8462, -0.0, -0.0, -4.3764, -0.0, -0.0, -3.922, -0.0, -0.0, -3.4845,
    -0.0, -0.0, -3.0658, -0.0, -0.0, -2.6676, -0.0, -0.0, -2.2915, -0.0, -0.0,
    -1.9392, -0.0, -0.0, -1.6124, -0.0, -0.0, -0.0, -10.5338, -0.0, -0.0,
    -10.0064, -0.0, -0.0, -9.479, -0.0, -0.0, -10.9516, -0.0, -0.0, -8.4206,
    -0.0, -0.0, -7.8912, -0.0, -0.0, -7.3653, -0.0, -0.0, -6.8444, -0.0, -0.0,
    -6.3304, -0.0, -0.0, -5.8249, -0.0, -0.0, -5.3296, -0.0, -0.0, -4.8462, -0.0,
    -0.0, -4.3764, -0.0, -0.0, -3.922, -0.0, -0.0, -3.4845, -0.0, -0.0, -3.0658,
    -0.0, -0.0, -2.6676, -0.0, -0.0, -2.2915, -0.0, -0.0, -1.9392, -0.0, -0.0,
    -1.6124, -9.8779, -0.0, -0.0, -9.3921, -0.0, -0.0, -8.9063, -0.0, -0.0,
    -8.4206, -0.0, -0.0, -9.9348, -0.0, -0.0, -7.4453, -0.0, -0.0, -6.9576, -0.0,
    -0.0, -6.4732, -0.0, -0.0, -5.994, -0.0, -0.0, -5.5216, -0.0, -0.0, -5.0577,
    -0.0, -0.0, -4.604, -0.0, -0.0, -4.1622, -0.0, -0.0, -3.734, -0.0, -0.0,
    -3.3212, -0.0, -0.0, -2.9254, -0.0, -0.0, -2.5483, -0.0, -0.0, -2.1916, -0.0,
    -0.0, -1.8571, -0.0, -0.0, -1.5464, -0.0, -0.0, -0.0, -9.8779, -0.0, -0.0,
    -9.3921, -0.0, -0.0, -8.9063, -0.0, -0.0, -8.4206, -0.0, -0.0, -9.9348, -0.0,
    -0.0, -7.4453, -0.0, -0.0, -6.9576, -0.0, -0.0, -6.4732, -0.0, -0.0, -5.994,
    -0.0, -0.0, -5.5216, -0.0, -0.0, -5.0577, -0.0, -0.0, -4.604, -0.0, -0.0,
    -4.1622, -0.0, -0.0, -3.734, -0.0, -0.0, -3.3212, -0.0, -0.0, -2.9254, -0.0,
    -0.0, -2.5483, -0.0, -0.0, -2.1916, -0.0, -0.0, -1.8571, -0.0, -0.0, -1.5464,
    -0.0, -0.0, -0.0, -9.8779, -0.0, -0.0, -9.3921, -0.0, -0.0, -8.9063, -0.0,
    -0.0, -8.4206, -0.0, -0.0, -9.9348, -0.0, -0.0, -7.4453, -0.0, -0.0, -6.9576,
    -0.0, -0.0, -6.4732, -0.0, -0.0, -5.994, -0.0, -0.0, -5.5216, -0.0, -0.0,
    -5.0577, -0.0, -0.0, -4.604, -0.0, -0.0, -4.1622, -0.0, -0.0, -3.734, -0.0,
    -0.0, -3.3212, -0.0, -0.0, -2.9254, -0.0, -0.0, -2.5483, -0.0, -0.0, -2.1916,
    -0.0, -0.0, -1.8571, -0.0, -0.0, -1.5464, -9.2289, -0.0, -0.0, -8.783, -0.0,
    -0.0, -8.3371, -0.0, -0.0, -7.8912, -0.0, -0.0, -7.4453, -0.0, -0.0, -8.9994,
    -0.0, -0.0, -6.5499, -0.0, -0.0, -6.102, -0.0, -0.0, -5.6576, -0.0, -0.0,
    -5.2183, -0.0, -0.0, -4.7858, -0.0, -0.0, -4.3618, -0.0, -0.0, -3.948, -0.0,
    -0.0, -3.5461, -0.0, -0.0, -3.1578, -0.0, -0.0, -2.7849, -0.0, -0.0, -2.429,
    -0.0, -0.0, -2.0918, -0.0, -0.0, -1.7751, -0.0, -0.0, -1.4805, -0.0, -0.0,
    -0.0, -9.2289, -0.0, -0.0, -8.783, -0.0, -0.0, -8.3371, -0.0, -0.0, -7.8912,
    -0.0, -0.0, -7.4453, -0.0, -0.0, -8.9994, -0.0, -0.0, -6.5499, -0.0, -0.0,
    -6.102, -0.0, -0.0, -5.6576, -0.0, -0.0, -5.2183, -0.0, -0.0, -4.7858, -0.0,
    -0.0, -4.3618, -0.0, -0.0, -3.948, -0.0, -0.0, -3.5461, -0.0, -0.0, -3.1578,
    -0.0, -0.0, -2.7849, -0.0, -0.0, -2.429, -0.0, -0.0, -2.0918, -0.0, -0.0,
    -1.7751, -0.0, -0.0, -1.4805, -0.0, -0.0, -0.0, -9.2289, -0.0, -0.0, -8.783,
    -0.0, -0.0, -8.3371, -0.0, -0.0, -7.8912, -0.0, -0.0, -7.4453, -0.0, -0.0,
    -8.9994, -0.0, -0.0, -6.5499, -0.0, -0.0, -6.102, -0.0, -0.0, -5.6576, -0.0,
    -0.0, -5.2183, -0.0, -0.0, -4.7858, -0.0, -0.0, -4.3618, -0.0, -0.0, -3.948,
    -0.0, -0.0, -3.5461, -0.0, -0.0, -3.1578, -0.0, -0.0, -2.7849, -0.0, -0.0,
    -2.429, -0.0, -0.0, -2.0918, -0.0, -0.0, -1.7751, -0.0, -0.0, -1.4805,
    -8.5883, -0.0, -0.0, -8.1806, -0.0, -0.0, -7.7729, -0.0, -0.0, -7.3653, -0.0,
    -0.0, -6.9576, -0.0, -0.0, -6.5499, -0.0, -0.0, -8.1422, -0.0, -0.0, -5.7309,
    -0.0, -0.0, -5.3212, -0.0, -0.0, -4.915, -0.0, -0.0, -4.5139, -0.0, -0.0,
    -4.1195, -0.0, -0.0, -3.7337, -0.0, -0.0, -3.3582, -0.0, -0.0, -2.9945, -0.0,
    -0.0, -2.6444, -0.0, -0.0, -2.3097, -0.0, -0.0, -1.992, -0.0, -0.0, -1.693,
    -0.0, -0.0, -1.4145, -0.0, -0.0, -0.0, -8.5883, -0.0, -0.0, -8.1806, -0.0,
    -0.0, -7.7729, -0.0, -0.0, -7.3653, -0.0, -0.0, -6.9576, -0.0, -0.0, -6.5499,
    -0.0, -0.0, -8.1422, -0.0, -0.0, -5.7309, -0.0, -0.0, -5.3212, -0.0, -0.0,
    -4.915, -0.0, -0.0, -4.5139, -0.0, -0.0, -4.1195, -0.0, -0.0, -3.7337, -0.0,
    -0.0, -3.3582, -0.0, -0.0, -2.9945, -0.0, -0.0, -2.6444, -0.0, -0.0, -2.3097,
    -0.0, -0.0, -1.992, -0.0, -0.0, -1.693, -0.0, -0.0, -1.4145, -0.0, -0.0,
    -0.0, -8.5883, -0.0, -0.0, -8.1806, -0.0, -0.0, -7.7729, -0.0, -0.0, -7.3653,
    -0.0, -0.0, -6.9576, -0.0, -0.0, -6.5499, -0.0, -0.0, -8.1422, -0.0, -0.0,
    -5.7309, -0.0, -0.0, -5.3212, -0.0, -0.0, -4.915, -0.0, -0.0, -4.5139, -0.0,
    -0.0, -4.1195, -0.0, -0.0, -3.7337, -0.0, -0.0, -3.3582, -0.0, -0.0, -2.9945,
    -0.0, -0.0, -2.6444, -0.0, -0.0, -2.3097, -0.0, -0.0, -1.992, -0.0, -0.0,
    -1.693, -0.0, -0.0, -1.4145, -7.958, -0.0, -0.0, -7.5868, -0.0, -0.0,
    -7.2156, -0.0, -0.0, -6.8444, -0.0, -0.0, -6.4732, -0.0, -0.0, -6.102, -0.0,
    -0.0, -5.7309, -0.0, -0.0, -7.3597, -0.0, -0.0, -4.9848, -0.0, -0.0, -4.6117,
    -0.0, -0.0, -4.242, -0.0, -0.0, -3.8773, -0.0, -0.0, -3.5195, -0.0, -0.0,
    -3.1702, -0.0, -0.0, -2.8311, -0.0, -0.0, -2.504, -0.0, -0.0, -2.1904, -0.0,
    -0.0, -1.8922, -0.0, -0.0, -1.611, -0.0, -0.0, -1.3485, -0.0, -0.0, -0.0,
    -7.958, -0.0, -0.0, -7.5868, -0.0, -0.0, -7.2156, -0.0, -0.0, -6.8444, -0.0,
    -0.0, -6.4732, -0.0, -0.0, -6.102, -0.0, -0.0, -5.7309, -0.0, -0.0, -7.3597,
    -0.0, -0.0, -4.9848, -0.0, -0.0, -4.6117, -0.0, -0.0, -4.242, -0.0, -0.0,
    -3.8773, -0.0, -0.0, -3.5195, -0.0, -0.0, -3.1702, -0.0, -0.0, -2.8311, -0.0,
    -0.0, -2.504, -0.0, -0.0, -2.1904, -0.0, -0.0, -1.8922, -0.0, -0.0, -1.611,
    -0.0, -0.0, -1.3485, -0.0, -0.0, -0.0, -7.958, -0.0, -0.0, -7.5868, -0.0,
    -0.0, -7.2156, -0.0, -0.0, -6.8444, -0.0, -0.0, -6.4732, -0.0, -0.0, -6.102,
    -0.0, -0.0, -5.7309, -0.0, -0.0, -7.3597, -0.0, -0.0, -4.9848, -0.0, -0.0,
    -4.6117, -0.0, -0.0, -4.242, -0.0, -0.0, -3.8773, -0.0, -0.0, -3.5195, -0.0,
    -0.0, -3.1702, -0.0, -0.0, -2.8311, -0.0, -0.0, -2.504, -0.0, -0.0, -2.1904,
    -0.0, -0.0, -1.8922, -0.0, -0.0, -1.611, -0.0, -0.0, -1.3485, -7.3395, -0.0,
    -0.0, -7.0031, -0.0, -0.0, -6.6668, -0.0, -0.0, -6.3304, -0.0, -0.0, -5.994,
    -0.0, -0.0, -5.6576, -0.0, -0.0, -5.3212, -0.0, -0.0, -4.9848, -0.0, -0.0,
    -6.6485, -0.0, -0.0, -4.3084, -0.0, -0.0, -3.9701, -0.0, -0.0, -3.6351, -0.0,
    -0.0, -3.3053, -0.0, -0.0, -2.9823, -0.0, -0.0, -2.6678, -0.0, -0.0, -2.3635,
    -0.0, -0.0, -2.0711, -0.0, -0.0, -1.7924, -0.0, -0.0, -1.5289, -0.0, -0.0,
    -1.2825, -0.0, -0.0, -0.0, -7.3395, -0.0, -0.0, -7.0031, -0.0, -0.0, -6.6668,
    -0.0, -0.0, -6.3304, -0.0, -0.0, -5.994, -0.0, -0.0, -5.6576, -0.0, -0.0,
    -5.3212, -0.0, -0.0, -4.9848, -0.0, -0.0, -6.6485, -0.0, -0.0, -4.3084, -0.0,
    -0.0, -3.9701, -0.0, -0.0, -3.6351, -0.0, -0.0, -3.3053, -0.0, -0.0, -2.9823,
    -0.0, -0.0, -2.6678, -0.0, -0.0, -2.3635, -0.0, -0.0, -2.0711, -0.0, -0.0,
    -1.7924, -0.0, -0.0, -1.5289, -0.0, -0.0, -1.2825, -0.0, -0.0, -0.0, -7.3395,
    -0.0, -0.0, -7.0031, -0.0, -0.0, -6.6668, -0.0, -0.0, -6.3304, -0.0, -0.0,
    -5.994, -0.0, -0.0, -5.6576, -0.0, -0.0, -5.3212, -0.0, -0.0, -4.9848, -0.0,
    -0.0, -6.6485, -0.0, -0.0, -4.3084, -0.0, -0.0, -3.9701, -0.0, -0.0, -3.6351,
    -0.0, -0.0, -3.3053, -0.0, -0.0, -2.9823, -0.0, -0.0, -2.6678, -0.0, -0.0,
    -2.3635, -0.0, -0.0, -2.0711, -0.0, -0.0, -1.7924, -0.0, -0.0, -1.5289, -0.0,
    -0.0, -1.2825, -6.7347, -0.0, -0.0, -6.4314, -0.0, -0.0, -6.1281, -0.0, -0.0,
    -5.8249, -0.0, -0.0, -5.5216, -0.0, -0.0, -5.2183, -0.0, -0.0, -4.915, -0.0,
    -0.0, -4.6117, -0.0, -0.0, -4.3084, -0.0, -0.0, -6.0051, -0.0, -0.0, -3.6982,
    -0.0, -0.0, -3.3929, -0.0, -0.0, -3.0911, -0.0, -0.0, -2.7944, -0.0, -0.0,
    -2.5044, -0.0, -0.0, -2.223, -0.0, -0.0, -1.9518, -0.0, -0.0, -1.6925, -0.0,
    -0.0, -1.4469, -0.0, -0.0, -1.2165, -0.0, -0.0, -0.0, -6.7347, -0.0, -0.0,
    -6.4314, -0.0, -0.0, -6.1281, -0.0, -0.0, -5.8249, -0.0, -0.0, -5.5216, -0.0,
    -0.0, -5.2183, -0.0, -0.0, -4.915, -0.0, -0.0, -4.6117, -0.0, -0.0, -4.3084,
    -0.0, -0.0, -6.0051, -0.0, -0.0, -3.6982, -0.0, -0.0, -3.3929, -0.0, -0.0,
    -3.0911, -0.0, -0.0, -2.7944, -0.0, -0.0, -2.5044, -0.0, -0.0, -2.223, -0.0,
    -0.0, -1.9518, -0.0, -0.0, -1.6925, -0.0, -0.0, -1.4469, -0.0, -0.0, -1.2165,
    -0.0, -0.0, -0.0, -6.7347, -0.0, -0.0, -6.4314, -0.0, -0.0, -6.1281, -0.0,
    -0.0, -5.8249, -0.0, -0.0, -5.5216, -0.0, -0.0, -5.2183, -0.0, -0.0, -4.915,
    -0.0, -0.0, -4.6117, -0.0, -0.0, -4.3084, -0.0, -0.0, -6.0051, -0.0, -0.0,
    -3.6982, -0.0, -0.0, -3.3929, -0.0, -0.0, -3.0911, -0.0, -0.0, -2.7944, -0.0,
    -0.0, -2.5044, -0.0, -0.0, -2.223, -0.0, -0.0, -1.9518, -0.0, -0.0, -1.6925,
    -0.0, -0.0, -1.4469, -0.0, -0.0, -1.2165, -6.1452, -0.0, -0.0, -5.8733, -0.0,
    -0.0, -5.6014, -0.0, -0.0, -5.3296, -0.0, -0.0, -5.0577, -0.0, -0.0, -4.7858,
    -0.0, -0.0, -4.5139, -0.0, -0.0, -4.242, -0.0, -0.0, -3.9701, -0.0, -0.0,
    -3.6982, -0.0, -0.0, -5.4263, -0.0, -0.0, -3.1507, -0.0, -0.0, -2.8769, -0.0,
    -0.0, -2.6064, -0.0, -0.0, -2.3411, -0.0, -0.0, -2.0825, -0.0, -0.0, -1.8325,
    -0.0, -0.0, -1.5927, -0.0, -0.0, -1.3648, -0.0, -0.0, -1.1506, -0.0, -0.0,
    -0.0, -6.1452, -0.0, -0.0, -5.8733, -0.0, -0.0, -5.6014, -0.0, -0.0, -5.3296,
    -0.0, -0.0, -5.0577, -0.0, -0.0, -4.7858, -0.0, -0.0, -4.5139, -0.0, -0.0,
    -4.242, -0.0, -0.0, -3.9701, -0.0, -0.0, -3.6982, -0.0, -0.0, -5.4263, -0.0,
    -0.0, -3.1507, -0.0, -0.0, -2.8769, -0.0, -0.0, -2.6064, -0.0, -0.0, -2.3411,
    -0.0, -0.0, -2.0825, -0.0, -0.0, -1.8325, -0.0, -0.0, -1.5927, -0.0, -0.0,
    -1.3648, -0.0, -0.0, -1.1506, -0.0, -0.0, -0.0, -6.1452, -0.0, -0.0, -5.8733,
    -0.0, -0.0, -5.6014, -0.0, -0.0, -5.3296, -0.0, -0.0, -5.0577, -0.0, -0.0,
    -4.7858, -0.0, -0.0, -4.5139, -0.0, -0.0, -4.242, -0.0, -0.0, -3.9701, -0.0,
    -0.0, -3.6982, -0.0, -0.0, -5.4263, -0.0, -0.0, -3.1507, -0.0, -0.0, -2.8769,
    -0.0, -0.0, -2.6064, -0.0, -0.0, -2.3411, -0.0, -0.0, -2.0825, -0.0, -0.0,
    -1.8325, -0.0, -0.0, -1.5927, -0.0, -0.0, -1.3648, -0.0, -0.0, -1.1506,
    -5.5728, -0.0, -0.0, -5.3306, -0.0, -0.0, -5.0884, -0.0, -0.0, -4.8462, -0.0,
    -0.0, -4.604, -0.0, -0.0, -4.3618, -0.0, -0.0, -4.1195, -0.0, -0.0, -3.8773,
    -0.0, -0.0, -3.6351, -0.0, -0.0, -3.3929, -0.0, -0.0, -3.1507, -0.0, -0.0,
    -4.9085, -0.0, -0.0, -2.6626, -0.0, -0.0, -2.4185, -0.0, -0.0, -2.1777, -0.0,
    -0.0, -1.9421, -0.0, -0.0, -1.7132, -0.0, -0.0, -1.4929, -0.0, -0.0, -1.2828,
    -0.0, -0.0, -1.0846, -0.0, -0.0, -0.0, -5.5728, -0.0, -0.0, -5.3306, -0.0,
    -0.0, -5.0884, -0.0, -0.0, -4.8462, -0.0, -0.0, -4.604, -0.0, -0.0, -4.3618,
    -0.0, -0.0, -4.1195, -0.0, -0.0, -3.8773, -0.0, -0.0, -3.6351, -0.0, -0.0,
    -3.3929, -0.0, -0.0, -3.1507, -0.0, -0.0, -4.9085, -0.0, -0.0, -2.6626, -0.0,
    -0.0, -2.4185, -0.0, -0.0, -2.1777, -0.0, -0.0, -1.9421, -0.0, -0.0, -1.7132,
    -0.0, -0.0, -1.4929, -0.0, -0.0, -1.2828, -0.0, -0.0, -1.0846, -0.0, -0.0,
    -0.0, -5.5728, -0.0, -0.0, -5.3306, -0.0, -0.0, -5.0884, -0.0, -0.0, -4.8462,
    -0.0, -0.0, -4.604, -0.0, -0.0, -4.3618, -0.0, -0.0, -4.1195, -0.0, -0.0,
    -3.8773, -0.0, -0.0, -3.6351, -0.0, -0.0, -3.3929, -0.0, -0.0, -3.1507, -0.0,
    -0.0, -4.9085, -0.0, -0.0, -2.6626, -0.0, -0.0, -2.4185, -0.0, -0.0, -2.1777,
    -0.0, -0.0, -1.9421, -0.0, -0.0, -1.7132, -0.0, -0.0, -1.4929, -0.0, -0.0,
    -1.2828, -0.0, -0.0, -1.0846, -5.0191, -0.0, -0.0, -4.8048, -0.0, -0.0,
    -4.5906, -0.0, -0.0, -4.3764, -0.0, -0.0, -4.1622, -0.0, -0.0, -3.948, -0.0,
    -0.0, -3.7337, -0.0, -0.0, -3.5195, -0.0, -0.0, -3.3053, -0.0, -0.0, -3.0911,
    -0.0, -0.0, -2.8769, -0.0, -0.0, -2.6626, -0.0, -0.0, -4.4484, -0.0, -0.0,
    -2.2305, -0.0, -0.0, -2.0144, -0.0, -0.0, -1.8016, -0.0, -0.0, -1.5939, -0.0,
    -0.0, -1.3931, -0.0, -0.0, -1.2007, -0.0, -0.0, -1.0186, -0.0, -0.0, -0.0,
    -5.0191, -0.0, -0.0, -4.8048, -0.0, -0.0, -4.5906, -0.0, -0.0, -4.3764, -0.0,
    -0.0, -4.1622, -0.0, -0.0, -3.948, -0.0, -0.0, -3.7337, -0.0, -0.0, -3.5195,
    -0.0, -0.0, -3.3053, -0.0, -0.0, -3.0911, -0.0, -0.0, -2.8769, -0.0, -0.0,
    -2.6626, -0.0, -0.0, -4.4484, -0.0, -0.0, -2.2305, -0.0, -0.0, -2.0144, -0.0,
    -0.0, -1.8016, -0.0, -0.0, -1.5939, -0.0, -0.0, -1.3931, -0.0, -0.0, -1.2007,
    -0.0, -0.0, -1.0186, -0.0, -0.0, -0.0, -5.0191, -0.0, -0.0, -4.8048, -0.0,
    -0.0, -4.5906, -0.0, -0.0, -4.3764, -0.0, -0.0, -4.1622, -0.0, -0.0, -3.948,
    -0.0, -0.0, -3.7337, -0.0, -0.0, -3.5195, -0.0, -0.0, -3.3053, -0.0, -0.0,
    -3.0911, -0.0, -0.0, -2.8769, -0.0, -0.0, -2.6626, -0.0, -0.0, -4.4484, -0.0,
    -0.0, -2.2305, -0.0, -0.0, -2.0144, -0.0, -0.0, -1.8016, -0.0, -0.0, -1.5939,
    -0.0, -0.0, -1.3931, -0.0, -0.0, -1.2007, -0.0, -0.0, -1.0186, -4.4858, -0.0,
    -0.0, -4.2978, -0.0, -0.0, -4.1099, -0.0, -0.0, -3.922, -0.0, -0.0, -3.734,
    -0.0, -0.0, -3.5461, -0.0, -0.0, -3.3582, -0.0, -0.0, -3.1702, -0.0, -0.0,
    -2.9823, -0.0, -0.0, -2.7944, -0.0, -0.0, -2.6064, -0.0, -0.0, -2.4185, -0.0,
    -0.0, -2.2305, -0.0, -0.0, -4.0426, -0.0, -0.0, -1.851, -0.0, -0.0, -1.6611,
    -0.0, -0.0, -1.4746, -0.0, -0.0, -1.2933, -0.0, -0.0, -1.1187, -0.0, -0.0,
    -0.95264, -0.0, -0.0, -0.0, -4.4858, -0.0, -0.0, -4.2978, -0.0, -0.0,
    -4.1099, -0.0, -0.0, -3.922, -0.0, -0.0, -3.734, -0.0, -0.0, -3.5461, -0.0,
    -0.0, -3.3582, -0.0, -0.0, -3.1702, -0.0, -0.0, -2.9823, -0.0, -0.0, -2.7944,
    -0.0, -0.0, -2.6064, -0.0, -0.0, -2.4185, -0.0, -0.0, -2.2305, -0.0, -0.0,
    -4.0426, -0.0, -0.0, -1.851, -0.0, -0.0, -1.6611, -0.0, -0.0, -1.4746, -0.0,
    -0.0, -1.2933, -0.0, -0.0, -1.1187, -0.0, -0.0, -0.95264, -0.0, -0.0, -0.0,
    -4.4858, -0.0, -0.0, -4.2978, -0.0, -0.0, -4.1099, -0.0, -0.0, -3.922, -0.0,
    -0.0, -3.734, -0.0, -0.0, -3.5461, -0.0, -0.0, -3.3582, -0.0, -0.0, -3.1702,
    -0.0, -0.0, -2.9823, -0.0, -0.0, -2.7944, -0.0, -0.0, -2.6064, -0.0, -0.0,
    -2.4185, -0.0, -0.0, -2.2305, -0.0, -0.0, -4.0426, -0.0, -0.0, -1.851, -0.0,
    -0.0, -1.6611, -0.0, -0.0, -1.4746, -0.0, -0.0, -1.2933, -0.0, -0.0, -1.1187,
    -0.0, -0.0, -0.95264, -3.9746, -0.0, -0.0, -3.8112, -0.0, -0.0, -3.6479,
    -0.0, -0.0, -3.4845, -0.0, -0.0, -3.3212, -0.0, -0.0, -3.1578, -0.0, -0.0,
    -2.9945, -0.0, -0.0, -2.8311, -0.0, -0.0, -2.6678, -0.0, -0.0, -2.5044, -0.0,
    -0.0, -2.3411, -0.0, -0.0, -2.1777, -0.0, -0.0, -2.0144, -0.0, -0.0, -1.851,
    -0.0, -0.0, -3.6877, -0.0, -0.0, -1.5207, -0.0, -0.0, -1.3553, -0.0, -0.0,
    -1.1934, -0.0, -0.0, -1.0366, -0.0, -0.0, -0.88666, -0.0, -0.0, -0.0,
    -3.9746, -0.0, -0.0, -3.8112, -0.0, -0.0, -3.6479, -0.0, -0.0, -3.4845, -0.0,
    -0.0, -3.3212, -0.0, -0.0, -3.1578, -0.0, -0.0, -2.9945, -0.0, -0.0, -2.8311,
    -0.0, -0.0, -2.6678, -0.0, -0.0, -2.5044, -0.0, -0.0, -2.3411, -0.0, -0.0,
    -2.1777, -0.0, -0.0, -2.0144, -0.0, -0.0, -1.851, -0.0, -0.0, -3.6877, -0.0,
    -0.0, -1.5207, -0.0, -0.0, -1.3553, -0.0, -0.0, -1.1934, -0.0, -0.0, -1.0366,
    -0.0, -0.0, -0.88666, -0.0, -0.0, -0.0, -3.9746, -0.0, -0.0, -3.8112, -0.0,
    -0.0, -3.6479, -0.0, -0.0, -3.4845, -0.0, -0.0, -3.3212, -0.0, -0.0, -3.1578,
    -0.0, -0.0, -2.9945, -0.0, -0.0, -2.8311, -0.0, -0.0, -2.6678, -0.0, -0.0,
    -2.5044, -0.0, -0.0, -2.3411, -0.0, -0.0, -2.1777, -0.0, -0.0, -2.0144, -0.0,
    -0.0, -1.851, -0.0, -0.0, -3.6877, -0.0, -0.0, -1.5207, -0.0, -0.0, -1.3553,
    -0.0, -0.0, -1.1934, -0.0, -0.0, -1.0366, -0.0, -0.0, -0.88666, -3.4873,
    -0.0, -0.0, -3.3468, -0.0, -0.0, -3.2063, -0.0, -0.0, -3.0658, -0.0, -0.0,
    -2.9254, -0.0, -0.0, -2.7849, -0.0, -0.0, -2.6444, -0.0, -0.0, -2.504, -0.0,
    -0.0, -2.3635, -0.0, -0.0, -2.223, -0.0, -0.0, -2.0825, -0.0, -0.0, -1.9421,
    -0.0, -0.0, -1.8016, -0.0, -0.0, -1.6611, -0.0, -0.0, -1.5207, -0.0, -0.0,
    -3.3802, -0.0, -0.0, -1.2361, -0.0, -0.0, -1.0936, -0.0, -0.0, -0.9546, -0.0,
    -0.0, -0.82068, -0.0, -0.0, -0.0, -3.4873, -0.0, -0.0, -3.3468, -0.0, -0.0,
    -3.2063, -0.0, -0.0, -3.0658, -0.0, -0.0, -2.9254, -0.0, -0.0, -2.7849, -0.0,
    -0.0, -2.6444, -0.0, -0.0, -2.504, -0.0, -0.0, -2.3635, -0.0, -0.0, -2.223,
    -0.0, -0.0, -2.0825, -0.0, -0.0, -1.9421, -0.0, -0.0, -1.8016, -0.0, -0.0,
    -1.6611, -0.0, -0.0, -1.5207, -0.0, -0.0, -3.3802, -0.0, -0.0, -1.2361, -0.0,
    -0.0, -1.0936, -0.0, -0.0, -0.9546, -0.0, -0.0, -0.82068, -0.0, -0.0, -0.0,
    -3.4873, -0.0, -0.0, -3.3468, -0.0, -0.0, -3.2063, -0.0, -0.0, -3.0658, -0.0,
    -0.0, -2.9254, -0.0, -0.0, -2.7849, -0.0, -0.0, -2.6444, -0.0, -0.0, -2.504,
    -0.0, -0.0, -2.3635, -0.0, -0.0, -2.223, -0.0, -0.0, -2.0825, -0.0, -0.0,
    -1.9421, -0.0, -0.0, -1.8016, -0.0, -0.0, -1.6611, -0.0, -0.0, -1.5207, -0.0,
    -0.0, -3.3802, -0.0, -0.0, -1.2361, -0.0, -0.0, -1.0936, -0.0, -0.0, -0.9546,
    -0.0, -0.0, -0.82068, -3.0255, -0.0, -0.0, -2.9062, -0.0, -0.0, -2.7869,
    -0.0, -0.0, -2.6676, -0.0, -0.0, -2.5483, -0.0, -0.0, -2.429, -0.0, -0.0,
    -2.3097, -0.0, -0.0, -2.1904, -0.0, -0.0, -2.0711, -0.0, -0.0, -1.9518, -0.0,
    -0.0, -1.8325, -0.0, -0.0, -1.7132, -0.0, -0.0, -1.5939, -0.0, -0.0, -1.4746,
    -0.0, -0.0, -1.3553, -0.0, -0.0, -1.2361, -0.0, -0.0, -3.1168, -0.0, -0.0,
    -0.99381, -0.0, -0.0, -0.87255, -0.0, -0.0, -0.75471, -0.0, -0.0, -0.0,
    -3.0255, -0.0, -0.0, -2.9062, -0.0, -0.0, -2.7869, -0.0, -0.0, -2.6676, -0.0,
    -0.0, -2.5483, -0.0, -0.0, -2.429, -0.0, -0.0, -2.3097, -0.0, -0.0, -2.1904,
    -0.0, -0.0, -2.0711, -0.0, -0.0, -1.9518, -0.0, -0.0, -1.8325, -0.0, -0.0,
    -1.7132, -0.0, -0.0, -1.5939, -0.0, -0.0, -1.4746, -0.0, -0.0, -1.3553, -0.0,
    -0.0, -1.2361, -0.0, -0.0, -3.1168, -0.0, -0.0, -0.99381, -0.0, -0.0,
    -0.87255, -0.0, -0.0, -0.75471, -0.0, -0.0, -0.0, -3.0255, -0.0, -0.0,
    -2.9062, -0.0, -0.0, -2.7869, -0.0, -0.0, -2.6676, -0.0, -0.0, -2.5483, -0.0,
    -0.0, -2.429, -0.0, -0.0, -2.3097, -0.0, -0.0, -2.1904, -0.0, -0.0, -2.0711,
    -0.0, -0.0, -1.9518, -0.0, -0.0, -1.8325, -0.0, -0.0, -1.7132, -0.0, -0.0,
    -1.5939, -0.0, -0.0, -1.4746, -0.0, -0.0, -1.3553, -0.0, -0.0, -1.2361, -0.0,
    -0.0, -3.1168, -0.0, -0.0, -0.99381, -0.0, -0.0, -0.87255, -0.0, -0.0,
    -0.75471, -2.5909, -0.0, -0.0, -2.4911, -0.0, -0.0, -2.3913, -0.0, -0.0,
    -2.2915, -0.0, -0.0, -2.1916, -0.0, -0.0, -2.0918, -0.0, -0.0, -1.992, -0.0,
    -0.0, -1.8922, -0.0, -0.0, -1.7924, -0.0, -0.0, -1.6925, -0.0, -0.0, -1.5927,
    -0.0, -0.0, -1.4929, -0.0, -0.0, -1.3931, -0.0, -0.0, -1.2933, -0.0, -0.0,
    -1.1934, -0.0, -0.0, -1.0936, -0.0, -0.0, -0.99381, -0.0, -0.0, -2.894, -0.0,
    -0.0, -0.79051, -0.0, -0.0, -0.68873, -0.0, -0.0, -0.0, -2.5909, -0.0, -0.0,
    -2.4911, -0.0, -0.0, -2.3913, -0.0, -0.0, -2.2915, -0.0, -0.0, -2.1916, -0.0,
    -0.0, -2.0918, -0.0, -0.0, -1.992, -0.0, -0.0, -1.8922, -0.0, -0.0, -1.7924,
    -0.0, -0.0, -1.6925, -0.0, -0.0, -1.5927, -0.0, -0.0, -1.4929, -0.0, -0.0,
    -1.3931, -0.0, -0.0, -1.2933, -0.0, -0.0, -1.1934, -0.0, -0.0, -1.0936, -0.0,
    -0.0, -0.99381, -0.0, -0.0, -2.894, -0.0, -0.0, -0.79051, -0.0, -0.0,
    -0.68873, -0.0, -0.0, -0.0, -2.5909, -0.0, -0.0, -2.4911, -0.0, -0.0,
    -2.3913, -0.0, -0.0, -2.2915, -0.0, -0.0, -2.1916, -0.0, -0.0, -2.0918, -0.0,
    -0.0, -1.992, -0.0, -0.0, -1.8922, -0.0, -0.0, -1.7924, -0.0, -0.0, -1.6925,
    -0.0, -0.0, -1.5927, -0.0, -0.0, -1.4929, -0.0, -0.0, -1.3931, -0.0, -0.0,
    -1.2933, -0.0, -0.0, -1.1934, -0.0, -0.0, -1.0936, -0.0, -0.0, -0.99381,
    -0.0, -0.0, -2.894, -0.0, -0.0, -0.79051, -0.0, -0.0, -0.68873, -2.1853,
    -0.0, -0.0, -2.1033, -0.0, -0.0, -2.0212, -0.0, -0.0, -1.9392, -0.0, -0.0,
    -1.8571, -0.0, -0.0, -1.7751, -0.0, -0.0, -1.693, -0.0, -0.0, -1.611, -0.0,
    -0.0, -1.5289, -0.0, -0.0, -1.4469, -0.0, -0.0, -1.3648, -0.0, -0.0, -1.2828,
    -0.0, -0.0, -1.2007, -0.0, -0.0, -1.1187, -0.0, -0.0, -1.0366, -0.0, -0.0,
    -0.9546, -0.0, -0.0, -0.87255, -0.0, -0.0, -0.79051, -0.0, -0.0, -2.7085,
    -0.0, -0.0, -0.62275, -0.0, -0.0, -0.0, -2.1853, -0.0, -0.0, -2.1033, -0.0,
    -0.0, -2.0212, -0.0, -0.0, -1.9392, -0.0, -0.0, -1.8571, -0.0, -0.0, -1.7751,
    -0.0, -0.0, -1.693, -0.0, -0.0, -1.611, -0.0, -0.0, -1.5289, -0.0, -0.0,
    -1.4469, -0.0, -0.0, -1.3648, -0.0, -0.0, -1.2828, -0.0, -0.0, -1.2007, -0.0,
    -0.0, -1.1187, -0.0, -0.0, -1.0366, -0.0, -0.0, -0.9546, -0.0, -0.0,
    -0.87255, -0.0, -0.0, -0.79051, -0.0, -0.0, -2.7085, -0.0, -0.0, -0.62275,
    -0.0, -0.0, -0.0, -2.1853, -0.0, -0.0, -2.1033, -0.0, -0.0, -2.0212, -0.0,
    -0.0, -1.9392, -0.0, -0.0, -1.8571, -0.0, -0.0, -1.7751, -0.0, -0.0, -1.693,
    -0.0, -0.0, -1.611, -0.0, -0.0, -1.5289, -0.0, -0.0, -1.4469, -0.0, -0.0,
    -1.3648, -0.0, -0.0, -1.2828, -0.0, -0.0, -1.2007, -0.0, -0.0, -1.1187, -0.0,
    -0.0, -1.0366, -0.0, -0.0, -0.9546, -0.0, -0.0, -0.87255, -0.0, -0.0,
    -0.79051, -0.0, -0.0, -2.7085, -0.0, -0.0, -0.62275, -1.8103, -0.0, -0.0,
    -1.7444, -0.0, -0.0, -1.6784, -0.0, -0.0, -1.6124, -0.0, -0.0, -1.5464, -0.0,
    -0.0, -1.4805, -0.0, -0.0, -1.4145, -0.0, -0.0, -1.3485, -0.0, -0.0, -1.2825,
    -0.0, -0.0, -1.2165, -0.0, -0.0, -1.1506, -0.0, -0.0, -1.0846, -0.0, -0.0,
    -1.0186, -0.0, -0.0, -0.95264, -0.0, -0.0, -0.88666, -0.0, -0.0, -0.82068,
    -0.0, -0.0, -0.75471, -0.0, -0.0, -0.68873, -0.0, -0.0, -0.62275, -0.0, -0.0,
    -2.5568, -0.0, -0.0, -0.0, -1.8103, -0.0, -0.0, -1.7444, -0.0, -0.0, -1.6784,
    -0.0, -0.0, -1.6124, -0.0, -0.0, -1.5464, -0.0, -0.0, -1.4805, -0.0, -0.0,
    -1.4145, -0.0, -0.0, -1.3485, -0.0, -0.0, -1.2825, -0.0, -0.0, -1.2165, -0.0,
    -0.0, -1.1506, -0.0, -0.0, -1.0846, -0.0, -0.0, -1.0186, -0.0, -0.0,
    -0.95264, -0.0, -0.0, -0.88666, -0.0, -0.0, -0.82068, -0.0, -0.0, -0.75471,
    -0.0, -0.0, -0.68873, -0.0, -0.0, -0.62275, -0.0, -0.0, -2.5568, -0.0, -0.0,
    -0.0, -1.8103, -0.0, -0.0, -1.7444, -0.0, -0.0, -1.6784, -0.0, -0.0, -1.6124,
    -0.0, -0.0, -1.5464, -0.0, -0.0, -1.4805, -0.0, -0.0, -1.4145, -0.0, -0.0,
    -1.3485, -0.0, -0.0, -1.2825, -0.0, -0.0, -1.2165, -0.0, -0.0, -1.1506, -0.0,
    -0.0, -1.0846, -0.0, -0.0, -1.0186, -0.0, -0.0, -0.95264, -0.0, -0.0,
    -0.88666, -0.0, -0.0, -0.82068, -0.0, -0.0, -0.75471, -0.0, -0.0, -0.68873,
    -0.0, -0.0, -0.62275, -0.0, -0.0, -2.5568 };

  static const double b[360] = { 25.3845, 0.0, 0.0, 246.3304, 0.0, 0.0, 0.0,
    25.3845, 0.0, 0.0, 246.3304, 0.0, 0.0, 0.0, 25.3845, 0.0, 0.0, 246.3304,
    23.5944, 0.0, 0.0, 233.1205, 0.0, 0.0, 0.0, 23.5944, 0.0, 0.0, 233.1205, 0.0,
    0.0, 0.0, 23.5944, 0.0, 0.0, 233.1205, 21.8695, 0.0, 0.0, 219.9595, 0.0, 0.0,
    0.0, 21.8695, 0.0, 0.0, 219.9595, 0.0, 0.0, 0.0, 21.8695, 0.0, 0.0, 219.9595,
    20.2098, 0.0, 0.0, 206.8801, 0.0, 0.0, 0.0, 20.2098, 0.0, 0.0, 206.8801, 0.0,
    0.0, 0.0, 20.2098, 0.0, 0.0, 206.8801, 18.6154, 0.0, 0.0, 193.9148, 0.0, 0.0,
    0.0, 18.6154, 0.0, 0.0, 193.9148, 0.0, 0.0, 0.0, 18.6154, 0.0, 0.0, 193.9148,
    17.0863, 0.0, 0.0, 181.0964, 0.0, 0.0, 0.0, 17.0863, 0.0, 0.0, 181.0964, 0.0,
    0.0, 0.0, 17.0863, 0.0, 0.0, 181.0964, 15.6224, 0.0, 0.0, 168.4573, 0.0, 0.0,
    0.0, 15.6224, 0.0, 0.0, 168.4573, 0.0, 0.0, 0.0, 15.6224, 0.0, 0.0, 168.4573,
    14.2237, 0.0, 0.0, 156.0303, 0.0, 0.0, 0.0, 14.2237, 0.0, 0.0, 156.0303, 0.0,
    0.0, 0.0, 14.2237, 0.0, 0.0, 156.0303, 12.8902, 0.0, 0.0, 143.8479, 0.0, 0.0,
    0.0, 12.8902, 0.0, 0.0, 143.8479, 0.0, 0.0, 0.0, 12.8902, 0.0, 0.0, 143.8479,
    11.622, 0.0, 0.0, 131.9428, 0.0, 0.0, 0.0, 11.622, 0.0, 0.0, 131.9428, 0.0,
    0.0, 0.0, 11.622, 0.0, 0.0, 131.9428, 10.4191, 0.0, 0.0, 120.3476, 0.0, 0.0,
    0.0, 10.4191, 0.0, 0.0, 120.3476, 0.0, 0.0, 0.0, 10.4191, 0.0, 0.0, 120.3476,
    9.2814, 0.0, 0.0, 109.0949, 0.0, 0.0, 0.0, 9.2814, 0.0, 0.0, 109.0949, 0.0,
    0.0, 0.0, 9.2814, 0.0, 0.0, 109.0949, 8.2089, 0.0, 0.0, 98.2173, 0.0, 0.0,
    0.0, 8.2089, 0.0, 0.0, 98.2173, 0.0, 0.0, 0.0, 8.2089, 0.0, 0.0, 98.2173,
    7.2016, 0.0, 0.0, 87.7475, 0.0, 0.0, 0.0, 7.2016, 0.0, 0.0, 87.7475, 0.0,
    0.0, 0.0, 7.2016, 0.0, 0.0, 87.7475, 6.2596, 0.0, 0.0, 77.7181, 0.0, 0.0,
    0.0, 6.2596, 0.0, 0.0, 77.7181, 0.0, 0.0, 0.0, 6.2596, 0.0, 0.0, 77.7181,
    5.3829, 0.0, 0.0, 68.1616, 0.0, 0.0, 0.0, 5.3829, 0.0, 0.0, 68.1616, 0.0,
    0.0, 0.0, 5.3829, 0.0, 0.0, 68.1616, 4.5713, 0.0, 0.0, 59.1108, 0.0, 0.0,
    0.0, 4.5713, 0.0, 0.0, 59.1108, 0.0, 0.0, 0.0, 4.5713, 0.0, 0.0, 59.1108,
    3.8251, 0.0, 0.0, 50.5981, 0.0, 0.0, 0.0, 3.8251, 0.0, 0.0, 50.5981, 0.0,
    0.0, 0.0, 3.8251, 0.0, 0.0, 50.5981, 3.144, 0.0, 0.0, 42.6564, 0.0, 0.0, 0.0,
    3.144, 0.0, 0.0, 42.6564, 0.0, 0.0, 0.0, 3.144, 0.0, 0.0, 42.6564, 2.5282,
    0.0, 0.0, 35.3181, 0.0, 0.0, 0.0, 2.5282, 0.0, 0.0, 35.3181, 0.0, 0.0, 0.0,
    2.5282, 0.0, 0.0, 35.3181 };

  static const signed char A[7200] = { -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
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

  static const signed char a[7200] = { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
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

  static const signed char c_a[7200] = { -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
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

  double IGA[7200];
  double A_data[3600];
  double y_tmp[3600];
  double b_del_lam[240];
  double result[240];
  double varargin_1_data[240];
  double del_lam[120];
  double esig[120];
  double ftol[120];
  double igr2[120];
  double ilam[120];
  double lam[120];
  double mesil[120];
  double nlamt[120];
  double H[60];
  double del_z[60];
  double r1[60];
  double x[60];
  double d;
  double mu;
  double mu_old;
  double ssq;
  int b_i;
  int exitflag;
  int i;
  int idx;
  int iy;
  int unusedU0;
  unsigned char ii_data[240];
  bool b_x[240];

  //   QP Solver
  //  Call QP Solver
  //  Solve quadratic programming problem using Wright's (1997) Method
  //  Minimise J(x) = 1/2x'Hx + f'x
  //  Subject to: Ax <= b
  //  Reference: S. J. Wright, "Applying New Optimization Algorithms to Model
  //  Predictive Control," in Chemical Process Control-V, CACHE, AIChE
  //  Symposium, 1997, pp. 147-155.
  // Number of decision variables
  //  p = 0;
  // Test for Cold Start
  // Warm Start
  for (i = 0; i < 60; i++) {
    d = 0.0;
    for (b_i = 0; b_i < 6; b_i++) {
      d += (2.0 * x0_vl[b_i]) * b[b_i + (6 * i)];
    }

    H[i] = d;
    x[i] = X_QP_vl[i];
  }

  // to tune
  // to tune
  // Default Values
  for (i = 0; i < 120; i++) {
    lam[i] = 100.0;
    ftol[i] = 100.0;
    esig[i] = 0.001;
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
    int i1;
    int ii_size_idx_0;
    int info;
    int j;
    bool exitg1;
    bool y;

    // Create common matrices
    for (i = 0; i < 120; i++) {
      d = lam[i];
      ssq = 1.0 / d;
      ilam[i] = ssq;
      nlamt[i] = (-d) / ftol[i];
      mesil[i] = (mu * esig[i]) * ssq;
    }

    // RHS
    for (i = 0; i < 60; i++) {
      for (b_i = 0; b_i < 120; b_i++) {
        idx = b_i + (120 * i);
        IGA[idx] = nlamt[b_i] * (static_cast<double>(A[idx]));
      }

      d = 0.0;
      for (b_i = 0; b_i < 60; b_i++) {
        d += b_a[i + (60 * b_i)] * x[b_i];
      }

      ssq = 0.0;
      for (b_i = 0; b_i < 120; b_i++) {
        ssq += (static_cast<double>(c_a[i + (60 * b_i)])) * lam[b_i];
      }

      r1[i] = (d - ssq) - H[i];
    }

    for (b_i = 0; b_i < 120; b_i++) {
      d = 0.0;
      for (i1 = 0; i1 < 60; i1++) {
        d += (static_cast<double>(a[b_i + (120 * i1)])) * x[i1];
      }

      igr2[b_i] = nlamt[b_i] * ((d + 0.15) - mesil[b_i]);
    }

    // Solve
    for (b_i = 0; b_i < 60; b_i++) {
      for (i1 = 0; i1 < 60; i1++) {
        d = 0.0;
        for (i = 0; i < 120; i++) {
          d += (static_cast<double>(c_a[b_i + (60 * i)])) * IGA[i + (120 * i1)];
        }

        y_tmp[b_i + (60 * i1)] = d;
      }
    }

    for (b_i = 0; b_i < 3600; b_i++) {
      A_data[b_i] = b_H[b_i] - y_tmp[b_i];
    }

    info = -1;
    j = 0;
    exitg1 = false;
    while ((!exitg1) && (j < 60)) {
      int idxA1j;
      int idxAjj;
      int ix;
      idxA1j = j * 60;
      idxAjj = idxA1j + j;
      ssq = 0.0;
      if (j >= 1) {
        ix = idxA1j;
        iy = idxA1j;
        for (i = 0; i < j; i++) {
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
          int idxAjjp1;
          i = idxA1j + 61;
          idxAjjp1 = idxAjj + 61;
          if (j != 0) {
            iy = idxAjj + 60;
            b_i = (idxA1j + (60 * (58 - j))) + 61;
            for (idx = i; idx <= b_i; idx += 60) {
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
          b_i = (idxAjj + (60 * (58 - j))) + 61;
          for (i = idxAjjp1; i <= b_i; i += 60) {
            A_data[i - 1] *= ssq;
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
      for (b_i = 0; b_i < 60; b_i++) {
        d = 0.0;
        for (i1 = 0; i1 < 120; i1++) {
          d += (static_cast<double>(c_a[b_i + (60 * i1)])) * igr2[i1];
        }

        del_z[b_i] = r1[b_i] - d;
      }

      for (b_i = 0; b_i < 3600; b_i++) {
        y_tmp[b_i] = b_H[b_i] - y_tmp[b_i];
      }

      mldivide(y_tmp, del_z);

      // old method (LU?)
      //   del_z = linsolve (R, linsolve (R, (r1-At*igr2), opUT), opU); %exploit matrix properties for solving 
    } else {
      // Not Positive Definite (problem? eg infeasible)
      for (b_i = 0; b_i < 60; b_i++) {
        d = 0.0;
        for (i1 = 0; i1 < 120; i1++) {
          d += (static_cast<double>(c_a[b_i + (60 * i1)])) * igr2[i1];
        }

        del_z[b_i] = r1[b_i] - d;
      }

      for (b_i = 0; b_i < 3600; b_i++) {
        y_tmp[b_i] = b_H[b_i] - y_tmp[b_i];
      }

      mldivide(y_tmp, del_z);

      // old method (LU?)
    }

    // Decide on suitable alpha (from Wright's paper)
    // Try Max Increment (alpha = 1)
    // Check lam and ftol > 0
    for (i = 0; i < 120; i++) {
      d = 0.0;
      for (b_i = 0; b_i < 60; b_i++) {
        d += IGA[i + (120 * b_i)] * del_z[b_i];
      }

      d = igr2[i] - d;
      del_lam[i] = d;
      ssq = ftol[i];
      mu_old = ((-ssq) + mesil[i]) - ((ilam[i] * ssq) * d);
      mesil[i] = mu_old;
      d += lam[i];
      nlamt[i] = d;
      ssq += mu_old;
      ilam[i] = ssq;
      b_x[i] = (d < 2.2204460492503131E-16);
      b_x[i + 120] = (ssq < 2.2204460492503131E-16);
    }

    y = false;
    i = 0;
    exitg1 = false;
    while ((!exitg1) && (i < 240)) {
      if (!b_x[i]) {
        i++;
      } else {
        y = true;
        exitg1 = true;
      }
    }

    if (!y) {
      // KKT met
      std::memcpy(&lam[0], &nlamt[0], 120U * (sizeof(double)));
      std::memcpy(&ftol[0], &ilam[0], 120U * (sizeof(double)));
      for (b_i = 0; b_i < 60; b_i++) {
        x[b_i] += del_z[b_i];
      }
    } else {
      // KKT failed - solve by finding minimum ratio
      for (b_i = 0; b_i < 120; b_i++) {
        result[b_i] = nlamt[b_i];
        result[b_i + 120] = ilam[b_i];
      }

      idx = 0;
      i = 0;
      exitg1 = false;
      while ((!exitg1) && (i < 240)) {
        if (result[i] < 2.2204460492503131E-16) {
          idx++;
          ii_data[idx - 1] = static_cast<unsigned char>(static_cast<int>(i + 1));
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
        ii_size_idx_0 = 0;
      } else {
        ii_size_idx_0 = idx;
      }

      // detects elements breaking KKT condition
      for (b_i = 0; b_i < 120; b_i++) {
        b_del_lam[b_i] = del_lam[b_i];
        b_del_lam[b_i + 120] = mesil[b_i];
      }

      for (b_i = 0; b_i < ii_size_idx_0; b_i++) {
        i = (static_cast<int>(ii_data[b_i])) - 1;
        varargin_1_data[b_i] = 1.0 - (result[i] / b_del_lam[i]);
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
          i = 2;
          exitg1 = false;
          while ((!exitg1) && (i <= ii_size_idx_0)) {
            if (!rtIsNaN(varargin_1_data[i - 1])) {
              idx = i;
              exitg1 = true;
            } else {
              i++;
            }
          }
        }

        if (idx == 0) {
          ssq = varargin_1_data[0];
        } else {
          ssq = varargin_1_data[idx - 1];
          b_i = idx + 1;
          for (i = b_i; i <= ii_size_idx_0; i++) {
            d = varargin_1_data[i - 1];
            if (ssq > d) {
              ssq = d;
            }
          }
        }
      }

      ssq *= 0.995;

      // solves for min ratio (max value of alpha allowed)
      // Increment
      for (b_i = 0; b_i < 120; b_i++) {
        lam[b_i] += ssq * del_lam[b_i];
        ftol[b_i] += ssq * mesil[b_i];
      }

      for (b_i = 0; b_i < 60; b_i++) {
        x[b_i] += ssq * del_z[b_i];
      }
    }

    // Complimentary Gap
    mu_old = mu;
    ssq = 0.0;
    for (b_i = 0; b_i < 120; b_i++) {
      ssq += ftol[b_i] * lam[b_i];
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

      for (i = 0; i < 120; i++) {
        esig[i] = ssq;
      }
    }

    unusedU0++;
  }

  // Check for failure
  //  if(iter == maxiter)
  //      exitflag = -1;
  //  end
  Fx = x[0];
  Fy = x[1];
  Fz = x[2];
  
  for (i = 0; i < 60; i++) {
    X_QP_vl[i]=x[i];
  }



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
void CoordinatorBase<T>::MPC()
{
 static const double b_H[3600] = { 14.5253, 0.0, 0.0, 11.8592, 0.0, 0.0,
    11.1948, 0.0, 0.0, 10.5338, 0.0, 0.0, 9.8779, 0.0, 0.0, 9.2289, 0.0, 0.0,
    8.5883, 0.0, 0.0, 7.958, 0.0, 0.0, 7.3395, 0.0, 0.0, 6.7347, 0.0, 0.0,
    6.1452, 0.0, 0.0, 5.5728, 0.0, 0.0, 5.0191, 0.0, 0.0, 4.4858, 0.0, 0.0,
    3.9746, 0.0, 0.0, 3.4873, 0.0, 0.0, 3.0255, 0.0, 0.0, 2.5909, 0.0, 0.0,
    2.1853, 0.0, 0.0, 1.8103, 0.0, 0.0, 0.0, 14.5253, 0.0, 0.0, 11.8592, 0.0,
    0.0, 11.1948, 0.0, 0.0, 10.5338, 0.0, 0.0, 9.8779, 0.0, 0.0, 9.2289, 0.0,
    0.0, 8.5883, 0.0, 0.0, 7.958, 0.0, 0.0, 7.3395, 0.0, 0.0, 6.7347, 0.0, 0.0,
    6.1452, 0.0, 0.0, 5.5728, 0.0, 0.0, 5.0191, 0.0, 0.0, 4.4858, 0.0, 0.0,
    3.9746, 0.0, 0.0, 3.4873, 0.0, 0.0, 3.0255, 0.0, 0.0, 2.5909, 0.0, 0.0,
    2.1853, 0.0, 0.0, 1.8103, 0.0, 0.0, 0.0, 14.5253, 0.0, 0.0, 11.8592, 0.0,
    0.0, 11.1948, 0.0, 0.0, 10.5338, 0.0, 0.0, 9.8779, 0.0, 0.0, 9.2289, 0.0,
    0.0, 8.5883, 0.0, 0.0, 7.958, 0.0, 0.0, 7.3395, 0.0, 0.0, 6.7347, 0.0, 0.0,
    6.1452, 0.0, 0.0, 5.5728, 0.0, 0.0, 5.0191, 0.0, 0.0, 4.4858, 0.0, 0.0,
    3.9746, 0.0, 0.0, 3.4873, 0.0, 0.0, 3.0255, 0.0, 0.0, 2.5909, 0.0, 0.0,
    2.1853, 0.0, 0.0, 1.8103, 11.8592, 0.0, 0.0, 13.2435, 0.0, 0.0, 10.6241, 0.0,
    0.0, 10.0064, 0.0, 0.0, 9.3921, 0.0, 0.0, 8.783, 0.0, 0.0, 8.1806, 0.0, 0.0,
    7.5868, 0.0, 0.0, 7.0031, 0.0, 0.0, 6.4314, 0.0, 0.0, 5.8733, 0.0, 0.0,
    5.3306, 0.0, 0.0, 4.8048, 0.0, 0.0, 4.2978, 0.0, 0.0, 3.8112, 0.0, 0.0,
    3.3468, 0.0, 0.0, 2.9062, 0.0, 0.0, 2.4911, 0.0, 0.0, 2.1033, 0.0, 0.0,
    1.7444, 0.0, 0.0, 0.0, 11.8592, 0.0, 0.0, 13.2435, 0.0, 0.0, 10.6241, 0.0,
    0.0, 10.0064, 0.0, 0.0, 9.3921, 0.0, 0.0, 8.783, 0.0, 0.0, 8.1806, 0.0, 0.0,
    7.5868, 0.0, 0.0, 7.0031, 0.0, 0.0, 6.4314, 0.0, 0.0, 5.8733, 0.0, 0.0,
    5.3306, 0.0, 0.0, 4.8048, 0.0, 0.0, 4.2978, 0.0, 0.0, 3.8112, 0.0, 0.0,
    3.3468, 0.0, 0.0, 2.9062, 0.0, 0.0, 2.4911, 0.0, 0.0, 2.1033, 0.0, 0.0,
    1.7444, 0.0, 0.0, 0.0, 11.8592, 0.0, 0.0, 13.2435, 0.0, 0.0, 10.6241, 0.0,
    0.0, 10.0064, 0.0, 0.0, 9.3921, 0.0, 0.0, 8.783, 0.0, 0.0, 8.1806, 0.0, 0.0,
    7.5868, 0.0, 0.0, 7.0031, 0.0, 0.0, 6.4314, 0.0, 0.0, 5.8733, 0.0, 0.0,
    5.3306, 0.0, 0.0, 4.8048, 0.0, 0.0, 4.2978, 0.0, 0.0, 3.8112, 0.0, 0.0,
    3.3468, 0.0, 0.0, 2.9062, 0.0, 0.0, 2.4911, 0.0, 0.0, 2.1033, 0.0, 0.0,
    1.7444, 11.1948, 0.0, 0.0, 10.6241, 0.0, 0.0, 12.0534, 0.0, 0.0, 9.479, 0.0,
    0.0, 8.9063, 0.0, 0.0, 8.3371, 0.0, 0.0, 7.7729, 0.0, 0.0, 7.2156, 0.0, 0.0,
    6.6668, 0.0, 0.0, 6.1281, 0.0, 0.0, 5.6014, 0.0, 0.0, 5.0884, 0.0, 0.0,
    4.5906, 0.0, 0.0, 4.1099, 0.0, 0.0, 3.6479, 0.0, 0.0, 3.2063, 0.0, 0.0,
    2.7869, 0.0, 0.0, 2.3913, 0.0, 0.0, 2.0212, 0.0, 0.0, 1.6784, 0.0, 0.0, 0.0,
    11.1948, 0.0, 0.0, 10.6241, 0.0, 0.0, 12.0534, 0.0, 0.0, 9.479, 0.0, 0.0,
    8.9063, 0.0, 0.0, 8.3371, 0.0, 0.0, 7.7729, 0.0, 0.0, 7.2156, 0.0, 0.0,
    6.6668, 0.0, 0.0, 6.1281, 0.0, 0.0, 5.6014, 0.0, 0.0, 5.0884, 0.0, 0.0,
    4.5906, 0.0, 0.0, 4.1099, 0.0, 0.0, 3.6479, 0.0, 0.0, 3.2063, 0.0, 0.0,
    2.7869, 0.0, 0.0, 2.3913, 0.0, 0.0, 2.0212, 0.0, 0.0, 1.6784, 0.0, 0.0, 0.0,
    11.1948, 0.0, 0.0, 10.6241, 0.0, 0.0, 12.0534, 0.0, 0.0, 9.479, 0.0, 0.0,
    8.9063, 0.0, 0.0, 8.3371, 0.0, 0.0, 7.7729, 0.0, 0.0, 7.2156, 0.0, 0.0,
    6.6668, 0.0, 0.0, 6.1281, 0.0, 0.0, 5.6014, 0.0, 0.0, 5.0884, 0.0, 0.0,
    4.5906, 0.0, 0.0, 4.1099, 0.0, 0.0, 3.6479, 0.0, 0.0, 3.2063, 0.0, 0.0,
    2.7869, 0.0, 0.0, 2.3913, 0.0, 0.0, 2.0212, 0.0, 0.0, 1.6784, 10.5338, 0.0,
    0.0, 10.0064, 0.0, 0.0, 9.479, 0.0, 0.0, 10.9516, 0.0, 0.0, 8.4206, 0.0, 0.0,
    7.8912, 0.0, 0.0, 7.3653, 0.0, 0.0, 6.8444, 0.0, 0.0, 6.3304, 0.0, 0.0,
    5.8249, 0.0, 0.0, 5.3296, 0.0, 0.0, 4.8462, 0.0, 0.0, 4.3764, 0.0, 0.0,
    3.922, 0.0, 0.0, 3.4845, 0.0, 0.0, 3.0658, 0.0, 0.0, 2.6676, 0.0, 0.0,
    2.2915, 0.0, 0.0, 1.9392, 0.0, 0.0, 1.6124, 0.0, 0.0, 0.0, 10.5338, 0.0, 0.0,
    10.0064, 0.0, 0.0, 9.479, 0.0, 0.0, 10.9516, 0.0, 0.0, 8.4206, 0.0, 0.0,
    7.8912, 0.0, 0.0, 7.3653, 0.0, 0.0, 6.8444, 0.0, 0.0, 6.3304, 0.0, 0.0,
    5.8249, 0.0, 0.0, 5.3296, 0.0, 0.0, 4.8462, 0.0, 0.0, 4.3764, 0.0, 0.0,
    3.922, 0.0, 0.0, 3.4845, 0.0, 0.0, 3.0658, 0.0, 0.0, 2.6676, 0.0, 0.0,
    2.2915, 0.0, 0.0, 1.9392, 0.0, 0.0, 1.6124, 0.0, 0.0, 0.0, 10.5338, 0.0, 0.0,
    10.0064, 0.0, 0.0, 9.479, 0.0, 0.0, 10.9516, 0.0, 0.0, 8.4206, 0.0, 0.0,
    7.8912, 0.0, 0.0, 7.3653, 0.0, 0.0, 6.8444, 0.0, 0.0, 6.3304, 0.0, 0.0,
    5.8249, 0.0, 0.0, 5.3296, 0.0, 0.0, 4.8462, 0.0, 0.0, 4.3764, 0.0, 0.0,
    3.922, 0.0, 0.0, 3.4845, 0.0, 0.0, 3.0658, 0.0, 0.0, 2.6676, 0.0, 0.0,
    2.2915, 0.0, 0.0, 1.9392, 0.0, 0.0, 1.6124, 9.8779, 0.0, 0.0, 9.3921, 0.0,
    0.0, 8.9063, 0.0, 0.0, 8.4206, 0.0, 0.0, 9.9348, 0.0, 0.0, 7.4453, 0.0, 0.0,
    6.9576, 0.0, 0.0, 6.4732, 0.0, 0.0, 5.994, 0.0, 0.0, 5.5216, 0.0, 0.0,
    5.0577, 0.0, 0.0, 4.604, 0.0, 0.0, 4.1622, 0.0, 0.0, 3.734, 0.0, 0.0, 3.3212,
    0.0, 0.0, 2.9254, 0.0, 0.0, 2.5483, 0.0, 0.0, 2.1916, 0.0, 0.0, 1.8571, 0.0,
    0.0, 1.5464, 0.0, 0.0, 0.0, 9.8779, 0.0, 0.0, 9.3921, 0.0, 0.0, 8.9063, 0.0,
    0.0, 8.4206, 0.0, 0.0, 9.9348, 0.0, 0.0, 7.4453, 0.0, 0.0, 6.9576, 0.0, 0.0,
    6.4732, 0.0, 0.0, 5.994, 0.0, 0.0, 5.5216, 0.0, 0.0, 5.0577, 0.0, 0.0, 4.604,
    0.0, 0.0, 4.1622, 0.0, 0.0, 3.734, 0.0, 0.0, 3.3212, 0.0, 0.0, 2.9254, 0.0,
    0.0, 2.5483, 0.0, 0.0, 2.1916, 0.0, 0.0, 1.8571, 0.0, 0.0, 1.5464, 0.0, 0.0,
    0.0, 9.8779, 0.0, 0.0, 9.3921, 0.0, 0.0, 8.9063, 0.0, 0.0, 8.4206, 0.0, 0.0,
    9.9348, 0.0, 0.0, 7.4453, 0.0, 0.0, 6.9576, 0.0, 0.0, 6.4732, 0.0, 0.0,
    5.994, 0.0, 0.0, 5.5216, 0.0, 0.0, 5.0577, 0.0, 0.0, 4.604, 0.0, 0.0, 4.1622,
    0.0, 0.0, 3.734, 0.0, 0.0, 3.3212, 0.0, 0.0, 2.9254, 0.0, 0.0, 2.5483, 0.0,
    0.0, 2.1916, 0.0, 0.0, 1.8571, 0.0, 0.0, 1.5464, 9.2289, 0.0, 0.0, 8.783,
    0.0, 0.0, 8.3371, 0.0, 0.0, 7.8912, 0.0, 0.0, 7.4453, 0.0, 0.0, 8.9994, 0.0,
    0.0, 6.5499, 0.0, 0.0, 6.102, 0.0, 0.0, 5.6576, 0.0, 0.0, 5.2183, 0.0, 0.0,
    4.7858, 0.0, 0.0, 4.3618, 0.0, 0.0, 3.948, 0.0, 0.0, 3.5461, 0.0, 0.0,
    3.1578, 0.0, 0.0, 2.7849, 0.0, 0.0, 2.429, 0.0, 0.0, 2.0918, 0.0, 0.0,
    1.7751, 0.0, 0.0, 1.4805, 0.0, 0.0, 0.0, 9.2289, 0.0, 0.0, 8.783, 0.0, 0.0,
    8.3371, 0.0, 0.0, 7.8912, 0.0, 0.0, 7.4453, 0.0, 0.0, 8.9994, 0.0, 0.0,
    6.5499, 0.0, 0.0, 6.102, 0.0, 0.0, 5.6576, 0.0, 0.0, 5.2183, 0.0, 0.0,
    4.7858, 0.0, 0.0, 4.3618, 0.0, 0.0, 3.948, 0.0, 0.0, 3.5461, 0.0, 0.0,
    3.1578, 0.0, 0.0, 2.7849, 0.0, 0.0, 2.429, 0.0, 0.0, 2.0918, 0.0, 0.0,
    1.7751, 0.0, 0.0, 1.4805, 0.0, 0.0, 0.0, 9.2289, 0.0, 0.0, 8.783, 0.0, 0.0,
    8.3371, 0.0, 0.0, 7.8912, 0.0, 0.0, 7.4453, 0.0, 0.0, 8.9994, 0.0, 0.0,
    6.5499, 0.0, 0.0, 6.102, 0.0, 0.0, 5.6576, 0.0, 0.0, 5.2183, 0.0, 0.0,
    4.7858, 0.0, 0.0, 4.3618, 0.0, 0.0, 3.948, 0.0, 0.0, 3.5461, 0.0, 0.0,
    3.1578, 0.0, 0.0, 2.7849, 0.0, 0.0, 2.429, 0.0, 0.0, 2.0918, 0.0, 0.0,
    1.7751, 0.0, 0.0, 1.4805, 8.5883, 0.0, 0.0, 8.1806, 0.0, 0.0, 7.7729, 0.0,
    0.0, 7.3653, 0.0, 0.0, 6.9576, 0.0, 0.0, 6.5499, 0.0, 0.0, 8.1422, 0.0, 0.0,
    5.7309, 0.0, 0.0, 5.3212, 0.0, 0.0, 4.915, 0.0, 0.0, 4.5139, 0.0, 0.0,
    4.1195, 0.0, 0.0, 3.7337, 0.0, 0.0, 3.3582, 0.0, 0.0, 2.9945, 0.0, 0.0,
    2.6444, 0.0, 0.0, 2.3097, 0.0, 0.0, 1.992, 0.0, 0.0, 1.693, 0.0, 0.0, 1.4145,
    0.0, 0.0, 0.0, 8.5883, 0.0, 0.0, 8.1806, 0.0, 0.0, 7.7729, 0.0, 0.0, 7.3653,
    0.0, 0.0, 6.9576, 0.0, 0.0, 6.5499, 0.0, 0.0, 8.1422, 0.0, 0.0, 5.7309, 0.0,
    0.0, 5.3212, 0.0, 0.0, 4.915, 0.0, 0.0, 4.5139, 0.0, 0.0, 4.1195, 0.0, 0.0,
    3.7337, 0.0, 0.0, 3.3582, 0.0, 0.0, 2.9945, 0.0, 0.0, 2.6444, 0.0, 0.0,
    2.3097, 0.0, 0.0, 1.992, 0.0, 0.0, 1.693, 0.0, 0.0, 1.4145, 0.0, 0.0, 0.0,
    8.5883, 0.0, 0.0, 8.1806, 0.0, 0.0, 7.7729, 0.0, 0.0, 7.3653, 0.0, 0.0,
    6.9576, 0.0, 0.0, 6.5499, 0.0, 0.0, 8.1422, 0.0, 0.0, 5.7309, 0.0, 0.0,
    5.3212, 0.0, 0.0, 4.915, 0.0, 0.0, 4.5139, 0.0, 0.0, 4.1195, 0.0, 0.0,
    3.7337, 0.0, 0.0, 3.3582, 0.0, 0.0, 2.9945, 0.0, 0.0, 2.6444, 0.0, 0.0,
    2.3097, 0.0, 0.0, 1.992, 0.0, 0.0, 1.693, 0.0, 0.0, 1.4145, 7.958, 0.0, 0.0,
    7.5868, 0.0, 0.0, 7.2156, 0.0, 0.0, 6.8444, 0.0, 0.0, 6.4732, 0.0, 0.0,
    6.102, 0.0, 0.0, 5.7309, 0.0, 0.0, 7.3597, 0.0, 0.0, 4.9848, 0.0, 0.0,
    4.6117, 0.0, 0.0, 4.242, 0.0, 0.0, 3.8773, 0.0, 0.0, 3.5195, 0.0, 0.0,
    3.1702, 0.0, 0.0, 2.8311, 0.0, 0.0, 2.504, 0.0, 0.0, 2.1904, 0.0, 0.0,
    1.8922, 0.0, 0.0, 1.611, 0.0, 0.0, 1.3485, 0.0, 0.0, 0.0, 7.958, 0.0, 0.0,
    7.5868, 0.0, 0.0, 7.2156, 0.0, 0.0, 6.8444, 0.0, 0.0, 6.4732, 0.0, 0.0,
    6.102, 0.0, 0.0, 5.7309, 0.0, 0.0, 7.3597, 0.0, 0.0, 4.9848, 0.0, 0.0,
    4.6117, 0.0, 0.0, 4.242, 0.0, 0.0, 3.8773, 0.0, 0.0, 3.5195, 0.0, 0.0,
    3.1702, 0.0, 0.0, 2.8311, 0.0, 0.0, 2.504, 0.0, 0.0, 2.1904, 0.0, 0.0,
    1.8922, 0.0, 0.0, 1.611, 0.0, 0.0, 1.3485, 0.0, 0.0, 0.0, 7.958, 0.0, 0.0,
    7.5868, 0.0, 0.0, 7.2156, 0.0, 0.0, 6.8444, 0.0, 0.0, 6.4732, 0.0, 0.0,
    6.102, 0.0, 0.0, 5.7309, 0.0, 0.0, 7.3597, 0.0, 0.0, 4.9848, 0.0, 0.0,
    4.6117, 0.0, 0.0, 4.242, 0.0, 0.0, 3.8773, 0.0, 0.0, 3.5195, 0.0, 0.0,
    3.1702, 0.0, 0.0, 2.8311, 0.0, 0.0, 2.504, 0.0, 0.0, 2.1904, 0.0, 0.0,
    1.8922, 0.0, 0.0, 1.611, 0.0, 0.0, 1.3485, 7.3395, 0.0, 0.0, 7.0031, 0.0,
    0.0, 6.6668, 0.0, 0.0, 6.3304, 0.0, 0.0, 5.994, 0.0, 0.0, 5.6576, 0.0, 0.0,
    5.3212, 0.0, 0.0, 4.9848, 0.0, 0.0, 6.6485, 0.0, 0.0, 4.3084, 0.0, 0.0,
    3.9701, 0.0, 0.0, 3.6351, 0.0, 0.0, 3.3053, 0.0, 0.0, 2.9823, 0.0, 0.0,
    2.6678, 0.0, 0.0, 2.3635, 0.0, 0.0, 2.0711, 0.0, 0.0, 1.7924, 0.0, 0.0,
    1.5289, 0.0, 0.0, 1.2825, 0.0, 0.0, 0.0, 7.3395, 0.0, 0.0, 7.0031, 0.0, 0.0,
    6.6668, 0.0, 0.0, 6.3304, 0.0, 0.0, 5.994, 0.0, 0.0, 5.6576, 0.0, 0.0,
    5.3212, 0.0, 0.0, 4.9848, 0.0, 0.0, 6.6485, 0.0, 0.0, 4.3084, 0.0, 0.0,
    3.9701, 0.0, 0.0, 3.6351, 0.0, 0.0, 3.3053, 0.0, 0.0, 2.9823, 0.0, 0.0,
    2.6678, 0.0, 0.0, 2.3635, 0.0, 0.0, 2.0711, 0.0, 0.0, 1.7924, 0.0, 0.0,
    1.5289, 0.0, 0.0, 1.2825, 0.0, 0.0, 0.0, 7.3395, 0.0, 0.0, 7.0031, 0.0, 0.0,
    6.6668, 0.0, 0.0, 6.3304, 0.0, 0.0, 5.994, 0.0, 0.0, 5.6576, 0.0, 0.0,
    5.3212, 0.0, 0.0, 4.9848, 0.0, 0.0, 6.6485, 0.0, 0.0, 4.3084, 0.0, 0.0,
    3.9701, 0.0, 0.0, 3.6351, 0.0, 0.0, 3.3053, 0.0, 0.0, 2.9823, 0.0, 0.0,
    2.6678, 0.0, 0.0, 2.3635, 0.0, 0.0, 2.0711, 0.0, 0.0, 1.7924, 0.0, 0.0,
    1.5289, 0.0, 0.0, 1.2825, 6.7347, 0.0, 0.0, 6.4314, 0.0, 0.0, 6.1281, 0.0,
    0.0, 5.8249, 0.0, 0.0, 5.5216, 0.0, 0.0, 5.2183, 0.0, 0.0, 4.915, 0.0, 0.0,
    4.6117, 0.0, 0.0, 4.3084, 0.0, 0.0, 6.0051, 0.0, 0.0, 3.6982, 0.0, 0.0,
    3.3929, 0.0, 0.0, 3.0911, 0.0, 0.0, 2.7944, 0.0, 0.0, 2.5044, 0.0, 0.0,
    2.223, 0.0, 0.0, 1.9518, 0.0, 0.0, 1.6925, 0.0, 0.0, 1.4469, 0.0, 0.0,
    1.2165, 0.0, 0.0, 0.0, 6.7347, 0.0, 0.0, 6.4314, 0.0, 0.0, 6.1281, 0.0, 0.0,
    5.8249, 0.0, 0.0, 5.5216, 0.0, 0.0, 5.2183, 0.0, 0.0, 4.915, 0.0, 0.0,
    4.6117, 0.0, 0.0, 4.3084, 0.0, 0.0, 6.0051, 0.0, 0.0, 3.6982, 0.0, 0.0,
    3.3929, 0.0, 0.0, 3.0911, 0.0, 0.0, 2.7944, 0.0, 0.0, 2.5044, 0.0, 0.0,
    2.223, 0.0, 0.0, 1.9518, 0.0, 0.0, 1.6925, 0.0, 0.0, 1.4469, 0.0, 0.0,
    1.2165, 0.0, 0.0, 0.0, 6.7347, 0.0, 0.0, 6.4314, 0.0, 0.0, 6.1281, 0.0, 0.0,
    5.8249, 0.0, 0.0, 5.5216, 0.0, 0.0, 5.2183, 0.0, 0.0, 4.915, 0.0, 0.0,
    4.6117, 0.0, 0.0, 4.3084, 0.0, 0.0, 6.0051, 0.0, 0.0, 3.6982, 0.0, 0.0,
    3.3929, 0.0, 0.0, 3.0911, 0.0, 0.0, 2.7944, 0.0, 0.0, 2.5044, 0.0, 0.0,
    2.223, 0.0, 0.0, 1.9518, 0.0, 0.0, 1.6925, 0.0, 0.0, 1.4469, 0.0, 0.0,
    1.2165, 6.1452, 0.0, 0.0, 5.8733, 0.0, 0.0, 5.6014, 0.0, 0.0, 5.3296, 0.0,
    0.0, 5.0577, 0.0, 0.0, 4.7858, 0.0, 0.0, 4.5139, 0.0, 0.0, 4.242, 0.0, 0.0,
    3.9701, 0.0, 0.0, 3.6982, 0.0, 0.0, 5.4263, 0.0, 0.0, 3.1507, 0.0, 0.0,
    2.8769, 0.0, 0.0, 2.6064, 0.0, 0.0, 2.3411, 0.0, 0.0, 2.0825, 0.0, 0.0,
    1.8325, 0.0, 0.0, 1.5927, 0.0, 0.0, 1.3648, 0.0, 0.0, 1.1506, 0.0, 0.0, 0.0,
    6.1452, 0.0, 0.0, 5.8733, 0.0, 0.0, 5.6014, 0.0, 0.0, 5.3296, 0.0, 0.0,
    5.0577, 0.0, 0.0, 4.7858, 0.0, 0.0, 4.5139, 0.0, 0.0, 4.242, 0.0, 0.0,
    3.9701, 0.0, 0.0, 3.6982, 0.0, 0.0, 5.4263, 0.0, 0.0, 3.1507, 0.0, 0.0,
    2.8769, 0.0, 0.0, 2.6064, 0.0, 0.0, 2.3411, 0.0, 0.0, 2.0825, 0.0, 0.0,
    1.8325, 0.0, 0.0, 1.5927, 0.0, 0.0, 1.3648, 0.0, 0.0, 1.1506, 0.0, 0.0, 0.0,
    6.1452, 0.0, 0.0, 5.8733, 0.0, 0.0, 5.6014, 0.0, 0.0, 5.3296, 0.0, 0.0,
    5.0577, 0.0, 0.0, 4.7858, 0.0, 0.0, 4.5139, 0.0, 0.0, 4.242, 0.0, 0.0,
    3.9701, 0.0, 0.0, 3.6982, 0.0, 0.0, 5.4263, 0.0, 0.0, 3.1507, 0.0, 0.0,
    2.8769, 0.0, 0.0, 2.6064, 0.0, 0.0, 2.3411, 0.0, 0.0, 2.0825, 0.0, 0.0,
    1.8325, 0.0, 0.0, 1.5927, 0.0, 0.0, 1.3648, 0.0, 0.0, 1.1506, 5.5728, 0.0,
    0.0, 5.3306, 0.0, 0.0, 5.0884, 0.0, 0.0, 4.8462, 0.0, 0.0, 4.604, 0.0, 0.0,
    4.3618, 0.0, 0.0, 4.1195, 0.0, 0.0, 3.8773, 0.0, 0.0, 3.6351, 0.0, 0.0,
    3.3929, 0.0, 0.0, 3.1507, 0.0, 0.0, 4.9085, 0.0, 0.0, 2.6626, 0.0, 0.0,
    2.4185, 0.0, 0.0, 2.1777, 0.0, 0.0, 1.9421, 0.0, 0.0, 1.7132, 0.0, 0.0,
    1.4929, 0.0, 0.0, 1.2828, 0.0, 0.0, 1.0846, 0.0, 0.0, 0.0, 5.5728, 0.0, 0.0,
    5.3306, 0.0, 0.0, 5.0884, 0.0, 0.0, 4.8462, 0.0, 0.0, 4.604, 0.0, 0.0,
    4.3618, 0.0, 0.0, 4.1195, 0.0, 0.0, 3.8773, 0.0, 0.0, 3.6351, 0.0, 0.0,
    3.3929, 0.0, 0.0, 3.1507, 0.0, 0.0, 4.9085, 0.0, 0.0, 2.6626, 0.0, 0.0,
    2.4185, 0.0, 0.0, 2.1777, 0.0, 0.0, 1.9421, 0.0, 0.0, 1.7132, 0.0, 0.0,
    1.4929, 0.0, 0.0, 1.2828, 0.0, 0.0, 1.0846, 0.0, 0.0, 0.0, 5.5728, 0.0, 0.0,
    5.3306, 0.0, 0.0, 5.0884, 0.0, 0.0, 4.8462, 0.0, 0.0, 4.604, 0.0, 0.0,
    4.3618, 0.0, 0.0, 4.1195, 0.0, 0.0, 3.8773, 0.0, 0.0, 3.6351, 0.0, 0.0,
    3.3929, 0.0, 0.0, 3.1507, 0.0, 0.0, 4.9085, 0.0, 0.0, 2.6626, 0.0, 0.0,
    2.4185, 0.0, 0.0, 2.1777, 0.0, 0.0, 1.9421, 0.0, 0.0, 1.7132, 0.0, 0.0,
    1.4929, 0.0, 0.0, 1.2828, 0.0, 0.0, 1.0846, 5.0191, 0.0, 0.0, 4.8048, 0.0,
    0.0, 4.5906, 0.0, 0.0, 4.3764, 0.0, 0.0, 4.1622, 0.0, 0.0, 3.948, 0.0, 0.0,
    3.7337, 0.0, 0.0, 3.5195, 0.0, 0.0, 3.3053, 0.0, 0.0, 3.0911, 0.0, 0.0,
    2.8769, 0.0, 0.0, 2.6626, 0.0, 0.0, 4.4484, 0.0, 0.0, 2.2305, 0.0, 0.0,
    2.0144, 0.0, 0.0, 1.8016, 0.0, 0.0, 1.5939, 0.0, 0.0, 1.3931, 0.0, 0.0,
    1.2007, 0.0, 0.0, 1.0186, 0.0, 0.0, 0.0, 5.0191, 0.0, 0.0, 4.8048, 0.0, 0.0,
    4.5906, 0.0, 0.0, 4.3764, 0.0, 0.0, 4.1622, 0.0, 0.0, 3.948, 0.0, 0.0,
    3.7337, 0.0, 0.0, 3.5195, 0.0, 0.0, 3.3053, 0.0, 0.0, 3.0911, 0.0, 0.0,
    2.8769, 0.0, 0.0, 2.6626, 0.0, 0.0, 4.4484, 0.0, 0.0, 2.2305, 0.0, 0.0,
    2.0144, 0.0, 0.0, 1.8016, 0.0, 0.0, 1.5939, 0.0, 0.0, 1.3931, 0.0, 0.0,
    1.2007, 0.0, 0.0, 1.0186, 0.0, 0.0, 0.0, 5.0191, 0.0, 0.0, 4.8048, 0.0, 0.0,
    4.5906, 0.0, 0.0, 4.3764, 0.0, 0.0, 4.1622, 0.0, 0.0, 3.948, 0.0, 0.0,
    3.7337, 0.0, 0.0, 3.5195, 0.0, 0.0, 3.3053, 0.0, 0.0, 3.0911, 0.0, 0.0,
    2.8769, 0.0, 0.0, 2.6626, 0.0, 0.0, 4.4484, 0.0, 0.0, 2.2305, 0.0, 0.0,
    2.0144, 0.0, 0.0, 1.8016, 0.0, 0.0, 1.5939, 0.0, 0.0, 1.3931, 0.0, 0.0,
    1.2007, 0.0, 0.0, 1.0186, 4.4858, 0.0, 0.0, 4.2978, 0.0, 0.0, 4.1099, 0.0,
    0.0, 3.922, 0.0, 0.0, 3.734, 0.0, 0.0, 3.5461, 0.0, 0.0, 3.3582, 0.0, 0.0,
    3.1702, 0.0, 0.0, 2.9823, 0.0, 0.0, 2.7944, 0.0, 0.0, 2.6064, 0.0, 0.0,
    2.4185, 0.0, 0.0, 2.2305, 0.0, 0.0, 4.0426, 0.0, 0.0, 1.851, 0.0, 0.0,
    1.6611, 0.0, 0.0, 1.4746, 0.0, 0.0, 1.2933, 0.0, 0.0, 1.1187, 0.0, 0.0,
    0.95264, 0.0, 0.0, 0.0, 4.4858, 0.0, 0.0, 4.2978, 0.0, 0.0, 4.1099, 0.0, 0.0,
    3.922, 0.0, 0.0, 3.734, 0.0, 0.0, 3.5461, 0.0, 0.0, 3.3582, 0.0, 0.0, 3.1702,
    0.0, 0.0, 2.9823, 0.0, 0.0, 2.7944, 0.0, 0.0, 2.6064, 0.0, 0.0, 2.4185, 0.0,
    0.0, 2.2305, 0.0, 0.0, 4.0426, 0.0, 0.0, 1.851, 0.0, 0.0, 1.6611, 0.0, 0.0,
    1.4746, 0.0, 0.0, 1.2933, 0.0, 0.0, 1.1187, 0.0, 0.0, 0.95264, 0.0, 0.0, 0.0,
    4.4858, 0.0, 0.0, 4.2978, 0.0, 0.0, 4.1099, 0.0, 0.0, 3.922, 0.0, 0.0, 3.734,
    0.0, 0.0, 3.5461, 0.0, 0.0, 3.3582, 0.0, 0.0, 3.1702, 0.0, 0.0, 2.9823, 0.0,
    0.0, 2.7944, 0.0, 0.0, 2.6064, 0.0, 0.0, 2.4185, 0.0, 0.0, 2.2305, 0.0, 0.0,
    4.0426, 0.0, 0.0, 1.851, 0.0, 0.0, 1.6611, 0.0, 0.0, 1.4746, 0.0, 0.0,
    1.2933, 0.0, 0.0, 1.1187, 0.0, 0.0, 0.95264, 3.9746, 0.0, 0.0, 3.8112, 0.0,
    0.0, 3.6479, 0.0, 0.0, 3.4845, 0.0, 0.0, 3.3212, 0.0, 0.0, 3.1578, 0.0, 0.0,
    2.9945, 0.0, 0.0, 2.8311, 0.0, 0.0, 2.6678, 0.0, 0.0, 2.5044, 0.0, 0.0,
    2.3411, 0.0, 0.0, 2.1777, 0.0, 0.0, 2.0144, 0.0, 0.0, 1.851, 0.0, 0.0,
    3.6877, 0.0, 0.0, 1.5207, 0.0, 0.0, 1.3553, 0.0, 0.0, 1.1934, 0.0, 0.0,
    1.0366, 0.0, 0.0, 0.88666, 0.0, 0.0, 0.0, 3.9746, 0.0, 0.0, 3.8112, 0.0, 0.0,
    3.6479, 0.0, 0.0, 3.4845, 0.0, 0.0, 3.3212, 0.0, 0.0, 3.1578, 0.0, 0.0,
    2.9945, 0.0, 0.0, 2.8311, 0.0, 0.0, 2.6678, 0.0, 0.0, 2.5044, 0.0, 0.0,
    2.3411, 0.0, 0.0, 2.1777, 0.0, 0.0, 2.0144, 0.0, 0.0, 1.851, 0.0, 0.0,
    3.6877, 0.0, 0.0, 1.5207, 0.0, 0.0, 1.3553, 0.0, 0.0, 1.1934, 0.0, 0.0,
    1.0366, 0.0, 0.0, 0.88666, 0.0, 0.0, 0.0, 3.9746, 0.0, 0.0, 3.8112, 0.0, 0.0,
    3.6479, 0.0, 0.0, 3.4845, 0.0, 0.0, 3.3212, 0.0, 0.0, 3.1578, 0.0, 0.0,
    2.9945, 0.0, 0.0, 2.8311, 0.0, 0.0, 2.6678, 0.0, 0.0, 2.5044, 0.0, 0.0,
    2.3411, 0.0, 0.0, 2.1777, 0.0, 0.0, 2.0144, 0.0, 0.0, 1.851, 0.0, 0.0,
    3.6877, 0.0, 0.0, 1.5207, 0.0, 0.0, 1.3553, 0.0, 0.0, 1.1934, 0.0, 0.0,
    1.0366, 0.0, 0.0, 0.88666, 3.4873, 0.0, 0.0, 3.3468, 0.0, 0.0, 3.2063, 0.0,
    0.0, 3.0658, 0.0, 0.0, 2.9254, 0.0, 0.0, 2.7849, 0.0, 0.0, 2.6444, 0.0, 0.0,
    2.504, 0.0, 0.0, 2.3635, 0.0, 0.0, 2.223, 0.0, 0.0, 2.0825, 0.0, 0.0, 1.9421,
    0.0, 0.0, 1.8016, 0.0, 0.0, 1.6611, 0.0, 0.0, 1.5207, 0.0, 0.0, 3.3802, 0.0,
    0.0, 1.2361, 0.0, 0.0, 1.0936, 0.0, 0.0, 0.9546, 0.0, 0.0, 0.82068, 0.0, 0.0,
    0.0, 3.4873, 0.0, 0.0, 3.3468, 0.0, 0.0, 3.2063, 0.0, 0.0, 3.0658, 0.0, 0.0,
    2.9254, 0.0, 0.0, 2.7849, 0.0, 0.0, 2.6444, 0.0, 0.0, 2.504, 0.0, 0.0,
    2.3635, 0.0, 0.0, 2.223, 0.0, 0.0, 2.0825, 0.0, 0.0, 1.9421, 0.0, 0.0,
    1.8016, 0.0, 0.0, 1.6611, 0.0, 0.0, 1.5207, 0.0, 0.0, 3.3802, 0.0, 0.0,
    1.2361, 0.0, 0.0, 1.0936, 0.0, 0.0, 0.9546, 0.0, 0.0, 0.82068, 0.0, 0.0, 0.0,
    3.4873, 0.0, 0.0, 3.3468, 0.0, 0.0, 3.2063, 0.0, 0.0, 3.0658, 0.0, 0.0,
    2.9254, 0.0, 0.0, 2.7849, 0.0, 0.0, 2.6444, 0.0, 0.0, 2.504, 0.0, 0.0,
    2.3635, 0.0, 0.0, 2.223, 0.0, 0.0, 2.0825, 0.0, 0.0, 1.9421, 0.0, 0.0,
    1.8016, 0.0, 0.0, 1.6611, 0.0, 0.0, 1.5207, 0.0, 0.0, 3.3802, 0.0, 0.0,
    1.2361, 0.0, 0.0, 1.0936, 0.0, 0.0, 0.9546, 0.0, 0.0, 0.82068, 3.0255, 0.0,
    0.0, 2.9062, 0.0, 0.0, 2.7869, 0.0, 0.0, 2.6676, 0.0, 0.0, 2.5483, 0.0, 0.0,
    2.429, 0.0, 0.0, 2.3097, 0.0, 0.0, 2.1904, 0.0, 0.0, 2.0711, 0.0, 0.0,
    1.9518, 0.0, 0.0, 1.8325, 0.0, 0.0, 1.7132, 0.0, 0.0, 1.5939, 0.0, 0.0,
    1.4746, 0.0, 0.0, 1.3553, 0.0, 0.0, 1.2361, 0.0, 0.0, 3.1168, 0.0, 0.0,
    0.99381, 0.0, 0.0, 0.87255, 0.0, 0.0, 0.75471, 0.0, 0.0, 0.0, 3.0255, 0.0,
    0.0, 2.9062, 0.0, 0.0, 2.7869, 0.0, 0.0, 2.6676, 0.0, 0.0, 2.5483, 0.0, 0.0,
    2.429, 0.0, 0.0, 2.3097, 0.0, 0.0, 2.1904, 0.0, 0.0, 2.0711, 0.0, 0.0,
    1.9518, 0.0, 0.0, 1.8325, 0.0, 0.0, 1.7132, 0.0, 0.0, 1.5939, 0.0, 0.0,
    1.4746, 0.0, 0.0, 1.3553, 0.0, 0.0, 1.2361, 0.0, 0.0, 3.1168, 0.0, 0.0,
    0.99381, 0.0, 0.0, 0.87255, 0.0, 0.0, 0.75471, 0.0, 0.0, 0.0, 3.0255, 0.0,
    0.0, 2.9062, 0.0, 0.0, 2.7869, 0.0, 0.0, 2.6676, 0.0, 0.0, 2.5483, 0.0, 0.0,
    2.429, 0.0, 0.0, 2.3097, 0.0, 0.0, 2.1904, 0.0, 0.0, 2.0711, 0.0, 0.0,
    1.9518, 0.0, 0.0, 1.8325, 0.0, 0.0, 1.7132, 0.0, 0.0, 1.5939, 0.0, 0.0,
    1.4746, 0.0, 0.0, 1.3553, 0.0, 0.0, 1.2361, 0.0, 0.0, 3.1168, 0.0, 0.0,
    0.99381, 0.0, 0.0, 0.87255, 0.0, 0.0, 0.75471, 2.5909, 0.0, 0.0, 2.4911, 0.0,
    0.0, 2.3913, 0.0, 0.0, 2.2915, 0.0, 0.0, 2.1916, 0.0, 0.0, 2.0918, 0.0, 0.0,
    1.992, 0.0, 0.0, 1.8922, 0.0, 0.0, 1.7924, 0.0, 0.0, 1.6925, 0.0, 0.0,
    1.5927, 0.0, 0.0, 1.4929, 0.0, 0.0, 1.3931, 0.0, 0.0, 1.2933, 0.0, 0.0,
    1.1934, 0.0, 0.0, 1.0936, 0.0, 0.0, 0.99381, 0.0, 0.0, 2.894, 0.0, 0.0,
    0.79051, 0.0, 0.0, 0.68873, 0.0, 0.0, 0.0, 2.5909, 0.0, 0.0, 2.4911, 0.0,
    0.0, 2.3913, 0.0, 0.0, 2.2915, 0.0, 0.0, 2.1916, 0.0, 0.0, 2.0918, 0.0, 0.0,
    1.992, 0.0, 0.0, 1.8922, 0.0, 0.0, 1.7924, 0.0, 0.0, 1.6925, 0.0, 0.0,
    1.5927, 0.0, 0.0, 1.4929, 0.0, 0.0, 1.3931, 0.0, 0.0, 1.2933, 0.0, 0.0,
    1.1934, 0.0, 0.0, 1.0936, 0.0, 0.0, 0.99381, 0.0, 0.0, 2.894, 0.0, 0.0,
    0.79051, 0.0, 0.0, 0.68873, 0.0, 0.0, 0.0, 2.5909, 0.0, 0.0, 2.4911, 0.0,
    0.0, 2.3913, 0.0, 0.0, 2.2915, 0.0, 0.0, 2.1916, 0.0, 0.0, 2.0918, 0.0, 0.0,
    1.992, 0.0, 0.0, 1.8922, 0.0, 0.0, 1.7924, 0.0, 0.0, 1.6925, 0.0, 0.0,
    1.5927, 0.0, 0.0, 1.4929, 0.0, 0.0, 1.3931, 0.0, 0.0, 1.2933, 0.0, 0.0,
    1.1934, 0.0, 0.0, 1.0936, 0.0, 0.0, 0.99381, 0.0, 0.0, 2.894, 0.0, 0.0,
    0.79051, 0.0, 0.0, 0.68873, 2.1853, 0.0, 0.0, 2.1033, 0.0, 0.0, 2.0212, 0.0,
    0.0, 1.9392, 0.0, 0.0, 1.8571, 0.0, 0.0, 1.7751, 0.0, 0.0, 1.693, 0.0, 0.0,
    1.611, 0.0, 0.0, 1.5289, 0.0, 0.0, 1.4469, 0.0, 0.0, 1.3648, 0.0, 0.0,
    1.2828, 0.0, 0.0, 1.2007, 0.0, 0.0, 1.1187, 0.0, 0.0, 1.0366, 0.0, 0.0,
    0.9546, 0.0, 0.0, 0.87255, 0.0, 0.0, 0.79051, 0.0, 0.0, 2.7085, 0.0, 0.0,
    0.62275, 0.0, 0.0, 0.0, 2.1853, 0.0, 0.0, 2.1033, 0.0, 0.0, 2.0212, 0.0, 0.0,
    1.9392, 0.0, 0.0, 1.8571, 0.0, 0.0, 1.7751, 0.0, 0.0, 1.693, 0.0, 0.0, 1.611,
    0.0, 0.0, 1.5289, 0.0, 0.0, 1.4469, 0.0, 0.0, 1.3648, 0.0, 0.0, 1.2828, 0.0,
    0.0, 1.2007, 0.0, 0.0, 1.1187, 0.0, 0.0, 1.0366, 0.0, 0.0, 0.9546, 0.0, 0.0,
    0.87255, 0.0, 0.0, 0.79051, 0.0, 0.0, 2.7085, 0.0, 0.0, 0.62275, 0.0, 0.0,
    0.0, 2.1853, 0.0, 0.0, 2.1033, 0.0, 0.0, 2.0212, 0.0, 0.0, 1.9392, 0.0, 0.0,
    1.8571, 0.0, 0.0, 1.7751, 0.0, 0.0, 1.693, 0.0, 0.0, 1.611, 0.0, 0.0, 1.5289,
    0.0, 0.0, 1.4469, 0.0, 0.0, 1.3648, 0.0, 0.0, 1.2828, 0.0, 0.0, 1.2007, 0.0,
    0.0, 1.1187, 0.0, 0.0, 1.0366, 0.0, 0.0, 0.9546, 0.0, 0.0, 0.87255, 0.0, 0.0,
    0.79051, 0.0, 0.0, 2.7085, 0.0, 0.0, 0.62275, 1.8103, 0.0, 0.0, 1.7444, 0.0,
    0.0, 1.6784, 0.0, 0.0, 1.6124, 0.0, 0.0, 1.5464, 0.0, 0.0, 1.4805, 0.0, 0.0,
    1.4145, 0.0, 0.0, 1.3485, 0.0, 0.0, 1.2825, 0.0, 0.0, 1.2165, 0.0, 0.0,
    1.1506, 0.0, 0.0, 1.0846, 0.0, 0.0, 1.0186, 0.0, 0.0, 0.95264, 0.0, 0.0,
    0.88666, 0.0, 0.0, 0.82068, 0.0, 0.0, 0.75471, 0.0, 0.0, 0.68873, 0.0, 0.0,
    0.62275, 0.0, 0.0, 2.5568, 0.0, 0.0, 0.0, 1.8103, 0.0, 0.0, 1.7444, 0.0, 0.0,
    1.6784, 0.0, 0.0, 1.6124, 0.0, 0.0, 1.5464, 0.0, 0.0, 1.4805, 0.0, 0.0,
    1.4145, 0.0, 0.0, 1.3485, 0.0, 0.0, 1.2825, 0.0, 0.0, 1.2165, 0.0, 0.0,
    1.1506, 0.0, 0.0, 1.0846, 0.0, 0.0, 1.0186, 0.0, 0.0, 0.95264, 0.0, 0.0,
    0.88666, 0.0, 0.0, 0.82068, 0.0, 0.0, 0.75471, 0.0, 0.0, 0.68873, 0.0, 0.0,
    0.62275, 0.0, 0.0, 2.5568, 0.0, 0.0, 0.0, 1.8103, 0.0, 0.0, 1.7444, 0.0, 0.0,
    1.6784, 0.0, 0.0, 1.6124, 0.0, 0.0, 1.5464, 0.0, 0.0, 1.4805, 0.0, 0.0,
    1.4145, 0.0, 0.0, 1.3485, 0.0, 0.0, 1.2825, 0.0, 0.0, 1.2165, 0.0, 0.0,
    1.1506, 0.0, 0.0, 1.0846, 0.0, 0.0, 1.0186, 0.0, 0.0, 0.95264, 0.0, 0.0,
    0.88666, 0.0, 0.0, 0.82068, 0.0, 0.0, 0.75471, 0.0, 0.0, 0.68873, 0.0, 0.0,
    0.62275, 0.0, 0.0, 2.5568 };

  static const double b_a[3600] = { -14.5253, -0.0, -0.0, -11.8592, -0.0, -0.0,
    -11.1948, -0.0, -0.0, -10.5338, -0.0, -0.0, -9.8779, -0.0, -0.0, -9.2289,
    -0.0, -0.0, -8.5883, -0.0, -0.0, -7.958, -0.0, -0.0, -7.3395, -0.0, -0.0,
    -6.7347, -0.0, -0.0, -6.1452, -0.0, -0.0, -5.5728, -0.0, -0.0, -5.0191, -0.0,
    -0.0, -4.4858, -0.0, -0.0, -3.9746, -0.0, -0.0, -3.4873, -0.0, -0.0, -3.0255,
    -0.0, -0.0, -2.5909, -0.0, -0.0, -2.1853, -0.0, -0.0, -1.8103, -0.0, -0.0,
    -0.0, -14.5253, -0.0, -0.0, -11.8592, -0.0, -0.0, -11.1948, -0.0, -0.0,
    -10.5338, -0.0, -0.0, -9.8779, -0.0, -0.0, -9.2289, -0.0, -0.0, -8.5883,
    -0.0, -0.0, -7.958, -0.0, -0.0, -7.3395, -0.0, -0.0, -6.7347, -0.0, -0.0,
    -6.1452, -0.0, -0.0, -5.5728, -0.0, -0.0, -5.0191, -0.0, -0.0, -4.4858, -0.0,
    -0.0, -3.9746, -0.0, -0.0, -3.4873, -0.0, -0.0, -3.0255, -0.0, -0.0, -2.5909,
    -0.0, -0.0, -2.1853, -0.0, -0.0, -1.8103, -0.0, -0.0, -0.0, -14.5253, -0.0,
    -0.0, -11.8592, -0.0, -0.0, -11.1948, -0.0, -0.0, -10.5338, -0.0, -0.0,
    -9.8779, -0.0, -0.0, -9.2289, -0.0, -0.0, -8.5883, -0.0, -0.0, -7.958, -0.0,
    -0.0, -7.3395, -0.0, -0.0, -6.7347, -0.0, -0.0, -6.1452, -0.0, -0.0, -5.5728,
    -0.0, -0.0, -5.0191, -0.0, -0.0, -4.4858, -0.0, -0.0, -3.9746, -0.0, -0.0,
    -3.4873, -0.0, -0.0, -3.0255, -0.0, -0.0, -2.5909, -0.0, -0.0, -2.1853, -0.0,
    -0.0, -1.8103, -11.8592, -0.0, -0.0, -13.2435, -0.0, -0.0, -10.6241, -0.0,
    -0.0, -10.0064, -0.0, -0.0, -9.3921, -0.0, -0.0, -8.783, -0.0, -0.0, -8.1806,
    -0.0, -0.0, -7.5868, -0.0, -0.0, -7.0031, -0.0, -0.0, -6.4314, -0.0, -0.0,
    -5.8733, -0.0, -0.0, -5.3306, -0.0, -0.0, -4.8048, -0.0, -0.0, -4.2978, -0.0,
    -0.0, -3.8112, -0.0, -0.0, -3.3468, -0.0, -0.0, -2.9062, -0.0, -0.0, -2.4911,
    -0.0, -0.0, -2.1033, -0.0, -0.0, -1.7444, -0.0, -0.0, -0.0, -11.8592, -0.0,
    -0.0, -13.2435, -0.0, -0.0, -10.6241, -0.0, -0.0, -10.0064, -0.0, -0.0,
    -9.3921, -0.0, -0.0, -8.783, -0.0, -0.0, -8.1806, -0.0, -0.0, -7.5868, -0.0,
    -0.0, -7.0031, -0.0, -0.0, -6.4314, -0.0, -0.0, -5.8733, -0.0, -0.0, -5.3306,
    -0.0, -0.0, -4.8048, -0.0, -0.0, -4.2978, -0.0, -0.0, -3.8112, -0.0, -0.0,
    -3.3468, -0.0, -0.0, -2.9062, -0.0, -0.0, -2.4911, -0.0, -0.0, -2.1033, -0.0,
    -0.0, -1.7444, -0.0, -0.0, -0.0, -11.8592, -0.0, -0.0, -13.2435, -0.0, -0.0,
    -10.6241, -0.0, -0.0, -10.0064, -0.0, -0.0, -9.3921, -0.0, -0.0, -8.783,
    -0.0, -0.0, -8.1806, -0.0, -0.0, -7.5868, -0.0, -0.0, -7.0031, -0.0, -0.0,
    -6.4314, -0.0, -0.0, -5.8733, -0.0, -0.0, -5.3306, -0.0, -0.0, -4.8048, -0.0,
    -0.0, -4.2978, -0.0, -0.0, -3.8112, -0.0, -0.0, -3.3468, -0.0, -0.0, -2.9062,
    -0.0, -0.0, -2.4911, -0.0, -0.0, -2.1033, -0.0, -0.0, -1.7444, -11.1948,
    -0.0, -0.0, -10.6241, -0.0, -0.0, -12.0534, -0.0, -0.0, -9.479, -0.0, -0.0,
    -8.9063, -0.0, -0.0, -8.3371, -0.0, -0.0, -7.7729, -0.0, -0.0, -7.2156, -0.0,
    -0.0, -6.6668, -0.0, -0.0, -6.1281, -0.0, -0.0, -5.6014, -0.0, -0.0, -5.0884,
    -0.0, -0.0, -4.5906, -0.0, -0.0, -4.1099, -0.0, -0.0, -3.6479, -0.0, -0.0,
    -3.2063, -0.0, -0.0, -2.7869, -0.0, -0.0, -2.3913, -0.0, -0.0, -2.0212, -0.0,
    -0.0, -1.6784, -0.0, -0.0, -0.0, -11.1948, -0.0, -0.0, -10.6241, -0.0, -0.0,
    -12.0534, -0.0, -0.0, -9.479, -0.0, -0.0, -8.9063, -0.0, -0.0, -8.3371, -0.0,
    -0.0, -7.7729, -0.0, -0.0, -7.2156, -0.0, -0.0, -6.6668, -0.0, -0.0, -6.1281,
    -0.0, -0.0, -5.6014, -0.0, -0.0, -5.0884, -0.0, -0.0, -4.5906, -0.0, -0.0,
    -4.1099, -0.0, -0.0, -3.6479, -0.0, -0.0, -3.2063, -0.0, -0.0, -2.7869, -0.0,
    -0.0, -2.3913, -0.0, -0.0, -2.0212, -0.0, -0.0, -1.6784, -0.0, -0.0, -0.0,
    -11.1948, -0.0, -0.0, -10.6241, -0.0, -0.0, -12.0534, -0.0, -0.0, -9.479,
    -0.0, -0.0, -8.9063, -0.0, -0.0, -8.3371, -0.0, -0.0, -7.7729, -0.0, -0.0,
    -7.2156, -0.0, -0.0, -6.6668, -0.0, -0.0, -6.1281, -0.0, -0.0, -5.6014, -0.0,
    -0.0, -5.0884, -0.0, -0.0, -4.5906, -0.0, -0.0, -4.1099, -0.0, -0.0, -3.6479,
    -0.0, -0.0, -3.2063, -0.0, -0.0, -2.7869, -0.0, -0.0, -2.3913, -0.0, -0.0,
    -2.0212, -0.0, -0.0, -1.6784, -10.5338, -0.0, -0.0, -10.0064, -0.0, -0.0,
    -9.479, -0.0, -0.0, -10.9516, -0.0, -0.0, -8.4206, -0.0, -0.0, -7.8912, -0.0,
    -0.0, -7.3653, -0.0, -0.0, -6.8444, -0.0, -0.0, -6.3304, -0.0, -0.0, -5.8249,
    -0.0, -0.0, -5.3296, -0.0, -0.0, -4.8462, -0.0, -0.0, -4.3764, -0.0, -0.0,
    -3.922, -0.0, -0.0, -3.4845, -0.0, -0.0, -3.0658, -0.0, -0.0, -2.6676, -0.0,
    -0.0, -2.2915, -0.0, -0.0, -1.9392, -0.0, -0.0, -1.6124, -0.0, -0.0, -0.0,
    -10.5338, -0.0, -0.0, -10.0064, -0.0, -0.0, -9.479, -0.0, -0.0, -10.9516,
    -0.0, -0.0, -8.4206, -0.0, -0.0, -7.8912, -0.0, -0.0, -7.3653, -0.0, -0.0,
    -6.8444, -0.0, -0.0, -6.3304, -0.0, -0.0, -5.8249, -0.0, -0.0, -5.3296, -0.0,
    -0.0, -4.8462, -0.0, -0.0, -4.3764, -0.0, -0.0, -3.922, -0.0, -0.0, -3.4845,
    -0.0, -0.0, -3.0658, -0.0, -0.0, -2.6676, -0.0, -0.0, -2.2915, -0.0, -0.0,
    -1.9392, -0.0, -0.0, -1.6124, -0.0, -0.0, -0.0, -10.5338, -0.0, -0.0,
    -10.0064, -0.0, -0.0, -9.479, -0.0, -0.0, -10.9516, -0.0, -0.0, -8.4206,
    -0.0, -0.0, -7.8912, -0.0, -0.0, -7.3653, -0.0, -0.0, -6.8444, -0.0, -0.0,
    -6.3304, -0.0, -0.0, -5.8249, -0.0, -0.0, -5.3296, -0.0, -0.0, -4.8462, -0.0,
    -0.0, -4.3764, -0.0, -0.0, -3.922, -0.0, -0.0, -3.4845, -0.0, -0.0, -3.0658,
    -0.0, -0.0, -2.6676, -0.0, -0.0, -2.2915, -0.0, -0.0, -1.9392, -0.0, -0.0,
    -1.6124, -9.8779, -0.0, -0.0, -9.3921, -0.0, -0.0, -8.9063, -0.0, -0.0,
    -8.4206, -0.0, -0.0, -9.9348, -0.0, -0.0, -7.4453, -0.0, -0.0, -6.9576, -0.0,
    -0.0, -6.4732, -0.0, -0.0, -5.994, -0.0, -0.0, -5.5216, -0.0, -0.0, -5.0577,
    -0.0, -0.0, -4.604, -0.0, -0.0, -4.1622, -0.0, -0.0, -3.734, -0.0, -0.0,
    -3.3212, -0.0, -0.0, -2.9254, -0.0, -0.0, -2.5483, -0.0, -0.0, -2.1916, -0.0,
    -0.0, -1.8571, -0.0, -0.0, -1.5464, -0.0, -0.0, -0.0, -9.8779, -0.0, -0.0,
    -9.3921, -0.0, -0.0, -8.9063, -0.0, -0.0, -8.4206, -0.0, -0.0, -9.9348, -0.0,
    -0.0, -7.4453, -0.0, -0.0, -6.9576, -0.0, -0.0, -6.4732, -0.0, -0.0, -5.994,
    -0.0, -0.0, -5.5216, -0.0, -0.0, -5.0577, -0.0, -0.0, -4.604, -0.0, -0.0,
    -4.1622, -0.0, -0.0, -3.734, -0.0, -0.0, -3.3212, -0.0, -0.0, -2.9254, -0.0,
    -0.0, -2.5483, -0.0, -0.0, -2.1916, -0.0, -0.0, -1.8571, -0.0, -0.0, -1.5464,
    -0.0, -0.0, -0.0, -9.8779, -0.0, -0.0, -9.3921, -0.0, -0.0, -8.9063, -0.0,
    -0.0, -8.4206, -0.0, -0.0, -9.9348, -0.0, -0.0, -7.4453, -0.0, -0.0, -6.9576,
    -0.0, -0.0, -6.4732, -0.0, -0.0, -5.994, -0.0, -0.0, -5.5216, -0.0, -0.0,
    -5.0577, -0.0, -0.0, -4.604, -0.0, -0.0, -4.1622, -0.0, -0.0, -3.734, -0.0,
    -0.0, -3.3212, -0.0, -0.0, -2.9254, -0.0, -0.0, -2.5483, -0.0, -0.0, -2.1916,
    -0.0, -0.0, -1.8571, -0.0, -0.0, -1.5464, -9.2289, -0.0, -0.0, -8.783, -0.0,
    -0.0, -8.3371, -0.0, -0.0, -7.8912, -0.0, -0.0, -7.4453, -0.0, -0.0, -8.9994,
    -0.0, -0.0, -6.5499, -0.0, -0.0, -6.102, -0.0, -0.0, -5.6576, -0.0, -0.0,
    -5.2183, -0.0, -0.0, -4.7858, -0.0, -0.0, -4.3618, -0.0, -0.0, -3.948, -0.0,
    -0.0, -3.5461, -0.0, -0.0, -3.1578, -0.0, -0.0, -2.7849, -0.0, -0.0, -2.429,
    -0.0, -0.0, -2.0918, -0.0, -0.0, -1.7751, -0.0, -0.0, -1.4805, -0.0, -0.0,
    -0.0, -9.2289, -0.0, -0.0, -8.783, -0.0, -0.0, -8.3371, -0.0, -0.0, -7.8912,
    -0.0, -0.0, -7.4453, -0.0, -0.0, -8.9994, -0.0, -0.0, -6.5499, -0.0, -0.0,
    -6.102, -0.0, -0.0, -5.6576, -0.0, -0.0, -5.2183, -0.0, -0.0, -4.7858, -0.0,
    -0.0, -4.3618, -0.0, -0.0, -3.948, -0.0, -0.0, -3.5461, -0.0, -0.0, -3.1578,
    -0.0, -0.0, -2.7849, -0.0, -0.0, -2.429, -0.0, -0.0, -2.0918, -0.0, -0.0,
    -1.7751, -0.0, -0.0, -1.4805, -0.0, -0.0, -0.0, -9.2289, -0.0, -0.0, -8.783,
    -0.0, -0.0, -8.3371, -0.0, -0.0, -7.8912, -0.0, -0.0, -7.4453, -0.0, -0.0,
    -8.9994, -0.0, -0.0, -6.5499, -0.0, -0.0, -6.102, -0.0, -0.0, -5.6576, -0.0,
    -0.0, -5.2183, -0.0, -0.0, -4.7858, -0.0, -0.0, -4.3618, -0.0, -0.0, -3.948,
    -0.0, -0.0, -3.5461, -0.0, -0.0, -3.1578, -0.0, -0.0, -2.7849, -0.0, -0.0,
    -2.429, -0.0, -0.0, -2.0918, -0.0, -0.0, -1.7751, -0.0, -0.0, -1.4805,
    -8.5883, -0.0, -0.0, -8.1806, -0.0, -0.0, -7.7729, -0.0, -0.0, -7.3653, -0.0,
    -0.0, -6.9576, -0.0, -0.0, -6.5499, -0.0, -0.0, -8.1422, -0.0, -0.0, -5.7309,
    -0.0, -0.0, -5.3212, -0.0, -0.0, -4.915, -0.0, -0.0, -4.5139, -0.0, -0.0,
    -4.1195, -0.0, -0.0, -3.7337, -0.0, -0.0, -3.3582, -0.0, -0.0, -2.9945, -0.0,
    -0.0, -2.6444, -0.0, -0.0, -2.3097, -0.0, -0.0, -1.992, -0.0, -0.0, -1.693,
    -0.0, -0.0, -1.4145, -0.0, -0.0, -0.0, -8.5883, -0.0, -0.0, -8.1806, -0.0,
    -0.0, -7.7729, -0.0, -0.0, -7.3653, -0.0, -0.0, -6.9576, -0.0, -0.0, -6.5499,
    -0.0, -0.0, -8.1422, -0.0, -0.0, -5.7309, -0.0, -0.0, -5.3212, -0.0, -0.0,
    -4.915, -0.0, -0.0, -4.5139, -0.0, -0.0, -4.1195, -0.0, -0.0, -3.7337, -0.0,
    -0.0, -3.3582, -0.0, -0.0, -2.9945, -0.0, -0.0, -2.6444, -0.0, -0.0, -2.3097,
    -0.0, -0.0, -1.992, -0.0, -0.0, -1.693, -0.0, -0.0, -1.4145, -0.0, -0.0,
    -0.0, -8.5883, -0.0, -0.0, -8.1806, -0.0, -0.0, -7.7729, -0.0, -0.0, -7.3653,
    -0.0, -0.0, -6.9576, -0.0, -0.0, -6.5499, -0.0, -0.0, -8.1422, -0.0, -0.0,
    -5.7309, -0.0, -0.0, -5.3212, -0.0, -0.0, -4.915, -0.0, -0.0, -4.5139, -0.0,
    -0.0, -4.1195, -0.0, -0.0, -3.7337, -0.0, -0.0, -3.3582, -0.0, -0.0, -2.9945,
    -0.0, -0.0, -2.6444, -0.0, -0.0, -2.3097, -0.0, -0.0, -1.992, -0.0, -0.0,
    -1.693, -0.0, -0.0, -1.4145, -7.958, -0.0, -0.0, -7.5868, -0.0, -0.0,
    -7.2156, -0.0, -0.0, -6.8444, -0.0, -0.0, -6.4732, -0.0, -0.0, -6.102, -0.0,
    -0.0, -5.7309, -0.0, -0.0, -7.3597, -0.0, -0.0, -4.9848, -0.0, -0.0, -4.6117,
    -0.0, -0.0, -4.242, -0.0, -0.0, -3.8773, -0.0, -0.0, -3.5195, -0.0, -0.0,
    -3.1702, -0.0, -0.0, -2.8311, -0.0, -0.0, -2.504, -0.0, -0.0, -2.1904, -0.0,
    -0.0, -1.8922, -0.0, -0.0, -1.611, -0.0, -0.0, -1.3485, -0.0, -0.0, -0.0,
    -7.958, -0.0, -0.0, -7.5868, -0.0, -0.0, -7.2156, -0.0, -0.0, -6.8444, -0.0,
    -0.0, -6.4732, -0.0, -0.0, -6.102, -0.0, -0.0, -5.7309, -0.0, -0.0, -7.3597,
    -0.0, -0.0, -4.9848, -0.0, -0.0, -4.6117, -0.0, -0.0, -4.242, -0.0, -0.0,
    -3.8773, -0.0, -0.0, -3.5195, -0.0, -0.0, -3.1702, -0.0, -0.0, -2.8311, -0.0,
    -0.0, -2.504, -0.0, -0.0, -2.1904, -0.0, -0.0, -1.8922, -0.0, -0.0, -1.611,
    -0.0, -0.0, -1.3485, -0.0, -0.0, -0.0, -7.958, -0.0, -0.0, -7.5868, -0.0,
    -0.0, -7.2156, -0.0, -0.0, -6.8444, -0.0, -0.0, -6.4732, -0.0, -0.0, -6.102,
    -0.0, -0.0, -5.7309, -0.0, -0.0, -7.3597, -0.0, -0.0, -4.9848, -0.0, -0.0,
    -4.6117, -0.0, -0.0, -4.242, -0.0, -0.0, -3.8773, -0.0, -0.0, -3.5195, -0.0,
    -0.0, -3.1702, -0.0, -0.0, -2.8311, -0.0, -0.0, -2.504, -0.0, -0.0, -2.1904,
    -0.0, -0.0, -1.8922, -0.0, -0.0, -1.611, -0.0, -0.0, -1.3485, -7.3395, -0.0,
    -0.0, -7.0031, -0.0, -0.0, -6.6668, -0.0, -0.0, -6.3304, -0.0, -0.0, -5.994,
    -0.0, -0.0, -5.6576, -0.0, -0.0, -5.3212, -0.0, -0.0, -4.9848, -0.0, -0.0,
    -6.6485, -0.0, -0.0, -4.3084, -0.0, -0.0, -3.9701, -0.0, -0.0, -3.6351, -0.0,
    -0.0, -3.3053, -0.0, -0.0, -2.9823, -0.0, -0.0, -2.6678, -0.0, -0.0, -2.3635,
    -0.0, -0.0, -2.0711, -0.0, -0.0, -1.7924, -0.0, -0.0, -1.5289, -0.0, -0.0,
    -1.2825, -0.0, -0.0, -0.0, -7.3395, -0.0, -0.0, -7.0031, -0.0, -0.0, -6.6668,
    -0.0, -0.0, -6.3304, -0.0, -0.0, -5.994, -0.0, -0.0, -5.6576, -0.0, -0.0,
    -5.3212, -0.0, -0.0, -4.9848, -0.0, -0.0, -6.6485, -0.0, -0.0, -4.3084, -0.0,
    -0.0, -3.9701, -0.0, -0.0, -3.6351, -0.0, -0.0, -3.3053, -0.0, -0.0, -2.9823,
    -0.0, -0.0, -2.6678, -0.0, -0.0, -2.3635, -0.0, -0.0, -2.0711, -0.0, -0.0,
    -1.7924, -0.0, -0.0, -1.5289, -0.0, -0.0, -1.2825, -0.0, -0.0, -0.0, -7.3395,
    -0.0, -0.0, -7.0031, -0.0, -0.0, -6.6668, -0.0, -0.0, -6.3304, -0.0, -0.0,
    -5.994, -0.0, -0.0, -5.6576, -0.0, -0.0, -5.3212, -0.0, -0.0, -4.9848, -0.0,
    -0.0, -6.6485, -0.0, -0.0, -4.3084, -0.0, -0.0, -3.9701, -0.0, -0.0, -3.6351,
    -0.0, -0.0, -3.3053, -0.0, -0.0, -2.9823, -0.0, -0.0, -2.6678, -0.0, -0.0,
    -2.3635, -0.0, -0.0, -2.0711, -0.0, -0.0, -1.7924, -0.0, -0.0, -1.5289, -0.0,
    -0.0, -1.2825, -6.7347, -0.0, -0.0, -6.4314, -0.0, -0.0, -6.1281, -0.0, -0.0,
    -5.8249, -0.0, -0.0, -5.5216, -0.0, -0.0, -5.2183, -0.0, -0.0, -4.915, -0.0,
    -0.0, -4.6117, -0.0, -0.0, -4.3084, -0.0, -0.0, -6.0051, -0.0, -0.0, -3.6982,
    -0.0, -0.0, -3.3929, -0.0, -0.0, -3.0911, -0.0, -0.0, -2.7944, -0.0, -0.0,
    -2.5044, -0.0, -0.0, -2.223, -0.0, -0.0, -1.9518, -0.0, -0.0, -1.6925, -0.0,
    -0.0, -1.4469, -0.0, -0.0, -1.2165, -0.0, -0.0, -0.0, -6.7347, -0.0, -0.0,
    -6.4314, -0.0, -0.0, -6.1281, -0.0, -0.0, -5.8249, -0.0, -0.0, -5.5216, -0.0,
    -0.0, -5.2183, -0.0, -0.0, -4.915, -0.0, -0.0, -4.6117, -0.0, -0.0, -4.3084,
    -0.0, -0.0, -6.0051, -0.0, -0.0, -3.6982, -0.0, -0.0, -3.3929, -0.0, -0.0,
    -3.0911, -0.0, -0.0, -2.7944, -0.0, -0.0, -2.5044, -0.0, -0.0, -2.223, -0.0,
    -0.0, -1.9518, -0.0, -0.0, -1.6925, -0.0, -0.0, -1.4469, -0.0, -0.0, -1.2165,
    -0.0, -0.0, -0.0, -6.7347, -0.0, -0.0, -6.4314, -0.0, -0.0, -6.1281, -0.0,
    -0.0, -5.8249, -0.0, -0.0, -5.5216, -0.0, -0.0, -5.2183, -0.0, -0.0, -4.915,
    -0.0, -0.0, -4.6117, -0.0, -0.0, -4.3084, -0.0, -0.0, -6.0051, -0.0, -0.0,
    -3.6982, -0.0, -0.0, -3.3929, -0.0, -0.0, -3.0911, -0.0, -0.0, -2.7944, -0.0,
    -0.0, -2.5044, -0.0, -0.0, -2.223, -0.0, -0.0, -1.9518, -0.0, -0.0, -1.6925,
    -0.0, -0.0, -1.4469, -0.0, -0.0, -1.2165, -6.1452, -0.0, -0.0, -5.8733, -0.0,
    -0.0, -5.6014, -0.0, -0.0, -5.3296, -0.0, -0.0, -5.0577, -0.0, -0.0, -4.7858,
    -0.0, -0.0, -4.5139, -0.0, -0.0, -4.242, -0.0, -0.0, -3.9701, -0.0, -0.0,
    -3.6982, -0.0, -0.0, -5.4263, -0.0, -0.0, -3.1507, -0.0, -0.0, -2.8769, -0.0,
    -0.0, -2.6064, -0.0, -0.0, -2.3411, -0.0, -0.0, -2.0825, -0.0, -0.0, -1.8325,
    -0.0, -0.0, -1.5927, -0.0, -0.0, -1.3648, -0.0, -0.0, -1.1506, -0.0, -0.0,
    -0.0, -6.1452, -0.0, -0.0, -5.8733, -0.0, -0.0, -5.6014, -0.0, -0.0, -5.3296,
    -0.0, -0.0, -5.0577, -0.0, -0.0, -4.7858, -0.0, -0.0, -4.5139, -0.0, -0.0,
    -4.242, -0.0, -0.0, -3.9701, -0.0, -0.0, -3.6982, -0.0, -0.0, -5.4263, -0.0,
    -0.0, -3.1507, -0.0, -0.0, -2.8769, -0.0, -0.0, -2.6064, -0.0, -0.0, -2.3411,
    -0.0, -0.0, -2.0825, -0.0, -0.0, -1.8325, -0.0, -0.0, -1.5927, -0.0, -0.0,
    -1.3648, -0.0, -0.0, -1.1506, -0.0, -0.0, -0.0, -6.1452, -0.0, -0.0, -5.8733,
    -0.0, -0.0, -5.6014, -0.0, -0.0, -5.3296, -0.0, -0.0, -5.0577, -0.0, -0.0,
    -4.7858, -0.0, -0.0, -4.5139, -0.0, -0.0, -4.242, -0.0, -0.0, -3.9701, -0.0,
    -0.0, -3.6982, -0.0, -0.0, -5.4263, -0.0, -0.0, -3.1507, -0.0, -0.0, -2.8769,
    -0.0, -0.0, -2.6064, -0.0, -0.0, -2.3411, -0.0, -0.0, -2.0825, -0.0, -0.0,
    -1.8325, -0.0, -0.0, -1.5927, -0.0, -0.0, -1.3648, -0.0, -0.0, -1.1506,
    -5.5728, -0.0, -0.0, -5.3306, -0.0, -0.0, -5.0884, -0.0, -0.0, -4.8462, -0.0,
    -0.0, -4.604, -0.0, -0.0, -4.3618, -0.0, -0.0, -4.1195, -0.0, -0.0, -3.8773,
    -0.0, -0.0, -3.6351, -0.0, -0.0, -3.3929, -0.0, -0.0, -3.1507, -0.0, -0.0,
    -4.9085, -0.0, -0.0, -2.6626, -0.0, -0.0, -2.4185, -0.0, -0.0, -2.1777, -0.0,
    -0.0, -1.9421, -0.0, -0.0, -1.7132, -0.0, -0.0, -1.4929, -0.0, -0.0, -1.2828,
    -0.0, -0.0, -1.0846, -0.0, -0.0, -0.0, -5.5728, -0.0, -0.0, -5.3306, -0.0,
    -0.0, -5.0884, -0.0, -0.0, -4.8462, -0.0, -0.0, -4.604, -0.0, -0.0, -4.3618,
    -0.0, -0.0, -4.1195, -0.0, -0.0, -3.8773, -0.0, -0.0, -3.6351, -0.0, -0.0,
    -3.3929, -0.0, -0.0, -3.1507, -0.0, -0.0, -4.9085, -0.0, -0.0, -2.6626, -0.0,
    -0.0, -2.4185, -0.0, -0.0, -2.1777, -0.0, -0.0, -1.9421, -0.0, -0.0, -1.7132,
    -0.0, -0.0, -1.4929, -0.0, -0.0, -1.2828, -0.0, -0.0, -1.0846, -0.0, -0.0,
    -0.0, -5.5728, -0.0, -0.0, -5.3306, -0.0, -0.0, -5.0884, -0.0, -0.0, -4.8462,
    -0.0, -0.0, -4.604, -0.0, -0.0, -4.3618, -0.0, -0.0, -4.1195, -0.0, -0.0,
    -3.8773, -0.0, -0.0, -3.6351, -0.0, -0.0, -3.3929, -0.0, -0.0, -3.1507, -0.0,
    -0.0, -4.9085, -0.0, -0.0, -2.6626, -0.0, -0.0, -2.4185, -0.0, -0.0, -2.1777,
    -0.0, -0.0, -1.9421, -0.0, -0.0, -1.7132, -0.0, -0.0, -1.4929, -0.0, -0.0,
    -1.2828, -0.0, -0.0, -1.0846, -5.0191, -0.0, -0.0, -4.8048, -0.0, -0.0,
    -4.5906, -0.0, -0.0, -4.3764, -0.0, -0.0, -4.1622, -0.0, -0.0, -3.948, -0.0,
    -0.0, -3.7337, -0.0, -0.0, -3.5195, -0.0, -0.0, -3.3053, -0.0, -0.0, -3.0911,
    -0.0, -0.0, -2.8769, -0.0, -0.0, -2.6626, -0.0, -0.0, -4.4484, -0.0, -0.0,
    -2.2305, -0.0, -0.0, -2.0144, -0.0, -0.0, -1.8016, -0.0, -0.0, -1.5939, -0.0,
    -0.0, -1.3931, -0.0, -0.0, -1.2007, -0.0, -0.0, -1.0186, -0.0, -0.0, -0.0,
    -5.0191, -0.0, -0.0, -4.8048, -0.0, -0.0, -4.5906, -0.0, -0.0, -4.3764, -0.0,
    -0.0, -4.1622, -0.0, -0.0, -3.948, -0.0, -0.0, -3.7337, -0.0, -0.0, -3.5195,
    -0.0, -0.0, -3.3053, -0.0, -0.0, -3.0911, -0.0, -0.0, -2.8769, -0.0, -0.0,
    -2.6626, -0.0, -0.0, -4.4484, -0.0, -0.0, -2.2305, -0.0, -0.0, -2.0144, -0.0,
    -0.0, -1.8016, -0.0, -0.0, -1.5939, -0.0, -0.0, -1.3931, -0.0, -0.0, -1.2007,
    -0.0, -0.0, -1.0186, -0.0, -0.0, -0.0, -5.0191, -0.0, -0.0, -4.8048, -0.0,
    -0.0, -4.5906, -0.0, -0.0, -4.3764, -0.0, -0.0, -4.1622, -0.0, -0.0, -3.948,
    -0.0, -0.0, -3.7337, -0.0, -0.0, -3.5195, -0.0, -0.0, -3.3053, -0.0, -0.0,
    -3.0911, -0.0, -0.0, -2.8769, -0.0, -0.0, -2.6626, -0.0, -0.0, -4.4484, -0.0,
    -0.0, -2.2305, -0.0, -0.0, -2.0144, -0.0, -0.0, -1.8016, -0.0, -0.0, -1.5939,
    -0.0, -0.0, -1.3931, -0.0, -0.0, -1.2007, -0.0, -0.0, -1.0186, -4.4858, -0.0,
    -0.0, -4.2978, -0.0, -0.0, -4.1099, -0.0, -0.0, -3.922, -0.0, -0.0, -3.734,
    -0.0, -0.0, -3.5461, -0.0, -0.0, -3.3582, -0.0, -0.0, -3.1702, -0.0, -0.0,
    -2.9823, -0.0, -0.0, -2.7944, -0.0, -0.0, -2.6064, -0.0, -0.0, -2.4185, -0.0,
    -0.0, -2.2305, -0.0, -0.0, -4.0426, -0.0, -0.0, -1.851, -0.0, -0.0, -1.6611,
    -0.0, -0.0, -1.4746, -0.0, -0.0, -1.2933, -0.0, -0.0, -1.1187, -0.0, -0.0,
    -0.95264, -0.0, -0.0, -0.0, -4.4858, -0.0, -0.0, -4.2978, -0.0, -0.0,
    -4.1099, -0.0, -0.0, -3.922, -0.0, -0.0, -3.734, -0.0, -0.0, -3.5461, -0.0,
    -0.0, -3.3582, -0.0, -0.0, -3.1702, -0.0, -0.0, -2.9823, -0.0, -0.0, -2.7944,
    -0.0, -0.0, -2.6064, -0.0, -0.0, -2.4185, -0.0, -0.0, -2.2305, -0.0, -0.0,
    -4.0426, -0.0, -0.0, -1.851, -0.0, -0.0, -1.6611, -0.0, -0.0, -1.4746, -0.0,
    -0.0, -1.2933, -0.0, -0.0, -1.1187, -0.0, -0.0, -0.95264, -0.0, -0.0, -0.0,
    -4.4858, -0.0, -0.0, -4.2978, -0.0, -0.0, -4.1099, -0.0, -0.0, -3.922, -0.0,
    -0.0, -3.734, -0.0, -0.0, -3.5461, -0.0, -0.0, -3.3582, -0.0, -0.0, -3.1702,
    -0.0, -0.0, -2.9823, -0.0, -0.0, -2.7944, -0.0, -0.0, -2.6064, -0.0, -0.0,
    -2.4185, -0.0, -0.0, -2.2305, -0.0, -0.0, -4.0426, -0.0, -0.0, -1.851, -0.0,
    -0.0, -1.6611, -0.0, -0.0, -1.4746, -0.0, -0.0, -1.2933, -0.0, -0.0, -1.1187,
    -0.0, -0.0, -0.95264, -3.9746, -0.0, -0.0, -3.8112, -0.0, -0.0, -3.6479,
    -0.0, -0.0, -3.4845, -0.0, -0.0, -3.3212, -0.0, -0.0, -3.1578, -0.0, -0.0,
    -2.9945, -0.0, -0.0, -2.8311, -0.0, -0.0, -2.6678, -0.0, -0.0, -2.5044, -0.0,
    -0.0, -2.3411, -0.0, -0.0, -2.1777, -0.0, -0.0, -2.0144, -0.0, -0.0, -1.851,
    -0.0, -0.0, -3.6877, -0.0, -0.0, -1.5207, -0.0, -0.0, -1.3553, -0.0, -0.0,
    -1.1934, -0.0, -0.0, -1.0366, -0.0, -0.0, -0.88666, -0.0, -0.0, -0.0,
    -3.9746, -0.0, -0.0, -3.8112, -0.0, -0.0, -3.6479, -0.0, -0.0, -3.4845, -0.0,
    -0.0, -3.3212, -0.0, -0.0, -3.1578, -0.0, -0.0, -2.9945, -0.0, -0.0, -2.8311,
    -0.0, -0.0, -2.6678, -0.0, -0.0, -2.5044, -0.0, -0.0, -2.3411, -0.0, -0.0,
    -2.1777, -0.0, -0.0, -2.0144, -0.0, -0.0, -1.851, -0.0, -0.0, -3.6877, -0.0,
    -0.0, -1.5207, -0.0, -0.0, -1.3553, -0.0, -0.0, -1.1934, -0.0, -0.0, -1.0366,
    -0.0, -0.0, -0.88666, -0.0, -0.0, -0.0, -3.9746, -0.0, -0.0, -3.8112, -0.0,
    -0.0, -3.6479, -0.0, -0.0, -3.4845, -0.0, -0.0, -3.3212, -0.0, -0.0, -3.1578,
    -0.0, -0.0, -2.9945, -0.0, -0.0, -2.8311, -0.0, -0.0, -2.6678, -0.0, -0.0,
    -2.5044, -0.0, -0.0, -2.3411, -0.0, -0.0, -2.1777, -0.0, -0.0, -2.0144, -0.0,
    -0.0, -1.851, -0.0, -0.0, -3.6877, -0.0, -0.0, -1.5207, -0.0, -0.0, -1.3553,
    -0.0, -0.0, -1.1934, -0.0, -0.0, -1.0366, -0.0, -0.0, -0.88666, -3.4873,
    -0.0, -0.0, -3.3468, -0.0, -0.0, -3.2063, -0.0, -0.0, -3.0658, -0.0, -0.0,
    -2.9254, -0.0, -0.0, -2.7849, -0.0, -0.0, -2.6444, -0.0, -0.0, -2.504, -0.0,
    -0.0, -2.3635, -0.0, -0.0, -2.223, -0.0, -0.0, -2.0825, -0.0, -0.0, -1.9421,
    -0.0, -0.0, -1.8016, -0.0, -0.0, -1.6611, -0.0, -0.0, -1.5207, -0.0, -0.0,
    -3.3802, -0.0, -0.0, -1.2361, -0.0, -0.0, -1.0936, -0.0, -0.0, -0.9546, -0.0,
    -0.0, -0.82068, -0.0, -0.0, -0.0, -3.4873, -0.0, -0.0, -3.3468, -0.0, -0.0,
    -3.2063, -0.0, -0.0, -3.0658, -0.0, -0.0, -2.9254, -0.0, -0.0, -2.7849, -0.0,
    -0.0, -2.6444, -0.0, -0.0, -2.504, -0.0, -0.0, -2.3635, -0.0, -0.0, -2.223,
    -0.0, -0.0, -2.0825, -0.0, -0.0, -1.9421, -0.0, -0.0, -1.8016, -0.0, -0.0,
    -1.6611, -0.0, -0.0, -1.5207, -0.0, -0.0, -3.3802, -0.0, -0.0, -1.2361, -0.0,
    -0.0, -1.0936, -0.0, -0.0, -0.9546, -0.0, -0.0, -0.82068, -0.0, -0.0, -0.0,
    -3.4873, -0.0, -0.0, -3.3468, -0.0, -0.0, -3.2063, -0.0, -0.0, -3.0658, -0.0,
    -0.0, -2.9254, -0.0, -0.0, -2.7849, -0.0, -0.0, -2.6444, -0.0, -0.0, -2.504,
    -0.0, -0.0, -2.3635, -0.0, -0.0, -2.223, -0.0, -0.0, -2.0825, -0.0, -0.0,
    -1.9421, -0.0, -0.0, -1.8016, -0.0, -0.0, -1.6611, -0.0, -0.0, -1.5207, -0.0,
    -0.0, -3.3802, -0.0, -0.0, -1.2361, -0.0, -0.0, -1.0936, -0.0, -0.0, -0.9546,
    -0.0, -0.0, -0.82068, -3.0255, -0.0, -0.0, -2.9062, -0.0, -0.0, -2.7869,
    -0.0, -0.0, -2.6676, -0.0, -0.0, -2.5483, -0.0, -0.0, -2.429, -0.0, -0.0,
    -2.3097, -0.0, -0.0, -2.1904, -0.0, -0.0, -2.0711, -0.0, -0.0, -1.9518, -0.0,
    -0.0, -1.8325, -0.0, -0.0, -1.7132, -0.0, -0.0, -1.5939, -0.0, -0.0, -1.4746,
    -0.0, -0.0, -1.3553, -0.0, -0.0, -1.2361, -0.0, -0.0, -3.1168, -0.0, -0.0,
    -0.99381, -0.0, -0.0, -0.87255, -0.0, -0.0, -0.75471, -0.0, -0.0, -0.0,
    -3.0255, -0.0, -0.0, -2.9062, -0.0, -0.0, -2.7869, -0.0, -0.0, -2.6676, -0.0,
    -0.0, -2.5483, -0.0, -0.0, -2.429, -0.0, -0.0, -2.3097, -0.0, -0.0, -2.1904,
    -0.0, -0.0, -2.0711, -0.0, -0.0, -1.9518, -0.0, -0.0, -1.8325, -0.0, -0.0,
    -1.7132, -0.0, -0.0, -1.5939, -0.0, -0.0, -1.4746, -0.0, -0.0, -1.3553, -0.0,
    -0.0, -1.2361, -0.0, -0.0, -3.1168, -0.0, -0.0, -0.99381, -0.0, -0.0,
    -0.87255, -0.0, -0.0, -0.75471, -0.0, -0.0, -0.0, -3.0255, -0.0, -0.0,
    -2.9062, -0.0, -0.0, -2.7869, -0.0, -0.0, -2.6676, -0.0, -0.0, -2.5483, -0.0,
    -0.0, -2.429, -0.0, -0.0, -2.3097, -0.0, -0.0, -2.1904, -0.0, -0.0, -2.0711,
    -0.0, -0.0, -1.9518, -0.0, -0.0, -1.8325, -0.0, -0.0, -1.7132, -0.0, -0.0,
    -1.5939, -0.0, -0.0, -1.4746, -0.0, -0.0, -1.3553, -0.0, -0.0, -1.2361, -0.0,
    -0.0, -3.1168, -0.0, -0.0, -0.99381, -0.0, -0.0, -0.87255, -0.0, -0.0,
    -0.75471, -2.5909, -0.0, -0.0, -2.4911, -0.0, -0.0, -2.3913, -0.0, -0.0,
    -2.2915, -0.0, -0.0, -2.1916, -0.0, -0.0, -2.0918, -0.0, -0.0, -1.992, -0.0,
    -0.0, -1.8922, -0.0, -0.0, -1.7924, -0.0, -0.0, -1.6925, -0.0, -0.0, -1.5927,
    -0.0, -0.0, -1.4929, -0.0, -0.0, -1.3931, -0.0, -0.0, -1.2933, -0.0, -0.0,
    -1.1934, -0.0, -0.0, -1.0936, -0.0, -0.0, -0.99381, -0.0, -0.0, -2.894, -0.0,
    -0.0, -0.79051, -0.0, -0.0, -0.68873, -0.0, -0.0, -0.0, -2.5909, -0.0, -0.0,
    -2.4911, -0.0, -0.0, -2.3913, -0.0, -0.0, -2.2915, -0.0, -0.0, -2.1916, -0.0,
    -0.0, -2.0918, -0.0, -0.0, -1.992, -0.0, -0.0, -1.8922, -0.0, -0.0, -1.7924,
    -0.0, -0.0, -1.6925, -0.0, -0.0, -1.5927, -0.0, -0.0, -1.4929, -0.0, -0.0,
    -1.3931, -0.0, -0.0, -1.2933, -0.0, -0.0, -1.1934, -0.0, -0.0, -1.0936, -0.0,
    -0.0, -0.99381, -0.0, -0.0, -2.894, -0.0, -0.0, -0.79051, -0.0, -0.0,
    -0.68873, -0.0, -0.0, -0.0, -2.5909, -0.0, -0.0, -2.4911, -0.0, -0.0,
    -2.3913, -0.0, -0.0, -2.2915, -0.0, -0.0, -2.1916, -0.0, -0.0, -2.0918, -0.0,
    -0.0, -1.992, -0.0, -0.0, -1.8922, -0.0, -0.0, -1.7924, -0.0, -0.0, -1.6925,
    -0.0, -0.0, -1.5927, -0.0, -0.0, -1.4929, -0.0, -0.0, -1.3931, -0.0, -0.0,
    -1.2933, -0.0, -0.0, -1.1934, -0.0, -0.0, -1.0936, -0.0, -0.0, -0.99381,
    -0.0, -0.0, -2.894, -0.0, -0.0, -0.79051, -0.0, -0.0, -0.68873, -2.1853,
    -0.0, -0.0, -2.1033, -0.0, -0.0, -2.0212, -0.0, -0.0, -1.9392, -0.0, -0.0,
    -1.8571, -0.0, -0.0, -1.7751, -0.0, -0.0, -1.693, -0.0, -0.0, -1.611, -0.0,
    -0.0, -1.5289, -0.0, -0.0, -1.4469, -0.0, -0.0, -1.3648, -0.0, -0.0, -1.2828,
    -0.0, -0.0, -1.2007, -0.0, -0.0, -1.1187, -0.0, -0.0, -1.0366, -0.0, -0.0,
    -0.9546, -0.0, -0.0, -0.87255, -0.0, -0.0, -0.79051, -0.0, -0.0, -2.7085,
    -0.0, -0.0, -0.62275, -0.0, -0.0, -0.0, -2.1853, -0.0, -0.0, -2.1033, -0.0,
    -0.0, -2.0212, -0.0, -0.0, -1.9392, -0.0, -0.0, -1.8571, -0.0, -0.0, -1.7751,
    -0.0, -0.0, -1.693, -0.0, -0.0, -1.611, -0.0, -0.0, -1.5289, -0.0, -0.0,
    -1.4469, -0.0, -0.0, -1.3648, -0.0, -0.0, -1.2828, -0.0, -0.0, -1.2007, -0.0,
    -0.0, -1.1187, -0.0, -0.0, -1.0366, -0.0, -0.0, -0.9546, -0.0, -0.0,
    -0.87255, -0.0, -0.0, -0.79051, -0.0, -0.0, -2.7085, -0.0, -0.0, -0.62275,
    -0.0, -0.0, -0.0, -2.1853, -0.0, -0.0, -2.1033, -0.0, -0.0, -2.0212, -0.0,
    -0.0, -1.9392, -0.0, -0.0, -1.8571, -0.0, -0.0, -1.7751, -0.0, -0.0, -1.693,
    -0.0, -0.0, -1.611, -0.0, -0.0, -1.5289, -0.0, -0.0, -1.4469, -0.0, -0.0,
    -1.3648, -0.0, -0.0, -1.2828, -0.0, -0.0, -1.2007, -0.0, -0.0, -1.1187, -0.0,
    -0.0, -1.0366, -0.0, -0.0, -0.9546, -0.0, -0.0, -0.87255, -0.0, -0.0,
    -0.79051, -0.0, -0.0, -2.7085, -0.0, -0.0, -0.62275, -1.8103, -0.0, -0.0,
    -1.7444, -0.0, -0.0, -1.6784, -0.0, -0.0, -1.6124, -0.0, -0.0, -1.5464, -0.0,
    -0.0, -1.4805, -0.0, -0.0, -1.4145, -0.0, -0.0, -1.3485, -0.0, -0.0, -1.2825,
    -0.0, -0.0, -1.2165, -0.0, -0.0, -1.1506, -0.0, -0.0, -1.0846, -0.0, -0.0,
    -1.0186, -0.0, -0.0, -0.95264, -0.0, -0.0, -0.88666, -0.0, -0.0, -0.82068,
    -0.0, -0.0, -0.75471, -0.0, -0.0, -0.68873, -0.0, -0.0, -0.62275, -0.0, -0.0,
    -2.5568, -0.0, -0.0, -0.0, -1.8103, -0.0, -0.0, -1.7444, -0.0, -0.0, -1.6784,
    -0.0, -0.0, -1.6124, -0.0, -0.0, -1.5464, -0.0, -0.0, -1.4805, -0.0, -0.0,
    -1.4145, -0.0, -0.0, -1.3485, -0.0, -0.0, -1.2825, -0.0, -0.0, -1.2165, -0.0,
    -0.0, -1.1506, -0.0, -0.0, -1.0846, -0.0, -0.0, -1.0186, -0.0, -0.0,
    -0.95264, -0.0, -0.0, -0.88666, -0.0, -0.0, -0.82068, -0.0, -0.0, -0.75471,
    -0.0, -0.0, -0.68873, -0.0, -0.0, -0.62275, -0.0, -0.0, -2.5568, -0.0, -0.0,
    -0.0, -1.8103, -0.0, -0.0, -1.7444, -0.0, -0.0, -1.6784, -0.0, -0.0, -1.6124,
    -0.0, -0.0, -1.5464, -0.0, -0.0, -1.4805, -0.0, -0.0, -1.4145, -0.0, -0.0,
    -1.3485, -0.0, -0.0, -1.2825, -0.0, -0.0, -1.2165, -0.0, -0.0, -1.1506, -0.0,
    -0.0, -1.0846, -0.0, -0.0, -1.0186, -0.0, -0.0, -0.95264, -0.0, -0.0,
    -0.88666, -0.0, -0.0, -0.82068, -0.0, -0.0, -0.75471, -0.0, -0.0, -0.68873,
    -0.0, -0.0, -0.62275, -0.0, -0.0, -2.5568 };

  static const double b[360] = { 25.3845, 0.0, 0.0, 246.3304, 0.0, 0.0, 0.0,
    25.3845, 0.0, 0.0, 246.3304, 0.0, 0.0, 0.0, 25.3845, 0.0, 0.0, 246.3304,
    23.5944, 0.0, 0.0, 233.1205, 0.0, 0.0, 0.0, 23.5944, 0.0, 0.0, 233.1205, 0.0,
    0.0, 0.0, 23.5944, 0.0, 0.0, 233.1205, 21.8695, 0.0, 0.0, 219.9595, 0.0, 0.0,
    0.0, 21.8695, 0.0, 0.0, 219.9595, 0.0, 0.0, 0.0, 21.8695, 0.0, 0.0, 219.9595,
    20.2098, 0.0, 0.0, 206.8801, 0.0, 0.0, 0.0, 20.2098, 0.0, 0.0, 206.8801, 0.0,
    0.0, 0.0, 20.2098, 0.0, 0.0, 206.8801, 18.6154, 0.0, 0.0, 193.9148, 0.0, 0.0,
    0.0, 18.6154, 0.0, 0.0, 193.9148, 0.0, 0.0, 0.0, 18.6154, 0.0, 0.0, 193.9148,
    17.0863, 0.0, 0.0, 181.0964, 0.0, 0.0, 0.0, 17.0863, 0.0, 0.0, 181.0964, 0.0,
    0.0, 0.0, 17.0863, 0.0, 0.0, 181.0964, 15.6224, 0.0, 0.0, 168.4573, 0.0, 0.0,
    0.0, 15.6224, 0.0, 0.0, 168.4573, 0.0, 0.0, 0.0, 15.6224, 0.0, 0.0, 168.4573,
    14.2237, 0.0, 0.0, 156.0303, 0.0, 0.0, 0.0, 14.2237, 0.0, 0.0, 156.0303, 0.0,
    0.0, 0.0, 14.2237, 0.0, 0.0, 156.0303, 12.8902, 0.0, 0.0, 143.8479, 0.0, 0.0,
    0.0, 12.8902, 0.0, 0.0, 143.8479, 0.0, 0.0, 0.0, 12.8902, 0.0, 0.0, 143.8479,
    11.622, 0.0, 0.0, 131.9428, 0.0, 0.0, 0.0, 11.622, 0.0, 0.0, 131.9428, 0.0,
    0.0, 0.0, 11.622, 0.0, 0.0, 131.9428, 10.4191, 0.0, 0.0, 120.3476, 0.0, 0.0,
    0.0, 10.4191, 0.0, 0.0, 120.3476, 0.0, 0.0, 0.0, 10.4191, 0.0, 0.0, 120.3476,
    9.2814, 0.0, 0.0, 109.0949, 0.0, 0.0, 0.0, 9.2814, 0.0, 0.0, 109.0949, 0.0,
    0.0, 0.0, 9.2814, 0.0, 0.0, 109.0949, 8.2089, 0.0, 0.0, 98.2173, 0.0, 0.0,
    0.0, 8.2089, 0.0, 0.0, 98.2173, 0.0, 0.0, 0.0, 8.2089, 0.0, 0.0, 98.2173,
    7.2016, 0.0, 0.0, 87.7475, 0.0, 0.0, 0.0, 7.2016, 0.0, 0.0, 87.7475, 0.0,
    0.0, 0.0, 7.2016, 0.0, 0.0, 87.7475, 6.2596, 0.0, 0.0, 77.7181, 0.0, 0.0,
    0.0, 6.2596, 0.0, 0.0, 77.7181, 0.0, 0.0, 0.0, 6.2596, 0.0, 0.0, 77.7181,
    5.3829, 0.0, 0.0, 68.1616, 0.0, 0.0, 0.0, 5.3829, 0.0, 0.0, 68.1616, 0.0,
    0.0, 0.0, 5.3829, 0.0, 0.0, 68.1616, 4.5713, 0.0, 0.0, 59.1108, 0.0, 0.0,
    0.0, 4.5713, 0.0, 0.0, 59.1108, 0.0, 0.0, 0.0, 4.5713, 0.0, 0.0, 59.1108,
    3.8251, 0.0, 0.0, 50.5981, 0.0, 0.0, 0.0, 3.8251, 0.0, 0.0, 50.5981, 0.0,
    0.0, 0.0, 3.8251, 0.0, 0.0, 50.5981, 3.144, 0.0, 0.0, 42.6564, 0.0, 0.0, 0.0,
    3.144, 0.0, 0.0, 42.6564, 0.0, 0.0, 0.0, 3.144, 0.0, 0.0, 42.6564, 2.5282,
    0.0, 0.0, 35.3181, 0.0, 0.0, 0.0, 2.5282, 0.0, 0.0, 35.3181, 0.0, 0.0, 0.0,
    2.5282, 0.0, 0.0, 35.3181 };

  static const signed char A[7200] = { -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
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

  static const signed char a[7200] = { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
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

  static const signed char c_a[7200] = { -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
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

  double IGA[7200];
  double A_data[3600];
  double y_tmp[3600];
  double b_del_lam[240];
  double result[240];
  double varargin_1_data[240];
  double del_lam[120];
  double esig[120];
  double ftol[120];
  double igr2[120];
  double ilam[120];
  double lam[120];
  double mesil[120];
  double nlamt[120];
  double H[60];
  double del_z[60];
  double r1[60];
  double x[60];
  double d;
  double mu;
  double mu_old;
  double ssq;
  int b_i;
  int exitflag;
  int i;
  int idx;
  int iy;
  int unusedU0;
  unsigned char ii_data[240];
  bool b_x[240];

  //   QP Solver
  //  Call QP Solver
  //  Solve quadratic programming problem using Wright's (1997) Method
  //  Minimise J(x) = 1/2x'Hx + f'x
  //  Subject to: Ax <= b
  //  Reference: S. J. Wright, "Applying New Optimization Algorithms to Model
  //  Predictive Control," in Chemical Process Control-V, CACHE, AIChE
  //  Symposium, 1997, pp. 147-155.
  // Number of decision variables
  //  p = 0;
  // Test for Cold Start
  // Warm Start
  for (i = 0; i < 60; i++) {
    d = 0.0;
    for (b_i = 0; b_i < 6; b_i++) {
      d += (2.0 * x0[b_i]) * b[b_i + (6 * i)];
    }

    H[i] = d;
    x[i] = X_QP[i];
  }

  // to tune
  // to tune
  // Default Values
  for (i = 0; i < 120; i++) {
    lam[i] = 100.0;
    ftol[i] = 100.0;
    esig[i] = 0.001;
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
    int i1;
    int ii_size_idx_0;
    int info;
    int j;
    bool exitg1;
    bool y;

    // Create common matrices
    for (i = 0; i < 120; i++) {
      d = lam[i];
      ssq = 1.0 / d;
      ilam[i] = ssq;
      nlamt[i] = (-d) / ftol[i];
      mesil[i] = (mu * esig[i]) * ssq;
    }

    // RHS
    for (i = 0; i < 60; i++) {
      for (b_i = 0; b_i < 120; b_i++) {
        idx = b_i + (120 * i);
        IGA[idx] = nlamt[b_i] * (static_cast<double>(A[idx]));
      }

      d = 0.0;
      for (b_i = 0; b_i < 60; b_i++) {
        d += b_a[i + (60 * b_i)] * x[b_i];
      }

      ssq = 0.0;
      for (b_i = 0; b_i < 120; b_i++) {
        ssq += (static_cast<double>(c_a[i + (60 * b_i)])) * lam[b_i];
      }

      r1[i] = (d - ssq) - H[i];
    }

    for (b_i = 0; b_i < 120; b_i++) {
      d = 0.0;
      for (i1 = 0; i1 < 60; i1++) {
        d += (static_cast<double>(a[b_i + (120 * i1)])) * x[i1];
      }

      igr2[b_i] = nlamt[b_i] * ((d + 0.15) - mesil[b_i]);
    }

    // Solve
    for (b_i = 0; b_i < 60; b_i++) {
      for (i1 = 0; i1 < 60; i1++) {
        d = 0.0;
        for (i = 0; i < 120; i++) {
          d += (static_cast<double>(c_a[b_i + (60 * i)])) * IGA[i + (120 * i1)];
        }

        y_tmp[b_i + (60 * i1)] = d;
      }
    }

    for (b_i = 0; b_i < 3600; b_i++) {
      A_data[b_i] = b_H[b_i] - y_tmp[b_i];
    }

    info = -1;
    j = 0;
    exitg1 = false;
    while ((!exitg1) && (j < 60)) {
      int idxA1j;
      int idxAjj;
      int ix;
      idxA1j = j * 60;
      idxAjj = idxA1j + j;
      ssq = 0.0;
      if (j >= 1) {
        ix = idxA1j;
        iy = idxA1j;
        for (i = 0; i < j; i++) {
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
          int idxAjjp1;
          i = idxA1j + 61;
          idxAjjp1 = idxAjj + 61;
          if (j != 0) {
            iy = idxAjj + 60;
            b_i = (idxA1j + (60 * (58 - j))) + 61;
            for (idx = i; idx <= b_i; idx += 60) {
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
          b_i = (idxAjj + (60 * (58 - j))) + 61;
          for (i = idxAjjp1; i <= b_i; i += 60) {
            A_data[i - 1] *= ssq;
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
      for (b_i = 0; b_i < 60; b_i++) {
        d = 0.0;
        for (i1 = 0; i1 < 120; i1++) {
          d += (static_cast<double>(c_a[b_i + (60 * i1)])) * igr2[i1];
        }

        del_z[b_i] = r1[b_i] - d;
      }

      for (b_i = 0; b_i < 3600; b_i++) {
        y_tmp[b_i] = b_H[b_i] - y_tmp[b_i];
      }

      mldivide(y_tmp, del_z);

      // old method (LU?)
      //   del_z = linsolve (R, linsolve (R, (r1-At*igr2), opUT), opU); %exploit matrix properties for solving 
    } else {
      // Not Positive Definite (problem? eg infeasible)
      for (b_i = 0; b_i < 60; b_i++) {
        d = 0.0;
        for (i1 = 0; i1 < 120; i1++) {
          d += (static_cast<double>(c_a[b_i + (60 * i1)])) * igr2[i1];
        }

        del_z[b_i] = r1[b_i] - d;
      }

      for (b_i = 0; b_i < 3600; b_i++) {
        y_tmp[b_i] = b_H[b_i] - y_tmp[b_i];
      }

      mldivide(y_tmp, del_z);

      // old method (LU?)
    }

    // Decide on suitable alpha (from Wright's paper)
    // Try Max Increment (alpha = 1)
    // Check lam and ftol > 0
    for (i = 0; i < 120; i++) {
      d = 0.0;
      for (b_i = 0; b_i < 60; b_i++) {
        d += IGA[i + (120 * b_i)] * del_z[b_i];
      }

      d = igr2[i] - d;
      del_lam[i] = d;
      ssq = ftol[i];
      mu_old = ((-ssq) + mesil[i]) - ((ilam[i] * ssq) * d);
      mesil[i] = mu_old;
      d += lam[i];
      nlamt[i] = d;
      ssq += mu_old;
      ilam[i] = ssq;
      b_x[i] = (d < 2.2204460492503131E-16);
      b_x[i + 120] = (ssq < 2.2204460492503131E-16);
    }

    y = false;
    i = 0;
    exitg1 = false;
    while ((!exitg1) && (i < 240)) {
      if (!b_x[i]) {
        i++;
      } else {
        y = true;
        exitg1 = true;
      }
    }

    if (!y) {
      // KKT met
      std::memcpy(&lam[0], &nlamt[0], 120U * (sizeof(double)));
      std::memcpy(&ftol[0], &ilam[0], 120U * (sizeof(double)));
      for (b_i = 0; b_i < 60; b_i++) {
        x[b_i] += del_z[b_i];
      }
    } else {
      // KKT failed - solve by finding minimum ratio
      for (b_i = 0; b_i < 120; b_i++) {
        result[b_i] = nlamt[b_i];
        result[b_i + 120] = ilam[b_i];
      }

      idx = 0;
      i = 0;
      exitg1 = false;
      while ((!exitg1) && (i < 240)) {
        if (result[i] < 2.2204460492503131E-16) {
          idx++;
          ii_data[idx - 1] = static_cast<unsigned char>(static_cast<int>(i + 1));
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
        ii_size_idx_0 = 0;
      } else {
        ii_size_idx_0 = idx;
      }

      // detects elements breaking KKT condition
      for (b_i = 0; b_i < 120; b_i++) {
        b_del_lam[b_i] = del_lam[b_i];
        b_del_lam[b_i + 120] = mesil[b_i];
      }

      for (b_i = 0; b_i < ii_size_idx_0; b_i++) {
        i = (static_cast<int>(ii_data[b_i])) - 1;
        varargin_1_data[b_i] = 1.0 - (result[i] / b_del_lam[i]);
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
          i = 2;
          exitg1 = false;
          while ((!exitg1) && (i <= ii_size_idx_0)) {
            if (!rtIsNaN(varargin_1_data[i - 1])) {
              idx = i;
              exitg1 = true;
            } else {
              i++;
            }
          }
        }

        if (idx == 0) {
          ssq = varargin_1_data[0];
        } else {
          ssq = varargin_1_data[idx - 1];
          b_i = idx + 1;
          for (i = b_i; i <= ii_size_idx_0; i++) {
            d = varargin_1_data[i - 1];
            if (ssq > d) {
              ssq = d;
            }
          }
        }
      }

      ssq *= 0.995;

      // solves for min ratio (max value of alpha allowed)
      // Increment
      for (b_i = 0; b_i < 120; b_i++) {
        lam[b_i] += ssq * del_lam[b_i];
        ftol[b_i] += ssq * mesil[b_i];
      }

      for (b_i = 0; b_i < 60; b_i++) {
        x[b_i] += ssq * del_z[b_i];
      }
    }

    // Complimentary Gap
    mu_old = mu;
    ssq = 0.0;
    for (b_i = 0; b_i < 120; b_i++) {
      ssq += ftol[b_i] * lam[b_i];
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

      for (i = 0; i < 120; i++) {
        esig[i] = ssq;
      }
    }

    unusedU0++;
  }

  // Check for failure
  //  if(iter == maxiter)
  //      exitflag = -1;
  //  end
  Fx = x[0];
  Fy = x[1];
  Fz = x[2];
  
  for (i = 0; i < 60; i++) {
    X_QP[i]=x[i];
  }



}

template<typename T>
void CoordinatorBase<T>::MPC_rev_1()
{

  static const double b_H[3600] = { 204.966, 0.0, 0.0, 4.9108, 0.0, 0.0, 4.8555,
    0.0, 0.0, 4.8003, 0.0, 0.0, 4.745, 0.0, 0.0, 4.6898, 0.0, 0.0, 4.6345, 0.0,
    0.0, 4.5793, 0.0, 0.0, 4.524, 0.0, 0.0, 4.4688, 0.0, 0.0, 4.4136, 0.0, 0.0,
    4.3584, 0.0, 0.0, 4.3032, 0.0, 0.0, 4.2481, 0.0, 0.0, 4.1929, 0.0, 0.0,
    4.1378, 0.0, 0.0, 4.0826, 0.0, 0.0, 4.0275, 0.0, 0.0, 3.9724, 0.0, 0.0,
    3.9174, 0.0, 0.0, 0.0, 204.966, 0.0, 0.0, 4.9108, 0.0, 0.0, 4.8555, 0.0, 0.0,
    4.8003, 0.0, 0.0, 4.745, 0.0, 0.0, 4.6898, 0.0, 0.0, 4.6345, 0.0, 0.0,
    4.5793, 0.0, 0.0, 4.524, 0.0, 0.0, 4.4688, 0.0, 0.0, 4.4136, 0.0, 0.0,
    4.3584, 0.0, 0.0, 4.3032, 0.0, 0.0, 4.2481, 0.0, 0.0, 4.1929, 0.0, 0.0,
    4.1378, 0.0, 0.0, 4.0826, 0.0, 0.0, 4.0275, 0.0, 0.0, 3.9724, 0.0, 0.0,
    3.9174, 0.0, 0.0, 0.0, 204.966, 0.0, 0.0, 4.9108, 0.0, 0.0, 4.8555, 0.0, 0.0,
    4.8003, 0.0, 0.0, 4.745, 0.0, 0.0, 4.6898, 0.0, 0.0, 4.6345, 0.0, 0.0,
    4.5793, 0.0, 0.0, 4.524, 0.0, 0.0, 4.4688, 0.0, 0.0, 4.4136, 0.0, 0.0,
    4.3584, 0.0, 0.0, 4.3032, 0.0, 0.0, 4.2481, 0.0, 0.0, 4.1929, 0.0, 0.0,
    4.1378, 0.0, 0.0, 4.0826, 0.0, 0.0, 4.0275, 0.0, 0.0, 3.9724, 0.0, 0.0,
    3.9174, 4.9108, 0.0, 0.0, 204.8907, 0.0, 0.0, 4.8357, 0.0, 0.0, 4.7807, 0.0,
    0.0, 4.7258, 0.0, 0.0, 4.6708, 0.0, 0.0, 4.6158, 0.0, 0.0, 4.5609, 0.0, 0.0,
    4.5059, 0.0, 0.0, 4.451, 0.0, 0.0, 4.3961, 0.0, 0.0, 4.3412, 0.0, 0.0,
    4.2862, 0.0, 0.0, 4.2314, 0.0, 0.0, 4.1765, 0.0, 0.0, 4.1216, 0.0, 0.0,
    4.0668, 0.0, 0.0, 4.0119, 0.0, 0.0, 3.9571, 0.0, 0.0, 3.9023, 0.0, 0.0, 0.0,
    4.9108, 0.0, 0.0, 204.8907, 0.0, 0.0, 4.8357, 0.0, 0.0, 4.7807, 0.0, 0.0,
    4.7258, 0.0, 0.0, 4.6708, 0.0, 0.0, 4.6158, 0.0, 0.0, 4.5609, 0.0, 0.0,
    4.5059, 0.0, 0.0, 4.451, 0.0, 0.0, 4.3961, 0.0, 0.0, 4.3412, 0.0, 0.0,
    4.2862, 0.0, 0.0, 4.2314, 0.0, 0.0, 4.1765, 0.0, 0.0, 4.1216, 0.0, 0.0,
    4.0668, 0.0, 0.0, 4.0119, 0.0, 0.0, 3.9571, 0.0, 0.0, 3.9023, 0.0, 0.0, 0.0,
    4.9108, 0.0, 0.0, 204.8907, 0.0, 0.0, 4.8357, 0.0, 0.0, 4.7807, 0.0, 0.0,
    4.7258, 0.0, 0.0, 4.6708, 0.0, 0.0, 4.6158, 0.0, 0.0, 4.5609, 0.0, 0.0,
    4.5059, 0.0, 0.0, 4.451, 0.0, 0.0, 4.3961, 0.0, 0.0, 4.3412, 0.0, 0.0,
    4.2862, 0.0, 0.0, 4.2314, 0.0, 0.0, 4.1765, 0.0, 0.0, 4.1216, 0.0, 0.0,
    4.0668, 0.0, 0.0, 4.0119, 0.0, 0.0, 3.9571, 0.0, 0.0, 3.9023, 4.8555, 0.0,
    0.0, 4.8357, 0.0, 0.0, 204.8159, 0.0, 0.0, 4.7612, 0.0, 0.0, 4.7065, 0.0,
    0.0, 4.6518, 0.0, 0.0, 4.5972, 0.0, 0.0, 4.5425, 0.0, 0.0, 4.4878, 0.0, 0.0,
    4.4332, 0.0, 0.0, 4.3785, 0.0, 0.0, 4.3239, 0.0, 0.0, 4.2693, 0.0, 0.0,
    4.2146, 0.0, 0.0, 4.16, 0.0, 0.0, 4.1054, 0.0, 0.0, 4.0509, 0.0, 0.0, 3.9963,
    0.0, 0.0, 3.9418, 0.0, 0.0, 3.8872, 0.0, 0.0, 0.0, 4.8555, 0.0, 0.0, 4.8357,
    0.0, 0.0, 204.8159, 0.0, 0.0, 4.7612, 0.0, 0.0, 4.7065, 0.0, 0.0, 4.6518,
    0.0, 0.0, 4.5972, 0.0, 0.0, 4.5425, 0.0, 0.0, 4.4878, 0.0, 0.0, 4.4332, 0.0,
    0.0, 4.3785, 0.0, 0.0, 4.3239, 0.0, 0.0, 4.2693, 0.0, 0.0, 4.2146, 0.0, 0.0,
    4.16, 0.0, 0.0, 4.1054, 0.0, 0.0, 4.0509, 0.0, 0.0, 3.9963, 0.0, 0.0, 3.9418,
    0.0, 0.0, 3.8872, 0.0, 0.0, 0.0, 4.8555, 0.0, 0.0, 4.8357, 0.0, 0.0,
    204.8159, 0.0, 0.0, 4.7612, 0.0, 0.0, 4.7065, 0.0, 0.0, 4.6518, 0.0, 0.0,
    4.5972, 0.0, 0.0, 4.5425, 0.0, 0.0, 4.4878, 0.0, 0.0, 4.4332, 0.0, 0.0,
    4.3785, 0.0, 0.0, 4.3239, 0.0, 0.0, 4.2693, 0.0, 0.0, 4.2146, 0.0, 0.0, 4.16,
    0.0, 0.0, 4.1054, 0.0, 0.0, 4.0509, 0.0, 0.0, 3.9963, 0.0, 0.0, 3.9418, 0.0,
    0.0, 3.8872, 4.8003, 0.0, 0.0, 4.7807, 0.0, 0.0, 4.7612, 0.0, 0.0, 204.7417,
    0.0, 0.0, 4.6873, 0.0, 0.0, 4.6329, 0.0, 0.0, 4.5785, 0.0, 0.0, 4.5241, 0.0,
    0.0, 4.4697, 0.0, 0.0, 4.4154, 0.0, 0.0, 4.361, 0.0, 0.0, 4.3066, 0.0, 0.0,
    4.2523, 0.0, 0.0, 4.1979, 0.0, 0.0, 4.1436, 0.0, 0.0, 4.0893, 0.0, 0.0,
    4.035, 0.0, 0.0, 3.9807, 0.0, 0.0, 3.9264, 0.0, 0.0, 3.8722, 0.0, 0.0, 0.0,
    4.8003, 0.0, 0.0, 4.7807, 0.0, 0.0, 4.7612, 0.0, 0.0, 204.7417, 0.0, 0.0,
    4.6873, 0.0, 0.0, 4.6329, 0.0, 0.0, 4.5785, 0.0, 0.0, 4.5241, 0.0, 0.0,
    4.4697, 0.0, 0.0, 4.4154, 0.0, 0.0, 4.361, 0.0, 0.0, 4.3066, 0.0, 0.0,
    4.2523, 0.0, 0.0, 4.1979, 0.0, 0.0, 4.1436, 0.0, 0.0, 4.0893, 0.0, 0.0,
    4.035, 0.0, 0.0, 3.9807, 0.0, 0.0, 3.9264, 0.0, 0.0, 3.8722, 0.0, 0.0, 0.0,
    4.8003, 0.0, 0.0, 4.7807, 0.0, 0.0, 4.7612, 0.0, 0.0, 204.7417, 0.0, 0.0,
    4.6873, 0.0, 0.0, 4.6329, 0.0, 0.0, 4.5785, 0.0, 0.0, 4.5241, 0.0, 0.0,
    4.4697, 0.0, 0.0, 4.4154, 0.0, 0.0, 4.361, 0.0, 0.0, 4.3066, 0.0, 0.0,
    4.2523, 0.0, 0.0, 4.1979, 0.0, 0.0, 4.1436, 0.0, 0.0, 4.0893, 0.0, 0.0,
    4.035, 0.0, 0.0, 3.9807, 0.0, 0.0, 3.9264, 0.0, 0.0, 3.8722, 4.745, 0.0, 0.0,
    4.7258, 0.0, 0.0, 4.7065, 0.0, 0.0, 4.6873, 0.0, 0.0, 204.668, 0.0, 0.0,
    4.6139, 0.0, 0.0, 4.5598, 0.0, 0.0, 4.5057, 0.0, 0.0, 4.4516, 0.0, 0.0,
    4.3975, 0.0, 0.0, 4.3434, 0.0, 0.0, 4.2894, 0.0, 0.0, 4.2353, 0.0, 0.0,
    4.1812, 0.0, 0.0, 4.1272, 0.0, 0.0, 4.0731, 0.0, 0.0, 4.0191, 0.0, 0.0,
    3.9651, 0.0, 0.0, 3.9111, 0.0, 0.0, 3.8571, 0.0, 0.0, 0.0, 4.745, 0.0, 0.0,
    4.7258, 0.0, 0.0, 4.7065, 0.0, 0.0, 4.6873, 0.0, 0.0, 204.668, 0.0, 0.0,
    4.6139, 0.0, 0.0, 4.5598, 0.0, 0.0, 4.5057, 0.0, 0.0, 4.4516, 0.0, 0.0,
    4.3975, 0.0, 0.0, 4.3434, 0.0, 0.0, 4.2894, 0.0, 0.0, 4.2353, 0.0, 0.0,
    4.1812, 0.0, 0.0, 4.1272, 0.0, 0.0, 4.0731, 0.0, 0.0, 4.0191, 0.0, 0.0,
    3.9651, 0.0, 0.0, 3.9111, 0.0, 0.0, 3.8571, 0.0, 0.0, 0.0, 4.745, 0.0, 0.0,
    4.7258, 0.0, 0.0, 4.7065, 0.0, 0.0, 4.6873, 0.0, 0.0, 204.668, 0.0, 0.0,
    4.6139, 0.0, 0.0, 4.5598, 0.0, 0.0, 4.5057, 0.0, 0.0, 4.4516, 0.0, 0.0,
    4.3975, 0.0, 0.0, 4.3434, 0.0, 0.0, 4.2894, 0.0, 0.0, 4.2353, 0.0, 0.0,
    4.1812, 0.0, 0.0, 4.1272, 0.0, 0.0, 4.0731, 0.0, 0.0, 4.0191, 0.0, 0.0,
    3.9651, 0.0, 0.0, 3.9111, 0.0, 0.0, 3.8571, 4.6898, 0.0, 0.0, 4.6708, 0.0,
    0.0, 4.6518, 0.0, 0.0, 4.6329, 0.0, 0.0, 4.6139, 0.0, 0.0, 204.595, 0.0, 0.0,
    4.5412, 0.0, 0.0, 4.4873, 0.0, 0.0, 4.4335, 0.0, 0.0, 4.3797, 0.0, 0.0,
    4.3259, 0.0, 0.0, 4.2721, 0.0, 0.0, 4.2183, 0.0, 0.0, 4.1645, 0.0, 0.0,
    4.1107, 0.0, 0.0, 4.057, 0.0, 0.0, 4.0032, 0.0, 0.0, 3.9495, 0.0, 0.0,
    3.8957, 0.0, 0.0, 3.842, 0.0, 0.0, 0.0, 4.6898, 0.0, 0.0, 4.6708, 0.0, 0.0,
    4.6518, 0.0, 0.0, 4.6329, 0.0, 0.0, 4.6139, 0.0, 0.0, 204.595, 0.0, 0.0,
    4.5412, 0.0, 0.0, 4.4873, 0.0, 0.0, 4.4335, 0.0, 0.0, 4.3797, 0.0, 0.0,
    4.3259, 0.0, 0.0, 4.2721, 0.0, 0.0, 4.2183, 0.0, 0.0, 4.1645, 0.0, 0.0,
    4.1107, 0.0, 0.0, 4.057, 0.0, 0.0, 4.0032, 0.0, 0.0, 3.9495, 0.0, 0.0,
    3.8957, 0.0, 0.0, 3.842, 0.0, 0.0, 0.0, 4.6898, 0.0, 0.0, 4.6708, 0.0, 0.0,
    4.6518, 0.0, 0.0, 4.6329, 0.0, 0.0, 4.6139, 0.0, 0.0, 204.595, 0.0, 0.0,
    4.5412, 0.0, 0.0, 4.4873, 0.0, 0.0, 4.4335, 0.0, 0.0, 4.3797, 0.0, 0.0,
    4.3259, 0.0, 0.0, 4.2721, 0.0, 0.0, 4.2183, 0.0, 0.0, 4.1645, 0.0, 0.0,
    4.1107, 0.0, 0.0, 4.057, 0.0, 0.0, 4.0032, 0.0, 0.0, 3.9495, 0.0, 0.0,
    3.8957, 0.0, 0.0, 3.842, 4.6345, 0.0, 0.0, 4.6158, 0.0, 0.0, 4.5972, 0.0,
    0.0, 4.5785, 0.0, 0.0, 4.5598, 0.0, 0.0, 4.5412, 0.0, 0.0, 204.5225, 0.0,
    0.0, 4.469, 0.0, 0.0, 4.4154, 0.0, 0.0, 4.3619, 0.0, 0.0, 4.3084, 0.0, 0.0,
    4.2548, 0.0, 0.0, 4.2013, 0.0, 0.0, 4.1478, 0.0, 0.0, 4.0943, 0.0, 0.0,
    4.0408, 0.0, 0.0, 3.9873, 0.0, 0.0, 3.9338, 0.0, 0.0, 3.8804, 0.0, 0.0,
    3.8269, 0.0, 0.0, 0.0, 4.6345, 0.0, 0.0, 4.6158, 0.0, 0.0, 4.5972, 0.0, 0.0,
    4.5785, 0.0, 0.0, 4.5598, 0.0, 0.0, 4.5412, 0.0, 0.0, 204.5225, 0.0, 0.0,
    4.469, 0.0, 0.0, 4.4154, 0.0, 0.0, 4.3619, 0.0, 0.0, 4.3084, 0.0, 0.0,
    4.2548, 0.0, 0.0, 4.2013, 0.0, 0.0, 4.1478, 0.0, 0.0, 4.0943, 0.0, 0.0,
    4.0408, 0.0, 0.0, 3.9873, 0.0, 0.0, 3.9338, 0.0, 0.0, 3.8804, 0.0, 0.0,
    3.8269, 0.0, 0.0, 0.0, 4.6345, 0.0, 0.0, 4.6158, 0.0, 0.0, 4.5972, 0.0, 0.0,
    4.5785, 0.0, 0.0, 4.5598, 0.0, 0.0, 4.5412, 0.0, 0.0, 204.5225, 0.0, 0.0,
    4.469, 0.0, 0.0, 4.4154, 0.0, 0.0, 4.3619, 0.0, 0.0, 4.3084, 0.0, 0.0,
    4.2548, 0.0, 0.0, 4.2013, 0.0, 0.0, 4.1478, 0.0, 0.0, 4.0943, 0.0, 0.0,
    4.0408, 0.0, 0.0, 3.9873, 0.0, 0.0, 3.9338, 0.0, 0.0, 3.8804, 0.0, 0.0,
    3.8269, 4.5793, 0.0, 0.0, 4.5609, 0.0, 0.0, 4.5425, 0.0, 0.0, 4.5241, 0.0,
    0.0, 4.5057, 0.0, 0.0, 4.4873, 0.0, 0.0, 4.469, 0.0, 0.0, 204.4506, 0.0, 0.0,
    4.3973, 0.0, 0.0, 4.3441, 0.0, 0.0, 4.2908, 0.0, 0.0, 4.2376, 0.0, 0.0,
    4.1843, 0.0, 0.0, 4.1311, 0.0, 0.0, 4.0779, 0.0, 0.0, 4.0246, 0.0, 0.0,
    3.9714, 0.0, 0.0, 3.9182, 0.0, 0.0, 3.865, 0.0, 0.0, 3.8119, 0.0, 0.0, 0.0,
    4.5793, 0.0, 0.0, 4.5609, 0.0, 0.0, 4.5425, 0.0, 0.0, 4.5241, 0.0, 0.0,
    4.5057, 0.0, 0.0, 4.4873, 0.0, 0.0, 4.469, 0.0, 0.0, 204.4506, 0.0, 0.0,
    4.3973, 0.0, 0.0, 4.3441, 0.0, 0.0, 4.2908, 0.0, 0.0, 4.2376, 0.0, 0.0,
    4.1843, 0.0, 0.0, 4.1311, 0.0, 0.0, 4.0779, 0.0, 0.0, 4.0246, 0.0, 0.0,
    3.9714, 0.0, 0.0, 3.9182, 0.0, 0.0, 3.865, 0.0, 0.0, 3.8119, 0.0, 0.0, 0.0,
    4.5793, 0.0, 0.0, 4.5609, 0.0, 0.0, 4.5425, 0.0, 0.0, 4.5241, 0.0, 0.0,
    4.5057, 0.0, 0.0, 4.4873, 0.0, 0.0, 4.469, 0.0, 0.0, 204.4506, 0.0, 0.0,
    4.3973, 0.0, 0.0, 4.3441, 0.0, 0.0, 4.2908, 0.0, 0.0, 4.2376, 0.0, 0.0,
    4.1843, 0.0, 0.0, 4.1311, 0.0, 0.0, 4.0779, 0.0, 0.0, 4.0246, 0.0, 0.0,
    3.9714, 0.0, 0.0, 3.9182, 0.0, 0.0, 3.865, 0.0, 0.0, 3.8119, 4.524, 0.0, 0.0,
    4.5059, 0.0, 0.0, 4.4878, 0.0, 0.0, 4.4697, 0.0, 0.0, 4.4516, 0.0, 0.0,
    4.4335, 0.0, 0.0, 4.4154, 0.0, 0.0, 4.3973, 0.0, 0.0, 204.3792, 0.0, 0.0,
    4.3262, 0.0, 0.0, 4.2733, 0.0, 0.0, 4.2203, 0.0, 0.0, 4.1673, 0.0, 0.0,
    4.1144, 0.0, 0.0, 4.0614, 0.0, 0.0, 4.0085, 0.0, 0.0, 3.9555, 0.0, 0.0,
    3.9026, 0.0, 0.0, 3.8497, 0.0, 0.0, 3.7968, 0.0, 0.0, 0.0, 4.524, 0.0, 0.0,
    4.5059, 0.0, 0.0, 4.4878, 0.0, 0.0, 4.4697, 0.0, 0.0, 4.4516, 0.0, 0.0,
    4.4335, 0.0, 0.0, 4.4154, 0.0, 0.0, 4.3973, 0.0, 0.0, 204.3792, 0.0, 0.0,
    4.3262, 0.0, 0.0, 4.2733, 0.0, 0.0, 4.2203, 0.0, 0.0, 4.1673, 0.0, 0.0,
    4.1144, 0.0, 0.0, 4.0614, 0.0, 0.0, 4.0085, 0.0, 0.0, 3.9555, 0.0, 0.0,
    3.9026, 0.0, 0.0, 3.8497, 0.0, 0.0, 3.7968, 0.0, 0.0, 0.0, 4.524, 0.0, 0.0,
    4.5059, 0.0, 0.0, 4.4878, 0.0, 0.0, 4.4697, 0.0, 0.0, 4.4516, 0.0, 0.0,
    4.4335, 0.0, 0.0, 4.4154, 0.0, 0.0, 4.3973, 0.0, 0.0, 204.3792, 0.0, 0.0,
    4.3262, 0.0, 0.0, 4.2733, 0.0, 0.0, 4.2203, 0.0, 0.0, 4.1673, 0.0, 0.0,
    4.1144, 0.0, 0.0, 4.0614, 0.0, 0.0, 4.0085, 0.0, 0.0, 3.9555, 0.0, 0.0,
    3.9026, 0.0, 0.0, 3.8497, 0.0, 0.0, 3.7968, 4.4688, 0.0, 0.0, 4.451, 0.0,
    0.0, 4.4332, 0.0, 0.0, 4.4154, 0.0, 0.0, 4.3975, 0.0, 0.0, 4.3797, 0.0, 0.0,
    4.3619, 0.0, 0.0, 4.3441, 0.0, 0.0, 4.3262, 0.0, 0.0, 204.3084, 0.0, 0.0,
    4.2557, 0.0, 0.0, 4.203, 0.0, 0.0, 4.1504, 0.0, 0.0, 4.0977, 0.0, 0.0, 4.045,
    0.0, 0.0, 3.9923, 0.0, 0.0, 3.9397, 0.0, 0.0, 3.887, 0.0, 0.0, 3.8344, 0.0,
    0.0, 3.7817, 0.0, 0.0, 0.0, 4.4688, 0.0, 0.0, 4.451, 0.0, 0.0, 4.4332, 0.0,
    0.0, 4.4154, 0.0, 0.0, 4.3975, 0.0, 0.0, 4.3797, 0.0, 0.0, 4.3619, 0.0, 0.0,
    4.3441, 0.0, 0.0, 4.3262, 0.0, 0.0, 204.3084, 0.0, 0.0, 4.2557, 0.0, 0.0,
    4.203, 0.0, 0.0, 4.1504, 0.0, 0.0, 4.0977, 0.0, 0.0, 4.045, 0.0, 0.0, 3.9923,
    0.0, 0.0, 3.9397, 0.0, 0.0, 3.887, 0.0, 0.0, 3.8344, 0.0, 0.0, 3.7817, 0.0,
    0.0, 0.0, 4.4688, 0.0, 0.0, 4.451, 0.0, 0.0, 4.4332, 0.0, 0.0, 4.4154, 0.0,
    0.0, 4.3975, 0.0, 0.0, 4.3797, 0.0, 0.0, 4.3619, 0.0, 0.0, 4.3441, 0.0, 0.0,
    4.3262, 0.0, 0.0, 204.3084, 0.0, 0.0, 4.2557, 0.0, 0.0, 4.203, 0.0, 0.0,
    4.1504, 0.0, 0.0, 4.0977, 0.0, 0.0, 4.045, 0.0, 0.0, 3.9923, 0.0, 0.0,
    3.9397, 0.0, 0.0, 3.887, 0.0, 0.0, 3.8344, 0.0, 0.0, 3.7817, 4.4136, 0.0,
    0.0, 4.3961, 0.0, 0.0, 4.3785, 0.0, 0.0, 4.361, 0.0, 0.0, 4.3434, 0.0, 0.0,
    4.3259, 0.0, 0.0, 4.3084, 0.0, 0.0, 4.2908, 0.0, 0.0, 4.2733, 0.0, 0.0,
    4.2557, 0.0, 0.0, 204.2382, 0.0, 0.0, 4.1858, 0.0, 0.0, 4.1334, 0.0, 0.0,
    4.081, 0.0, 0.0, 4.0286, 0.0, 0.0, 3.9762, 0.0, 0.0, 3.9238, 0.0, 0.0,
    3.8714, 0.0, 0.0, 3.819, 0.0, 0.0, 3.7666, 0.0, 0.0, 0.0, 4.4136, 0.0, 0.0,
    4.3961, 0.0, 0.0, 4.3785, 0.0, 0.0, 4.361, 0.0, 0.0, 4.3434, 0.0, 0.0,
    4.3259, 0.0, 0.0, 4.3084, 0.0, 0.0, 4.2908, 0.0, 0.0, 4.2733, 0.0, 0.0,
    4.2557, 0.0, 0.0, 204.2382, 0.0, 0.0, 4.1858, 0.0, 0.0, 4.1334, 0.0, 0.0,
    4.081, 0.0, 0.0, 4.0286, 0.0, 0.0, 3.9762, 0.0, 0.0, 3.9238, 0.0, 0.0,
    3.8714, 0.0, 0.0, 3.819, 0.0, 0.0, 3.7666, 0.0, 0.0, 0.0, 4.4136, 0.0, 0.0,
    4.3961, 0.0, 0.0, 4.3785, 0.0, 0.0, 4.361, 0.0, 0.0, 4.3434, 0.0, 0.0,
    4.3259, 0.0, 0.0, 4.3084, 0.0, 0.0, 4.2908, 0.0, 0.0, 4.2733, 0.0, 0.0,
    4.2557, 0.0, 0.0, 204.2382, 0.0, 0.0, 4.1858, 0.0, 0.0, 4.1334, 0.0, 0.0,
    4.081, 0.0, 0.0, 4.0286, 0.0, 0.0, 3.9762, 0.0, 0.0, 3.9238, 0.0, 0.0,
    3.8714, 0.0, 0.0, 3.819, 0.0, 0.0, 3.7666, 4.3584, 0.0, 0.0, 4.3412, 0.0,
    0.0, 4.3239, 0.0, 0.0, 4.3066, 0.0, 0.0, 4.2894, 0.0, 0.0, 4.2721, 0.0, 0.0,
    4.2548, 0.0, 0.0, 4.2376, 0.0, 0.0, 4.2203, 0.0, 0.0, 4.203, 0.0, 0.0,
    4.1858, 0.0, 0.0, 204.1685, 0.0, 0.0, 4.1164, 0.0, 0.0, 4.0643, 0.0, 0.0,
    4.0121, 0.0, 0.0, 3.96, 0.0, 0.0, 3.9079, 0.0, 0.0, 3.8558, 0.0, 0.0, 3.8037,
    0.0, 0.0, 3.7516, 0.0, 0.0, 0.0, 4.3584, 0.0, 0.0, 4.3412, 0.0, 0.0, 4.3239,
    0.0, 0.0, 4.3066, 0.0, 0.0, 4.2894, 0.0, 0.0, 4.2721, 0.0, 0.0, 4.2548, 0.0,
    0.0, 4.2376, 0.0, 0.0, 4.2203, 0.0, 0.0, 4.203, 0.0, 0.0, 4.1858, 0.0, 0.0,
    204.1685, 0.0, 0.0, 4.1164, 0.0, 0.0, 4.0643, 0.0, 0.0, 4.0121, 0.0, 0.0,
    3.96, 0.0, 0.0, 3.9079, 0.0, 0.0, 3.8558, 0.0, 0.0, 3.8037, 0.0, 0.0, 3.7516,
    0.0, 0.0, 0.0, 4.3584, 0.0, 0.0, 4.3412, 0.0, 0.0, 4.3239, 0.0, 0.0, 4.3066,
    0.0, 0.0, 4.2894, 0.0, 0.0, 4.2721, 0.0, 0.0, 4.2548, 0.0, 0.0, 4.2376, 0.0,
    0.0, 4.2203, 0.0, 0.0, 4.203, 0.0, 0.0, 4.1858, 0.0, 0.0, 204.1685, 0.0, 0.0,
    4.1164, 0.0, 0.0, 4.0643, 0.0, 0.0, 4.0121, 0.0, 0.0, 3.96, 0.0, 0.0, 3.9079,
    0.0, 0.0, 3.8558, 0.0, 0.0, 3.8037, 0.0, 0.0, 3.7516, 4.3032, 0.0, 0.0,
    4.2862, 0.0, 0.0, 4.2693, 0.0, 0.0, 4.2523, 0.0, 0.0, 4.2353, 0.0, 0.0,
    4.2183, 0.0, 0.0, 4.2013, 0.0, 0.0, 4.1843, 0.0, 0.0, 4.1673, 0.0, 0.0,
    4.1504, 0.0, 0.0, 4.1334, 0.0, 0.0, 4.1164, 0.0, 0.0, 204.0994, 0.0, 0.0,
    4.0476, 0.0, 0.0, 3.9957, 0.0, 0.0, 3.9438, 0.0, 0.0, 3.892, 0.0, 0.0,
    3.8402, 0.0, 0.0, 3.7883, 0.0, 0.0, 3.7365, 0.0, 0.0, 0.0, 4.3032, 0.0, 0.0,
    4.2862, 0.0, 0.0, 4.2693, 0.0, 0.0, 4.2523, 0.0, 0.0, 4.2353, 0.0, 0.0,
    4.2183, 0.0, 0.0, 4.2013, 0.0, 0.0, 4.1843, 0.0, 0.0, 4.1673, 0.0, 0.0,
    4.1504, 0.0, 0.0, 4.1334, 0.0, 0.0, 4.1164, 0.0, 0.0, 204.0994, 0.0, 0.0,
    4.0476, 0.0, 0.0, 3.9957, 0.0, 0.0, 3.9438, 0.0, 0.0, 3.892, 0.0, 0.0,
    3.8402, 0.0, 0.0, 3.7883, 0.0, 0.0, 3.7365, 0.0, 0.0, 0.0, 4.3032, 0.0, 0.0,
    4.2862, 0.0, 0.0, 4.2693, 0.0, 0.0, 4.2523, 0.0, 0.0, 4.2353, 0.0, 0.0,
    4.2183, 0.0, 0.0, 4.2013, 0.0, 0.0, 4.1843, 0.0, 0.0, 4.1673, 0.0, 0.0,
    4.1504, 0.0, 0.0, 4.1334, 0.0, 0.0, 4.1164, 0.0, 0.0, 204.0994, 0.0, 0.0,
    4.0476, 0.0, 0.0, 3.9957, 0.0, 0.0, 3.9438, 0.0, 0.0, 3.892, 0.0, 0.0,
    3.8402, 0.0, 0.0, 3.7883, 0.0, 0.0, 3.7365, 4.2481, 0.0, 0.0, 4.2314, 0.0,
    0.0, 4.2146, 0.0, 0.0, 4.1979, 0.0, 0.0, 4.1812, 0.0, 0.0, 4.1645, 0.0, 0.0,
    4.1478, 0.0, 0.0, 4.1311, 0.0, 0.0, 4.1144, 0.0, 0.0, 4.0977, 0.0, 0.0,
    4.081, 0.0, 0.0, 4.0643, 0.0, 0.0, 4.0476, 0.0, 0.0, 204.0308, 0.0, 0.0,
    3.9793, 0.0, 0.0, 3.9277, 0.0, 0.0, 3.8761, 0.0, 0.0, 3.8245, 0.0, 0.0,
    3.773, 0.0, 0.0, 3.7214, 0.0, 0.0, 0.0, 4.2481, 0.0, 0.0, 4.2314, 0.0, 0.0,
    4.2146, 0.0, 0.0, 4.1979, 0.0, 0.0, 4.1812, 0.0, 0.0, 4.1645, 0.0, 0.0,
    4.1478, 0.0, 0.0, 4.1311, 0.0, 0.0, 4.1144, 0.0, 0.0, 4.0977, 0.0, 0.0,
    4.081, 0.0, 0.0, 4.0643, 0.0, 0.0, 4.0476, 0.0, 0.0, 204.0308, 0.0, 0.0,
    3.9793, 0.0, 0.0, 3.9277, 0.0, 0.0, 3.8761, 0.0, 0.0, 3.8245, 0.0, 0.0,
    3.773, 0.0, 0.0, 3.7214, 0.0, 0.0, 0.0, 4.2481, 0.0, 0.0, 4.2314, 0.0, 0.0,
    4.2146, 0.0, 0.0, 4.1979, 0.0, 0.0, 4.1812, 0.0, 0.0, 4.1645, 0.0, 0.0,
    4.1478, 0.0, 0.0, 4.1311, 0.0, 0.0, 4.1144, 0.0, 0.0, 4.0977, 0.0, 0.0,
    4.081, 0.0, 0.0, 4.0643, 0.0, 0.0, 4.0476, 0.0, 0.0, 204.0308, 0.0, 0.0,
    3.9793, 0.0, 0.0, 3.9277, 0.0, 0.0, 3.8761, 0.0, 0.0, 3.8245, 0.0, 0.0,
    3.773, 0.0, 0.0, 3.7214, 4.1929, 0.0, 0.0, 4.1765, 0.0, 0.0, 4.16, 0.0, 0.0,
    4.1436, 0.0, 0.0, 4.1272, 0.0, 0.0, 4.1107, 0.0, 0.0, 4.0943, 0.0, 0.0,
    4.0779, 0.0, 0.0, 4.0614, 0.0, 0.0, 4.045, 0.0, 0.0, 4.0286, 0.0, 0.0,
    4.0121, 0.0, 0.0, 3.9957, 0.0, 0.0, 3.9793, 0.0, 0.0, 203.9628, 0.0, 0.0,
    3.9115, 0.0, 0.0, 3.8602, 0.0, 0.0, 3.8089, 0.0, 0.0, 3.7576, 0.0, 0.0,
    3.7063, 0.0, 0.0, 0.0, 4.1929, 0.0, 0.0, 4.1765, 0.0, 0.0, 4.16, 0.0, 0.0,
    4.1436, 0.0, 0.0, 4.1272, 0.0, 0.0, 4.1107, 0.0, 0.0, 4.0943, 0.0, 0.0,
    4.0779, 0.0, 0.0, 4.0614, 0.0, 0.0, 4.045, 0.0, 0.0, 4.0286, 0.0, 0.0,
    4.0121, 0.0, 0.0, 3.9957, 0.0, 0.0, 3.9793, 0.0, 0.0, 203.9628, 0.0, 0.0,
    3.9115, 0.0, 0.0, 3.8602, 0.0, 0.0, 3.8089, 0.0, 0.0, 3.7576, 0.0, 0.0,
    3.7063, 0.0, 0.0, 0.0, 4.1929, 0.0, 0.0, 4.1765, 0.0, 0.0, 4.16, 0.0, 0.0,
    4.1436, 0.0, 0.0, 4.1272, 0.0, 0.0, 4.1107, 0.0, 0.0, 4.0943, 0.0, 0.0,
    4.0779, 0.0, 0.0, 4.0614, 0.0, 0.0, 4.045, 0.0, 0.0, 4.0286, 0.0, 0.0,
    4.0121, 0.0, 0.0, 3.9957, 0.0, 0.0, 3.9793, 0.0, 0.0, 203.9628, 0.0, 0.0,
    3.9115, 0.0, 0.0, 3.8602, 0.0, 0.0, 3.8089, 0.0, 0.0, 3.7576, 0.0, 0.0,
    3.7063, 4.1378, 0.0, 0.0, 4.1216, 0.0, 0.0, 4.1054, 0.0, 0.0, 4.0893, 0.0,
    0.0, 4.0731, 0.0, 0.0, 4.057, 0.0, 0.0, 4.0408, 0.0, 0.0, 4.0246, 0.0, 0.0,
    4.0085, 0.0, 0.0, 3.9923, 0.0, 0.0, 3.9762, 0.0, 0.0, 3.96, 0.0, 0.0, 3.9438,
    0.0, 0.0, 3.9277, 0.0, 0.0, 3.9115, 0.0, 0.0, 203.8954, 0.0, 0.0, 3.8443,
    0.0, 0.0, 3.7933, 0.0, 0.0, 3.7423, 0.0, 0.0, 3.6913, 0.0, 0.0, 0.0, 4.1378,
    0.0, 0.0, 4.1216, 0.0, 0.0, 4.1054, 0.0, 0.0, 4.0893, 0.0, 0.0, 4.0731, 0.0,
    0.0, 4.057, 0.0, 0.0, 4.0408, 0.0, 0.0, 4.0246, 0.0, 0.0, 4.0085, 0.0, 0.0,
    3.9923, 0.0, 0.0, 3.9762, 0.0, 0.0, 3.96, 0.0, 0.0, 3.9438, 0.0, 0.0, 3.9277,
    0.0, 0.0, 3.9115, 0.0, 0.0, 203.8954, 0.0, 0.0, 3.8443, 0.0, 0.0, 3.7933,
    0.0, 0.0, 3.7423, 0.0, 0.0, 3.6913, 0.0, 0.0, 0.0, 4.1378, 0.0, 0.0, 4.1216,
    0.0, 0.0, 4.1054, 0.0, 0.0, 4.0893, 0.0, 0.0, 4.0731, 0.0, 0.0, 4.057, 0.0,
    0.0, 4.0408, 0.0, 0.0, 4.0246, 0.0, 0.0, 4.0085, 0.0, 0.0, 3.9923, 0.0, 0.0,
    3.9762, 0.0, 0.0, 3.96, 0.0, 0.0, 3.9438, 0.0, 0.0, 3.9277, 0.0, 0.0, 3.9115,
    0.0, 0.0, 203.8954, 0.0, 0.0, 3.8443, 0.0, 0.0, 3.7933, 0.0, 0.0, 3.7423,
    0.0, 0.0, 3.6913, 4.0826, 0.0, 0.0, 4.0668, 0.0, 0.0, 4.0509, 0.0, 0.0,
    4.035, 0.0, 0.0, 4.0191, 0.0, 0.0, 4.0032, 0.0, 0.0, 3.9873, 0.0, 0.0,
    3.9714, 0.0, 0.0, 3.9555, 0.0, 0.0, 3.9397, 0.0, 0.0, 3.9238, 0.0, 0.0,
    3.9079, 0.0, 0.0, 3.892, 0.0, 0.0, 3.8761, 0.0, 0.0, 3.8602, 0.0, 0.0,
    3.8443, 0.0, 0.0, 203.8285, 0.0, 0.0, 3.7777, 0.0, 0.0, 3.727, 0.0, 0.0,
    3.6762, 0.0, 0.0, 0.0, 4.0826, 0.0, 0.0, 4.0668, 0.0, 0.0, 4.0509, 0.0, 0.0,
    4.035, 0.0, 0.0, 4.0191, 0.0, 0.0, 4.0032, 0.0, 0.0, 3.9873, 0.0, 0.0,
    3.9714, 0.0, 0.0, 3.9555, 0.0, 0.0, 3.9397, 0.0, 0.0, 3.9238, 0.0, 0.0,
    3.9079, 0.0, 0.0, 3.892, 0.0, 0.0, 3.8761, 0.0, 0.0, 3.8602, 0.0, 0.0,
    3.8443, 0.0, 0.0, 203.8285, 0.0, 0.0, 3.7777, 0.0, 0.0, 3.727, 0.0, 0.0,
    3.6762, 0.0, 0.0, 0.0, 4.0826, 0.0, 0.0, 4.0668, 0.0, 0.0, 4.0509, 0.0, 0.0,
    4.035, 0.0, 0.0, 4.0191, 0.0, 0.0, 4.0032, 0.0, 0.0, 3.9873, 0.0, 0.0,
    3.9714, 0.0, 0.0, 3.9555, 0.0, 0.0, 3.9397, 0.0, 0.0, 3.9238, 0.0, 0.0,
    3.9079, 0.0, 0.0, 3.892, 0.0, 0.0, 3.8761, 0.0, 0.0, 3.8602, 0.0, 0.0,
    3.8443, 0.0, 0.0, 203.8285, 0.0, 0.0, 3.7777, 0.0, 0.0, 3.727, 0.0, 0.0,
    3.6762, 4.0275, 0.0, 0.0, 4.0119, 0.0, 0.0, 3.9963, 0.0, 0.0, 3.9807, 0.0,
    0.0, 3.9651, 0.0, 0.0, 3.9495, 0.0, 0.0, 3.9338, 0.0, 0.0, 3.9182, 0.0, 0.0,
    3.9026, 0.0, 0.0, 3.887, 0.0, 0.0, 3.8714, 0.0, 0.0, 3.8558, 0.0, 0.0,
    3.8402, 0.0, 0.0, 3.8245, 0.0, 0.0, 3.8089, 0.0, 0.0, 3.7933, 0.0, 0.0,
    3.7777, 0.0, 0.0, 203.7621, 0.0, 0.0, 3.7116, 0.0, 0.0, 3.6611, 0.0, 0.0,
    0.0, 4.0275, 0.0, 0.0, 4.0119, 0.0, 0.0, 3.9963, 0.0, 0.0, 3.9807, 0.0, 0.0,
    3.9651, 0.0, 0.0, 3.9495, 0.0, 0.0, 3.9338, 0.0, 0.0, 3.9182, 0.0, 0.0,
    3.9026, 0.0, 0.0, 3.887, 0.0, 0.0, 3.8714, 0.0, 0.0, 3.8558, 0.0, 0.0,
    3.8402, 0.0, 0.0, 3.8245, 0.0, 0.0, 3.8089, 0.0, 0.0, 3.7933, 0.0, 0.0,
    3.7777, 0.0, 0.0, 203.7621, 0.0, 0.0, 3.7116, 0.0, 0.0, 3.6611, 0.0, 0.0,
    0.0, 4.0275, 0.0, 0.0, 4.0119, 0.0, 0.0, 3.9963, 0.0, 0.0, 3.9807, 0.0, 0.0,
    3.9651, 0.0, 0.0, 3.9495, 0.0, 0.0, 3.9338, 0.0, 0.0, 3.9182, 0.0, 0.0,
    3.9026, 0.0, 0.0, 3.887, 0.0, 0.0, 3.8714, 0.0, 0.0, 3.8558, 0.0, 0.0,
    3.8402, 0.0, 0.0, 3.8245, 0.0, 0.0, 3.8089, 0.0, 0.0, 3.7933, 0.0, 0.0,
    3.7777, 0.0, 0.0, 203.7621, 0.0, 0.0, 3.7116, 0.0, 0.0, 3.6611, 3.9724, 0.0,
    0.0, 3.9571, 0.0, 0.0, 3.9418, 0.0, 0.0, 3.9264, 0.0, 0.0, 3.9111, 0.0, 0.0,
    3.8957, 0.0, 0.0, 3.8804, 0.0, 0.0, 3.865, 0.0, 0.0, 3.8497, 0.0, 0.0,
    3.8344, 0.0, 0.0, 3.819, 0.0, 0.0, 3.8037, 0.0, 0.0, 3.7883, 0.0, 0.0, 3.773,
    0.0, 0.0, 3.7576, 0.0, 0.0, 3.7423, 0.0, 0.0, 3.727, 0.0, 0.0, 3.7116, 0.0,
    0.0, 203.6963, 0.0, 0.0, 3.6461, 0.0, 0.0, 0.0, 3.9724, 0.0, 0.0, 3.9571,
    0.0, 0.0, 3.9418, 0.0, 0.0, 3.9264, 0.0, 0.0, 3.9111, 0.0, 0.0, 3.8957, 0.0,
    0.0, 3.8804, 0.0, 0.0, 3.865, 0.0, 0.0, 3.8497, 0.0, 0.0, 3.8344, 0.0, 0.0,
    3.819, 0.0, 0.0, 3.8037, 0.0, 0.0, 3.7883, 0.0, 0.0, 3.773, 0.0, 0.0, 3.7576,
    0.0, 0.0, 3.7423, 0.0, 0.0, 3.727, 0.0, 0.0, 3.7116, 0.0, 0.0, 203.6963, 0.0,
    0.0, 3.6461, 0.0, 0.0, 0.0, 3.9724, 0.0, 0.0, 3.9571, 0.0, 0.0, 3.9418, 0.0,
    0.0, 3.9264, 0.0, 0.0, 3.9111, 0.0, 0.0, 3.8957, 0.0, 0.0, 3.8804, 0.0, 0.0,
    3.865, 0.0, 0.0, 3.8497, 0.0, 0.0, 3.8344, 0.0, 0.0, 3.819, 0.0, 0.0, 3.8037,
    0.0, 0.0, 3.7883, 0.0, 0.0, 3.773, 0.0, 0.0, 3.7576, 0.0, 0.0, 3.7423, 0.0,
    0.0, 3.727, 0.0, 0.0, 3.7116, 0.0, 0.0, 203.6963, 0.0, 0.0, 3.6461, 3.9174,
    0.0, 0.0, 3.9023, 0.0, 0.0, 3.8872, 0.0, 0.0, 3.8722, 0.0, 0.0, 3.8571, 0.0,
    0.0, 3.842, 0.0, 0.0, 3.8269, 0.0, 0.0, 3.8119, 0.0, 0.0, 3.7968, 0.0, 0.0,
    3.7817, 0.0, 0.0, 3.7666, 0.0, 0.0, 3.7516, 0.0, 0.0, 3.7365, 0.0, 0.0,
    3.7214, 0.0, 0.0, 3.7063, 0.0, 0.0, 3.6913, 0.0, 0.0, 3.6762, 0.0, 0.0,
    3.6611, 0.0, 0.0, 3.6461, 0.0, 0.0, 203.631, 0.0, 0.0, 0.0, 3.9174, 0.0, 0.0,
    3.9023, 0.0, 0.0, 3.8872, 0.0, 0.0, 3.8722, 0.0, 0.0, 3.8571, 0.0, 0.0,
    3.842, 0.0, 0.0, 3.8269, 0.0, 0.0, 3.8119, 0.0, 0.0, 3.7968, 0.0, 0.0,
    3.7817, 0.0, 0.0, 3.7666, 0.0, 0.0, 3.7516, 0.0, 0.0, 3.7365, 0.0, 0.0,
    3.7214, 0.0, 0.0, 3.7063, 0.0, 0.0, 3.6913, 0.0, 0.0, 3.6762, 0.0, 0.0,
    3.6611, 0.0, 0.0, 3.6461, 0.0, 0.0, 203.631, 0.0, 0.0, 0.0, 3.9174, 0.0, 0.0,
    3.9023, 0.0, 0.0, 3.8872, 0.0, 0.0, 3.8722, 0.0, 0.0, 3.8571, 0.0, 0.0,
    3.842, 0.0, 0.0, 3.8269, 0.0, 0.0, 3.8119, 0.0, 0.0, 3.7968, 0.0, 0.0,
    3.7817, 0.0, 0.0, 3.7666, 0.0, 0.0, 3.7516, 0.0, 0.0, 3.7365, 0.0, 0.0,
    3.7214, 0.0, 0.0, 3.7063, 0.0, 0.0, 3.6913, 0.0, 0.0, 3.6762, 0.0, 0.0,
    3.6611, 0.0, 0.0, 3.6461, 0.0, 0.0, 203.631 };

  static const double b_a[3600] = { -204.966, -0.0, -0.0, -4.9108, -0.0, -0.0,
    -4.8555, -0.0, -0.0, -4.8003, -0.0, -0.0, -4.745, -0.0, -0.0, -4.6898, -0.0,
    -0.0, -4.6345, -0.0, -0.0, -4.5793, -0.0, -0.0, -4.524, -0.0, -0.0, -4.4688,
    -0.0, -0.0, -4.4136, -0.0, -0.0, -4.3584, -0.0, -0.0, -4.3032, -0.0, -0.0,
    -4.2481, -0.0, -0.0, -4.1929, -0.0, -0.0, -4.1378, -0.0, -0.0, -4.0826, -0.0,
    -0.0, -4.0275, -0.0, -0.0, -3.9724, -0.0, -0.0, -3.9174, -0.0, -0.0, -0.0,
    -204.966, -0.0, -0.0, -4.9108, -0.0, -0.0, -4.8555, -0.0, -0.0, -4.8003,
    -0.0, -0.0, -4.745, -0.0, -0.0, -4.6898, -0.0, -0.0, -4.6345, -0.0, -0.0,
    -4.5793, -0.0, -0.0, -4.524, -0.0, -0.0, -4.4688, -0.0, -0.0, -4.4136, -0.0,
    -0.0, -4.3584, -0.0, -0.0, -4.3032, -0.0, -0.0, -4.2481, -0.0, -0.0, -4.1929,
    -0.0, -0.0, -4.1378, -0.0, -0.0, -4.0826, -0.0, -0.0, -4.0275, -0.0, -0.0,
    -3.9724, -0.0, -0.0, -3.9174, -0.0, -0.0, -0.0, -204.966, -0.0, -0.0,
    -4.9108, -0.0, -0.0, -4.8555, -0.0, -0.0, -4.8003, -0.0, -0.0, -4.745, -0.0,
    -0.0, -4.6898, -0.0, -0.0, -4.6345, -0.0, -0.0, -4.5793, -0.0, -0.0, -4.524,
    -0.0, -0.0, -4.4688, -0.0, -0.0, -4.4136, -0.0, -0.0, -4.3584, -0.0, -0.0,
    -4.3032, -0.0, -0.0, -4.2481, -0.0, -0.0, -4.1929, -0.0, -0.0, -4.1378, -0.0,
    -0.0, -4.0826, -0.0, -0.0, -4.0275, -0.0, -0.0, -3.9724, -0.0, -0.0, -3.9174,
    -4.9108, -0.0, -0.0, -204.8907, -0.0, -0.0, -4.8357, -0.0, -0.0, -4.7807,
    -0.0, -0.0, -4.7258, -0.0, -0.0, -4.6708, -0.0, -0.0, -4.6158, -0.0, -0.0,
    -4.5609, -0.0, -0.0, -4.5059, -0.0, -0.0, -4.451, -0.0, -0.0, -4.3961, -0.0,
    -0.0, -4.3412, -0.0, -0.0, -4.2862, -0.0, -0.0, -4.2314, -0.0, -0.0, -4.1765,
    -0.0, -0.0, -4.1216, -0.0, -0.0, -4.0668, -0.0, -0.0, -4.0119, -0.0, -0.0,
    -3.9571, -0.0, -0.0, -3.9023, -0.0, -0.0, -0.0, -4.9108, -0.0, -0.0,
    -204.8907, -0.0, -0.0, -4.8357, -0.0, -0.0, -4.7807, -0.0, -0.0, -4.7258,
    -0.0, -0.0, -4.6708, -0.0, -0.0, -4.6158, -0.0, -0.0, -4.5609, -0.0, -0.0,
    -4.5059, -0.0, -0.0, -4.451, -0.0, -0.0, -4.3961, -0.0, -0.0, -4.3412, -0.0,
    -0.0, -4.2862, -0.0, -0.0, -4.2314, -0.0, -0.0, -4.1765, -0.0, -0.0, -4.1216,
    -0.0, -0.0, -4.0668, -0.0, -0.0, -4.0119, -0.0, -0.0, -3.9571, -0.0, -0.0,
    -3.9023, -0.0, -0.0, -0.0, -4.9108, -0.0, -0.0, -204.8907, -0.0, -0.0,
    -4.8357, -0.0, -0.0, -4.7807, -0.0, -0.0, -4.7258, -0.0, -0.0, -4.6708, -0.0,
    -0.0, -4.6158, -0.0, -0.0, -4.5609, -0.0, -0.0, -4.5059, -0.0, -0.0, -4.451,
    -0.0, -0.0, -4.3961, -0.0, -0.0, -4.3412, -0.0, -0.0, -4.2862, -0.0, -0.0,
    -4.2314, -0.0, -0.0, -4.1765, -0.0, -0.0, -4.1216, -0.0, -0.0, -4.0668, -0.0,
    -0.0, -4.0119, -0.0, -0.0, -3.9571, -0.0, -0.0, -3.9023, -4.8555, -0.0, -0.0,
    -4.8357, -0.0, -0.0, -204.8159, -0.0, -0.0, -4.7612, -0.0, -0.0, -4.7065,
    -0.0, -0.0, -4.6518, -0.0, -0.0, -4.5972, -0.0, -0.0, -4.5425, -0.0, -0.0,
    -4.4878, -0.0, -0.0, -4.4332, -0.0, -0.0, -4.3785, -0.0, -0.0, -4.3239, -0.0,
    -0.0, -4.2693, -0.0, -0.0, -4.2146, -0.0, -0.0, -4.16, -0.0, -0.0, -4.1054,
    -0.0, -0.0, -4.0509, -0.0, -0.0, -3.9963, -0.0, -0.0, -3.9418, -0.0, -0.0,
    -3.8872, -0.0, -0.0, -0.0, -4.8555, -0.0, -0.0, -4.8357, -0.0, -0.0,
    -204.8159, -0.0, -0.0, -4.7612, -0.0, -0.0, -4.7065, -0.0, -0.0, -4.6518,
    -0.0, -0.0, -4.5972, -0.0, -0.0, -4.5425, -0.0, -0.0, -4.4878, -0.0, -0.0,
    -4.4332, -0.0, -0.0, -4.3785, -0.0, -0.0, -4.3239, -0.0, -0.0, -4.2693, -0.0,
    -0.0, -4.2146, -0.0, -0.0, -4.16, -0.0, -0.0, -4.1054, -0.0, -0.0, -4.0509,
    -0.0, -0.0, -3.9963, -0.0, -0.0, -3.9418, -0.0, -0.0, -3.8872, -0.0, -0.0,
    -0.0, -4.8555, -0.0, -0.0, -4.8357, -0.0, -0.0, -204.8159, -0.0, -0.0,
    -4.7612, -0.0, -0.0, -4.7065, -0.0, -0.0, -4.6518, -0.0, -0.0, -4.5972, -0.0,
    -0.0, -4.5425, -0.0, -0.0, -4.4878, -0.0, -0.0, -4.4332, -0.0, -0.0, -4.3785,
    -0.0, -0.0, -4.3239, -0.0, -0.0, -4.2693, -0.0, -0.0, -4.2146, -0.0, -0.0,
    -4.16, -0.0, -0.0, -4.1054, -0.0, -0.0, -4.0509, -0.0, -0.0, -3.9963, -0.0,
    -0.0, -3.9418, -0.0, -0.0, -3.8872, -4.8003, -0.0, -0.0, -4.7807, -0.0, -0.0,
    -4.7612, -0.0, -0.0, -204.7417, -0.0, -0.0, -4.6873, -0.0, -0.0, -4.6329,
    -0.0, -0.0, -4.5785, -0.0, -0.0, -4.5241, -0.0, -0.0, -4.4697, -0.0, -0.0,
    -4.4154, -0.0, -0.0, -4.361, -0.0, -0.0, -4.3066, -0.0, -0.0, -4.2523, -0.0,
    -0.0, -4.1979, -0.0, -0.0, -4.1436, -0.0, -0.0, -4.0893, -0.0, -0.0, -4.035,
    -0.0, -0.0, -3.9807, -0.0, -0.0, -3.9264, -0.0, -0.0, -3.8722, -0.0, -0.0,
    -0.0, -4.8003, -0.0, -0.0, -4.7807, -0.0, -0.0, -4.7612, -0.0, -0.0,
    -204.7417, -0.0, -0.0, -4.6873, -0.0, -0.0, -4.6329, -0.0, -0.0, -4.5785,
    -0.0, -0.0, -4.5241, -0.0, -0.0, -4.4697, -0.0, -0.0, -4.4154, -0.0, -0.0,
    -4.361, -0.0, -0.0, -4.3066, -0.0, -0.0, -4.2523, -0.0, -0.0, -4.1979, -0.0,
    -0.0, -4.1436, -0.0, -0.0, -4.0893, -0.0, -0.0, -4.035, -0.0, -0.0, -3.9807,
    -0.0, -0.0, -3.9264, -0.0, -0.0, -3.8722, -0.0, -0.0, -0.0, -4.8003, -0.0,
    -0.0, -4.7807, -0.0, -0.0, -4.7612, -0.0, -0.0, -204.7417, -0.0, -0.0,
    -4.6873, -0.0, -0.0, -4.6329, -0.0, -0.0, -4.5785, -0.0, -0.0, -4.5241, -0.0,
    -0.0, -4.4697, -0.0, -0.0, -4.4154, -0.0, -0.0, -4.361, -0.0, -0.0, -4.3066,
    -0.0, -0.0, -4.2523, -0.0, -0.0, -4.1979, -0.0, -0.0, -4.1436, -0.0, -0.0,
    -4.0893, -0.0, -0.0, -4.035, -0.0, -0.0, -3.9807, -0.0, -0.0, -3.9264, -0.0,
    -0.0, -3.8722, -4.745, -0.0, -0.0, -4.7258, -0.0, -0.0, -4.7065, -0.0, -0.0,
    -4.6873, -0.0, -0.0, -204.668, -0.0, -0.0, -4.6139, -0.0, -0.0, -4.5598,
    -0.0, -0.0, -4.5057, -0.0, -0.0, -4.4516, -0.0, -0.0, -4.3975, -0.0, -0.0,
    -4.3434, -0.0, -0.0, -4.2894, -0.0, -0.0, -4.2353, -0.0, -0.0, -4.1812, -0.0,
    -0.0, -4.1272, -0.0, -0.0, -4.0731, -0.0, -0.0, -4.0191, -0.0, -0.0, -3.9651,
    -0.0, -0.0, -3.9111, -0.0, -0.0, -3.8571, -0.0, -0.0, -0.0, -4.745, -0.0,
    -0.0, -4.7258, -0.0, -0.0, -4.7065, -0.0, -0.0, -4.6873, -0.0, -0.0,
    -204.668, -0.0, -0.0, -4.6139, -0.0, -0.0, -4.5598, -0.0, -0.0, -4.5057,
    -0.0, -0.0, -4.4516, -0.0, -0.0, -4.3975, -0.0, -0.0, -4.3434, -0.0, -0.0,
    -4.2894, -0.0, -0.0, -4.2353, -0.0, -0.0, -4.1812, -0.0, -0.0, -4.1272, -0.0,
    -0.0, -4.0731, -0.0, -0.0, -4.0191, -0.0, -0.0, -3.9651, -0.0, -0.0, -3.9111,
    -0.0, -0.0, -3.8571, -0.0, -0.0, -0.0, -4.745, -0.0, -0.0, -4.7258, -0.0,
    -0.0, -4.7065, -0.0, -0.0, -4.6873, -0.0, -0.0, -204.668, -0.0, -0.0,
    -4.6139, -0.0, -0.0, -4.5598, -0.0, -0.0, -4.5057, -0.0, -0.0, -4.4516, -0.0,
    -0.0, -4.3975, -0.0, -0.0, -4.3434, -0.0, -0.0, -4.2894, -0.0, -0.0, -4.2353,
    -0.0, -0.0, -4.1812, -0.0, -0.0, -4.1272, -0.0, -0.0, -4.0731, -0.0, -0.0,
    -4.0191, -0.0, -0.0, -3.9651, -0.0, -0.0, -3.9111, -0.0, -0.0, -3.8571,
    -4.6898, -0.0, -0.0, -4.6708, -0.0, -0.0, -4.6518, -0.0, -0.0, -4.6329, -0.0,
    -0.0, -4.6139, -0.0, -0.0, -204.595, -0.0, -0.0, -4.5412, -0.0, -0.0,
    -4.4873, -0.0, -0.0, -4.4335, -0.0, -0.0, -4.3797, -0.0, -0.0, -4.3259, -0.0,
    -0.0, -4.2721, -0.0, -0.0, -4.2183, -0.0, -0.0, -4.1645, -0.0, -0.0, -4.1107,
    -0.0, -0.0, -4.057, -0.0, -0.0, -4.0032, -0.0, -0.0, -3.9495, -0.0, -0.0,
    -3.8957, -0.0, -0.0, -3.842, -0.0, -0.0, -0.0, -4.6898, -0.0, -0.0, -4.6708,
    -0.0, -0.0, -4.6518, -0.0, -0.0, -4.6329, -0.0, -0.0, -4.6139, -0.0, -0.0,
    -204.595, -0.0, -0.0, -4.5412, -0.0, -0.0, -4.4873, -0.0, -0.0, -4.4335,
    -0.0, -0.0, -4.3797, -0.0, -0.0, -4.3259, -0.0, -0.0, -4.2721, -0.0, -0.0,
    -4.2183, -0.0, -0.0, -4.1645, -0.0, -0.0, -4.1107, -0.0, -0.0, -4.057, -0.0,
    -0.0, -4.0032, -0.0, -0.0, -3.9495, -0.0, -0.0, -3.8957, -0.0, -0.0, -3.842,
    -0.0, -0.0, -0.0, -4.6898, -0.0, -0.0, -4.6708, -0.0, -0.0, -4.6518, -0.0,
    -0.0, -4.6329, -0.0, -0.0, -4.6139, -0.0, -0.0, -204.595, -0.0, -0.0,
    -4.5412, -0.0, -0.0, -4.4873, -0.0, -0.0, -4.4335, -0.0, -0.0, -4.3797, -0.0,
    -0.0, -4.3259, -0.0, -0.0, -4.2721, -0.0, -0.0, -4.2183, -0.0, -0.0, -4.1645,
    -0.0, -0.0, -4.1107, -0.0, -0.0, -4.057, -0.0, -0.0, -4.0032, -0.0, -0.0,
    -3.9495, -0.0, -0.0, -3.8957, -0.0, -0.0, -3.842, -4.6345, -0.0, -0.0,
    -4.6158, -0.0, -0.0, -4.5972, -0.0, -0.0, -4.5785, -0.0, -0.0, -4.5598, -0.0,
    -0.0, -4.5412, -0.0, -0.0, -204.5225, -0.0, -0.0, -4.469, -0.0, -0.0,
    -4.4154, -0.0, -0.0, -4.3619, -0.0, -0.0, -4.3084, -0.0, -0.0, -4.2548, -0.0,
    -0.0, -4.2013, -0.0, -0.0, -4.1478, -0.0, -0.0, -4.0943, -0.0, -0.0, -4.0408,
    -0.0, -0.0, -3.9873, -0.0, -0.0, -3.9338, -0.0, -0.0, -3.8804, -0.0, -0.0,
    -3.8269, -0.0, -0.0, -0.0, -4.6345, -0.0, -0.0, -4.6158, -0.0, -0.0, -4.5972,
    -0.0, -0.0, -4.5785, -0.0, -0.0, -4.5598, -0.0, -0.0, -4.5412, -0.0, -0.0,
    -204.5225, -0.0, -0.0, -4.469, -0.0, -0.0, -4.4154, -0.0, -0.0, -4.3619,
    -0.0, -0.0, -4.3084, -0.0, -0.0, -4.2548, -0.0, -0.0, -4.2013, -0.0, -0.0,
    -4.1478, -0.0, -0.0, -4.0943, -0.0, -0.0, -4.0408, -0.0, -0.0, -3.9873, -0.0,
    -0.0, -3.9338, -0.0, -0.0, -3.8804, -0.0, -0.0, -3.8269, -0.0, -0.0, -0.0,
    -4.6345, -0.0, -0.0, -4.6158, -0.0, -0.0, -4.5972, -0.0, -0.0, -4.5785, -0.0,
    -0.0, -4.5598, -0.0, -0.0, -4.5412, -0.0, -0.0, -204.5225, -0.0, -0.0,
    -4.469, -0.0, -0.0, -4.4154, -0.0, -0.0, -4.3619, -0.0, -0.0, -4.3084, -0.0,
    -0.0, -4.2548, -0.0, -0.0, -4.2013, -0.0, -0.0, -4.1478, -0.0, -0.0, -4.0943,
    -0.0, -0.0, -4.0408, -0.0, -0.0, -3.9873, -0.0, -0.0, -3.9338, -0.0, -0.0,
    -3.8804, -0.0, -0.0, -3.8269, -4.5793, -0.0, -0.0, -4.5609, -0.0, -0.0,
    -4.5425, -0.0, -0.0, -4.5241, -0.0, -0.0, -4.5057, -0.0, -0.0, -4.4873, -0.0,
    -0.0, -4.469, -0.0, -0.0, -204.4506, -0.0, -0.0, -4.3973, -0.0, -0.0,
    -4.3441, -0.0, -0.0, -4.2908, -0.0, -0.0, -4.2376, -0.0, -0.0, -4.1843, -0.0,
    -0.0, -4.1311, -0.0, -0.0, -4.0779, -0.0, -0.0, -4.0246, -0.0, -0.0, -3.9714,
    -0.0, -0.0, -3.9182, -0.0, -0.0, -3.865, -0.0, -0.0, -3.8119, -0.0, -0.0,
    -0.0, -4.5793, -0.0, -0.0, -4.5609, -0.0, -0.0, -4.5425, -0.0, -0.0, -4.5241,
    -0.0, -0.0, -4.5057, -0.0, -0.0, -4.4873, -0.0, -0.0, -4.469, -0.0, -0.0,
    -204.4506, -0.0, -0.0, -4.3973, -0.0, -0.0, -4.3441, -0.0, -0.0, -4.2908,
    -0.0, -0.0, -4.2376, -0.0, -0.0, -4.1843, -0.0, -0.0, -4.1311, -0.0, -0.0,
    -4.0779, -0.0, -0.0, -4.0246, -0.0, -0.0, -3.9714, -0.0, -0.0, -3.9182, -0.0,
    -0.0, -3.865, -0.0, -0.0, -3.8119, -0.0, -0.0, -0.0, -4.5793, -0.0, -0.0,
    -4.5609, -0.0, -0.0, -4.5425, -0.0, -0.0, -4.5241, -0.0, -0.0, -4.5057, -0.0,
    -0.0, -4.4873, -0.0, -0.0, -4.469, -0.0, -0.0, -204.4506, -0.0, -0.0,
    -4.3973, -0.0, -0.0, -4.3441, -0.0, -0.0, -4.2908, -0.0, -0.0, -4.2376, -0.0,
    -0.0, -4.1843, -0.0, -0.0, -4.1311, -0.0, -0.0, -4.0779, -0.0, -0.0, -4.0246,
    -0.0, -0.0, -3.9714, -0.0, -0.0, -3.9182, -0.0, -0.0, -3.865, -0.0, -0.0,
    -3.8119, -4.524, -0.0, -0.0, -4.5059, -0.0, -0.0, -4.4878, -0.0, -0.0,
    -4.4697, -0.0, -0.0, -4.4516, -0.0, -0.0, -4.4335, -0.0, -0.0, -4.4154, -0.0,
    -0.0, -4.3973, -0.0, -0.0, -204.3792, -0.0, -0.0, -4.3262, -0.0, -0.0,
    -4.2733, -0.0, -0.0, -4.2203, -0.0, -0.0, -4.1673, -0.0, -0.0, -4.1144, -0.0,
    -0.0, -4.0614, -0.0, -0.0, -4.0085, -0.0, -0.0, -3.9555, -0.0, -0.0, -3.9026,
    -0.0, -0.0, -3.8497, -0.0, -0.0, -3.7968, -0.0, -0.0, -0.0, -4.524, -0.0,
    -0.0, -4.5059, -0.0, -0.0, -4.4878, -0.0, -0.0, -4.4697, -0.0, -0.0, -4.4516,
    -0.0, -0.0, -4.4335, -0.0, -0.0, -4.4154, -0.0, -0.0, -4.3973, -0.0, -0.0,
    -204.3792, -0.0, -0.0, -4.3262, -0.0, -0.0, -4.2733, -0.0, -0.0, -4.2203,
    -0.0, -0.0, -4.1673, -0.0, -0.0, -4.1144, -0.0, -0.0, -4.0614, -0.0, -0.0,
    -4.0085, -0.0, -0.0, -3.9555, -0.0, -0.0, -3.9026, -0.0, -0.0, -3.8497, -0.0,
    -0.0, -3.7968, -0.0, -0.0, -0.0, -4.524, -0.0, -0.0, -4.5059, -0.0, -0.0,
    -4.4878, -0.0, -0.0, -4.4697, -0.0, -0.0, -4.4516, -0.0, -0.0, -4.4335, -0.0,
    -0.0, -4.4154, -0.0, -0.0, -4.3973, -0.0, -0.0, -204.3792, -0.0, -0.0,
    -4.3262, -0.0, -0.0, -4.2733, -0.0, -0.0, -4.2203, -0.0, -0.0, -4.1673, -0.0,
    -0.0, -4.1144, -0.0, -0.0, -4.0614, -0.0, -0.0, -4.0085, -0.0, -0.0, -3.9555,
    -0.0, -0.0, -3.9026, -0.0, -0.0, -3.8497, -0.0, -0.0, -3.7968, -4.4688, -0.0,
    -0.0, -4.451, -0.0, -0.0, -4.4332, -0.0, -0.0, -4.4154, -0.0, -0.0, -4.3975,
    -0.0, -0.0, -4.3797, -0.0, -0.0, -4.3619, -0.0, -0.0, -4.3441, -0.0, -0.0,
    -4.3262, -0.0, -0.0, -204.3084, -0.0, -0.0, -4.2557, -0.0, -0.0, -4.203,
    -0.0, -0.0, -4.1504, -0.0, -0.0, -4.0977, -0.0, -0.0, -4.045, -0.0, -0.0,
    -3.9923, -0.0, -0.0, -3.9397, -0.0, -0.0, -3.887, -0.0, -0.0, -3.8344, -0.0,
    -0.0, -3.7817, -0.0, -0.0, -0.0, -4.4688, -0.0, -0.0, -4.451, -0.0, -0.0,
    -4.4332, -0.0, -0.0, -4.4154, -0.0, -0.0, -4.3975, -0.0, -0.0, -4.3797, -0.0,
    -0.0, -4.3619, -0.0, -0.0, -4.3441, -0.0, -0.0, -4.3262, -0.0, -0.0,
    -204.3084, -0.0, -0.0, -4.2557, -0.0, -0.0, -4.203, -0.0, -0.0, -4.1504,
    -0.0, -0.0, -4.0977, -0.0, -0.0, -4.045, -0.0, -0.0, -3.9923, -0.0, -0.0,
    -3.9397, -0.0, -0.0, -3.887, -0.0, -0.0, -3.8344, -0.0, -0.0, -3.7817, -0.0,
    -0.0, -0.0, -4.4688, -0.0, -0.0, -4.451, -0.0, -0.0, -4.4332, -0.0, -0.0,
    -4.4154, -0.0, -0.0, -4.3975, -0.0, -0.0, -4.3797, -0.0, -0.0, -4.3619, -0.0,
    -0.0, -4.3441, -0.0, -0.0, -4.3262, -0.0, -0.0, -204.3084, -0.0, -0.0,
    -4.2557, -0.0, -0.0, -4.203, -0.0, -0.0, -4.1504, -0.0, -0.0, -4.0977, -0.0,
    -0.0, -4.045, -0.0, -0.0, -3.9923, -0.0, -0.0, -3.9397, -0.0, -0.0, -3.887,
    -0.0, -0.0, -3.8344, -0.0, -0.0, -3.7817, -4.4136, -0.0, -0.0, -4.3961, -0.0,
    -0.0, -4.3785, -0.0, -0.0, -4.361, -0.0, -0.0, -4.3434, -0.0, -0.0, -4.3259,
    -0.0, -0.0, -4.3084, -0.0, -0.0, -4.2908, -0.0, -0.0, -4.2733, -0.0, -0.0,
    -4.2557, -0.0, -0.0, -204.2382, -0.0, -0.0, -4.1858, -0.0, -0.0, -4.1334,
    -0.0, -0.0, -4.081, -0.0, -0.0, -4.0286, -0.0, -0.0, -3.9762, -0.0, -0.0,
    -3.9238, -0.0, -0.0, -3.8714, -0.0, -0.0, -3.819, -0.0, -0.0, -3.7666, -0.0,
    -0.0, -0.0, -4.4136, -0.0, -0.0, -4.3961, -0.0, -0.0, -4.3785, -0.0, -0.0,
    -4.361, -0.0, -0.0, -4.3434, -0.0, -0.0, -4.3259, -0.0, -0.0, -4.3084, -0.0,
    -0.0, -4.2908, -0.0, -0.0, -4.2733, -0.0, -0.0, -4.2557, -0.0, -0.0,
    -204.2382, -0.0, -0.0, -4.1858, -0.0, -0.0, -4.1334, -0.0, -0.0, -4.081,
    -0.0, -0.0, -4.0286, -0.0, -0.0, -3.9762, -0.0, -0.0, -3.9238, -0.0, -0.0,
    -3.8714, -0.0, -0.0, -3.819, -0.0, -0.0, -3.7666, -0.0, -0.0, -0.0, -4.4136,
    -0.0, -0.0, -4.3961, -0.0, -0.0, -4.3785, -0.0, -0.0, -4.361, -0.0, -0.0,
    -4.3434, -0.0, -0.0, -4.3259, -0.0, -0.0, -4.3084, -0.0, -0.0, -4.2908, -0.0,
    -0.0, -4.2733, -0.0, -0.0, -4.2557, -0.0, -0.0, -204.2382, -0.0, -0.0,
    -4.1858, -0.0, -0.0, -4.1334, -0.0, -0.0, -4.081, -0.0, -0.0, -4.0286, -0.0,
    -0.0, -3.9762, -0.0, -0.0, -3.9238, -0.0, -0.0, -3.8714, -0.0, -0.0, -3.819,
    -0.0, -0.0, -3.7666, -4.3584, -0.0, -0.0, -4.3412, -0.0, -0.0, -4.3239, -0.0,
    -0.0, -4.3066, -0.0, -0.0, -4.2894, -0.0, -0.0, -4.2721, -0.0, -0.0, -4.2548,
    -0.0, -0.0, -4.2376, -0.0, -0.0, -4.2203, -0.0, -0.0, -4.203, -0.0, -0.0,
    -4.1858, -0.0, -0.0, -204.1685, -0.0, -0.0, -4.1164, -0.0, -0.0, -4.0643,
    -0.0, -0.0, -4.0121, -0.0, -0.0, -3.96, -0.0, -0.0, -3.9079, -0.0, -0.0,
    -3.8558, -0.0, -0.0, -3.8037, -0.0, -0.0, -3.7516, -0.0, -0.0, -0.0, -4.3584,
    -0.0, -0.0, -4.3412, -0.0, -0.0, -4.3239, -0.0, -0.0, -4.3066, -0.0, -0.0,
    -4.2894, -0.0, -0.0, -4.2721, -0.0, -0.0, -4.2548, -0.0, -0.0, -4.2376, -0.0,
    -0.0, -4.2203, -0.0, -0.0, -4.203, -0.0, -0.0, -4.1858, -0.0, -0.0,
    -204.1685, -0.0, -0.0, -4.1164, -0.0, -0.0, -4.0643, -0.0, -0.0, -4.0121,
    -0.0, -0.0, -3.96, -0.0, -0.0, -3.9079, -0.0, -0.0, -3.8558, -0.0, -0.0,
    -3.8037, -0.0, -0.0, -3.7516, -0.0, -0.0, -0.0, -4.3584, -0.0, -0.0, -4.3412,
    -0.0, -0.0, -4.3239, -0.0, -0.0, -4.3066, -0.0, -0.0, -4.2894, -0.0, -0.0,
    -4.2721, -0.0, -0.0, -4.2548, -0.0, -0.0, -4.2376, -0.0, -0.0, -4.2203, -0.0,
    -0.0, -4.203, -0.0, -0.0, -4.1858, -0.0, -0.0, -204.1685, -0.0, -0.0,
    -4.1164, -0.0, -0.0, -4.0643, -0.0, -0.0, -4.0121, -0.0, -0.0, -3.96, -0.0,
    -0.0, -3.9079, -0.0, -0.0, -3.8558, -0.0, -0.0, -3.8037, -0.0, -0.0, -3.7516,
    -4.3032, -0.0, -0.0, -4.2862, -0.0, -0.0, -4.2693, -0.0, -0.0, -4.2523, -0.0,
    -0.0, -4.2353, -0.0, -0.0, -4.2183, -0.0, -0.0, -4.2013, -0.0, -0.0, -4.1843,
    -0.0, -0.0, -4.1673, -0.0, -0.0, -4.1504, -0.0, -0.0, -4.1334, -0.0, -0.0,
    -4.1164, -0.0, -0.0, -204.0994, -0.0, -0.0, -4.0476, -0.0, -0.0, -3.9957,
    -0.0, -0.0, -3.9438, -0.0, -0.0, -3.892, -0.0, -0.0, -3.8402, -0.0, -0.0,
    -3.7883, -0.0, -0.0, -3.7365, -0.0, -0.0, -0.0, -4.3032, -0.0, -0.0, -4.2862,
    -0.0, -0.0, -4.2693, -0.0, -0.0, -4.2523, -0.0, -0.0, -4.2353, -0.0, -0.0,
    -4.2183, -0.0, -0.0, -4.2013, -0.0, -0.0, -4.1843, -0.0, -0.0, -4.1673, -0.0,
    -0.0, -4.1504, -0.0, -0.0, -4.1334, -0.0, -0.0, -4.1164, -0.0, -0.0,
    -204.0994, -0.0, -0.0, -4.0476, -0.0, -0.0, -3.9957, -0.0, -0.0, -3.9438,
    -0.0, -0.0, -3.892, -0.0, -0.0, -3.8402, -0.0, -0.0, -3.7883, -0.0, -0.0,
    -3.7365, -0.0, -0.0, -0.0, -4.3032, -0.0, -0.0, -4.2862, -0.0, -0.0, -4.2693,
    -0.0, -0.0, -4.2523, -0.0, -0.0, -4.2353, -0.0, -0.0, -4.2183, -0.0, -0.0,
    -4.2013, -0.0, -0.0, -4.1843, -0.0, -0.0, -4.1673, -0.0, -0.0, -4.1504, -0.0,
    -0.0, -4.1334, -0.0, -0.0, -4.1164, -0.0, -0.0, -204.0994, -0.0, -0.0,
    -4.0476, -0.0, -0.0, -3.9957, -0.0, -0.0, -3.9438, -0.0, -0.0, -3.892, -0.0,
    -0.0, -3.8402, -0.0, -0.0, -3.7883, -0.0, -0.0, -3.7365, -4.2481, -0.0, -0.0,
    -4.2314, -0.0, -0.0, -4.2146, -0.0, -0.0, -4.1979, -0.0, -0.0, -4.1812, -0.0,
    -0.0, -4.1645, -0.0, -0.0, -4.1478, -0.0, -0.0, -4.1311, -0.0, -0.0, -4.1144,
    -0.0, -0.0, -4.0977, -0.0, -0.0, -4.081, -0.0, -0.0, -4.0643, -0.0, -0.0,
    -4.0476, -0.0, -0.0, -204.0308, -0.0, -0.0, -3.9793, -0.0, -0.0, -3.9277,
    -0.0, -0.0, -3.8761, -0.0, -0.0, -3.8245, -0.0, -0.0, -3.773, -0.0, -0.0,
    -3.7214, -0.0, -0.0, -0.0, -4.2481, -0.0, -0.0, -4.2314, -0.0, -0.0, -4.2146,
    -0.0, -0.0, -4.1979, -0.0, -0.0, -4.1812, -0.0, -0.0, -4.1645, -0.0, -0.0,
    -4.1478, -0.0, -0.0, -4.1311, -0.0, -0.0, -4.1144, -0.0, -0.0, -4.0977, -0.0,
    -0.0, -4.081, -0.0, -0.0, -4.0643, -0.0, -0.0, -4.0476, -0.0, -0.0,
    -204.0308, -0.0, -0.0, -3.9793, -0.0, -0.0, -3.9277, -0.0, -0.0, -3.8761,
    -0.0, -0.0, -3.8245, -0.0, -0.0, -3.773, -0.0, -0.0, -3.7214, -0.0, -0.0,
    -0.0, -4.2481, -0.0, -0.0, -4.2314, -0.0, -0.0, -4.2146, -0.0, -0.0, -4.1979,
    -0.0, -0.0, -4.1812, -0.0, -0.0, -4.1645, -0.0, -0.0, -4.1478, -0.0, -0.0,
    -4.1311, -0.0, -0.0, -4.1144, -0.0, -0.0, -4.0977, -0.0, -0.0, -4.081, -0.0,
    -0.0, -4.0643, -0.0, -0.0, -4.0476, -0.0, -0.0, -204.0308, -0.0, -0.0,
    -3.9793, -0.0, -0.0, -3.9277, -0.0, -0.0, -3.8761, -0.0, -0.0, -3.8245, -0.0,
    -0.0, -3.773, -0.0, -0.0, -3.7214, -4.1929, -0.0, -0.0, -4.1765, -0.0, -0.0,
    -4.16, -0.0, -0.0, -4.1436, -0.0, -0.0, -4.1272, -0.0, -0.0, -4.1107, -0.0,
    -0.0, -4.0943, -0.0, -0.0, -4.0779, -0.0, -0.0, -4.0614, -0.0, -0.0, -4.045,
    -0.0, -0.0, -4.0286, -0.0, -0.0, -4.0121, -0.0, -0.0, -3.9957, -0.0, -0.0,
    -3.9793, -0.0, -0.0, -203.9628, -0.0, -0.0, -3.9115, -0.0, -0.0, -3.8602,
    -0.0, -0.0, -3.8089, -0.0, -0.0, -3.7576, -0.0, -0.0, -3.7063, -0.0, -0.0,
    -0.0, -4.1929, -0.0, -0.0, -4.1765, -0.0, -0.0, -4.16, -0.0, -0.0, -4.1436,
    -0.0, -0.0, -4.1272, -0.0, -0.0, -4.1107, -0.0, -0.0, -4.0943, -0.0, -0.0,
    -4.0779, -0.0, -0.0, -4.0614, -0.0, -0.0, -4.045, -0.0, -0.0, -4.0286, -0.0,
    -0.0, -4.0121, -0.0, -0.0, -3.9957, -0.0, -0.0, -3.9793, -0.0, -0.0,
    -203.9628, -0.0, -0.0, -3.9115, -0.0, -0.0, -3.8602, -0.0, -0.0, -3.8089,
    -0.0, -0.0, -3.7576, -0.0, -0.0, -3.7063, -0.0, -0.0, -0.0, -4.1929, -0.0,
    -0.0, -4.1765, -0.0, -0.0, -4.16, -0.0, -0.0, -4.1436, -0.0, -0.0, -4.1272,
    -0.0, -0.0, -4.1107, -0.0, -0.0, -4.0943, -0.0, -0.0, -4.0779, -0.0, -0.0,
    -4.0614, -0.0, -0.0, -4.045, -0.0, -0.0, -4.0286, -0.0, -0.0, -4.0121, -0.0,
    -0.0, -3.9957, -0.0, -0.0, -3.9793, -0.0, -0.0, -203.9628, -0.0, -0.0,
    -3.9115, -0.0, -0.0, -3.8602, -0.0, -0.0, -3.8089, -0.0, -0.0, -3.7576, -0.0,
    -0.0, -3.7063, -4.1378, -0.0, -0.0, -4.1216, -0.0, -0.0, -4.1054, -0.0, -0.0,
    -4.0893, -0.0, -0.0, -4.0731, -0.0, -0.0, -4.057, -0.0, -0.0, -4.0408, -0.0,
    -0.0, -4.0246, -0.0, -0.0, -4.0085, -0.0, -0.0, -3.9923, -0.0, -0.0, -3.9762,
    -0.0, -0.0, -3.96, -0.0, -0.0, -3.9438, -0.0, -0.0, -3.9277, -0.0, -0.0,
    -3.9115, -0.0, -0.0, -203.8954, -0.0, -0.0, -3.8443, -0.0, -0.0, -3.7933,
    -0.0, -0.0, -3.7423, -0.0, -0.0, -3.6913, -0.0, -0.0, -0.0, -4.1378, -0.0,
    -0.0, -4.1216, -0.0, -0.0, -4.1054, -0.0, -0.0, -4.0893, -0.0, -0.0, -4.0731,
    -0.0, -0.0, -4.057, -0.0, -0.0, -4.0408, -0.0, -0.0, -4.0246, -0.0, -0.0,
    -4.0085, -0.0, -0.0, -3.9923, -0.0, -0.0, -3.9762, -0.0, -0.0, -3.96, -0.0,
    -0.0, -3.9438, -0.0, -0.0, -3.9277, -0.0, -0.0, -3.9115, -0.0, -0.0,
    -203.8954, -0.0, -0.0, -3.8443, -0.0, -0.0, -3.7933, -0.0, -0.0, -3.7423,
    -0.0, -0.0, -3.6913, -0.0, -0.0, -0.0, -4.1378, -0.0, -0.0, -4.1216, -0.0,
    -0.0, -4.1054, -0.0, -0.0, -4.0893, -0.0, -0.0, -4.0731, -0.0, -0.0, -4.057,
    -0.0, -0.0, -4.0408, -0.0, -0.0, -4.0246, -0.0, -0.0, -4.0085, -0.0, -0.0,
    -3.9923, -0.0, -0.0, -3.9762, -0.0, -0.0, -3.96, -0.0, -0.0, -3.9438, -0.0,
    -0.0, -3.9277, -0.0, -0.0, -3.9115, -0.0, -0.0, -203.8954, -0.0, -0.0,
    -3.8443, -0.0, -0.0, -3.7933, -0.0, -0.0, -3.7423, -0.0, -0.0, -3.6913,
    -4.0826, -0.0, -0.0, -4.0668, -0.0, -0.0, -4.0509, -0.0, -0.0, -4.035, -0.0,
    -0.0, -4.0191, -0.0, -0.0, -4.0032, -0.0, -0.0, -3.9873, -0.0, -0.0, -3.9714,
    -0.0, -0.0, -3.9555, -0.0, -0.0, -3.9397, -0.0, -0.0, -3.9238, -0.0, -0.0,
    -3.9079, -0.0, -0.0, -3.892, -0.0, -0.0, -3.8761, -0.0, -0.0, -3.8602, -0.0,
    -0.0, -3.8443, -0.0, -0.0, -203.8285, -0.0, -0.0, -3.7777, -0.0, -0.0,
    -3.727, -0.0, -0.0, -3.6762, -0.0, -0.0, -0.0, -4.0826, -0.0, -0.0, -4.0668,
    -0.0, -0.0, -4.0509, -0.0, -0.0, -4.035, -0.0, -0.0, -4.0191, -0.0, -0.0,
    -4.0032, -0.0, -0.0, -3.9873, -0.0, -0.0, -3.9714, -0.0, -0.0, -3.9555, -0.0,
    -0.0, -3.9397, -0.0, -0.0, -3.9238, -0.0, -0.0, -3.9079, -0.0, -0.0, -3.892,
    -0.0, -0.0, -3.8761, -0.0, -0.0, -3.8602, -0.0, -0.0, -3.8443, -0.0, -0.0,
    -203.8285, -0.0, -0.0, -3.7777, -0.0, -0.0, -3.727, -0.0, -0.0, -3.6762,
    -0.0, -0.0, -0.0, -4.0826, -0.0, -0.0, -4.0668, -0.0, -0.0, -4.0509, -0.0,
    -0.0, -4.035, -0.0, -0.0, -4.0191, -0.0, -0.0, -4.0032, -0.0, -0.0, -3.9873,
    -0.0, -0.0, -3.9714, -0.0, -0.0, -3.9555, -0.0, -0.0, -3.9397, -0.0, -0.0,
    -3.9238, -0.0, -0.0, -3.9079, -0.0, -0.0, -3.892, -0.0, -0.0, -3.8761, -0.0,
    -0.0, -3.8602, -0.0, -0.0, -3.8443, -0.0, -0.0, -203.8285, -0.0, -0.0,
    -3.7777, -0.0, -0.0, -3.727, -0.0, -0.0, -3.6762, -4.0275, -0.0, -0.0,
    -4.0119, -0.0, -0.0, -3.9963, -0.0, -0.0, -3.9807, -0.0, -0.0, -3.9651, -0.0,
    -0.0, -3.9495, -0.0, -0.0, -3.9338, -0.0, -0.0, -3.9182, -0.0, -0.0, -3.9026,
    -0.0, -0.0, -3.887, -0.0, -0.0, -3.8714, -0.0, -0.0, -3.8558, -0.0, -0.0,
    -3.8402, -0.0, -0.0, -3.8245, -0.0, -0.0, -3.8089, -0.0, -0.0, -3.7933, -0.0,
    -0.0, -3.7777, -0.0, -0.0, -203.7621, -0.0, -0.0, -3.7116, -0.0, -0.0,
    -3.6611, -0.0, -0.0, -0.0, -4.0275, -0.0, -0.0, -4.0119, -0.0, -0.0, -3.9963,
    -0.0, -0.0, -3.9807, -0.0, -0.0, -3.9651, -0.0, -0.0, -3.9495, -0.0, -0.0,
    -3.9338, -0.0, -0.0, -3.9182, -0.0, -0.0, -3.9026, -0.0, -0.0, -3.887, -0.0,
    -0.0, -3.8714, -0.0, -0.0, -3.8558, -0.0, -0.0, -3.8402, -0.0, -0.0, -3.8245,
    -0.0, -0.0, -3.8089, -0.0, -0.0, -3.7933, -0.0, -0.0, -3.7777, -0.0, -0.0,
    -203.7621, -0.0, -0.0, -3.7116, -0.0, -0.0, -3.6611, -0.0, -0.0, -0.0,
    -4.0275, -0.0, -0.0, -4.0119, -0.0, -0.0, -3.9963, -0.0, -0.0, -3.9807, -0.0,
    -0.0, -3.9651, -0.0, -0.0, -3.9495, -0.0, -0.0, -3.9338, -0.0, -0.0, -3.9182,
    -0.0, -0.0, -3.9026, -0.0, -0.0, -3.887, -0.0, -0.0, -3.8714, -0.0, -0.0,
    -3.8558, -0.0, -0.0, -3.8402, -0.0, -0.0, -3.8245, -0.0, -0.0, -3.8089, -0.0,
    -0.0, -3.7933, -0.0, -0.0, -3.7777, -0.0, -0.0, -203.7621, -0.0, -0.0,
    -3.7116, -0.0, -0.0, -3.6611, -3.9724, -0.0, -0.0, -3.9571, -0.0, -0.0,
    -3.9418, -0.0, -0.0, -3.9264, -0.0, -0.0, -3.9111, -0.0, -0.0, -3.8957, -0.0,
    -0.0, -3.8804, -0.0, -0.0, -3.865, -0.0, -0.0, -3.8497, -0.0, -0.0, -3.8344,
    -0.0, -0.0, -3.819, -0.0, -0.0, -3.8037, -0.0, -0.0, -3.7883, -0.0, -0.0,
    -3.773, -0.0, -0.0, -3.7576, -0.0, -0.0, -3.7423, -0.0, -0.0, -3.727, -0.0,
    -0.0, -3.7116, -0.0, -0.0, -203.6963, -0.0, -0.0, -3.6461, -0.0, -0.0, -0.0,
    -3.9724, -0.0, -0.0, -3.9571, -0.0, -0.0, -3.9418, -0.0, -0.0, -3.9264, -0.0,
    -0.0, -3.9111, -0.0, -0.0, -3.8957, -0.0, -0.0, -3.8804, -0.0, -0.0, -3.865,
    -0.0, -0.0, -3.8497, -0.0, -0.0, -3.8344, -0.0, -0.0, -3.819, -0.0, -0.0,
    -3.8037, -0.0, -0.0, -3.7883, -0.0, -0.0, -3.773, -0.0, -0.0, -3.7576, -0.0,
    -0.0, -3.7423, -0.0, -0.0, -3.727, -0.0, -0.0, -3.7116, -0.0, -0.0,
    -203.6963, -0.0, -0.0, -3.6461, -0.0, -0.0, -0.0, -3.9724, -0.0, -0.0,
    -3.9571, -0.0, -0.0, -3.9418, -0.0, -0.0, -3.9264, -0.0, -0.0, -3.9111, -0.0,
    -0.0, -3.8957, -0.0, -0.0, -3.8804, -0.0, -0.0, -3.865, -0.0, -0.0, -3.8497,
    -0.0, -0.0, -3.8344, -0.0, -0.0, -3.819, -0.0, -0.0, -3.8037, -0.0, -0.0,
    -3.7883, -0.0, -0.0, -3.773, -0.0, -0.0, -3.7576, -0.0, -0.0, -3.7423, -0.0,
    -0.0, -3.727, -0.0, -0.0, -3.7116, -0.0, -0.0, -203.6963, -0.0, -0.0,
    -3.6461, -3.9174, -0.0, -0.0, -3.9023, -0.0, -0.0, -3.8872, -0.0, -0.0,
    -3.8722, -0.0, -0.0, -3.8571, -0.0, -0.0, -3.842, -0.0, -0.0, -3.8269, -0.0,
    -0.0, -3.8119, -0.0, -0.0, -3.7968, -0.0, -0.0, -3.7817, -0.0, -0.0, -3.7666,
    -0.0, -0.0, -3.7516, -0.0, -0.0, -3.7365, -0.0, -0.0, -3.7214, -0.0, -0.0,
    -3.7063, -0.0, -0.0, -3.6913, -0.0, -0.0, -3.6762, -0.0, -0.0, -3.6611, -0.0,
    -0.0, -3.6461, -0.0, -0.0, -203.631, -0.0, -0.0, -0.0, -3.9174, -0.0, -0.0,
    -3.9023, -0.0, -0.0, -3.8872, -0.0, -0.0, -3.8722, -0.0, -0.0, -3.8571, -0.0,
    -0.0, -3.842, -0.0, -0.0, -3.8269, -0.0, -0.0, -3.8119, -0.0, -0.0, -3.7968,
    -0.0, -0.0, -3.7817, -0.0, -0.0, -3.7666, -0.0, -0.0, -3.7516, -0.0, -0.0,
    -3.7365, -0.0, -0.0, -3.7214, -0.0, -0.0, -3.7063, -0.0, -0.0, -3.6913, -0.0,
    -0.0, -3.6762, -0.0, -0.0, -3.6611, -0.0, -0.0, -3.6461, -0.0, -0.0,
    -203.631, -0.0, -0.0, -0.0, -3.9174, -0.0, -0.0, -3.9023, -0.0, -0.0,
    -3.8872, -0.0, -0.0, -3.8722, -0.0, -0.0, -3.8571, -0.0, -0.0, -3.842, -0.0,
    -0.0, -3.8269, -0.0, -0.0, -3.8119, -0.0, -0.0, -3.7968, -0.0, -0.0, -3.7817,
    -0.0, -0.0, -3.7666, -0.0, -0.0, -3.7516, -0.0, -0.0, -3.7365, -0.0, -0.0,
    -3.7214, -0.0, -0.0, -3.7063, -0.0, -0.0, -3.6913, -0.0, -0.0, -3.6762, -0.0,
    -0.0, -3.6611, -0.0, -0.0, -3.6461, -0.0, -0.0, -203.631 };

  static const double b[360] = { 30.5242, 0.0, 0.0, 595.9023, 0.0, 0.0, 0.0,
    30.5242, 0.0, 0.0, 595.9023, 0.0, 0.0, 0.0, 30.5242, 0.0, 0.0, 595.9023,
    30.0901, 0.0, 0.0, 589.2677, 0.0, 0.0, 0.0, 30.0901, 0.0, 0.0, 589.2677, 0.0,
    0.0, 0.0, 30.0901, 0.0, 0.0, 589.2677, 29.6577, 0.0, 0.0, 582.6333, 0.0, 0.0,
    0.0, 29.6577, 0.0, 0.0, 582.6333, 0.0, 0.0, 0.0, 29.6577, 0.0, 0.0, 582.6333,
    29.227, 0.0, 0.0, 575.9992, 0.0, 0.0, 0.0, 29.227, 0.0, 0.0, 575.9992, 0.0,
    0.0, 0.0, 29.227, 0.0, 0.0, 575.9992, 28.798, 0.0, 0.0, 569.3655, 0.0, 0.0,
    0.0, 28.798, 0.0, 0.0, 569.3655, 0.0, 0.0, 0.0, 28.798, 0.0, 0.0, 569.3655,
    28.3706, 0.0, 0.0, 562.7325, 0.0, 0.0, 0.0, 28.3706, 0.0, 0.0, 562.7325, 0.0,
    0.0, 0.0, 28.3706, 0.0, 0.0, 562.7325, 27.9449, 0.0, 0.0, 556.1002, 0.0, 0.0,
    0.0, 27.9449, 0.0, 0.0, 556.1002, 0.0, 0.0, 0.0, 27.9449, 0.0, 0.0, 556.1002,
    27.5209, 0.0, 0.0, 549.4688, 0.0, 0.0, 0.0, 27.5209, 0.0, 0.0, 549.4688, 0.0,
    0.0, 0.0, 27.5209, 0.0, 0.0, 549.4688, 27.0985, 0.0, 0.0, 542.8384, 0.0, 0.0,
    0.0, 27.0985, 0.0, 0.0, 542.8384, 0.0, 0.0, 0.0, 27.0985, 0.0, 0.0, 542.8384,
    26.6778, 0.0, 0.0, 536.2091, 0.0, 0.0, 0.0, 26.6778, 0.0, 0.0, 536.2091, 0.0,
    0.0, 0.0, 26.6778, 0.0, 0.0, 536.2091, 26.2588, 0.0, 0.0, 529.5811, 0.0, 0.0,
    0.0, 26.2588, 0.0, 0.0, 529.5811, 0.0, 0.0, 0.0, 26.2588, 0.0, 0.0, 529.5811,
    25.8414, 0.0, 0.0, 522.9545, 0.0, 0.0, 0.0, 25.8414, 0.0, 0.0, 522.9545, 0.0,
    0.0, 0.0, 25.8414, 0.0, 0.0, 522.9545, 25.4258, 0.0, 0.0, 516.3294, 0.0, 0.0,
    0.0, 25.4258, 0.0, 0.0, 516.3294, 0.0, 0.0, 0.0, 25.4258, 0.0, 0.0, 516.3294,
    25.0117, 0.0, 0.0, 509.706, 0.0, 0.0, 0.0, 25.0117, 0.0, 0.0, 509.706, 0.0,
    0.0, 0.0, 25.0117, 0.0, 0.0, 509.706, 24.5994, 0.0, 0.0, 503.0844, 0.0, 0.0,
    0.0, 24.5994, 0.0, 0.0, 503.0844, 0.0, 0.0, 0.0, 24.5994, 0.0, 0.0, 503.0844,
    24.1887, 0.0, 0.0, 496.4647, 0.0, 0.0, 0.0, 24.1887, 0.0, 0.0, 496.4647, 0.0,
    0.0, 0.0, 24.1887, 0.0, 0.0, 496.4647, 23.7797, 0.0, 0.0, 489.8471, 0.0, 0.0,
    0.0, 23.7797, 0.0, 0.0, 489.8471, 0.0, 0.0, 0.0, 23.7797, 0.0, 0.0, 489.8471,
    23.3724, 0.0, 0.0, 483.2318, 0.0, 0.0, 0.0, 23.3724, 0.0, 0.0, 483.2318, 0.0,
    0.0, 0.0, 23.3724, 0.0, 0.0, 483.2318, 22.9667, 0.0, 0.0, 476.6187, 0.0, 0.0,
    0.0, 22.9667, 0.0, 0.0, 476.6187, 0.0, 0.0, 0.0, 22.9667, 0.0, 0.0, 476.6187,
    22.5627, 0.0, 0.0, 470.0081, 0.0, 0.0, 0.0, 22.5627, 0.0, 0.0, 470.0081, 0.0,
    0.0, 0.0, 22.5627, 0.0, 0.0, 470.0081 };

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

  double IGA[7200];
  double A_data[3600];
  double y_tmp[3600];
  double b_del_lam[240];
  double result[240];
  double varargin_1_data[240];
  double del_lam[120];
  double esig[120];
  double ftol[120];
  double igr2[120];
  double ilam[120];
  double lam[120];
  double mesil[120];
  double nlamt[120];
  double H[60];
  double del_z[60];
  double r1[60];
  double x[60];
  double d;
  double mu;
  double mu_old;
  double ssq;
int  b_i;
int  exitflag;
int  i;
int  idx;
int  iy;
int  unusedU0;
unsigned char ii_data[240];
bool b_x[240];

  //   QP Solver
  //  Call QP Solver
  //  Solve quadratic programming problem using Wright's (1997) Method
  //  Minimise J(x) = 1/2x'Hx + f'x
  //  Subject to: Ax <= b
  //  Reference: S. J. Wright, "Applying New Optimization Algorithms to Model
  //  Predictive Control," in Chemical Process Control-V, CACHE, AIChE
  //  Symposium, 1997, pp. 147-155.
  // Number of decision variables
  //  p = 0;
  // Test for Cold Start
  // Warm Start
  for (i = 0; i < 60; i++) {
    d = 0.0;
    for (b_i = 0; b_i < 6; b_i++) {
      d += (2.0 * x0[b_i]) * b[b_i + (6 * i)];
    }

    H[i] = d;
    x[i] = X_QP[i];
  }

  // to tune
  // to tune
  // Default Values
  for (i = 0; i < 120; i++) {
    lam[i] = 100.0;
    ftol[i] = 100.0;
    esig[i] = 0.001;
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
  int  i1;
  int  ii_size_idx_0;
  int  info;
  int  j;
  bool exitg1;
  bool y;

    // Create common matrices
    for (i = 0; i < 120; i++) {
      d = lam[i];
      ssq = 1.0 / d;
      ilam[i] = ssq;
      nlamt[i] = (-d) / ftol[i];
      mesil[i] = (mu * esig[i]) * ssq;
    }

    // RHS
    for (i = 0; i < 60; i++) {
      for (b_i = 0; b_i < 120; b_i++) {
        idx = b_i + (120 * i);
        IGA[idx] = nlamt[b_i] * (static_cast<double>(A[idx]));
      }

      d = 0.0;
      for (b_i = 0; b_i < 60; b_i++) {
        d += b_a[i + (60 * b_i)] * x[b_i];
      }

      ssq = 0.0;
      for (b_i = 0; b_i < 120; b_i++) {
        ssq += (static_cast<double>(c_a[i + (60 * b_i)])) * lam[b_i];
      }

      r1[i] = (d - ssq) - H[i];
    }

    for (b_i = 0; b_i < 120; b_i++) {
      d = 0.0;
      for (i1 = 0; i1 < 60; i1++) {
        d += (static_cast<double>(a[b_i + (120 * i1)])) * x[i1];
      }

      igr2[b_i] = nlamt[b_i] * ((d + 0.15) - mesil[b_i]);
    }

    // Solve
    for (b_i = 0; b_i < 60; b_i++) {
      for (i1 = 0; i1 < 60; i1++) {
        d = 0.0;
        for (i = 0; i < 120; i++) {
          d += (static_cast<double>(c_a[b_i + (60 * i)])) * IGA[i + (120 * i1)];
        }

        y_tmp[b_i + (60 * i1)] = d;
      }
    }

    for (b_i = 0; b_i < 3600; b_i++) {
      A_data[b_i] = b_H[b_i] - y_tmp[b_i];
    }

    info = -1;
    j = 0;
    exitg1 = false;
    while ((!exitg1) && (j < 60)) {
    int  idxA1j;
    int  idxAjj;
    int  ix;
      idxA1j = j * 60;
      idxAjj = idxA1j + j;
      ssq = 0.0;
      if (j >= 1) {
        ix = idxA1j;
        iy = idxA1j;
        for (i = 0; i < j; i++) {
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
        int  idxAjjp1;
          i = idxA1j + 61;
          idxAjjp1 = idxAjj + 61;
          if (j != 0) {
            iy = idxAjj + 60;
            b_i = (idxA1j + (60 * (58 - j))) + 61;
            for (idx = i; idx <= b_i; idx += 60) {
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
          b_i = (idxAjj + (60 * (58 - j))) + 61;
          for (i = idxAjjp1; i <= b_i; i += 60) {
            A_data[i - 1] *= ssq;
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
      for (b_i = 0; b_i < 60; b_i++) {
        d = 0.0;
        for (i1 = 0; i1 < 120; i1++) {
          d += (static_cast<double>(c_a[b_i + (60 * i1)])) * igr2[i1];
        }

        del_z[b_i] = r1[b_i] - d;
      }

      for (b_i = 0; b_i < 3600; b_i++) {
        y_tmp[b_i] = b_H[b_i] - y_tmp[b_i];
      }

      mldivide(y_tmp, del_z);

      // old method (LU?)
      //   del_z = linsolve (R, linsolve (R, (r1-At*igr2), opUT), opU); %exploit matrix properties for solving 
    } else {
      // Not Positive Definite (problem? eg infeasible)
      for (b_i = 0; b_i < 60; b_i++) {
        d = 0.0;
        for (i1 = 0; i1 < 120; i1++) {
          d += (static_cast<double>(c_a[b_i + (60 * i1)])) * igr2[i1];
        }

        del_z[b_i] = r1[b_i] - d;
      }

      for (b_i = 0; b_i < 3600; b_i++) {
        y_tmp[b_i] = b_H[b_i] - y_tmp[b_i];
      }

      mldivide(y_tmp, del_z);

      // old method (LU?)
    }

    // Decide on suitable alpha (from Wright's paper)
    // Try Max Increment (alpha = 1)
    // Check lam and ftol > 0
    for (i = 0; i < 120; i++) {
      d = 0.0;
      for (b_i = 0; b_i < 60; b_i++) {
        d += IGA[i + (120 * b_i)] * del_z[b_i];
      }

      d = igr2[i] - d;
      del_lam[i] = d;
      ssq = ftol[i];
      mu_old = ((-ssq) + mesil[i]) - ((ilam[i] * ssq) * d);
      mesil[i] = mu_old;
      d += lam[i];
      nlamt[i] = d;
      ssq += mu_old;
      ilam[i] = ssq;
      b_x[i] = (d < 2.2204460492503131E-16);
      b_x[i + 120] = (ssq < 2.2204460492503131E-16);
    }

    y = false;
    i = 0;
    exitg1 = false;
    while ((!exitg1) && (i < 240)) {
      if (!b_x[i]) {
        i++;
      } else {
        y = true;
        exitg1 = true;
      }
    }

    if (!y) {
      // KKT met
      std::memcpy(&lam[0], &nlamt[0], 120U * (sizeof(double)));
      std::memcpy(&ftol[0], &ilam[0], 120U * (sizeof(double)));
      for (b_i = 0; b_i < 60; b_i++) {
        x[b_i] += del_z[b_i];
      }
    } else {
      // KKT failed - solve by finding minimum ratio
      for (b_i = 0; b_i < 120; b_i++) {
        result[b_i] = nlamt[b_i];
        result[b_i + 120] = ilam[b_i];
      }

      idx = 0;
      i = 0;
      exitg1 = false;
      while ((!exitg1) && (i < 240)) {
        if (result[i] < 2.2204460492503131E-16) {
          idx++;
          ii_data[idx - 1] = static_cast<uint8_T>(static_cast<int32_T>(i + 1));
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
        ii_size_idx_0 = 0;
      } else {
        ii_size_idx_0 = idx;
      }

      // detects elements breaking KKT condition
      for (b_i = 0; b_i < 120; b_i++) {
        b_del_lam[b_i] = del_lam[b_i];
        b_del_lam[b_i + 120] = mesil[b_i];
      }

      for (b_i = 0; b_i < ii_size_idx_0; b_i++) {
        i = (static_cast<int32_T>(ii_data[b_i])) - 1;
        varargin_1_data[b_i] = 1.0 - (result[i] / b_del_lam[i]);
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
          i = 2;
          exitg1 = false;
          while ((!exitg1) && (i <= ii_size_idx_0)) {
            if (!rtIsNaN(varargin_1_data[i - 1])) {
              idx = i;
              exitg1 = true;
            } else {
              i++;
            }
          }
        }

        if (idx == 0) {
          ssq = varargin_1_data[0];
        } else {
          ssq = varargin_1_data[idx - 1];
          b_i = idx + 1;
          for (i = b_i; i <= ii_size_idx_0; i++) {
            d = varargin_1_data[i - 1];
            if (ssq > d) {
              ssq = d;
            }
          }
        }
      }

      ssq *= 0.995;

      // solves for min ratio (max value of alpha allowed)
      // Increment
      for (b_i = 0; b_i < 120; b_i++) {
        lam[b_i] += ssq * del_lam[b_i];
        ftol[b_i] += ssq * mesil[b_i];
      }

      for (b_i = 0; b_i < 60; b_i++) {
        x[b_i] += ssq * del_z[b_i];
      }
    }

    // Complimentary Gap
    mu_old = mu;
    ssq = 0.0;
    for (b_i = 0; b_i < 120; b_i++) {
      ssq += ftol[b_i] * lam[b_i];
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

      for (i = 0; i < 120; i++) {
        esig[i] = ssq;
      }
    }

    unusedU0++;
  }

  // Check for failure
  //  if(iter == maxiter)
  //      exitflag = -1;
  //  end
  Fx = x[0];
  Fy = x[1];
  Fz = x[2];

  for (i = 0; i < 60; i++) 
  {
    X_QP[i]=x[i];
  }
 
} 

template<typename T>
bool  CoordinatorBase<T>::rtIsNaN(double value){
  return ((value!=value) ? true : false);
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

template<typename T>
void CoordinatorBase<T>::LSS_primary(double u_primary[3], double Tau[3])
{
  double M_R[15];
  double y_tmp[15];
  double x[9];
  double y[9];
  double b_u_primary[5];
  double absx11;
  double absx21;
  double absx31;
  int p1;
  int p2;
  int p3;
  M_R[0] = 1.0;
  M_R[1] = 0.0;
  M_R[2] = 0.0;
  M_R[5] = 0.0;
  M_R[6] = 1.0;
  M_R[7] = 0.0;
  M_R[10] = 0.0;
  M_R[11] = 0.0;
  M_R[12] = 1.0;
  M_R[3] = 0.0;
  M_R[8] = 0.0;
  M_R[13] = -L;
  M_R[4] = 0.0;
  M_R[9] = L;
  M_R[14] = 0.0;
  for (p1 = 0; p1 < 5; p1++) {
    y_tmp[3 * p1] = M_R[p1];
    y_tmp[(3 * p1) + 1] = M_R[p1 + 5];
    y_tmp[(3 * p1) + 2] = M_R[p1 + 10];
  }

  for (p1 = 0; p1 < 3; p1++) {
    for (p2 = 0; p2 < 3; p2++) {
      absx11 = 0.0;
      for (p3 = 0; p3 < 5; p3++) {
        absx11 += y_tmp[p1 + (3 * p3)] * M_R[p3 + (5 * p2)];
      }

      y[p1 + (3 * p2)] = absx11;
    }
  }

  std::memcpy(&x[0], &y[0], 9U * (sizeof(double)));
  p1 = 0;
  p2 = 3;
  p3 = 6;
  absx11 = std::abs(y[0]);
  absx21 = std::abs(y[1]);
  absx31 = std::abs(y[2]);
  if ((absx21 > absx11) && (absx21 > absx31)) {
    p1 = 3;
    p2 = 0;
    x[0] = y[1];
    x[1] = y[0];
    x[3] = y[4];
    x[4] = y[3];
    x[6] = y[7];
    x[7] = y[6];
  } else {
    if (absx31 > absx11) {
      p1 = 6;
      p3 = 0;
      x[0] = y[2];
      x[2] = y[0];
      x[3] = y[5];
      x[5] = y[3];
      x[6] = y[8];
      x[8] = y[6];
    }
  }

  x[1] /= x[0];
  x[2] /= x[0];
  x[4] -= x[1] * x[3];
  x[5] -= x[2] * x[3];
  x[7] -= x[1] * x[6];
  x[8] -= x[2] * x[6];
  if (std::abs(x[5]) > std::abs(x[4])) {
    int itmp;
    itmp = p2;
    p2 = p3;
    p3 = itmp;
    absx11 = x[1];
    x[1] = x[2];
    x[2] = absx11;
    absx11 = x[4];
    x[4] = x[5];
    x[5] = absx11;
    absx11 = x[7];
    x[7] = x[8];
    x[8] = absx11;
  }

  x[5] /= x[4];
  x[8] -= x[5] * x[7];
  absx11 = ((x[5] * x[1]) - x[2]) / x[8];
  absx21 = (-(x[1] + (x[7] * absx11))) / x[4];
  y[p1] = ((1.0 - (x[3] * absx21)) - (x[6] * absx11)) / x[0];
  y[p1 + 1] = absx21;
  y[p1 + 2] = absx11;
  absx11 = (-x[5]) / x[8];
  absx21 = (1.0 - (x[7] * absx11)) / x[4];
  y[p2] = (-((x[3] * absx21) + (x[6] * absx11))) / x[0];
  y[p2 + 1] = absx21;
  y[p2 + 2] = absx11;
  absx11 = 1.0 / x[8];
  absx21 = ((-x[7]) * absx11) / x[4];
  y[p3] = (-((x[3] * absx21) + (x[6] * absx11))) / x[0];
  y[p3 + 1] = absx21;
  y[p3 + 2] = absx11;
  b_u_primary[0] = u_primary[0];
  b_u_primary[1] = u_primary[1];
  b_u_primary[2] = u_primary[2];
  b_u_primary[3] = Tau[1] / 2.0;
  b_u_primary[4] = Tau[2] / 2.0;
  for (p1 = 0; p1 < 3; p1++) {
    double d;
    absx11 = 0.0;
    absx21 = y[p1];
    absx31 = y[p1 + 3];
    d = y[p1 + 6];
    for (p2 = 0; p2 < 5; p2++) {
      absx11 += (((absx21 * y_tmp[3 * p2]) + (absx31 * y_tmp[(3 * p2) + 1])) +
                 (d * y_tmp[(3 * p2) + 2])) * b_u_primary[p2];
    }

    F_r[p1] = absx11;
  }
}


template<typename T>
void CoordinatorBase<T>::LSS_secondary(double u_secondary[3], double Tau[3])
{
  double M_L[15];
  double y_tmp[15];
  double x[9];
  double y[9];
  double b_u_secondary[5];
  double absx11;
  double absx21;
  double absx31;
  int p1;
  int p2;
  int p3;
  M_L[0] = -1.0;
  M_L[1] = 0.0;
  M_L[2] = 0.0;
  M_L[5] = 0.0;
  M_L[6] = -1.0;
  M_L[7] = 0.0;
  M_L[10] = 0.0;
  M_L[11] = 0.0;
  M_L[12] = 1.0;
  M_L[3] = 0.0;
  M_L[8] = 0.0;
  M_L[13] = L;
  M_L[4] = 0.0;
  M_L[9] = L;
  M_L[14] = 0.0;
  for (p1 = 0; p1 < 5; p1++) {
    y_tmp[3 * p1] = M_L[p1];
    y_tmp[(3 * p1) + 1] = M_L[p1 + 5];
    y_tmp[(3 * p1) + 2] = M_L[p1 + 10];
  }

  for (p1 = 0; p1 < 3; p1++) {
    for (p2 = 0; p2 < 3; p2++) {
      absx11 = 0.0;
      for (p3 = 0; p3 < 5; p3++) {
        absx11 += y_tmp[p1 + (3 * p3)] * M_L[p3 + (5 * p2)];
      }

      y[p1 + (3 * p2)] = absx11;
    }
  }

  std::memcpy(&x[0], &y[0], 9U * (sizeof(double)));
  p1 = 0;
  p2 = 3;
  p3 = 6;
  absx11 = std::abs(y[0]);
  absx21 = std::abs(y[1]);
  absx31 = std::abs(y[2]);
  if ((absx21 > absx11) && (absx21 > absx31)) {
    p1 = 3;
    p2 = 0;
    x[0] = y[1];
    x[1] = y[0];
    x[3] = y[4];
    x[4] = y[3];
    x[6] = y[7];
    x[7] = y[6];
  } else {
    if (absx31 > absx11) {
      p1 = 6;
      p3 = 0;
      x[0] = y[2];
      x[2] = y[0];
      x[3] = y[5];
      x[5] = y[3];
      x[6] = y[8];
      x[8] = y[6];
    }
  }

  x[1] /= x[0];
  x[2] /= x[0];
  x[4] -= x[1] * x[3];
  x[5] -= x[2] * x[3];
  x[7] -= x[1] * x[6];
  x[8] -= x[2] * x[6];
  if (std::abs(x[5]) > std::abs(x[4])) {
    int itmp;
    itmp = p2;
    p2 = p3;
    p3 = itmp;
    absx11 = x[1];
    x[1] = x[2];
    x[2] = absx11;
    absx11 = x[4];
    x[4] = x[5];
    x[5] = absx11;
    absx11 = x[7];
    x[7] = x[8];
    x[8] = absx11;
  }

  x[5] /= x[4];
  x[8] -= x[5] * x[7];
  absx11 = ((x[5] * x[1]) - x[2]) / x[8];
  absx21 = (-(x[1] + (x[7] * absx11))) / x[4];
  y[p1] = ((1.0 - (x[3] * absx21)) - (x[6] * absx11)) / x[0];
  y[p1 + 1] = absx21;
  y[p1 + 2] = absx11;
  absx11 = (-x[5]) / x[8];
  absx21 = (1.0 - (x[7] * absx11)) / x[4];
  y[p2] = (-((x[3] * absx21) + (x[6] * absx11))) / x[0];
  y[p2 + 1] = absx21;
  y[p2 + 2] = absx11;
  absx11 = 1.0 / x[8];
  absx21 = ((-x[7]) * absx11) / x[4];
  y[p3] = (-((x[3] * absx21) + (x[6] * absx11))) / x[0];
  y[p3 + 1] = absx21;
  y[p3 + 2] = absx11;
  b_u_secondary[0] = u_secondary[0];
  b_u_secondary[1] = u_secondary[1];
  b_u_secondary[2] = u_secondary[2];
  b_u_secondary[3] = Tau[1] / 2.0;
  b_u_secondary[4] = Tau[2] / 2.0;
  for (p1 = 0; p1 < 3; p1++) {
    double d;
    absx11 = 0.0;
    absx21 = y[p1];
    absx31 = y[p1 + 3];
    d = y[p1 + 6];
    for (p2 = 0; p2 < 5; p2++) {
      absx11 += (((absx21 * y_tmp[3 * p2]) + (absx31 * y_tmp[(3 * p2) + 1])) +
                 (d * y_tmp[(3 * p2) + 2])) * b_u_secondary[p2];
    }

    F_l[p1] = absx11;
  }
}
