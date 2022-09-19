//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: controller1.h
//
// Code generated for Simulink model 'controller1'.
//
// Model version                  : 1.28
// Simulink Coder version         : 9.4 (R2020b) 29-Jul-2020
// C/C++ source code generated on : Sat Sep 17 17:31:41 2022
//
// Target selection: ert.tlc
// Embedded hardware selection: Intel->x86-64 (Windows64)
// Code generation objectives:
//    1. Execution efficiency
//    2. RAM efficiency
// Validation result: Not run
//
#ifndef RTW_HEADER_controller1_h_
#define RTW_HEADER_controller1_h_
#include "rtwtypes.h"

// Model Code Variants

// Macros for accessing real-time model data structure
#ifndef rtmGetErrorStatus
#define rtmGetErrorStatus(rtm)         ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
#define rtmSetErrorStatus(rtm, val)    ((rtm)->errorStatus = (val))
#endif

// Class declaration for model controller1
class controller1ModelClass {
  // public data and function members
 public:
  // Real-time Model Data Structure
  struct RT_MODEL_controller1_T {
    const char_T * volatile errorStatus;
  };

  // model initialize function
  void initialize();

  // model step function
  void step(real_T arg_x_e, real_T arg_y_e, real_T arg_z_e, real_T arg_vx,
            real_T arg_vy, real_T arg_vz, real_T arg_qx, real_T arg_qy, real_T
            arg_qz, real_T arg_qw, real_T arg_omegax, real_T arg_omegay, real_T
            arg_omegaz, real_T *arg_fx, real_T *arg_fy, real_T *arg_fz, real_T
            *arg_tau_x, real_T *arg_tau_y, real_T *arg_tau_z);

  // Constructor
  controller1ModelClass();

  // Destructor
  ~controller1ModelClass();

  // Real-Time Model get method
  controller1ModelClass::RT_MODEL_controller1_T * getRTM();

  // private data and function members
 private:
  // Real-Time Model
  RT_MODEL_controller1_T controller1_M;
};

//-
//  The generated code includes comments that allow you to trace directly
//  back to the appropriate location in the model.  The basic format
//  is <system>/block_name, where system is the system number (uniquely
//  assigned by Simulink) and block_name is the name of the block.
//
//  Use the MATLAB hilite_system command to trace the generated code back
//  to the model.  For example,
//
//  hilite_system('<S3>')    - opens system 3
//  hilite_system('<S3>/Kp') - opens and selects block Kp which resides in S3
//
//  Here is the system hierarchy for this model
//
//  '<Root>' : 'controller1'
//  '<S1>'   : 'controller1/Attitude Controller'
//  '<S2>'   : 'controller1/Position Controller'

#endif                                 // RTW_HEADER_controller1_h_

//
// File trailer for generated code.
//
// [EOF]
//
