//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: controller1.cpp
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
#include "controller1.h"

//===========*
//  Constants *
// ===========
#define RT_PI                          3.14159265358979323846
#define RT_PIF                         3.1415927F
#define RT_LN_10                       2.30258509299404568402
#define RT_LN_10F                      2.3025851F
#define RT_LOG10E                      0.43429448190325182765
#define RT_LOG10EF                     0.43429449F
#define RT_E                           2.7182818284590452354
#define RT_EF                          2.7182817F

//
//  UNUSED_PARAMETER(x)
//    Used to specify that a function parameter (argument) is required but not
//    accessed by the function body.

#ifndef UNUSED_PARAMETER
#if defined(__LCC__)
#define UNUSED_PARAMETER(x)                                      // do nothing
#else

//
//  This is the semi-ANSI standard way of indicating that an
//  unused function parameter is required.

#define UNUSED_PARAMETER(x)            (void) (x)
#endif
#endif

// Model step function
void controller1ModelClass::step(real_T arg_x_e, real_T arg_y_e, real_T arg_z_e,
  real_T arg_vx, real_T arg_vy, real_T arg_vz, real_T arg_qx, real_T arg_qy,
  real_T arg_qz, real_T arg_qw, real_T arg_omegax, real_T arg_omegay, real_T
  arg_omegaz, real_T *arg_fx, real_T *arg_fy, real_T *arg_fz, real_T *arg_tau_x,
  real_T *arg_tau_y, real_T *arg_tau_z)
{
  static const real_T a[18] = { -4.0, -0.0, -0.0, -0.0, -4.0, -0.0, -0.0, -0.0,
    -4.0, -2.828, -0.0, -0.0, -0.0, -2.828, -0.0, -0.0, -0.0, -2.828 };

  static const real_T a_0[18] = { -0.612, -0.0, -0.0, -0.0, -0.572, -0.0, -0.0,
    -0.0, -0.652, -0.43268399999999996, -0.0, -0.0, -0.0, -0.40440399999999993,
    -0.0, -0.0, -0.0, -0.460964 };

  real_T tmp[6];
  real_T u[3];
  int32_T i;
  int32_T i_0;
  UNUSED_PARAMETER(arg_qw);

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

  // Outport: '<Root>/fx' incorporates:
  //   MATLAB Function: '<Root>/Position Controller'

  *arg_fx = u[0];

  // Outport: '<Root>/fy' incorporates:
  //   MATLAB Function: '<Root>/Position Controller'

  *arg_fy = u[1];

  // Outport: '<Root>/fz' incorporates:
  //   MATLAB Function: '<Root>/Position Controller'

  *arg_fz = u[2];

  // MATLAB Function: '<Root>/Attitude Controller' incorporates:
  //   Inport: '<Root>/omegax'
  //   Inport: '<Root>/omegay'
  //   Inport: '<Root>/omegaz'
  //   Inport: '<Root>/qx'
  //   Inport: '<Root>/qy'
  //   Inport: '<Root>/qz'

  tmp[0] = arg_qx;
  tmp[1] = arg_qy;
  tmp[2] = arg_qz;
  tmp[3] = arg_omegax;
  tmp[4] = arg_omegay;
  tmp[5] = arg_omegaz;
  for (i = 0; i < 3; i++) {
    u[i] = 0.0;
    for (i_0 = 0; i_0 < 6; i_0++) {
      u[i] += a_0[3 * i_0 + i] * tmp[i_0];
    }
  }

  // Outport: '<Root>/tau_x' incorporates:
  //   MATLAB Function: '<Root>/Attitude Controller'

  *arg_tau_x = u[0];

  // Outport: '<Root>/tau_y' incorporates:
  //   MATLAB Function: '<Root>/Attitude Controller'

  *arg_tau_y = u[1];

  // Outport: '<Root>/tau_z' incorporates:
  //   MATLAB Function: '<Root>/Attitude Controller'

  *arg_tau_z = u[2];
}

// Model initialize function
void controller1ModelClass::initialize()
{
  // (no initialization code required)
}

// Constructor
controller1ModelClass::controller1ModelClass() :
  controller1_M()
{
  // Currently there is no constructor body generated.
}

// Destructor
controller1ModelClass::~controller1ModelClass()
{
  // Currently there is no destructor body generated.
}

// Real-Time Model get method
controller1ModelClass::RT_MODEL_controller1_T * controller1ModelClass::getRTM()
{
  return (&controller1_M);
}

//
// File trailer for generated code.
//
// [EOF]
//
