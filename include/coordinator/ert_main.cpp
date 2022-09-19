//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: ert_main.cpp
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
#include <stddef.h>
#include <stdio.h>                // This ert_main.c example uses printf/fflush
#include "controller1.h"               // Model's header file
#include "rtwtypes.h"

static controller1ModelClass controller1_Obj;// Instance of model class

// '<Root>/x_e'
static real_T arg_x_e = 0.0;

// '<Root>/y_e'
static real_T arg_y_e = 0.0;

// '<Root>/z_e'
static real_T arg_z_e = 0.0;

// '<Root>/vx'
static real_T arg_vx = 0.0;

// '<Root>/vy'
static real_T arg_vy = 0.0;

// '<Root>/vz'
static real_T arg_vz = 0.0;

// '<Root>/qx'
static real_T arg_qx = 0.0;

// '<Root>/qy'
static real_T arg_qy = 0.0;

// '<Root>/qz'
static real_T arg_qz = 0.0;

// '<Root>/qw'
static real_T arg_qw = 0.0;

// '<Root>/omegax'
static real_T arg_omegax = 0.0;

// '<Root>/omegay'
static real_T arg_omegay = 0.0;

// '<Root>/omegaz'
static real_T arg_omegaz = 0.0;

// '<Root>/fx'
static real_T arg_fx;

// '<Root>/fy'
static real_T arg_fy;

// '<Root>/fz'
static real_T arg_fz;

// '<Root>/tau_x'
static real_T arg_tau_x;

// '<Root>/tau_y'
static real_T arg_tau_y;

// '<Root>/tau_z'
static real_T arg_tau_z;

//
// Associating rt_OneStep with a real-time clock or interrupt service routine
// is what makes the generated code "real-time".  The function rt_OneStep is
// always associated with the base rate of the model.  Subrates are managed
// by the base rate from inside the generated code.  Enabling/disabling
// interrupts and floating point context switches are target specific.  This
// example code indicates where these should take place relative to executing
// the generated code step function.  Overrun behavior should be tailored to
// your application needs.  This example simply sets an error status in the
// real-time model and returns from rt_OneStep.
//
void rt_OneStep(void);
void rt_OneStep(void)
{
  static boolean_T OverrunFlag = false;

  // Disable interrupts here

  // Check for overrun
  if (OverrunFlag) {
    rtmSetErrorStatus(controller1_Obj.getRTM(), "Overrun");
    return;
  }

  OverrunFlag = true;

  // Save FPU context here (if necessary)
  // Re-enable timer or interrupt here
  // Set model inputs here

  // Step the model
  controller1_Obj.step(arg_x_e, arg_y_e, arg_z_e, arg_vx, arg_vy, arg_vz, arg_qx,
                       arg_qy, arg_qz, arg_qw, arg_omegax, arg_omegay,
                       arg_omegaz, &arg_fx, &arg_fy, &arg_fz, &arg_tau_x,
                       &arg_tau_y, &arg_tau_z);

  // Get model outputs here

  // Indicate task complete
  OverrunFlag = false;

  // Disable interrupts here
  // Restore FPU context here (if necessary)
  // Enable interrupts here
}

//
// The example "main" function illustrates what is required by your
// application code to initialize, execute, and terminate the generated code.
// Attaching rt_OneStep to a real-time clock is target specific.  This example
// illustrates how you do this relative to initializing the model.
//
int_T main(int_T argc, const char *argv[])
{
  // Unused arguments
  (void)(argc);
  (void)(argv);

  // Initialize model
  controller1_Obj.initialize();

  // Attach rt_OneStep to a timer or interrupt service routine with
  //  period 0.016 seconds (the model's base sample time) here.  The
  //  call syntax for rt_OneStep is
  //
  //   rt_OneStep();

  printf("Warning: The simulation will run forever. "
         "Generated ERT main won't simulate model step behavior. "
         "To change this behavior select the 'MAT-file logging' option.\n");
  fflush((NULL));
  while (rtmGetErrorStatus(controller1_Obj.getRTM()) == (NULL)) {
    //  Perform other application tasks here
  }

  // Disable rt_OneStep() here
  return 0;
}

//
// File trailer for generated code.
//
// [EOF]
//
