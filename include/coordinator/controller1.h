//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: controller1.h
//
// Code generated for Simulink model 'controller1'.
//
// Model version                  : 1.13
// Simulink Coder version         : 9.4 (R2020b) 29-Jul-2020
// C/C++ source code generated on : Tue Sep 13 20:07:33 2022
//
// Target selection: ert.tlc
// Embedded hardware selection: Intel->x86-64 (Linux 64)
// Code generation objective: Debugging
// Validation result: Not run
//
#ifndef RTW_HEADER_controller1_h_
#define RTW_HEADER_controller1_h_
#include <stddef.h>
#include "rtw_modelmap.h"
#include "rtwtypes.h"
#include <stddef.h>

// Model Code Variants

// Macros for accessing real-time model data structure
#ifndef rtmGetDataMapInfo
#define rtmGetDataMapInfo(rtm)         ((rtm)->DataMapInfo)
#endif

#ifndef rtmSetDataMapInfo
#define rtmSetDataMapInfo(rtm, val)    ((rtm)->DataMapInfo = (val))
#endif

#ifndef rtmGetErrorStatus
#define rtmGetErrorStatus(rtm)         ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
#define rtmSetErrorStatus(rtm, val)    ((rtm)->errorStatus = (val))
#endif

// Function to get C API Model Mapping Static Info
extern const rtwCAPI_ModelMappingStaticInfo*
  controller1_GetCAPIStaticMap(void);

// Class declaration for model controller1
class controller1ModelClass {
  // public data and function members
 public:
  // Real-time Model Data Structure
  struct RT_MODEL {
    const char_T * volatile errorStatus;

    //
    //  DataMapInfo:
    //  The following substructure contains information regarding
    //  structures generated in the model's C API.

    struct {
      rtwCAPI_ModelMappingInfo mmi;
    } DataMapInfo;
  };

  // model initialize function
  void initialize();

  // model step function
  void step(real_T &arg_x_e, real_T &arg_y_e, real_T &arg_z_e, real_T &arg_vx,
            real_T &arg_vy, real_T &arg_vz, real_T &arg_qx, real_T &arg_qy,
            real_T &arg_qz, real_T &arg_qw, real_T &arg_omegax, real_T
            &arg_omegay, real_T *arg_omegaz, real_T *arg_fx, real_T *arg_fy,
            real_T *arg_fz, real_T *arg_tau_x, real_T *arg_tau_y, real_T
            *arg_tau_z);

  // Constructor
  controller1ModelClass();

  // Destructor
  ~controller1ModelClass();

  // Real-Time Model get method
  controller1ModelClass::RT_MODEL * getRTM();

  // private data and function members
 private:
  // Real-Time Model
  RT_MODEL rtM;
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


//-
//  Requirements for '<Root>': controller1

#endif                                 // RTW_HEADER_controller1_h_

//
// File trailer for generated code.
//
// [EOF]
//
