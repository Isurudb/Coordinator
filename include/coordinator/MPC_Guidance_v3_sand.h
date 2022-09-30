//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: MPC_Guidance_v3_sand.h
//
// MATLAB Coder version            : 5.1
// C/C++ source code generated on  : 29-Sep-2022 00:21:21
//
#ifndef MPC_GUIDANCE_V3_SAND_H
#define MPC_GUIDANCE_V3_SAND_H

// Include Files
#include "rtwtypes.h"
#include <cstddef>
#include <cstdlib>

// Function Declarations
extern void MPC_Guidance_v3_sand( double x0[6], double Fx, double Fy,
  double Fz, double target_state[6], double dock_flag, double CollAvoid_flag,
  double dock_complete, double num_iter, double X_QP[60], double pt_sel,
  double dr);

#endif

//
// File trailer for MPC_Guidance_v3_sand.h
//
// [EOF]
//
