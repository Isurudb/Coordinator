//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: MPC_Guidance_v3_sand.h
//
// MATLAB Coder version            : 5.1
// C/C++ source code generated on  : 26-Sep-2022 17:35:01
//
#ifndef MPC_GUIDANCE_V3_SAND_H
#define MPC_GUIDANCE_V3_SAND_H

// Include Files
#include "MPC_Guidance_v3_sand_types.h"
#include "rtwtypes.h"
#include <cstddef>
#include <cstdlib>

// Function Declarations
extern void MPC_Guidance_v3_sand(const struct0_T *MPCParams, const double x0[6],
  const double xDeputy[6], double xfinal[6], const double X0[60], double
  prev_pt_sel, double *Fx, double *Fy, double *Fz, double target_state[6],
  double *ObstAvoid, double *AppCone, double *dock_flag, double *CollAvoid_flag,
  double *dock_complete, double *num_iter, double X_QP[60], double *pt_sel,
  double *dr);

#endif

//
// File trailer for MPC_Guidance_v3_sand.h
//
// [EOF]
//
