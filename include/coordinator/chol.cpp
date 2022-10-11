//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: chol.cpp
//
// MATLAB Coder version            : 5.1
// C/C++ source code generated on  : 10-Oct-2022 09:46:31
//

// Include Files
#include "chol.h"
#include "rt_nonfinite.h"
#include <cmath>
#include <cstring>

// Function Definitions
//
// Arguments    : double A_data[]
//                int A_size[2]
// Return Type  : int
//
namespace coder
{
  int cholesky(double A_data[], int A_size[2])
  {
    double b_A_data[3600];
    double c;
    double ssq;
    int b_info;
    int i;
    int i1;
    int iac;
    int idxA1j;
    int info;
    int iy;
    int j;
    int loop_ub;
    boolean_T exitg1;
    loop_ub = A_size[0] * A_size[1];
    if (0 <= loop_ub - 1) {
      std::memcpy(&b_A_data[0], &A_data[0], loop_ub * sizeof(double));
    }

    b_info = -1;
    j = 0;
    exitg1 = false;
    while ((!exitg1) && (j < 60)) {
      int idxAjj;
      int ix;
      idxA1j = j * 60;
      idxAjj = idxA1j + j;
      ssq = 0.0;
      if (j >= 1) {
        ix = idxA1j;
        iy = idxA1j;
        for (loop_ub = 0; loop_ub < j; loop_ub++) {
          ssq += b_A_data[ix] * b_A_data[iy];
          ix++;
          iy++;
        }
      }

      ssq = b_A_data[idxAjj] - ssq;
      if (ssq > 0.0) {
        ssq = std::sqrt(ssq);
        b_A_data[idxAjj] = ssq;
        if (j + 1 < 60) {
          int idxAjjp1;
          loop_ub = idxA1j + 61;
          idxAjjp1 = idxAjj + 61;
          if (j != 0) {
            iy = idxAjj + 60;
            i = (idxA1j + 60 * (58 - j)) + 61;
            for (iac = loop_ub; iac <= i; iac += 60) {
              ix = idxA1j;
              c = 0.0;
              i1 = (iac + j) - 1;
              for (int ia = iac; ia <= i1; ia++) {
                c += b_A_data[ia - 1] * b_A_data[ix];
                ix++;
              }

              b_A_data[iy] += -c;
              iy += 60;
            }
          }

          ssq = 1.0 / ssq;
          i = (idxAjj + 60 * (58 - j)) + 61;
          for (loop_ub = idxAjjp1; loop_ub <= i; loop_ub += 60) {
            b_A_data[loop_ub - 1] *= ssq;
          }
        }

        j++;
      } else {
        b_A_data[idxAjj] = ssq;
        b_info = j;
        exitg1 = true;
      }
    }

    A_size[0] = 60;
    A_size[1] = 60;
    for (i = 0; i < 60; i++) {
      for (i1 = 0; i1 < 60; i1++) {
        A_data[i1 + A_size[0] * i] = b_A_data[i1 + 60 * i];
      }
    }

    info = b_info + 1;
    if (b_info + 1 == 0) {
      idxA1j = 60;
    } else {
      idxA1j = b_info;
    }

    for (j = 0; j < idxA1j; j++) {
      i = j + 2;
      for (loop_ub = i; loop_ub <= idxA1j; loop_ub++) {
        A_data[(loop_ub + A_size[0] * j) - 1] = 0.0;
      }
    }

    if (1 > idxA1j) {
      loop_ub = 0;
      idxA1j = 0;
    } else {
      loop_ub = idxA1j;
    }

    for (i = 0; i < idxA1j; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_A_data[i1 + loop_ub * i] = A_data[i1 + A_size[0] * i];
      }
    }

    A_size[0] = loop_ub;
    A_size[1] = idxA1j;
    for (i = 0; i < idxA1j; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        A_data[i1 + A_size[0] * i] = b_A_data[i1 + loop_ub * i];
      }
    }

    return info;
  }
}

//
// File trailer for chol.cpp
//
// [EOF]
//
