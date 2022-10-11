//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: inv.cpp
//
// MATLAB Coder version            : 5.1
// C/C++ source code generated on  : 10-Oct-2022 09:46:31
//

// Include Files
#include "inv.h"
#include "rt_nonfinite.h"
#include <cmath>

// Function Definitions
//
// Arguments    : const double x[36]
//                double y[36]
// Return Type  : void
//
namespace coder
{
  void inv(const double x[36], double y[36])
  {
    double b_x[36];
    double smax;
    int b_i;
    int i;
    int i2;
    int ix;
    int iy;
    int j;
    int jA;
    int jp1j;
    int k;
    signed char ipiv[6];
    signed char p[6];
    for (i = 0; i < 36; i++) {
      y[i] = 0.0;
      b_x[i] = x[i];
    }

    for (i = 0; i < 6; i++) {
      ipiv[i] = static_cast<signed char>(i + 1);
    }

    for (j = 0; j < 5; j++) {
      int b_tmp;
      int mmj_tmp;
      mmj_tmp = 4 - j;
      b_tmp = j * 7;
      jp1j = b_tmp + 2;
      iy = 6 - j;
      jA = 0;
      ix = b_tmp;
      smax = std::abs(b_x[b_tmp]);
      for (k = 2; k <= iy; k++) {
        double s;
        ix++;
        s = std::abs(b_x[ix]);
        if (s > smax) {
          jA = k - 1;
          smax = s;
        }
      }

      if (b_x[b_tmp + jA] != 0.0) {
        if (jA != 0) {
          iy = j + jA;
          ipiv[j] = static_cast<signed char>(iy + 1);
          ix = j;
          for (k = 0; k < 6; k++) {
            smax = b_x[ix];
            b_x[ix] = b_x[iy];
            b_x[iy] = smax;
            ix += 6;
            iy += 6;
          }
        }

        i = (b_tmp - j) + 6;
        for (b_i = jp1j; b_i <= i; b_i++) {
          b_x[b_i - 1] /= b_x[b_tmp];
        }
      }

      iy = b_tmp + 6;
      jA = b_tmp;
      for (jp1j = 0; jp1j <= mmj_tmp; jp1j++) {
        smax = b_x[iy];
        if (b_x[iy] != 0.0) {
          ix = b_tmp + 1;
          i = jA + 8;
          i2 = (jA - j) + 12;
          for (b_i = i; b_i <= i2; b_i++) {
            b_x[b_i - 1] += b_x[ix] * -smax;
            ix++;
          }
        }

        iy += 6;
        jA += 6;
      }
    }

    for (i = 0; i < 6; i++) {
      p[i] = static_cast<signed char>(i + 1);
    }

    for (k = 0; k < 5; k++) {
      signed char i1;
      i1 = ipiv[k];
      if (i1 > k + 1) {
        iy = p[i1 - 1];
        p[i1 - 1] = p[k];
        p[k] = static_cast<signed char>(iy);
      }
    }

    for (k = 0; k < 6; k++) {
      jp1j = 6 * (p[k] - 1);
      y[k + jp1j] = 1.0;
      for (j = k + 1; j < 7; j++) {
        i = (j + jp1j) - 1;
        if (y[i] != 0.0) {
          i2 = j + 1;
          for (b_i = i2; b_i < 7; b_i++) {
            iy = (b_i + jp1j) - 1;
            y[iy] -= y[i] * b_x[(b_i + 6 * (j - 1)) - 1];
          }
        }
      }
    }

    for (j = 0; j < 6; j++) {
      iy = 6 * j;
      for (k = 5; k >= 0; k--) {
        jA = 6 * k;
        i = k + iy;
        smax = y[i];
        if (smax != 0.0) {
          y[i] = smax / b_x[k + jA];
          for (b_i = 0; b_i < k; b_i++) {
            jp1j = b_i + iy;
            y[jp1j] -= y[i] * b_x[b_i + jA];
          }
        }
      }
    }
  }
}

//
// File trailer for inv.cpp
//
// [EOF]
//
