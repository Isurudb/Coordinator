//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: minOrMax.cpp
//
// MATLAB Coder version            : 5.1
// C/C++ source code generated on  : 10-Oct-2022 09:46:31
//

// Include Files
#include "minOrMax.h"
#include "rt_nonfinite.h"
#include "rt_nonfinite.h"

// Function Definitions
//
// Arguments    : const double x_data[]
//                const int x_size[1]
// Return Type  : double
//
namespace coder
{
  namespace internal
  {
    double minimum(const double x_data[], const int x_size[1])
    {
      double ex;
      int n;
      n = x_size[0];
      if (x_size[0] <= 2) {
        if (x_size[0] == 1) {
          ex = x_data[0];
        } else if ((x_data[0] > x_data[1]) || (rtIsNaN(x_data[0]) && (!rtIsNaN
                     (x_data[1])))) {
          ex = x_data[1];
        } else {
          ex = x_data[0];
        }
      } else {
        int idx;
        int k;
        if (!rtIsNaN(x_data[0])) {
          idx = 1;
        } else {
          boolean_T exitg1;
          idx = 0;
          k = 2;
          exitg1 = false;
          while ((!exitg1) && (k <= x_size[0])) {
            if (!rtIsNaN(x_data[k - 1])) {
              idx = k;
              exitg1 = true;
            } else {
              k++;
            }
          }
        }

        if (idx == 0) {
          ex = x_data[0];
        } else {
          ex = x_data[idx - 1];
          idx++;
          for (k = idx; k <= n; k++) {
            double d;
            d = x_data[k - 1];
            if (ex > d) {
              ex = d;
            }
          }
        }
      }

      return ex;
    }
  }
}

//
// File trailer for minOrMax.cpp
//
// [EOF]
//
