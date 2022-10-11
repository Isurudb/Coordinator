//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: find.cpp
//
// MATLAB Coder version            : 5.1
// C/C++ source code generated on  : 10-Oct-2022 09:46:31
//

// Include Files
#include "find.h"
#include "rt_nonfinite.h"

// Function Definitions
//
// Arguments    : const boolean_T x[240]
//                int i_data[]
//                int i_size[1]
// Return Type  : void
//
namespace coder
{
  void eml_find(const boolean_T x[240], int i_data[], int i_size[1])
  {
    int idx;
    int ii;
    boolean_T exitg1;
    idx = 0;
    ii = 0;
    exitg1 = false;
    while ((!exitg1) && (ii < 240)) {
      if (x[ii]) {
        idx++;
        i_data[idx - 1] = ii + 1;
        if (idx >= 240) {
          exitg1 = true;
        } else {
          ii++;
        }
      } else {
        ii++;
      }
    }

    if (1 > idx) {
      i_size[0] = 0;
    } else {
      i_size[0] = idx;
    }
  }
}

//
// File trailer for find.cpp
//
// [EOF]
//
