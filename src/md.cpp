#include <cpp11.hpp>
#include "mdiag.h"
using namespace cpp11;

[[cpp11::register]]
data_frame C_murphydiag_expect(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const double& level) {
  md_expect md{ x, y, level };
  murphydiag mdiag{ md };
  return cpp11::writable::data_frame({
    "knot"_nm = mdiag.knots,
    "val_left"_nm = mdiag.vals_l,
    "val_right"_nm = mdiag.vals_r
  });
}


[[cpp11::register]]
data_frame C_murphydiag_median(
    const std::vector<double>& x,
    const std::vector<double>& y) {
  md_median md{ x, y };
  murphydiag mdiag{ md };
  return cpp11::writable::data_frame({
    "knot"_nm = mdiag.knots,
    "val_left"_nm = mdiag.vals_l,
    "val_right"_nm = mdiag.vals_r
  });
}

[[cpp11::register]]
data_frame C_murphydiag_prob(
    const std::vector<double>& x,
    const std::vector<double>& y) {
  md_prob md{ x, y };
  murphydiag mdiag{ md };
  return cpp11::writable::data_frame({
    "knot"_nm = mdiag.knots,
    "val_left"_nm = mdiag.vals_l,
    "val_right"_nm = mdiag.vals_r
  });
}

[[cpp11::register]]
data_frame C_murphydiag_prob_ref(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::vector<double>& ref) {
  md_prob md{ x, y, ref };
  murphydiag mdiag{ md };
  return cpp11::writable::data_frame({
    "knot"_nm = mdiag.knots,
      "val_left"_nm = mdiag.vals_l,
      "val_right"_nm = mdiag.vals_r
  });
}
