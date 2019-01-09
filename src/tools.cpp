#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */
  int n = estimations.size();
  VectorXd result(estimations[0].size());
  result.setZero();
  for (int i = 0; i < n; ++i) {
    VectorXd elem = (estimations[i] - ground_truth[i]);
    elem = elem.array() * elem.array();
    result = result + elem;
  }
  result = result.array() / estimations.size();
  result = result.array().sqrt();
  return result;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
}
