#include "kalman_filter.h"
#include <iostream>
#define PI 3.14159265

using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  x_ = F_ * x_;
  P_ = (F_ * P_ * F_.transpose()) + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  VectorXd y = z - H_ * x_;
  MatrixXd S = H_ * P_ * H_.transpose() + R_;
  MatrixXd K = P_ * H_.transpose() * S.inverse();

  x_ = x_ + K * y;
  P_ = (MatrixXd::Identity(P_.rows(), P_.cols()) - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  float px = x_(0), py = x_(1), vx = x_(2), vy = x_(3);
  VectorXd hx(3);
  hx << sqrt(px*px+py*py),
      atan2(py, px),
      (px*vx+py*vy)/sqrt(px*px+py*py);
  VectorXd y = z - hx;//H_*x_;
  if (y(1) < -PI) y(1) += 2*PI;
  if (y(1) > PI) y(1) -= 2*PI;
  MatrixXd S = H_ * P_ * H_.transpose() + R_;
  MatrixXd K = P_ * H_.transpose() * S.inverse();

  x_ = x_ + K * y;
  P_ = (MatrixXd::Identity(P_.rows(), P_.cols()) - K * H_) * P_;
}
