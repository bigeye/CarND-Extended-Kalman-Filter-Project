#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  Qv_.resize(2, 2);
  int noise_ax = 9;
  int noise_ay = 9;
  Qv_ << noise_ax, 0,
      0, noise_ay;

  H_.resize(2, 4);
  H_ << 1, 0, 0, 0,
        0, 1, 0, 0;
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.P_ = MatrixXd(4, 4);
    ekf_.P_ << 100, 0, 0, 0,
        0, 100, 0, 0,
        0, 0, 100, 0,
        0, 0, 0, 100;
    //ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      float rho = measurement_pack.raw_measurements_[0];
      float phi = measurement_pack.raw_measurements_[1];
      ekf_.x_ << rho * sin(phi),
          rho * cos(phi),
          0,
          0;
      previous_timestamp_ = measurement_pack.timestamp_;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
      previous_timestamp_ = measurement_pack.timestamp_;
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */

  float time_delta = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  MatrixXd F = MatrixXd::Identity(4, 4);
  F(0, 2) = F(1, 3) = time_delta;

  MatrixXd G(4, 2);
  G << (time_delta * time_delta / 2.0), 0,
      0, (time_delta * time_delta / 2.0),
      time_delta, 0,
      0, time_delta;

  MatrixXd Q = G * Qv_ * G.transpose();

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    ekf_.Init(ekf_.x_, ekf_.P_, F, H_, R_radar_, Q);
  } else {
    ekf_.Init(ekf_.x_, ekf_.P_, F, H_, R_laser_, Q);
  }

  ekf_.Predict();

  /**
   * Update
   */
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    float px = ekf_.x_(0), py = ekf_.x_(1), vx = ekf_.x_(2), vy = ekf_.x_(3);
    float px2py2 = px*px + py*py;
    float px2py2sqrt = sqrt(px2py2);
    MatrixXd H(3,4);
    H.setZero();
    if (px2py2 > 1e-7) {
      H(0, 0) = px / px2py2sqrt;
      H(0, 1) = py / px2py2sqrt;
      H(1, 0) = - py / px2py2;
      H(1, 1) = px / px2py2;
      H(2, 0) = py * (vx * py - vy * px) / (px2py2 * px2py2sqrt);
      H(2, 1) = px * (vy * px - vx * py) / (px2py2 * px2py2sqrt);
      H(2, 2) = px / px2py2sqrt;
      H(2, 3) = py / px2py2sqrt;
    }

    ekf_.H_ = H;
    VectorXd z(3);
    z << measurement_pack.raw_measurements_(0),
        measurement_pack.raw_measurements_(1),
        measurement_pack.raw_measurements_(2);
    ekf_.UpdateEKF(z);
  } else {
    VectorXd z(2);
    z << measurement_pack.raw_measurements_(0), measurement_pack.raw_measurements_(1);
    ekf_.Update(z);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
