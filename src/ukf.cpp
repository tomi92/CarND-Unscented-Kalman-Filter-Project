#include "ukf.h"
#include <iostream>
#include "Eigen/Dense"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF()
    : use_laser_(true),
      use_radar_(true),
      n_x_(5),
      n_aug_(n_x_ + 2),
      n_sig_(2 * n_aug_ + 1),
      lambda_(3 - n_aug_),  // tunable parameter
      std_a_(30),           // TODO tune
      std_yawdd_(30),       // TODO tune
      std_laspx_(0.15),
      std_laspy_(0.15),
      std_radr_(0.3),
      std_radphi_(0.03),
      std_radrd_(0.3),
      weights_(CalculateWeights(n_sig_, n_aug_, lambda_)),
      is_initialized_(false),
      x_(n_x_),        // value uninitialized
      P_(n_x_, n_x_),  // value uninitialized
      time_us_(0),
      Xsig_pred_(n_x_, n_sig_)  // value uninitialized
{}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
}

VectorXd UKF::CalculateWeights(int n_sig, int n_aug, double lambda) {
  VectorXd weights = VectorXd(n_sig);
  weights(0) = lambda / (lambda + n_aug);
  weights.segment(1, n_sig - 1).fill(0.5 / (n_aug + lambda));
  return weights;
}