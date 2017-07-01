#include "ukf.h"
#include <iostream>
#include "Eigen/Dense"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

const double UKF::kPi = 3.14159265358979323846;
const double UKF::kHalfPi = 0.5 * kPi;
const double UKF::kTwoPi = 2.0 * kPi;

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
      std_a_(1.5),
      std_yawdd_(kPi / 4),
      std_laspx_(0.15),
      std_laspy_(0.15),
      std_radr_(0.3),
      std_radphi_(0.03),
      std_radrd_(0.3),
      weights_(GenerateWeights(n_sig_, n_aug_, lambda_)),
      is_initialized_(false),
      x_(n_x_),        // value uninitialized
      P_(n_x_, n_x_),  // value uninitialized
      time_us_(0),
      nis_radar_(0),
      nis_lidar_(0),
      Xsig_aug_(n_aug_, n_sig_),  // value uninitialized
      Xsig_pred_(n_x_, n_sig_)    // value uninitialized
{}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  if (!is_initialized_) {
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      InitStateFromRadar(meas_package.raw_measurements_);
    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      InitStateFromLaser(meas_package.raw_measurements_);
    }

    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;

    return;
  }

  double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;

  // Prevent strange behaviour when changing datasets
  if (delta_t < 0 || delta_t > 1) {
    delta_t = 1;
  }

  time_us_ = meas_package.timestamp_;

  Prediction(delta_t);

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    if (use_radar_) {
      UpdateRadar(meas_package);
    }
  } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    if (use_laser_) {
      UpdateLidar(meas_package);
    }
  } else {
    cerr << "Unsupported sensor type: " << meas_package.sensor_type_ << endl;
  }

  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;
  // cout << nis_lidar_ << ';' << nis_radar_ << endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  GenerateAugmentedSigmaPoints();
  PredictSigmaPoints(delta_t);
  PredictMeanAndCovariance();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  You'll also need to calculate the lidar NIS.
  */
  const int n_z = 2;  // lidar
  const VectorXd z = meas_package.raw_measurements_;

  MatrixXd Zsig;
  VectorXd z_pred;
  MatrixXd S;
  PredictLidarMeasurement(Zsig, z_pred, S);
  MatrixXd S_inverse = S.inverse();

  MatrixXd Tc = MatrixXd::Zero(n_x_, n_z);  // cross correlation

  for (int i = 0; i < Xsig_pred_.cols(); i++) {
    VectorXd xdiff = Xsig_pred_.col(i) - x_;
    xdiff(3) = Normalize(xdiff(3));
    VectorXd zdiff = Zsig.col(i) - z_pred;
    Tc += weights_(i) * xdiff * zdiff.transpose();
  }

  MatrixXd K = Tc * S_inverse;  // Kalman gain

  VectorXd zdiff = z - z_pred;
  x_ += K * zdiff;
  P_ -= K * S * K.transpose();

  nis_lidar_ = zdiff.transpose() * S_inverse * zdiff;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  const int n_z = 3;  // radar
  const VectorXd z = meas_package.raw_measurements_;

  MatrixXd Zsig;
  VectorXd z_pred;
  MatrixXd S;
  PredictRadarMeasurement(Zsig, z_pred, S);
  MatrixXd S_inverse = S.inverse();

  MatrixXd Tc = MatrixXd::Zero(n_x_, n_z);  // cross correlation

  for (int i = 0; i < Xsig_pred_.cols(); i++) {
    VectorXd xdiff = Xsig_pred_.col(i) - x_;
    xdiff(3) = Normalize(xdiff(3));
    VectorXd zdiff = Zsig.col(i) - z_pred;
    zdiff(1) = Normalize(zdiff(1));
    Tc += weights_(i) * xdiff * zdiff.transpose();
  }

  MatrixXd K = Tc * S_inverse;  // Kalman gain

  VectorXd zdiff = z - z_pred;
  zdiff(1) = Normalize(zdiff(1));
  x_ += K * zdiff;
  P_ -= K * S * K.transpose();

  nis_radar_ = zdiff.transpose() * S_inverse * zdiff;
}

const VectorXd UKF::GenerateWeights(int n_sig, int n_aug, double lambda) {
  VectorXd weights = VectorXd(n_sig);
  weights(0) = lambda_ / (lambda + n_aug);
  for (int i = 1; i < 2 * n_aug + 1; i++) {
    double weight = 0.5 / (lambda + n_aug);
    weights(i) = weight;
  }
  return weights;
}

void UKF::InitStateFromRadar(const VectorXd& radar_data) {
  const double rho = radar_data(0);
  const double phi = radar_data(1);
  const double rho_dot = radar_data(2);

  const double px = rho * cos(phi);
  const double py = rho * sin(phi);
  double vel_abs = 0;
  double yaw_angle = 0;
  const double yaw_rate = 0;

  if (rho > 0.001) {  // rho is non-negative, so no need for abs
    vel_abs = abs(rho_dot);
    yaw_angle = rho_dot > 0 ? atan2(py, px) : Normalize(atan2(py, px) + kPi);
  }

  // clang-format off
  x_ << px, py, vel_abs, yaw_angle, yaw_rate;
  P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;
  // clang-format on
}

void UKF::InitStateFromLaser(const VectorXd& laser_data) {
  const double px = laser_data(0);
  const double py = laser_data(1);
  const double vel_abs = 0;
  const double yaw_angle = 0;
  const double yaw_rate = 0;

  // clang-format off
  x_ << px, py, vel_abs, yaw_angle, yaw_rate;
  P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;
  // clang-format on
}

const double UKF::Normalize(double rad) {
  while (rad <= -kPi) {
    rad += kTwoPi;
  }
  while (rad > kPi) {
    rad -= kTwoPi;
  }
  return rad;
}

const double UKF::RadToDeg(double rad) { return rad / kPi * 180.0; }

void UKF::GenerateAugmentedSigmaPoints() {
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug << x_, 0, 0;

  MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5, 5) = std_a_ * std_a_;
  P_aug(6, 6) = std_yawdd_ * std_yawdd_;

  // create square root matrix
  MatrixXd A = sqrt(lambda_ + n_aug_) * Sqrt(P_aug);
  Xsig_aug_ << x_aug, A.colwise() + x_aug, (-A).colwise() + x_aug;
}

const MatrixXd UKF::Sqrt(const MatrixXd& M) { return M.llt().matrixL(); }

void UKF::PredictSigmaPoints(double delta_t) {
  for (int i = 0; i < n_sig_; i++) {
    VectorXd xsig_aug = Xsig_aug_.col(i);
    VectorXd xsig_pred = VectorXd(n_x_);

    double x = xsig_aug(0);
    double y = xsig_aug(1);
    double v = xsig_aug(2);
    double ori = xsig_aug(3);
    double turnRate = xsig_aug(4);
    double noise_a = xsig_aug(5);
    double noise_rateChange = xsig_aug(6);

    double cos_ori = cos(ori);
    double sin_ori = sin(ori);
    double half_delta_t_2 = 0.5 * delta_t * delta_t;

    xsig_pred = xsig_aug.head(n_x_);

    // process model, adding integral + the effect of noise
    if (abs(turnRate) < 0.001) {
      // clang-format off
      xsig_pred(0) += v * cos_ori * delta_t + half_delta_t_2 * cos_ori * noise_a;
      xsig_pred(1) += v * sin_ori * delta_t + half_delta_t_2 * sin_ori * noise_a;
      xsig_pred(2) += delta_t * noise_a;
      xsig_pred(3) += half_delta_t_2 * noise_rateChange;
      xsig_pred(4) += delta_t * noise_rateChange;
      // clang-format on
    } else {
      // clang-format off
      const double ori2 = ori + turnRate*delta_t;
      xsig_pred(0) += (v / turnRate) * (sin(ori2) - sin_ori) + half_delta_t_2 * cos_ori * noise_a;
      xsig_pred(1) += (v / turnRate) * (-cos(ori2) + cos_ori) + half_delta_t_2 * sin_ori * noise_a;
      xsig_pred(2) += delta_t * noise_a;
      xsig_pred(3) += turnRate * delta_t + half_delta_t_2 * noise_rateChange; 
      xsig_pred(4) += delta_t * noise_rateChange;
      // clang-format on
    }

    Xsig_pred_.col(i) = xsig_pred;
  }
}

void UKF::PredictMeanAndCovariance() {
  x_.fill(0.0);
  // https://stackoverflow.com/questions/491738/how-do-you-calculate-the-average-of-a-set-of-circular-data
  double xx = 0;
  double yy = 0;
  for (int c = 0; c < Xsig_pred_.cols(); c++) {
    x_ += weights_(c) * Xsig_pred_.col(c);
    xx += weights_(c) * cos(Xsig_pred_(3, c));
    yy += weights_(c) * sin(Xsig_pred_(3, c));
  }
  x_(3) = sqrt(xx * xx + yy * yy) < 0.001 ? 0 : atan2(yy, xx);

  P_.fill(0.0);
  for (int i = 0; i < Xsig_pred_.cols(); i++) {
    VectorXd x_diff = (Xsig_pred_.col(i) - x_);
    x_diff(3) = Normalize(x_diff(3));
    P_ += weights_(i) * x_diff * x_diff.transpose();
  }
}

void UKF::PredictRadarMeasurement(MatrixXd& Zsig, VectorXd& z_pred,
                                  MatrixXd& S) {
  const int n_z = 3;  // radar

  Zsig = MatrixXd(n_z, n_sig_);
  z_pred = VectorXd(n_z);
  S = MatrixXd(n_z, n_z);  // measurement covariance matrix

  // transform sigma points into measurement space
  for (int i = 0; i < Zsig.cols(); i++) {
    double x = Xsig_pred_(0, i);
    double y = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double psi = Xsig_pred_(3, i);
    double rho = sqrt(x * x + y * y);
    double phi = (rho < 0.001) ? 0 : atan2(y, x);
    double rhod =
        (rho < 0.001) ? 0 : (x * cos(psi) * v + y * sin(psi) * v) / rho;
    Zsig(0, i) = rho;
    Zsig(1, i) = phi;
    Zsig(2, i) = rhod;
  }

  z_pred.fill(0.0);
  double xx = 0;
  double yy = 0;
  for (int i = 0; i < Zsig.cols(); i++) {
    z_pred += weights_(i) * Zsig.col(i);
    xx += weights_(i) * cos(Zsig(1, i));
    yy += weights_(i) * sin(Zsig(1, i));
  }
  z_pred(1) = sqrt(xx * xx + yy * yy) > 0.001 ? atan2(yy, xx) : 0;

  S.fill(0.0);
  for (int i = 0; i < Zsig.cols(); i++) {
    VectorXd z_diff = (Zsig.col(i) - z_pred);
    z_diff(1) = Normalize(z_diff(1));
    S += weights_(i) * z_diff * z_diff.transpose();
  }

  S(0, 0) += std_radr_ * std_radr_;
  S(1, 1) += std_radphi_ * std_radphi_;
  S(2, 2) += std_radrd_ * std_radrd_;
}

void UKF::PredictLidarMeasurement(MatrixXd& Zsig, VectorXd& z_pred,
                                  MatrixXd& S) {
  const int n_z = 2;  // lidar

  Zsig = MatrixXd(n_z, n_sig_);
  z_pred = VectorXd(n_z);
  S = MatrixXd(n_z, n_z);  // measurement covariance matrix

  // transform sigma points into measurement space
  for (int i = 0; i < Zsig.cols(); i++) {
    Zsig(0, i) = Xsig_pred_(0, i);
    Zsig(1, i) = Xsig_pred_(1, i);
  }

  z_pred.fill(0.0);
  for (int i = 0; i < Zsig.cols(); i++) {
    z_pred += weights_(i) * Zsig.col(i);
  }

  S.fill(0.0);
  for (int i = 0; i < Zsig.cols(); i++) {
    VectorXd tmp = (Zsig.col(i) - z_pred);
    S += weights_(i) * tmp * tmp.transpose();
  }

  S(0, 0) += std_laspx_ * std_laspx_;
  S(1, 1) += std_laspy_ * std_laspy_;
}
