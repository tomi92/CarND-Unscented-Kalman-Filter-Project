#ifndef UKF_H
#define UKF_H

#include <fstream>
#include <string>
#include <tuple>
#include <vector>
#include "Eigen/Dense"
#include "measurement_package.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
 public:
  static const double kPi;
  static const double kHalfPi;
  static const double kTwoPi;

  /** Settings **/

  ///* if this is false, laser measurements will be ignored (except for init)
  const bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  const bool use_radar_;

  /** Constants **/

  ///* State dimension
  const int n_x_;

  ///* Augmented state dimension
  const int n_aug_;

  ///* Number of sigma points
  const int n_sig_;

  ///* Sigma point spreading parameter
  const double lambda_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  const double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  const double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  const double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  const double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  const double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  const double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  const double std_radrd_;

  ///* Weights of sigma points
  const VectorXd weights_;

  /** State **/

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* time when the state is true, in us
  long long time_us_;

  /** Working variables **/

  ///* augmented sigma points matrix
  MatrixXd Xsig_aug_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);

 private:
  const double Normalize(double rad);
  const double RadToDeg(double rad);
  const Eigen::VectorXd GenerateWeights(int n_sig, int n_aug, double lambda);
  void InitStateFromRadar(const VectorXd& radar_data);
  void InitStateFromLaser(const Eigen::VectorXd& laser_data);
  void GenerateAugmentedSigmaPoints();
  const Eigen::MatrixXd Sqrt(const Eigen::MatrixXd& M);
  void PredictSigmaPoints(double delta_t);
  void PredictMeanAndCovariance();
  void PredictRadarMeasurement(Eigen::MatrixXd& Zsig, Eigen::VectorXd& z_pred,
                               Eigen::MatrixXd& S);
  void PredictLidarMeasurement(Eigen::MatrixXd& Zsig, Eigen::VectorXd& z_pred,
                               Eigen::MatrixXd& S);
};

#endif /* UKF_H */
