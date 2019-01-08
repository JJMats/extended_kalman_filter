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
  R_laser_ << 0.0225,      0,
                   0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09,      0,    0,
                 0, 0.0009,    0,
                 0,      0, 0.09;

  /**
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process and measurement noises
   */
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;
  
  Hj_ << 1, 1, 0, 0,
         1, 1, 0, 0,
         1, 1, 1, 1;
  
  //state transition matrix
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
             0, 1, 0, 1,
             0, 0, 1, 0,
             0, 0, 0, 1;
  
  //state covariance matrix
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1,    0,    0,    0,
             0,    1,    0,    0,
             0,    0, 1000,    0,
             0,    0,    0, 1000;

  noise_ax = 9.0;
  noise_ay = 9.0;
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
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    
    // **** Tune these values (index 2 and 3) below to improve RMSE from the beginning
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // TODO: Convert radar from polar to cartesian coordinates 
      //         and initialize state.
      
      std::cout << "Initializing radar state..." << std::endl;
      float r = measurement_pack.raw_measurements_[0];
      float theta = measurement_pack.raw_measurements_[1];
      float r_dot = measurement_pack.raw_measurements_[2];
      float px = r * cos(theta);
      float py = r * sin(theta);
      float vx = r_dot * cos(theta);
      float vy = r_dot * sin(theta);
      
      ekf_.x_ << px,
                 py,
                 vx,
                 vy;      
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // TODO: Initialize state.
      std::cout << "Initializing laser state..." << std::endl;

      ekf_.x_ << measurement_pack.raw_measurements_[0],
                 measurement_pack.raw_measurements_[1],
                 0,
                 0;      
    }
    
    previous_timestamp_ = measurement_pack.timestamp_;
    
    // done initializing, no need to predict or update
    is_initialized_ = true;

    return;
  }

  /**
   * Prediction
   */

  /**
   * TODO: Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * TODO: Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  
  // Compute the elapsed time between subsequent measurements
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  
  // Update the state transition matrix F to reflect the new elapsed time
  ekf_.F_(0,2) = dt;
  ekf_.F_(1,3) = dt; 
  
  // Update the process covariance matrix Q
  ekf_.Q_ = MatrixXd(4, 4);
  float t2 = dt * dt;
  float t3 = t2 * dt;
  float t4 = t3 * dt;
  float t4_ax = t4/4.0 * noise_ax;
  float t3_ax = t3/2.0 * noise_ax;
  float t2_ax = t2     * noise_ax;
  float t4_ay = t4/4.0 * noise_ay;
  float t3_ay = t3/2.0 * noise_ay;
  float t2_ay = t2     * noise_ay;
  
  ekf_.Q_ << t4_ax,     0, t3_ax,     0,
                 0, t4_ay,     0, t3_ay,
             t3_ax,     0, t2_ax,     0,
                 0, t3_ay,     0, t2_ay;
  
  ekf_.Predict();

  /**
   * Update
   */

  /**
   * TODO:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // TODO: Radar updates
    Tools tools;
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.R_ = R_radar_;    
	ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // TODO: Laser updates    
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);    
  }

  // print the output
  cout << endl;
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
