#include "kalman_filter.h"
#include "tools.h"
#include <iostream>
#include <math.h>

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
  /**
   * TODO: predict the state
   */
  std::cout << "In Predict()..." << std::endl;
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
  std::cout << "Leaving Predict()..." << std::endl;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */
  std::cout << "In Update()" << std::endl;
  
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;
  
  // New Estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
  
  std::cout << "Leaving Update()" << std::endl;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */
  
  std::cout << "In UpdateEKF()..." << std::endl;
  Tools tools;
  float px  = x_(0);
  float py  = x_(1);
  float vx = x_(2);
  float vy = x_(3);
  
  float rho = sqrt(px*px + py*py);
  float theta = atan2(py, px);
  float rho_dot = 0;
  if (fabs(rho) >= 0.0001) {
    rho_dot = (px*vx + py*vy) / rho;
  }
  VectorXd z_pred = VectorXd(3);
  z_pred << rho, theta, rho_dot;
  VectorXd y = z - z_pred;
  
  //Matrix Hj = Tools::CalculateJacobian(x_);
  
  MatrixXd Hjt = H_.transpose();
  //float rho = sqrt((x_[0]*x_[0] + x_[1]*x_[1]));
  //float phi = atan2(x_[1], x_[0]);
  //float rho_dot = 0.0;
  //if (fabs(rho) >= 0.0001){
  //	rho_dot = (x_[0]*x_[2] + x_[1]*x_[3])/rho;
  //}
  
  //VectorXd h_x_Prime = tools.Calculate_X_Prime(x_); 
  //VectorXd h_x_Prime(3);
  //h_x_Prime << rho, phi, rho_dot;
  //std::cout << "h_x_prime: " << h_x_Prime << std::endl;
  //VectorXd z_pred = Hj * x_; // Does this need to be x'?
  
  // *** View Lesson 25: Extended Kalman Filters: EKF Algorithm Generalization
  // I need the h function y = z - h(x'), using Hj instead of H
  // x' should be x' = f(x, u), and use Fj instead of F
  //    u can be set to zero (mean state variable)
  //VectorXd y = z - h_x_Prime;  // This should use the h function directly, not y = z - Hx'
  MatrixXd S = H_ * P_ * Hjt + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHjt = P_ * Hjt;
  MatrixXd K = PHjt * Si;
  
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
  std::cout << "Leaving UpdateEKF()..." << std::endl;
}
