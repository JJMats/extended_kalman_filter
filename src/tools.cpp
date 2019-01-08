#include "tools.h"
#include <iostream>
#include <math.h>

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
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;
  
  if(estimations.size() == 0 || estimations.size() != ground_truth.size()){
  	return rmse;
  }
  
  // Accumulate squared residuals
  for (unsigned int i=0; i < estimations.size(); ++i){
  	VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array() * residual.array();
    rmse += residual;
  }
  
  // Calculate the mean of the residuals
  rmse = rmse / estimations.size();
  
  // Calculate the square root of the mse value
  rmse = rmse.array().sqrt();

  return rmse;  
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
  MatrixXd Hj(3, 4);
  
  // Break state vector into its individual components
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);
  
  float denom = pow(px, 2) + pow(py, 2);
  // Check to verify there is no division by zero
  if(denom < 0.0001){
    std::cout << "CalculateJacobian () - Error - Cannot divide by zero!" << std::endl;
    return Hj;
  }
  
  float px_h1 = px / sqrt(denom);
  float py_h1 = py / sqrt(denom);
  float px_h2 = -1.0 * py / denom;
  float py_h2 = px / denom;
  float px_h3 = (py * ((vx*py) - (vy*px))) / pow(denom, 1.5);
  float py_h3 = (px * ((vy*px) - (vx*py))) / pow(denom, 1.5);
  
  // Compile the Jacobian Matrix
  Hj << px_h1, py_h1,     0,     0,
        px_h2, py_h2,     0,     0,
        px_h3, py_h3, px_h1, py_h1;
  
  return Hj;
}

float Tools::Normalize_phi(float phi){ 
  while (fabs(phi) > M_PI) {
    if (phi > M_PI){
      phi -= 2.0*M_PI;
    } else {
      phi += 2.0*M_PI;
    }
  } 

  return phi;
}