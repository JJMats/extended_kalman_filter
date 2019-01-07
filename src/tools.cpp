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
  
  //std::cout << "Estimations: " << std::endl;
  //std::cout << estimations << std::endl;
  //std::cout << "Ground Truth: " << std::endl << ground_truth << std::endl;
  
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
  //rmse = rmse.array() / estimations.size();
  rmse = rmse / estimations.size();
  
  // Calculate the square root of the mse value
  rmse = rmse.array().sqrt();
  std::cout << "RMSE: " << rmse << std::endl;
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
  float py_h1 = px / sqrt(denom);
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

VectorXd Tools::Calculate_X_Prime(const VectorXd& x){
  std::cout << "In Calculate_X_Prime()..." << std::endl;
  
  VectorXd h_x(3);
  float rho = sqrt((x[0] * x[0]) + (x[1] * x[1]));
  float phi = atan2(x[1], x[0]);
  //float phi = Normalize_phi(atan2(x[1], x[0]));
  std::cout << "Back from Normalize_phi()" << std::endl;
  float rho_dot = 0.0;
  std::cout << "x: " << x << std::endl;
  if(fabs(rho) > 0.0001){    
  	rho_dot = ((x[0]*x[2]) + (x[1]*x[3])) / rho;
  }
  std::cout << "rho_dot: " << rho_dot << std::endl;
  h_x << rho, phi, rho_dot;
  std::cout << "Leaving Calculate_X_Prime()..." << std::endl;
  return h_x;
}

float Tools::Normalize_phi(float phi){
  std::cout << "In Normalize_phi()..." << std::endl;
  std::cout << "Phi: " << phi << std::endl;
  
  float phi_norm = phi;
  while(abs(phi_norm) > M_PI){
  	phi_norm /= M_PI;
  }
  
  std::cout << "Leaving Normalize_phi()..." << std::endl;
  std::cout << "Phi_norm: " << phi_norm << std::endl;
  return phi_norm;
}