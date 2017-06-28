#include "kalman_filter.h"
#include "tools.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

Tools tools;

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
    MatrixXd Ft = F_.transpose();
    P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::calcKF(const VectorXd &y) {
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd PHt = P_ * Ht;
    MatrixXd K = PHt * Si;

    //new estimate
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
}

void KalmanFilter::Update(const VectorXd &z) {
    calcKF(z - H_ * x_);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  float pho = sqrt(x_(0) * x_(0) + x_(1) * x_(1));
  float chi = atan2(x_(1), x_(0));

  float phodot = 0;
  if (pho != 0) {
	// if pho == 0, means x_(0) and x_(1) == 0, so this fomular result is 0
  	phodot = (x_(0) * x_(2) + x_(1) * x_(3)) / pho;
  }

  VectorXd x = VectorXd(3); 
  x << pho, chi, phodot;
  
  VectorXd y = z - x;
  
  y(1) = atan2(sin(y(1)), cos(y(1)));
  calcKF(y);
}

