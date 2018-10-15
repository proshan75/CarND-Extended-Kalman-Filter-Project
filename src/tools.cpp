#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
    const vector<VectorXd> &ground_truth) {
    /**
      * Calculate the RMSE here.
    */
    VectorXd rmse(4);
    rmse << 0, 0, 0, 0;

    // check the validity of the following inputs:
    //  * the estimation vector size should not be zero
    //  * the estimation vector size should equal ground truth vector size
    // ... your code here
    if (estimations.size() == 0 || estimations.size() != ground_truth.size())
    {
        cout << "estimation size cannot be zero.";
    }

    //accumulate squared residuals
    int size = estimations.size();
    VectorXd residual(4), squared_residual(4), total_residual(4);
    residual << 0, 0, 0, 0;
    squared_residual << 0, 0, 0, 0;
    total_residual << 0, 0, 0, 0;
    for (int i = 0; i < size; ++i) {
        // ... your code here
        residual = estimations[i] - ground_truth[i];
        squared_residual = (residual.array() * residual.array());
        total_residual = total_residual + squared_residual;
    }

    //calculate the mean
    // ... your code here
    VectorXd mean_residual = total_residual / (1.0 * size);

    //calculate the squared root
    // ... your code here
    rmse = mean_residual.array().sqrt();

    //return the result
    return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
    /**
      * Calculate a Jacobian here.
    */
    MatrixXd Hj(3, 4);
    //recover state parameters
    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);

    float px2_plus_py2 = px * px + py * py;
    //check division by zero
    if (fabs(px2_plus_py2) < 0.0001)
    {
        cout << "px2_plus_py2 is zero" << endl;
        return Hj;
    }

    float sqrt_px2_plus_py2 = sqrt(px2_plus_py2);

    float Hj_0_0 = px / sqrt_px2_plus_py2;
    float Hj_0_1 = py / sqrt_px2_plus_py2;
    float Hj_1_0 = -py / px2_plus_py2;
    float Hj_1_1 = px / px2_plus_py2;
    float Hj_2_0 = py * (vx*py - vy * px) / pow(px2_plus_py2, 1.5);
    float Hj_2_1 = px * (vy*px - vx * py) / pow(px2_plus_py2, 1.5);
    float Hj_2_2 = px / sqrt_px2_plus_py2;
    float Hj_2_3 = py / sqrt_px2_plus_py2;

    //check division by zero

    //compute the Jacobian matrix
    Hj << Hj_0_0, Hj_0_1, 0, 0,
        Hj_1_0, Hj_1_1, 0, 0,
        Hj_2_0, Hj_2_1, Hj_2_2, Hj_2_3;

    return Hj;
}
