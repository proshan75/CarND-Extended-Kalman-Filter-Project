#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
    is_initialized_ = false;

    previous_timestamp_ = 0;

    // initializing matrices
    Hj_ = MatrixXd(3, 4);

    //measurement covariance matrix - laser
    R_laser_ = MatrixXd(2, 2);
    R_laser_ << 0.0225, 0,
        0, 0.0225;

    //measurement covariance matrix - radar
    R_radar_ = MatrixXd(3, 3);
    R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

    H_laser_ = MatrixXd(2, 4);
    H_laser_ << 1, 0, 0, 0,
        0, 1, 0, 0;

    /**
      * Finish initializing the FusionEKF.
      * Set the process and measurement noises
    */
    //ekf_.x_
    ekf_.P_ = MatrixXd(4, 4);
    ekf_.P_ << 1, 0, 1, 0,
        0, 1, 0, 1,
        0, 0, 1, 0,
        0, 0, 0, 1;
    ekf_.F_ = MatrixXd(4, 4);
    ekf_.F_ << 1, 0, 1, 0,
        0, 1, 0, 1,
        0, 0, 1, 0,
        0, 0, 0, 1;
    //ekf_.Q_
    //ekf_.H_
    ekf_.R_ = MatrixXd(2, 2);
    ekf_.R_ << 1, 0,
        0, 1;


}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


    /*****************************************************************************
     *  Initialization
     ****************************************************************************/
    if (!is_initialized_) {
        /**
        TODO:
          * Initialize the state ekf_.x_ with the first measurement.
          * Create the covariance matrix.
          * Remember: you'll need to convert radar from polar to cartesian coordinates.
        */
        // first measurement
        cout << "EKF: " << endl;
        ekf_.x_ = VectorXd(4);
        ekf_.x_ << 1, 1, 1, 1;

        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
            /**
            Convert radar from polar to cartesian coordinates and initialize state.
            */
            //cout << "Initializing radar package..." << endl;
            ekf_.x_(0) = measurement_pack.raw_measurements_[0] * cos(measurement_pack.raw_measurements_[1]);
            ekf_.x_(1) = measurement_pack.raw_measurements_[0] * sin(measurement_pack.raw_measurements_[1]);
            ekf_.x_(2) = 0;
            ekf_.x_(3) = 0;
        }
        else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
            /**
            Initialize state.
            */
            //cout << "Initializing laser package... " << endl;
            ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
        }

        //cout << " initialized ekf_.x_: " << ekf_.x_ << endl;

        ekf_.F_ << 1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1;

        previous_timestamp_ = measurement_pack.timestamp_;
        noise_ax = 9; // sigma_ax squared
        noise_ay = 9; // sigma_ay squared
        // done initializing, no need to predict or update
        is_initialized_ = true;
        return;
    }

    /*****************************************************************************
     *  Prediction
     ****************************************************************************/

     /**
        * Update the state transition matrix F according to the new elapsed time.
         - Time is measured in seconds.
        * Update the process noise covariance matrix.
        * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
      */
      //compute the time elapsed between the current and previous measurements
    float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
    previous_timestamp_ = measurement_pack.timestamp_;

    ekf_.F_ << 1, 0, dt, 0,
        0, 1, 0, dt,
        0, 0, 1, 0,
        0, 0, 0, 1;

    //cout << "F: " << ekf_.F_ << endl;

    ekf_.Q_ = MatrixXd(4, 4);

    //        --                                                                     --
    //    Q = |(Δt^4 * σax^2)/4  ​0                 (Δt^3 * σax^2)/2 0 ​                |
    //        |0                 (Δt^4 * σay^2)/4 ​ 0                (Δt^3 * σay^2)/2  |
    //        |​(Δt^3 * σax^2)/2 ​ 0                 (Δt^2 * σax^2)   0 ​                |
    //        |0                 (Δt^3 * σay^2)/2 ​ 0                Δt^2 * σay^2      |
    //        --                                                                     --
    // noise_ax = σax^2
    // noise_ay = σay^2
    ekf_.Q_ << pow(dt, 4)*noise_ax / 4, 0, pow(dt, 3)*noise_ax / 2, 0,
        0, pow(dt, 4)*noise_ay / 4, 0, pow(dt, 3)*noise_ay / 2,
        pow(dt, 3)*noise_ax / 2, 0, dt*dt*noise_ax, 0,
        0, pow(dt, 3)*noise_ay / 2, 0, dt*dt*noise_ay;

    //cout << "Covariance Matrix Q calculated as: " << ekf_.Q_ << endl;

    ekf_.Predict();

    /*****************************************************************************
     *  Update
     ****************************************************************************/

     /**
        * Use the sensor type to perform the update step.
        * Update the state and covariance matrices.
      */

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        // Radar updates
        //cout << "Updating radar measurement ..." << endl;
        Hj_ = tools.CalculateJacobian(ekf_.x_);
        ekf_.H_ = Hj_;
        ekf_.R_ = R_radar_;
        ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    }
    else {
        // Laser updates
        //cout << "Updating laser measurement ..." << endl;
        ekf_.H_ = H_laser_;
        ekf_.R_ = R_laser_;
        ekf_.Update(measurement_pack.raw_measurements_);

    }

    // print the output
    cout << "x_ = " << ekf_.x_ << endl;
    cout << "P_ = " << ekf_.P_ << endl;
}
