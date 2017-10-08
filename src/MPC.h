#ifndef MPC_H
#define MPC_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include <cppad/cppad.hpp>

using CppAD::AD;
using namespace std;

class MPC {
 public:
  //nested class for using ipopt optimization
  class FG_eval 
  {
  public:
      Eigen::VectorXd coeffs;
      FG_eval(Eigen::VectorXd coeffs, MPC& mpc) : mpc_(mpc) { this->coeffs = coeffs; }

      typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
      void operator()(MPC::FG_eval::ADvector& fg, const MPC::FG_eval::ADvector& vars);
  private:
      MPC& mpc_;
  };

  // This value assumes the model presented in the classroom is used.
  //
  // It was obtained by measuring the radius formed by running the vehicle in the
  // simulator around in a circle with a constant steering angle and velocity on a
  // flat terrain.
  //
  // Lf was tuned until the the radius formed by the simulating the model
  // presented in the classroom matched the previous radius.
  //
  // This is the length from front to CoG that has a similar radius.
  static constexpr double LF = 2.67;
  //N is number of timesteps, and 1/dt is evaluation frequency
  //T = N*dt is the prodiction horizon. dt should not be too large, otherwise
  //the discretization error will be big.
  //T should not be too big, otherwise we are predicting too far ahead into the future
  //and T should not be too small either, otherwise we are only looking at local condition
  //and error will be big in both cases
  MPC(size_t N, double dt, double desiredspeed);
  virtual ~MPC();

  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuatotions.
  vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);

private:
  size_t N_;
  double dt_;
  double desiredSpeed_;

  size_t x_start_;
  size_t y_start_;
  size_t psi_start_;
  size_t v_start_;
  size_t cte_start_;
  size_t epsi_start_;
  size_t delta_start_;
  size_t a_start_;
};

#endif /* MPC_H */
