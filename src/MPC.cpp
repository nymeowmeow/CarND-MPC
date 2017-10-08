#ifndef HAVE_STDDEF_H
#define HAVE_STDDEF_H
#include <coin/IpTNLP.hpp>
#undef HAVE_STDDEF_H
#endif

#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

MPC::MPC(size_t N, double dt, double desiredSpeed) : N_(N), dt_(dt), desiredSpeed_(desiredSpeed)
{
  x_start_     = 0;
  y_start_     = x_start_ + N_;
  psi_start_   = y_start_ + N_;
  v_start_     = psi_start_ + N_;
  cte_start_   = v_start_ + N_;
  epsi_start_  = cte_start_ + N_;
  delta_start_ = epsi_start_ + N_;
  a_start_     = delta_start_ + N_ - 1;
}

MPC::~MPC() {}

void
MPC::FG_eval::operator()(MPC::FG_eval::ADvector& fg, const MPC::FG_eval::ADvector& vars)
{
    //the cost is stored in its first element of fg, any additions to the cost
    //have to update fg[0]
    fg[0] = 0.0;

    //cost based on reference state
    for (int t = 0; t < mpc_.N_; ++t) {
        fg[0] += 1*CppAD::pow(vars[mpc_.cte_start_ + t], 2);
        fg[0] += 1500*CppAD::pow(vars[mpc_.epsi_start_ + t], 2);
        fg[0] += 1*CppAD::pow(vars[mpc_.v_start_ + t] - mpc_.desiredSpeed_, 2);
    }
    // Minimize the use of actuators.
    for (int t = 0; t < mpc_.N_ - 1; ++t) 
    {
        fg[0] += 20000*CppAD::pow(vars[mpc_.delta_start_ + t], 2);
        fg[0] += 1*CppAD::pow(vars[mpc_.a_start_ + t], 2);
    }
    // Minimize the value gap between sequential actuations.
    for (int t = 0; t < mpc_.N_ - 2; ++t) 
    {
        fg[0] += 500*(CppAD::pow(vars[mpc_.delta_start_ + t + 1] - vars[mpc_.delta_start_ + t], 2));
        fg[0] += 1*(CppAD::pow(vars[mpc_.a_start_ + t + 1] - vars[mpc_.a_start_ + t], 2));
    }
    // Setup Constraints
    //
    // NOTE: In this section you'll setup the model constraints.
    // Initial constraints
    //
    // We add 1 to each of the starting indices due to cost being located at
    // index 0 of `fg`.
    // This bumps up the position of all the other values.
    fg[1 + mpc_.x_start_]    = vars[mpc_.x_start_];
    fg[1 + mpc_.y_start_]    = vars[mpc_.y_start_];
    fg[1 + mpc_.psi_start_]  = vars[mpc_.psi_start_];
    fg[1 + mpc_.v_start_]    = vars[mpc_.v_start_];
    fg[1 + mpc_.cte_start_]  = vars[mpc_.cte_start_];
    fg[1 + mpc_.epsi_start_] = vars[mpc_.epsi_start_];
    // The rest of the constraints
    for (int t = 1; t < mpc_.N_; ++t) 
    {
      // The state at time t+1 .
      AD<double> x1    = vars[mpc_.x_start_ + t];
      AD<double> y1    = vars[mpc_.y_start_ + t];
      AD<double> psi1  = vars[mpc_.psi_start_ + t];
      AD<double> v1    = vars[mpc_.v_start_ + t];
      AD<double> cte1  = vars[mpc_.cte_start_ + t];
      AD<double> epsi1 = vars[mpc_.epsi_start_ + t];

      // The state at time t.
      AD<double> x0    = vars[mpc_.x_start_ + t - 1];
      AD<double> y0    = vars[mpc_.y_start_ + t - 1];
      AD<double> psi0  = vars[mpc_.psi_start_ + t - 1];
      AD<double> v0    = vars[mpc_.v_start_ + t - 1];
      AD<double> cte0  = vars[mpc_.cte_start_ + t - 1];
      AD<double> epsi0 = vars[mpc_.epsi_start_ + t - 1];

      // Only consider the actuation at time t.
      AD<double> delta0 = vars[mpc_.delta_start_ + t - 1];
      AD<double> a0     = vars[mpc_.a_start_ + t - 1];

      AD<double> f0 = coeffs[0] + coeffs[1] * x0 + coeffs[2]*pow(x0, 2) + coeffs[3]*pow(x0, 3);
      AD<double> psides0 = CppAD::atan(coeffs[1] + 2*coeffs[2]*x0 + 3*coeffs[3]*pow(x0, 2));

      // Here's `x` to get you started.
      // The idea here is to constraint this value to be 0.
      //
      // Recall the equations for the model:
      // x_[t+1] = x[t] + v[t] * cos(psi[t]) * dt
      // y_[t+1] = y[t] + v[t] * sin(psi[t]) * dt
      // psi_[t+1] = psi[t] + v[t] / Lf * delta[t] * dt
      // v_[t+1] = v[t] + a[t] * dt
      // cte[t+1] = f(x[t]) - y[t] + v[t] * sin(epsi[t]) * dt
      // epsi[t+1] = psi[t] - psides[t] + v[t] * delta[t] / Lf * dt
      fg[1 + mpc_.x_start_ + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * mpc_.dt_);
      fg[1 + mpc_.y_start_ + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * mpc_.dt_);
      fg[1 + mpc_.psi_start_ + t] = psi1 - (psi0 + v0 * delta0 / MPC::LF * mpc_.dt_);
      fg[1 + mpc_.v_start_ + t] = v1 - (v0 + a0 * mpc_.dt_);
      fg[1 + mpc_.cte_start_ + t] =
          cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * mpc_.dt_));
      fg[1 + mpc_.epsi_start_ + t] =
          epsi1 - ((psi0 - psides0) + v0 * delta0 / MPC::LF * mpc_.dt_);
    }
}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  double x    = state[0];
  double y    = state[1];
  double psi  = state[2];
  double v    = state[3];
  double cte  = state[4];
  double epsi = state[5];

  // Set the number of model variables (includes both states and inputs).
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:
  //
  // 4 * 10 + 2 * 9, N timesteps => N -1 actuations
  size_t n_vars = N_*6 + (N_-1)*2;
  // Set the number of constraints
  size_t n_constraints = N_*6;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; ++i) 
  {
    vars[i] = 0;
  }
  // Set the initial variable values
  vars[x_start_]    = x;
  vars[y_start_]    = y;
  vars[psi_start_]  = psi;
  vars[v_start_]    = v;
  vars[cte_start_]  = cte;
  vars[epsi_start_] = epsi;
  //Lower and Upper limits for x
  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  //Set lower and upper limits for variables.
  // Set all non-actuators upper and lowerlimits
  // to the max negative and positive values.
  for (int i = 0; i < delta_start_; ++i) 
  {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }
  // The upper and lower limits of delta are set to -25 and 25
  // degrees (values in radians).
  for (int i = delta_start_; i < a_start_; ++i) 
  {
    vars_lowerbound[i] = -0.436332;
    vars_upperbound[i] = 0.436332;
  }
  // Acceleration/decceleration upper and lower limits.
  for (int i = a_start_; i < n_vars; ++i) 
  {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }
  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < n_constraints; ++i) 
  {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }
  constraints_lowerbound[x_start_]    = x;
  constraints_lowerbound[y_start_]    = y;
  constraints_lowerbound[psi_start_]  = psi;
  constraints_lowerbound[v_start_]    = v;
  constraints_lowerbound[cte_start_]  = cte;
  constraints_lowerbound[epsi_start_] = epsi;

  constraints_upperbound[x_start_]    = x;
  constraints_upperbound[y_start_]    = y;
  constraints_upperbound[psi_start_]  = psi;
  constraints_upperbound[v_start_]    = v;
  constraints_upperbound[cte_start_]  = cte;
  constraints_upperbound[epsi_start_] = epsi;
  // object that computes objective and constraints
  MPC::FG_eval fg_eval(coeffs, *this);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, MPC::FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  bool ok = true;
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Status: " << ok << ", Cost " << cost << std::endl;

  vector<double> result = { solution.x[delta_start_], solution.x[a_start_], (double)N_};
  for (int i = 0; i < N_; ++i)
  {
      result.push_back(solution.x[i + x_start_]);
      result.push_back(solution.x[i + y_start_]);
  }

  return result;
}
