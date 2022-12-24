// Copyright 2019 Alexander Liniger

// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at

//     http://www.apache.org/licenses/LICENSE-2.0

// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

#include "constraints.h"
namespace mpcc{
Constraints::Constraints()
{   
    std::cout << "default constructor, not everything is initialized properly" << std::endl;
}

Constraints::Constraints(double Ts,const PathToJson &path) 
:model_(Ts,path),
param_(Param(path.param_path))
{
}


bool get_line_intersection(double p0_x, double p0_y, double p1_x, double p1_y,
                           double p2_x, double p2_y, double p3_x, double p3_y, double *i_x, double *i_y)
{
    double s1_x, s1_y, s2_x, s2_y;
    s1_x = p1_x - p0_x;     s1_y = p1_y - p0_y;
    s2_x = p3_x - p2_x;     s2_y = p3_y - p2_y;

    double s, t;
    s = (-s1_y * (p0_x - p2_x) + s1_x * (p0_y - p2_y)) / (-s2_x * s1_y + s1_x * s2_y);
    t = ( s2_x * (p0_y - p2_y) - s2_y * (p0_x - p2_x)) / (-s2_x * s1_y + s1_x * s2_y);

    if (s >= 0 && s <= 1 && t >= 0 && t <= 1)
    {
        // Collision detected
        if (i_x != nullptr)
            *i_x = p0_x + (t * s1_x);
        if (i_y != nullptr)
            *i_y = p0_y + (t * s1_y);
        return true;
    }

    return false; // No collision
}


OneDConstraint Constraints::getTrackConstraints(const ArcLengthSpline &track,const State &x) const
{
    // given arc length s and the track -> compute linearized track constraints
    const double s = x.s;

    // X-Y point of the center line
    const Eigen::Vector2d pos_center = track.getPostion(s);
    const Eigen::Vector2d d_center   = track.getDerivative(s);
    // Tangent of center line at s
    const Eigen::Vector2d tan_center = {-d_center(1),d_center(0)};

    // inner and outer track boundary given left and right width of track
    // TODO make R_out and R_in dependent on s
    Eigen::Vector2d pos_outer = pos_center + param_.r_out*tan_center;
    Eigen::Vector2d pos_inner = pos_center - param_.r_in*tan_center;

    Eigen::MatrixXd &X_obs = *track.X_obs;
    Eigen::MatrixXd &Y_obs = *track.Y_obs;

    double i_x = 0;
    double i_y = 0;
    double *pi_x = &i_x;
    double *pi_y = &i_y;

    Eigen::Vector2d pos_inner_vec = pos_inner;
    Eigen::Vector2d pos_outer_vec = pos_outer;

    for (uint32_t i = 0; i < track.X_obs->rows(); i++)
    {
        uint32_t line_segment_count = track.X_obs->cols() - 1;
        for (uint32_t j = 0; j < line_segment_count; j++)
        {
            // d^2 = (x2 - x1)^2 + (y2 - y1)^2
            auto minDistance = [&pos_center](const Eigen::Vector2d &v1, const Eigen::Vector2d &v2) {
                return pow((v1(0) - pos_center[0]), 2) +
                       pow((v1(1) - pos_center[1]) , 2)
                       < pow((v2(0) - pos_center[0]), 2) +
                         pow((v2(1) - pos_center[1]), 2);
            };

            if (get_line_intersection(X_obs(i, j), Y_obs(i, j),
                                      X_obs(i, j + 1), Y_obs(i, j + 1),
                                      pos_center[0], pos_center[1], pos_inner[0], pos_inner[1],
                                      pi_x, pi_y))
            {
                Eigen::Vector2d tmp(*pi_x, *pi_y);
                pos_inner_vec = std::min(pos_inner_vec, tmp, minDistance);
            }

            if (get_line_intersection(X_obs(i, j), Y_obs(i, j),
                                      X_obs(i, j + 1), Y_obs(i, j + 1),
                                      pos_center[0], pos_center[1], pos_outer[0], pos_outer[1],
                                      pi_x, pi_y))
            {
                Eigen::Vector2d tmp(*pi_x, *pi_y);
                pos_outer_vec = std::min(pos_outer, tmp, minDistance);
            }
        }
    }

    // set pos_inner and pos_outer
    // check if stuck inside obstacle
    if (pos_inner != pos_inner_vec && pos_outer != pos_outer_vec) {
        // pos_inner = pos_outer_vec; // pass left
        pos_outer = pos_inner_vec; // pass right
    } else {
        pos_inner = pos_inner_vec;
        pos_outer = pos_outer_vec;
    }
    // add buffer boundaries
    pos_inner = pos_inner + 0.025*tan_center;
    pos_outer = pos_outer - 0.025*tan_center;


    // Define track Jacobian as Perpendicular vector
    C_i_MPC C_track_constraint = C_i_MPC::Zero();
    C_track_constraint(0,0) = tan_center(0);
    C_track_constraint(0,1) = tan_center(1);
    // Compute bounds
    const double track_constraint_lower = tan_center(0)*pos_inner(0) + tan_center(1)*pos_inner(1);
    const double track_constraint_upper = tan_center(0)*pos_outer(0) + tan_center(1)*pos_outer(1);

    return {C_track_constraint,track_constraint_lower,track_constraint_upper, pos_outer, pos_inner};
}

OneDConstraint Constraints::getTireConstraintRear(const State &x) const
{
    // compute tire friction elipse constraints
    // (param_.E_long*Frx)^2 + Fry^2 <= (param_.E_eps*F_max)^2
    const StateVector x_vec = stateToVector(x);
    const TireForces f_rear = model_.getForceRear(x);
    const NormalForces f_normal = model_.getForceNormal(x);

    // compute tire friction constraint jacobean
    const C_i_MPC C_tire_constraint = getTireConstraintRearJac(x);

    // compute zero order term and max force
//    const double tireCon0 = std::sqrt(std::pow(param_.e_long*f_rear.F_x,2) + std::pow(f_rear.F_y,2)); //zero order term
    const double tireCon0 = std::pow(param_.e_long*f_rear.F_x/f_normal.F_N_rear,2) + std::pow(f_rear.F_y/f_normal.F_N_rear,2); //zero order term
    const double maxForce = std::pow(param_.e_eps*param_.Dr/f_normal.F_N_rear,2);//param_.e_eps*param_.Dr;// //max allowed force

    // set bounds given linearized constraint
    // 0 <= 'Jac TC' (x - x0) + TC(x0) <= F_max
    const double tire_constraint_lower = C_tire_constraint*x_vec-tireCon0;
    const double tire_constraint_upper = maxForce+C_tire_constraint*x_vec-tireCon0;

    return {C_tire_constraint,tire_constraint_lower,tire_constraint_upper};
}

C_i_MPC Constraints::getTireConstraintRearJac(const State &x) const
{
    // compute Jacobean of the tire constraints
    const TireForces f_rear = model_.getForceRear(x);
    const TireForcesDerivatives df_rear = model_.getForceRearDerivatives(x);
    const NormalForces f_normal = model_.getForceNormal(x);

//    const double TC = 2.0*std::sqrt(std::pow(param_.e_long*f_rear.F_x,2) + std::pow(f_rear.F_y,2));
//
//    // Tire constraint derivatives
//    // TC = (param_.E_long*Frx)^2 + Fry^2
//    const double dTC_dvx = (2.0*param_.e_long*f_rear.F_x*df_rear.dF_x_vx + 2.0*f_rear.F_y*df_rear.dF_y_vx)/TC;
//    const double dTC_dvy = (2.0*f_rear.F_y*df_rear.dF_y_vy)/TC;
//    const double dTC_dr  = (2.0*f_rear.F_y*df_rear.dF_y_r)/TC;
//    const double dTC_dD  = (2.0*param_.e_long*f_rear.F_x*df_rear.dF_x_D)/TC;

    const double TC = std::pow(param_.e_long*f_rear.F_x/f_normal.F_N_rear,2) + std::pow(f_rear.F_y/f_normal.F_N_rear,2);

    // Tire constraint derivatives
    // TC = (param_.e_long*Frx)^2 + Fry^2
    const double dTC_dvx = (2.0*param_.e_long*f_rear.F_x/f_normal.F_N_rear*param_.e_long*df_rear.dF_x_vx/f_normal.F_N_rear +
                            2.0*f_rear.F_y/f_normal.F_N_rear*df_rear.dF_y_vx/f_normal.F_N_rear);
    const double dTC_dvy = (2.0*f_rear.F_y/f_normal.F_N_rear*df_rear.dF_y_vy/f_normal.F_N_rear);
    const double dTC_dr  = (2.0*f_rear.F_y/f_normal.F_N_rear*df_rear.dF_y_r/f_normal.F_N_rear);
    const double dTC_dD  = (2.0*param_.e_long*f_rear.F_x/f_normal.F_N_rear*param_.e_long*df_rear.dF_x_D/f_normal.F_N_rear);

    // Copy partial derivatives in jacobean matrix
    C_i_MPC Jac_tireCon = C_i_MPC::Zero();
    Jac_tireCon(si_index.vx) = dTC_dvx;
    Jac_tireCon(si_index.vy) = dTC_dvy;
    Jac_tireCon(si_index.r)  = dTC_dr;
    Jac_tireCon(si_index.D)  = dTC_dD;

    return Jac_tireCon;
}

OneDConstraint Constraints::getAlphaConstraintFront(const State &x) const
{
    // compute linearized slip angle constraints
    // -alpha_max <= alpha_f <= alpha_max
    const StateVector x_vec = stateToVector(x);
    const double alpha_f = model_.getSlipAngleFront(x);
    // compute the jacobean of alpha_f
    const C_i_MPC C_alpha_constraint = getAlphaConstraintFrontJac(x);
    // compute the bounds given the Tylor series expansion
    const double alpha_constraint_lower = -param_.max_alpha-alpha_f+C_alpha_constraint*x_vec;
    const double alpha_constraint_upper =  param_.max_alpha-alpha_f+C_alpha_constraint*x_vec;

    return {C_alpha_constraint,alpha_constraint_lower,alpha_constraint_upper};
}

C_i_MPC Constraints::getAlphaConstraintFrontJac(const State &x) const
{
    // compute the alpha_f jacobian
    const double vx     = x.vx;
    const double vy     = x.vy;
    const double r      = x.r;
    const double delta  = x.delta;

    C_i_MPC Jac_alphaCon;
    // alpha_f = -atan(vy+r*param_.lf/vx) + delta;
    // compute partial derivatives
    const double dalpha_f_dvx    = (vy+r*param_.lf)/(std::pow(vy+r*param_.lf,2)+std::pow(vx,2));
    const double dalpha_f_dvy    = -vx/(std::pow(vy+r*param_.lf,2)+std::pow(vx,2));
    const double dalpha_f_dr     = -(vx*param_.lf)/(std::pow(vy+r*param_.lf,2)+std::pow(vx,2));
    const double dalpha_f_ddelta = 1.0;

    // Copy partial derivatives in jacobean matrix
    Jac_alphaCon.setZero();
    Jac_alphaCon(si_index.vx)   = dalpha_f_dvx;
    Jac_alphaCon(si_index.vy)   = dalpha_f_dvy;
    Jac_alphaCon(si_index.r)    = dalpha_f_dr;
    Jac_alphaCon(si_index.delta)= dalpha_f_ddelta;

    return Jac_alphaCon;

}

ConstrainsMatrix Constraints::getConstraints(const ArcLengthSpline &track,const State &x,const Input &u) const
{
    // compute all the polytopic state constraints
    // compute the three constraints

    ConstrainsMatrix constrains_matrix;
    const OneDConstraint track_constraints = getTrackConstraints(track,x);
    const OneDConstraint tire_constraints_rear = getTireConstraintRear(x);
    const OneDConstraint alpha_constraints_front = getAlphaConstraintFront(x);

    C_MPC C_constrains_matrix;
    d_MPC dl_constrains_matrix;
    d_MPC du_constrains_matrix;

    C_constrains_matrix.row(si_index.con_track) = track_constraints.C_i;
    dl_constrains_matrix(si_index.con_track) = track_constraints.dl_i;
    du_constrains_matrix(si_index.con_track) = track_constraints.du_i;

    C_constrains_matrix.row(si_index.con_tire) = tire_constraints_rear.C_i;
    dl_constrains_matrix(si_index.con_tire) = tire_constraints_rear.dl_i;
    du_constrains_matrix(si_index.con_tire) = tire_constraints_rear.du_i;

    C_constrains_matrix.row(si_index.con_alpha) = alpha_constraints_front.C_i;
    dl_constrains_matrix(si_index.con_alpha) = alpha_constraints_front.dl_i;
    du_constrains_matrix(si_index.con_alpha) = alpha_constraints_front.du_i;

    return {C_constrains_matrix,D_MPC::Zero(),dl_constrains_matrix,du_constrains_matrix,
            track_constraints.outer, track_constraints.inner};
}
}