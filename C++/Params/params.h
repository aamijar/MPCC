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

#ifndef MPCC_PARAMS_H
#define MPCC_PARAMS_H

#include <iostream>
#include <fstream>
#include <vector>
#include <nlohmann/json.hpp>
namespace mpcc{
//used namespace
using json = nlohmann::json;

class Param{
public:
    double Cm1;
    double Cm2;
    double Cr0;
    double Cr2;

    double Br;
    double Cr;
    double Dr;

    double Bf;
    double Cf;
    double Df;

    double m;
    double Iz;
    double lf;
    double lr;

    double r_in;
    double r_out;

    double e_long;
    double e_eps;

    double max_alpha;

    double initial_velocity;
    double sqp_mixing;
    double s_trust_region;

    Param();
};

class CostParam{
public:
    double q_c;
    double q_l;
    double q_vs;

    double q_r;

    double q_beta;

    double r_D;
    double r_delta;
    double r_vs;

    double r_dD;
    double r_dDelta;
    double r_dVs;

    double q_c_N_mult;
    double q_r_N_mult;

    double sc_quad_track;
    double sc_quad_tire;
    double sc_quad_alpha;

    double sc_lin_track;
    double sc_lin_tire;
    double sc_lin_alpha;

    CostParam();
};

class BoundsParam{
public:
    struct LowerStateBounds{
        double X_l;
        double Y_l;
        double phi_l;
        double vx_l;
        double vy_l;
        double r_l;
        double s_l;
        double D_l;
        double delta_l;
        double vs_l;
    };
    struct UpperStateBounds{
        double X_u;
        double Y_u;
        double phi_u;
        double vx_u;
        double vy_u;
        double r_u;
        double s_u;
        double D_u;
        double delta_u;
        double vs_u;
    };
    struct LowerInputBounds{
        double dD_l;
        double dDelta_l;
        double dVs_l;
    };
    struct UpperInputBounds{
        double dD_u;
        double dDelta_u;
        double dVs_u;
    };

    LowerStateBounds lower_state_bounds;
    UpperStateBounds upper_state_bounds;

    LowerInputBounds lower_input_bounds;
    UpperInputBounds upper_input_bounds;

    BoundsParam();
};
}
#endif //MPCC_PARAMS_H