//
// Created by muhammed on 8/13/20.
//
//#include "stdafx.h"
#include "CartilageGnRMaterialPoint.h"
#include <FECore/FEModel.h>
#include <iostream>
#include <cmath>
#include <vector>

using namespace std;



//-------------------Rate Parameters-------------------------------------------------------------//
// Rate parameters for living chondrocytes
double r_c1 = 0.045;
double r_c2 = 0.4;
double r_c3 = r_c1;
double r_c4 = 0.01;
double r_c5 = 0.1;

// Rate parameters for necrotic chondrocytes
double r_nc1 = 1.0;
// double r_nc1 = 0.045;

// Rate parameters for type II collagen
double r_co1 = 0.00064;
double r_co2 = 0.1;
double r_co3 = r_co1;
double r_co5 = 1.0;


// Rate parameters for type PG
double r_pg1 = 0.0023;
double r_pg2 = 0.1;
double r_pg3 = r_pg1;


// Rate parameters for concentraion of collagenase (MMPs)
double r_ca4 = 4.15;
double r_ca2 = 1.0;
double r_ca3 = 100.0;
double r_ca5 = 0.1 * r_ca4;
double r_ca6 = 600.0;
double r_ca1 = r_ca4 + r_ca5;


// Rate parameters for concentraion of aggrecanase (ADAMTS)
double r_ag4 = 138.62;
double r_ag2 = 1.0;
double r_ag3 = 100.0;
double r_ag5 = r_ca5;
double r_ag1 = r_ag4 + r_ag5;
double r_ag6 = 5000.0;

// Rate parameters for concentraion of TIMPs
double r_i3 = 138.62;
double r_i2 = 0.1;
double r_i5 = r_ca5;
double r_i6 = r_ag5;
double r_i1 = r_i3 + r_i5 + r_i6;
double r_i4 = 100.0;


// Rate parameters for concentraion of latent growth factors
double r_lbeta3 = 50.0;
double r_lbeta2 = 5000.0;
double r_lbeta1 = r_lbeta3;
double r_lbeta4 = 0.5;

// Rate parameters for concentraion of active growth factors
double r_beta1 = r_lbeta4;
double r_beta2 = 50.0;

// Rate parameters for concentraion of latent pro-inflammatory cytokines
double r_lp5 = 50.0;
double r_lp2 = 5000.0;
double r_lp3 = 1.0;
double r_lp1 = r_lp5;
double r_lp6 = 1.0;
double r_lp4 = 5.0;

// Rate parameters for concentraion of active pro-inflammatory cytokines
double r_p1 = r_lp6;
double r_p2 = 50.0;

// Rate parameters for living hypertrophic cells
double r_hc1 = r_c4;
double r_hc2 = 0.45;
double r_hc3 = 0.0;

// Rate parameters for suramin
double r_sm1 = 1.438;
//----------------------------------------------------------------------------------------------------//

double f(double *a) {


    return 3*a[1] + 2*a[0] - 1;
};

double g(double b, double c) {
    return 5*b + 6*c - 1;
};


//! Defining the previous values
// prev[0] = n_c_prev,              // prev[1] = n_nc_prev,
// prev[2] = n_hc_prev,             // prev[3] = m_co_fn_prev,
// prev[4] = m_co_dm_prev,          // prev[5] = m_co_prev,
// prev[6] = m_pg_prev,             // prev[7] = c_ca_prev,
// prev[8] = c_ag_prev,             // prev[9] = c_i_prev,
// prev[10] = c_lbeta_prev,         // prev[11] = c_beta_prev,
// prev[12] = c_lp_prev,            // prev[13] = c_p_prev,
// prev[14] = c_sm_prev,            // prev[15] = f_sigmash,
// prev[16] = f_sigma1

//------------------Evolution Functions and Derivatives for Newton-Raphson----------------//

double n_c_func(double N_C, double dt, double *prev){
    return prev[0] + dt*((r_c1 + r_c2*prev[11])*N_C - (r_c3 + r_c4*prev[11] +r_c5*prev[13])*N_C) - N_C;
}

double n_c_deriv(double N_C, double dt, double *prev){
    return dt*((r_c1 + r_c2*prev[11]) - (r_c3 + r_c4*prev[11] +r_c5*prev[13])) - 1.0;
}

double n_nc_func(double N_NC, double dt, double *prev){
    return prev[1] + dt*(- r_nc1*N_NC) - N_NC;
}

double n_nc_deriv(double N_NC, double dt, double *prev){
    return dt*(- r_nc1) - 1.0;
}

double m_co_fn_func(double M_CO_FUN, double dt, double *prev){
    return prev[3] + dt*((r_co1 + r_co2*prev[11])*prev[0] - r_co3*prev[7]*M_CO_FUN) - M_CO_FUN;
}

double m_co_fn_deriv(double M_CO_FUN, double dt, double *prev){
    return dt*( - r_co3*prev[7]) - 1.0;
}

double m_co_dm_func(double M_CO_DM, double dt, double *prev){
    return prev[4] + dt*( - r_co5*prev[7]*M_CO_DM) - M_CO_DM;
}

double m_co_dm_deriv(double M_CO_DM, double dt, double *prev){
    return dt*( - r_co5*prev[7]) - 1.0;
}

double m_pg_func(double M_PG, double dt, double *prev){
    return prev[6] + dt*((r_pg1 + r_pg2*prev[11] )*prev[0] - r_pg3*prev[8]*M_PG) - M_PG;
}

double m_pg_deriv(double M_PG, double dt, double *prev){
    return dt*( - r_pg3*prev[8]) - 1.0;
}

double c_ca_func(double C_CA, double dt, double *prev){
    return prev[7] + dt*((((r_ca1 + r_ca6*prev[13])/ (1.0 + r_ca2*prev[11])))*prev[0] + r_ca3*prev[2]
            - (r_ca4 + r_ca5*prev[9])*C_CA) - C_CA;
}

double c_ca_deriv(double C_CA, double dt, double *prev){
    return dt*( - (r_ca4 + r_ca5*prev[9])) - 1.0;
}

double c_ag_func(double C_AG, double dt, double *prev){
    return prev[8] + dt*((((r_ag1 + r_ag6*prev[13])/ (1.0 + r_ag2*prev[11])))*prev[0] + r_ag3*prev[2]
                         - (r_ag4 + r_ag5*prev[9])*C_AG) - C_AG;
}

double c_ag_deriv(double C_AG, double dt, double *prev){
    return dt*( - (r_ag4 + r_ag5*prev[9])) - 1.0;
}

double c_i_func(double C_I, double dt, double *prev){
    return prev[9] + dt*((r_i1 + r_i2*prev[11])*prev[0]
            - (r_i3*prev[0]/(1.0 + r_i4*prev[14]) + r_i5*prev[7] + r_i6*prev[8])*C_I) - C_I;
}

double c_i_deriv(double C_I, double dt, double *prev){
    return dt*( - (r_i3*prev[0]/(1.0 + r_i4*prev[14]) + r_i5*prev[7] + r_i6*prev[8])) - 1.0;
}

double c_sm_func(double C_SM, double dt, double *prev){
    return prev[14] + dt*( - r_sm1*C_SM) - C_SM;
}

double c_sm_deriv(double C_SM, double dt, double *prev){
    return dt*( - r_sm1) - 1.0;
}

double c_lp_func(double C_LP, double dt, double *prev){
    return prev[12] + dt*((((r_lp1+r_lp2*prev[13])/ (1+ r_lp3*prev[11])))*prev[0]*
                        (1.0 + r_lp4*prev[1])
                - (r_lp5 + r_lp6*prev[16])*C_LP) - C_LP;
}

double c_lp_deriv(double C_LP, double dt, double *prev){
    return dt*( - (r_lp5 + r_lp6*prev[16])) - 1.0;
}

double c_p_func(double C_P, double dt, double *prev){
    return prev[13] + dt*(r_p1*prev[16]*prev[12] - r_p2*C_P) - C_P;
}

double c_p_deriv(double C_P, double dt, double *prev){
    return dt*( - r_p2) - 1.0;
}

double c_lbeta_func(double C_LBETA, double dt, double *prev){
    return prev[10] + dt*( (r_lbeta1+r_lbeta2*prev[11])*prev[0]
            - (r_lbeta3 + r_lbeta4*prev[15])*C_LBETA) - C_LBETA;
}

double c_lbeta_deriv(double C_LBETA, double dt, double *prev){
    return dt*( - (r_lbeta3 + r_lbeta4*prev[15])) - 1.0;
}

double c_beta_func(double C_BETA, double dt, double *prev){
    return prev[11] + dt*(r_beta1*prev[15]*prev[10] - r_beta2*C_BETA) - C_BETA;
}

double c_beta_deriv(double C_BETA, double dt, double *prev){
    return dt*( - r_beta2) - 1.0;
}

double n_hc_func(double N_HC, double dt, double *prev){
    return prev[2] + dt*(r_hc1*prev[11]*prev[0] - r_hc2*N_HC) - N_HC;
}

double n_hc_deriv(double N_HC, double dt, double *prev){
    return dt*( - r_hc2) - 1.0;
}

//----------------------------------------------------------------------------------------//

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

#define MIN(a,b) ((a)<(b) ? (a) : (b))
#define MAX(a,b) ((a)>(b) ? (a) : (b))

//! Create a shallow copy of the material point data
FEMaterialPointData* CartilageGnRMaterialPoint::Copy()
{
    CartilageGnRMaterialPoint* pt = new CartilageGnRMaterialPoint(*this);
    if (m_pNext) pt->m_pNext = m_pNext->Copy();
    return pt;
}

//! Initializes material point data.
void CartilageGnRMaterialPoint::Init()
{
    FEMaterialPointData::Init();

    // intialize data to one or zero
    n_c = 1.0;
    n_nc = 0.0;
    n_hc = 0.0;
    m_co_fn = 1.0;
    m_co_dm = 0.0;
    m_co = 1.0;
    m_pg = 1.0;
//    m_pg = 0.5;
    c_ca = 1.0;
    c_ag = 1.0;
    c_i = 1.0;
    c_sm = 0.0;
    c_lp = 1.0;
    c_p = 0.0;
    c_lbeta = 1.0;
    c_beta = 0.0;
    nvolu = 1.0;

    f_sigmash = 0.0;
    f_sigma1 = 0.0;

    time_it = 0.0;
    time_check = 0.0;
    iterator_count = 10000;     // Assigning a larger value to avoid
                                // unwanted iteration

}

//! Update material point data.
void CartilageGnRMaterialPoint::Update(const FETimeInfo& timeInfo)
{
    // don't forget to call the base class
    FEMaterialPointData::Update(timeInfo);

}


void CartilageGnRMaterialPoint::Serialize(DumpStream& ar)
{
    if (ar.IsSaving())
    {
        ar << n_c << n_nc << n_hc << m_co_fn << m_co_dm << m_co << m_pg << c_ca << c_ag << c_i << c_sm << c_lp << c_p << c_lbeta << c_beta << nvolu;
    }
    else
    {
        ar >> n_c >> n_nc >> n_hc >> m_co_fn >> m_co_dm >> m_co >> m_pg >> c_ca >> c_ag >> c_i >> c_sm >> c_lp >> c_p >> c_lbeta >> c_beta >> nvolu;
    }

    FEMaterialPointData::Serialize(ar);

}

double CartilageGnRMaterialPoint::Volume_Return(FEMaterialPoint& pt, double v_f, float dt, double time_step)
{
    FEElasticMaterialPoint& ppt = *pt.ExtractData<FEElasticMaterialPoint>();
    CartilageGnRMaterialPoint& mpt = *pt.ExtractData<CartilageGnRMaterialPoint>();

    //<! only at a certain time point
    for (int i = 0; i < 251; ++i) {
//        if ( abs(dt - ((float)i*0.3 + 0.1) )<= 0.001) {
        if ( abs(dt - ((float)i*2.0 - 0.75) )<= 0.001) {
            if (mpt.iterator_count == i) {

                if (i == 104) {
                    int a = 1;
                }
                break;
            }
            mpt.iterator_count = i;

            time_step = 0.1;
            double n_c_prev = mpt.n_c;
            double n_nc_prev = mpt.n_nc;
            double n_hc_prev = mpt.n_hc;
            double m_co_fn_prev = mpt.m_co_fn;
            double m_co_dm_prev = mpt.m_co_dm;
            double m_co_prev = mpt.m_co;
            double m_pg_prev = mpt.m_pg;
            double c_ca_prev = mpt.c_ca;
            double c_ag_prev = mpt.c_ag;
            double c_i_prev = mpt.c_i;
            double c_sm_prev = mpt.c_sm;
            double c_lp_prev = mpt.c_lp;
            double c_p_prev = mpt.c_p;
            double c_lbeta_prev = mpt.c_lbeta;
            double c_beta_prev = mpt.c_beta;
            double f_sigmash_prev = mpt.f_sigmash;
            double f_sigma1_prev = mpt.f_sigma1;
            double nvolu_prev = mpt.nvolu;

            double prev[] = {
                          n_c_prev,
                          n_nc_prev,
                          n_hc_prev,
                          m_co_fn_prev,
                          m_co_dm_prev,
                          m_co_prev,
                          m_pg_prev,
                          c_ca_prev,
                          c_ag_prev,
                          c_i_prev,
                          c_lbeta_prev,
                          c_beta_prev,
                          c_lp_prev,
                          c_p_prev,
                          c_sm_prev,
                          f_sigmash_prev,
                          f_sigma1_prev};




            //! Defining the previous values
            // prev[0] = n_c_prev,              // prev[1] = n_nc_prev,
            // prev[2] = n_hc_prev,             // prev[3] = m_co_fn_prev,
            // prev[4] = m_co_dm_prev,          // prev[5] = m_co_prev,
            // prev[6] = m_pg_prev,             // prev[7] = c_ca_prev,
            // prev[8] = c_ag_prev,             // prev[9] = c_i_prev,
            // prev[10] = c_lbeta_prev,         // prev[11] = c_beta_prev,
            // prev[12] = c_lp_prev,            // prev[13] = c_p_prev,
            // prev[14] = c_sm_prev,            // prev[15] = f_sigmash,
            // prev[16] = f_sigma1


//            double tolerance =time_step*time_step;
            double tolerance = 1.0e-9;


            // Solution for normalized number of normal chondrocytes
            mpt.n_c = NewtonRaphsonManual(&n_c_func, &n_c_deriv, prev[0], time_step, tolerance, prev);;

            // Solution for normalized number of necrotic chondrocytes
            mpt.n_nc = NewtonRaphsonManual(&n_nc_func, &n_nc_deriv, prev[1], time_step, tolerance, prev);;

            // Solution for normalized quantity of collagen II
            mpt.m_co_fn = NewtonRaphsonManual(&m_co_fn_func, &m_co_fn_deriv, prev[3], time_step, tolerance, prev);

            mpt.m_co_dm = NewtonRaphsonManual(&m_co_dm_func, &m_co_dm_deriv, prev[4], time_step, tolerance, prev);

            mpt.m_co = mpt.m_co_fn + mpt.m_co_dm;

            // Solution for normalized quantity of proteoglycan
            mpt.m_pg = NewtonRaphsonManual(&m_pg_func, &m_pg_deriv, prev[6], time_step, tolerance, prev);

            // Solution for normalized concentration of collagenase (MMPs)
            mpt.c_ca = NewtonRaphsonManual(&c_ca_func, &c_ca_deriv, prev[7], time_step, tolerance, prev);

            // Solution for normalized concentration of aggrecanase (MMPs)
            mpt.c_ag = NewtonRaphsonManual(&c_ag_func, &c_ag_deriv, prev[8], time_step, tolerance, prev);

            // Solution for normalized concentration of TIMPs
            mpt.c_i = NewtonRaphsonManual(&c_i_func, &c_i_deriv, prev[9], time_step, tolerance, prev);

            // Solution for normalized concentration of Suramin
            mpt.c_sm = NewtonRaphsonManual(&c_sm_func, &c_sm_deriv, prev[14], time_step, tolerance, prev);

            if (mpt.iterator_count == 10)
            {
                int a = 0;
            }
            // Solution for normalized concentration of latent pro-inflammatory
            // cytokines (IL-1beta, TNF-alpha, etc.)
            mpt.c_lp = NewtonRaphsonManual(&c_lp_func, &c_lp_deriv, prev[12], time_step, tolerance, prev);

            // Solution for normalized concentration of active pro-inflammatory
            // cytokines (IL-1beta, TNF-alpha, etc.)
            mpt.c_p = NewtonRaphsonManual(&c_p_func, &c_p_deriv, prev[13], time_step, tolerance, prev);

            // Solution for normalized concentration of latent growth factors
            // (TGF-beta, IGF-1, etc.)
            mpt.c_lbeta = NewtonRaphsonManual(&c_lbeta_func, &c_lbeta_deriv, prev[10], time_step, tolerance, prev);

            // Solution for normalized concentration of active growth factors
            // (TGF-beta, IGF-1, etc.)
            mpt.c_beta = NewtonRaphsonManual(&c_beta_func, &c_beta_deriv, prev[11], time_step, tolerance, prev);

            // Solution for normalized number of hypertrophic cells
            mpt.n_hc = NewtonRaphsonManual(&n_hc_func, &n_hc_deriv, prev[2], time_step, tolerance, prev);


        }
    }



    mpt.nvolu = mpt.m_pg*(1-v_f) + (mpt.m_co_fn + mpt.m_co_dm)*v_f;
    return mpt.nvolu;
}

double CartilageGnRMaterialPoint::NewtonRaphsonManual(double (*func)(double, double, double *), double (*deriv)(double, double, double *), double initial, double dt, double tol, double *prev)
{
    ynew = initial;
    res = - func(ynew, dt, prev)/deriv(ynew, dt, prev);



    while (abs(res) > tol) {
        ynew = ynew + res;
        res = - func(ynew, dt, prev)/deriv(ynew, dt, prev);


    }

    ynew = ynew + res;

    return ynew ;


}

double CartilageGnRMaterialPoint::sigmash_stimuli(FEMaterialPoint& pt, float dt)
{
    FEElasticMaterialPoint& ppt = *pt.ExtractData<FEElasticMaterialPoint>();
    CartilageGnRMaterialPoint& mpt = *pt.ExtractData<CartilageGnRMaterialPoint>();

    double sigma_sh_stimuli;

    if (mpt.time_check == dt) {

        sigma_sh_stimuli = mpt.f_sigmash;

    }
    else {

        double homeostatic_setpoint_t0 = 1.88;
        double current_stimuli = 1.6;

        double simulation_time = 24;
        double steps_per_month = 30;
        double NSTEPS = simulation_time*steps_per_month;
        double time_history    = 3; // averaging period in months
        double time_history_steps = time_history*steps_per_month;
        double ths = time_history_steps;
        double time_delay      = 2.0; // months - this is a delay before homeostatic adaption
        double k=20;
        double width=0.25;

        double f_stimulus_value = 0.63;

        int total_time_steps = NSTEPS+time_history_steps;
        double time[total_time_steps];
        double cell_record_stimulus[total_time_steps];
        double homeostatic_a[total_time_steps];
        double f_stimulus[total_time_steps];


        for (int i = 0; i < total_time_steps; ++i) {

            time[i] = (i - time_history_steps)/steps_per_month;

            if (time[i] <= time_delay){
                cell_record_stimulus[i] = homeostatic_setpoint_t0;
            }
            else {
                cell_record_stimulus[i] = current_stimuli;
            }

            if (time[i] <= time_delay){
                homeostatic_a[i] = homeostatic_setpoint_t0;
            }
            else {

                double sum_cell_stimuli = 0.0;
                for (int j = (i-ths-1); j < (i-1); ++j) {
                    sum_cell_stimuli = sum_cell_stimuli + cell_record_stimulus[j];
                }

                homeostatic_a[i] = sum_cell_stimuli / time_history_steps;
            }

            f_stimulus[i] = (f_stimulus_value - 0.0)/(1.0
                                                      + exp(k*(current_stimuli - (homeostatic_a[i] - 0.5*width) )));

            if (f_stimulus[i] > 0.99){
                f_stimulus[i] = 1.0;
            }

            if (f_stimulus[i] < 0.01){
                f_stimulus[i] = 0.0;
            }

        }

        std::vector<double> xData;
        std::vector<double> yData;


//    vector<double> xData(begin(time), end(time));

//    std::copy(std::begin(time), std::end(time), std::back_inserter(xData));
//
        for (int i=0; i< total_time_steps; ++i) {
            xData.push_back(time[i]);
            yData.push_back(f_stimulus[i]);
        }

        time_it = time_it + 0.1;





        sigma_sh_stimuli = interpolate(xData, yData, time_it, false);

        mpt.time_check = dt;

    }




    return sigma_sh_stimuli;


}

double CartilageGnRMaterialPoint::interpolate(vector<double> &xData, vector<double> &yData, double x, bool extrapolate)
{
    int size = xData.size();

    int i = 0;                                                                  // find left end of interval for interpolation
    if ( x >= xData[size - 2] )                                                 // special case: beyond right end
    {
        i = size - 2;
    }
    else
    {
        while ( x > xData[i+1] ) i++;
    }
    double xL = xData[i], yL = yData[i], xR = xData[i+1], yR = yData[i+1];      // points on either side (unless beyond ends)
    if ( !extrapolate )                                                         // if beyond ends of array and not extrapolating
    {
        if ( x < xL ) yR = yL;
        if ( x > xR ) yL = yR;
    }

    double dydx = ( yR - yL ) / ( xR - xL );                                    // gradient

    return yL + dydx * ( x - xL );                                              // linear interpolation

}

double CartilageGnRMaterialPoint::Volume_Return_Test(FEMaterialPoint &pt, double v_f, float dt, double time_step)
{
    FEElasticMaterialPoint& ppt = *pt.ExtractData<FEElasticMaterialPoint>();
    CartilageGnRMaterialPoint& mpt = *pt.ExtractData<CartilageGnRMaterialPoint>();


    double time_const = 100;    //<! This time will control growth and remodeling
    //<! only at a certain time point



    if ( abs(dt - (float)1.5)<= 0.001) {

        // Solution for normalized quantity of collagen II
        mpt.m_co_fn = 0.95;
        mpt.m_co_dm = 0.0;
        mpt.m_co = mpt.m_co_fn + mpt.m_co_dm;

        // Solution for normalized quantity of proteoglycan
        mpt.m_pg = 0.95;

    }

    if ( abs(dt - (float)3.2)<= 0.001) {

        // Solution for normalized quantity of collagen II
        mpt.m_co_fn = 0.9;
        mpt.m_co_dm = 0.0;
        mpt.m_co = mpt.m_co_fn + mpt.m_co_dm;

        // Solution for normalized quantity of proteoglycan
        mpt.m_pg = 0.9;

    }


    if ( abs(dt - (float)5.2)<= 0.001) {

        // Solution for normalized quantity of collagen II
        mpt.m_co_fn = 0.85;
        mpt.m_co_dm = 0.0;
        mpt.m_co = mpt.m_co_fn + mpt.m_co_dm;

        // Solution for normalized quantity of proteoglycan
        mpt.m_pg = 0.85;

    }

    if ( abs(dt - (float)7.2)<= 0.001) {

        // Solution for normalized quantity of collagen II
        mpt.m_co_fn = 0.8;
        mpt.m_co_dm = 0.0;
        mpt.m_co = mpt.m_co_fn + mpt.m_co_dm;

        // Solution for normalized quantity of proteoglycan
        mpt.m_pg = 0.8;

    }



    mpt.nvolu = mpt.m_pg*(1-v_f) + (mpt.m_co_fn + mpt.m_co_dm)*v_f;
    return mpt.nvolu;
}