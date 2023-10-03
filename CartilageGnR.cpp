#include "stdafx.h"
#include <limits>
#include "CartilageGnR.h"
#include "FECore/FEModel.h"

#define MIN(a,b) ((a)<(b) ? (a) : (b))
#define MAX(a,b) ((a)>(b) ? (a) : (b))

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(CartilageGnR, FEElasticMaterial)
	ADD_PARAMETER(m_mu, FE_RANGE_GREATER(0.0), "mu");
    ADD_PARAMETER(m_K      , FE_RANGE_GREATER_OR_EQUAL(0.0), "kappa");
    ADD_PARAMETER(m_v, FE_RANGE_RIGHT_OPEN(-1.0, 0.5), "nu");
    ADD_PARAMETER(J_cp, FE_RANGE_GREATER   (0.0), "Jcp");
    ADD_PARAMETER(n_FIB      , FE_RANGE_GREATER_OR_EQUAL(0.0), "vf");
    ADD_PARAMETER(m_k1      , FE_RANGE_GREATER_OR_EQUAL(0.0), "k1");
    ADD_PARAMETER(m_k2      , FE_RANGE_GREATER_OR_EQUAL(0.0), "k2");
    ADD_PARAMETER(m_p      , "pressure_model");
    ADD_PARAMETER(vgnr      , "volumetric_growth");
// material properties
    ADD_PROPERTY(m_pFmat, "fibers", FEProperty::Optional);
//    ADD_PROPERTY(m_pBase, "solid");
    ADD_PROPERTY(m_pFDD, "distribution");
    ADD_PROPERTY(m_pFint, "scheme");

END_FECORE_CLASS();

//-----------------------------------------------------------------------------
CartilageGnR::CartilageGnR(FEModel* pfem) : FEElasticMaterial(pfem)
{
    m_pFmat = 0;
    m_pFDD = nullptr;
    m_pFint = nullptr;
    m_epsf = 0;
    m_beta = 2.0;
    n_FIB = 0.8;
    vgnr = 0;



}

//-----------------------------------------------------------------------------
CartilageGnR::~CartilageGnR() {}

//-----------------------------------------------------------------------------
bool CartilageGnR::Init()
{
    // initialize base class
    if (FEElasticMaterial::Init() == false) return false;

    // calculate the integrated fiber density
//    IntegrateFiberDensity();

    return true;
}

//-----------------------------------------------------------------------------
// returns a pointer to a new material point object
FEMaterialPointData* CartilageGnR::CreateMaterialPointData()
{
    FESolidMaterial* pme = GetBaseMaterial();
//    FEElasticMaterial* pme = GetBaseMaterial();
    FEMaterialPointData* ep = pme->CreateMaterialPointData();
//    FEMaterialPoint* mp = m_pFmat->CreateMaterialPointData();
//    mp->SetName(m_pFmat.GetName());

//    return new CartilageGnRMaterialPoint(mp);
    return new CartilageGnRMaterialPoint(ep);

//    return new CartilageGnRMaterialPoint(nullptr);
}

//-----------------------------------------------------------------------------
//! Calculate the stress
mat3ds CartilageGnR::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    CartilageGnRMaterialPoint& ppt = *mp.ExtractData<CartilageGnRMaterialPoint>();


	// get material parameters
	double c1 = m_mu/2.0;

//---------------------------------------------------------------------------

    mat3d F_prev = pt.m_F;
    double dt = GetFEModel()->GetTime().currentTime;
    double time_step = GetFEModel()->GetTime().timeIncrement;
//    double time_step = 0.1;

//    nvolume = ppt.Volume_Return(pt, n_FIB, dt, time_step);
    nvolume = ppt.Volume_Return(mp, n_FIB, dt, time_step);
    norm_col_fun = ppt.m_co_fn;
    norm_pg = ppt.m_pg;

    nuhat_col = (norm_col_fun*n_FIB)/nvolume;
    nuhat_pg = (norm_pg*(1.0-n_FIB))/nvolume;

    mat3dd I(1.0); // identity tensor

    double alpha = 1.0;
    double beta = nvolume - 1.0;

    vec3d n;
    n.x = 0; n.y = 0; n.z = 1;
    mat3d ndyadn = n & n;

    mat3d F_g = alpha*I + ndyadn*beta;

    mat3d F = F_prev * F_g.inverse();

    double J = F.det();

    double Jm23 = pow(J, -2.0/3.0);

    // calculate left Cauchy-Green tensor
    mat3ds B;
    //! Calculating isochoric Left Cauchy Green Tensor for anisotropic growth
    B.xx() = Jm23*(F[0][0]*F[0][0]+F[0][1]*F[0][1]+F[0][2]*F[0][2]); // = b[0][0]
    B.yy() = Jm23*(F[1][0]*F[1][0]+F[1][1]*F[1][1]+F[1][2]*F[1][2]); // = b[1][1]
    B.zz() = Jm23*(F[2][0]*F[2][0]+F[2][1]*F[2][1]+F[2][2]*F[2][2]); // = b[2][2]
    B.xy() = Jm23*(F[0][0]*F[1][0]+F[0][1]*F[1][1]+F[0][2]*F[1][2]); // = b[0][1]
    B.yz() = Jm23*(F[1][0]*F[2][0]+F[1][1]*F[2][1]+F[1][2]*F[2][2]); // = b[1][2]
    B.xz() = Jm23*(F[0][0]*F[2][0]+F[0][1]*F[2][1]+F[0][2]*F[2][2]); // = b[0][2]


//---------------------------------------------------------------------------

	double W1 = c1;

	mat3ds sigma_IM = B*(W1)*(2.0/J);

    mat3ds sigma_IM_dev = sigma_IM.dev();

    // Fiber Contribution Start

    // calculate stress
    mat3ds s_fiber; s_fiber.zero();


    // get the element's local coordinate system
    mat3d QT = GetLocalCS(mp).transpose();
    double IFD = IntegratedFiberDensity(mp);

    // obtain an integration point iterator
//    FEFiberIntegrationSchemeIterator* it = m_pFint->GetIterator(&pt);
    FEFiberIntegrationSchemeIterator* it = m_pFint->GetIterator(&mp);
    if (it->IsValid())
    {
        do
        {
            // get the global fiber direction
            vec3d& n0 = it->m_fiber;

            // convert to local coordinates
            vec3d n0a = QT*n0;

            // rotate to local configuration to evaluate ellipsoidally distributed material coefficients
            double R = m_pFDD->FiberDensity(mp, n0a) / IFD;

            // calculate the stress
            double wn = it->m_weight;
//            s_fiber += m_pFmat->FiberStress(pt, n0)*(R*wn);
            s_fiber += GnRDevFiberStress(F, n0)*(R*wn);
        }
        while (it->Next());
    }

    // don't forget to delete the iterator
    delete it;

    mat3ds sigma_FN_dev = s_fiber;

    // Fiber Contribution end


    //! pressure, i.e. first derivative of U(J)

    if (m_p == 1)
    {
        m_K = 2.0*m_mu*(1.0+m_v)  / (3.0*(1.0-2.0*m_v));
        UJ = m_K*(J-nvolume);
//        UJ = m_mu*(J-nvolume);
    } else
    {
        J_cp = J_cp * nvolume;
        double chi_cp = ((2.0*m_mu*m_v)/(1.0-2.0*m_v))                  //Lame parameter
                        *pow((1.0 + J_cp*(1.0 + pow(J_cp,2.0)/(1.0 - J_cp))),-1.0);


        UJ = chi_cp*(log(J)/J +
                               J_cp/J     + (1.0-J_cp)/(J_cp-2.0)*                                 // d_zeta/d_J
                                            (-1.0/(J_cp-J) - (J_cp-1.0)/(J*(J_cp-1.0)-J_cp))       // derivative
        );
//    - m_mu/J;

    }



    pt.m_F = F_prev ;

    mat3ds cauchy_stress = nuhat_pg*norm_pg*sigma_IM_dev + nuhat_col*norm_col_fun*sigma_FN_dev + mat3dd(UJ);

    if (vgnr == 1)
    {
        for (int i = 0; i < 251; ++i) {
            if (abs(dt - ((float) i * 2.0 + 20.6)) <= 0.001) {
                //!------------Calculating Principal Stresses-------------------------------//
                mat3ds dev_cauchy_stress = cauchy_stress.dev();
                double d[3];
                vec3d v[3];
                dev_cauchy_stress.eigen(d, v);

                // Finding maximum and minimum value TODO: Find an optimized way to do that
                double max_d = MAX(d[0], d[1]);
                max_d = MAX(max_d, d[2]);
                double min_d = MIN(d[0], d[1]);
                min_d = MIN(min_d, d[2]);

                double sigma_1 = max_d;
                double sigma_sh = 0.5 * abs(max_d - min_d);

//                ppt.f_sigmash = ppt.sigmash_stimuli(pt, dt);
                ppt.f_sigmash = ppt.sigmash_stimuli(mp, dt);
                ppt.f_sigma1 = ppt.f_sigmash;
            }
        }

    }

    if (vgnr == 2)
    {
        for (int i = 0; i < 251; ++i) {
            if (abs(dt - ((float) i * 2.0 + 20.5)) <= 0.001) {

                ppt.m_pg = 0.99;
                ppt.m_co_fn = 0.99;
            }
        }

    }





//	return 0.2*sigma_IM_dev + 0.8*sigma_FN_dev + mat3dd(UJ);
    return cauchy_stress;

}


//-----------------------------------------------------------------------------
//! Calculate the tangent
tens4ds CartilageGnR::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    CartilageGnRMaterialPoint& ppt = *mp.ExtractData<CartilageGnRMaterialPoint>();

	// get material parameters
	double c1 = m_mu/2.0;


//---------------------------------------------------------------------------

    mat3d F_prev = pt.m_F;
    double dt = GetFEModel()->GetTime().currentTime;
    double time_step = GetFEModel()->GetTime().timeIncrement;
//    double time_step = 0.1;

//    nvolume = ppt.Volume_Return_Test(pt,n_FIB,dt, time_step);
    nvolume = ppt.Volume_Return(mp, n_FIB, dt, time_step);
//    nvolume = ppt.Volume_Return(pt, n_FIB, dt, time_step);
    norm_col_fun = ppt.m_co_fn;
    norm_pg = ppt.m_pg;

    nuhat_col = (norm_col_fun*n_FIB)/nvolume;
    nuhat_pg = (norm_pg*(1.0-n_FIB))/nvolume;

    mat3dd I(1.0); // identity tensor


    double alpha = 1.0;
    double beta = nvolume - 1.0;

    vec3d n;
    n.x = 0; n.y = 0; n.z = 1;
    mat3d ndyadn = n & n;

    mat3d F_g = alpha*I + ndyadn*beta;

    mat3d F = F_prev * F_g.inverse();

    double J = F.det();

    double Jm23 = pow(J, -2.0/3.0);

    // calculate left Cauchy-Green tensor
    mat3ds B;
    //! Calculating Left Cauchy Green Tensor for anisotropic growth
    B.xx() = Jm23*(F[0][0]*F[0][0]+F[0][1]*F[0][1]+F[0][2]*F[0][2]); // = b[0][0]
    B.yy() = Jm23*(F[1][0]*F[1][0]+F[1][1]*F[1][1]+F[1][2]*F[1][2]); // = b[1][1]
    B.zz() = Jm23*(F[2][0]*F[2][0]+F[2][1]*F[2][1]+F[2][2]*F[2][2]); // = b[2][2]
    B.xy() = Jm23*(F[0][0]*F[1][0]+F[0][1]*F[1][1]+F[0][2]*F[1][2]); // = b[0][1]
    B.yz() = Jm23*(F[1][0]*F[2][0]+F[1][1]*F[2][1]+F[1][2]*F[2][2]); // = b[1][2]
    B.xz() = Jm23*(F[0][0]*F[2][0]+F[0][1]*F[2][1]+F[0][2]*F[2][2]); // = b[0][2]


//---------------------------------------------------------------------------

	double Ji = 1.0/J;

	// calculate deviatoric left Cauchy-Green tensor: B = F*Ft
//	mat3ds B = pt.DevLeftCauchyGreen();


	// Invariants of B (= invariants of C)
	double I1 = B.tr();

	// --- TODO: put strain energy derivatives here ---
	// Wi = dW/dIi
	double W1;
	W1 = c1;
	// ---

	// calculate dWdC:C
	double WC = W1*I1;


	mat3ds T = B*(W1);
    mat3ds T_dev = T.dev()*(2.0/J);

	// Identity tensor
//	mat3ds I(1,1,1,0,0,0);

	tens4ds IxI = dyad1s(I);
	tens4ds I4  = dyad4s(I);
	tens4ds BxB = dyad1s(B);
	tens4ds B4  = dyad4s(B);

    tens4ds c_IM_dev = dyad1s(T_dev, I)*(-4.0/3.0) + (I4 - IxI/3.0)*(4.0/3.0*Ji*WC);

//    tens4ds c_IM_dev = dyad1s(T_dev, I)*(-2.0/3.0) + (I4 - IxI/3.0)*(4.0/3.0*Ji*WC);

    if (m_p == 1)
    {
        m_K = 2.0*m_mu*(1.0+m_v)  / (3.0*(1.0-2.0*m_v));
        //! pressure, i.e. first derivative of U(J)
        UJ = m_K*(J-nvolume);
//        UJ = m_mu*(J-nvolume);

        //! second derivative of U(J)
        UJJ = m_K;
//        UJJ = m_mu;
    } else
    {
        J_cp = J_cp * nvolume;

        double chi_cp = ((2.0*m_mu*m_v)/(1.0-2.0*m_v))                  //Lame parameter
                        *pow((1.0 + J_cp*(1.0 + pow(J_cp,2.0)/(1.0 - J_cp))),-1.0);


        UJ = chi_cp*(log(J)/J +
                            J_cp/J     + (1.0-J_cp)/(J_cp-2.0)*                                 // d_zeta/d_J
                                         (-1.0/(J_cp-J) - (J_cp-1.0)/(J*(J_cp-1.0)-J_cp))       // derivative
        ) ;
//    - m_mu/J;


        UJJ = chi_cp*
                     ( pow(J,-2.0)*(1.0 - log(J) - J_cp)
                       + (1.0 - J_cp)/(J_cp - 2.0)*                          // d_2_zeta/d_J2
                         (- pow((J_cp - J),-2.0) +                   // 2nd
                          pow(((J_cp-1.0)/(J*(J_cp-1.0)-J_cp)),2.0)  // derivative
                         )
                     ) ;
//    + m_mu*pow(J,-2.0);
    }






    // Fiber Contribution start

    // get the element's local coordinate system
    mat3d QT = GetLocalCS(mp).transpose();
    double IFD = IntegratedFiberDensity(mp);

    // initialize stress tensor
    tens4ds tang_fiber;
    tang_fiber.zero();

//    FEFiberIntegrationSchemeIterator* it = m_pFint->GetIterator(&pt);
    FEFiberIntegrationSchemeIterator* it = m_pFint->GetIterator(&mp);
    if (it->IsValid())
    {
        do
        {
            // get the global fiber direction
            vec3d& n0e = it->m_fiber;

            // convert to local
            vec3d n0a = QT*n0e;

            // rotate to local configuration to evaluate ellipsoidally distributed material coefficients
            double R = m_pFDD->FiberDensity(mp, n0a) / IFD;

            // calculate the tangent
//            tang_fiber += m_pFmat->FiberTangent(mp, n0e)*(R*it->m_weight);
            tang_fiber += GnRDevFiberTangent(F, n0e)*(R*it->m_weight);

        }
        while (it->Next());
    }

    // don't forget to delete the iterator
    delete it;

    tens4ds c_FN_dev = tang_fiber;

//    tens4ds c = 0.2*c_IM_dev +
//                0.8*c_FN_dev +
//                (IxI - I4*2)*UJ + IxI*(UJJ*J);  // c_pressure + c_k
//

    tens4ds c = nuhat_pg*norm_pg*c_IM_dev +
                nuhat_col*norm_col_fun*c_FN_dev +
                (IxI - I4*2)*UJ + IxI*(UJJ*J);  // c_pressure + c_k

//    tens4ds c = dyad1s(T_dev, I)*(-2.0/3.0) + (I4 - IxI/3.0)*(4.0/3.0*Ji*WC)
//                + (IxI - I4*2)*UJ + IxI*(UJJ*J);

// MMR: I am confused whether it should be (-4.0/3.0) or (-2.0/3.0)
// in the first term. Using formulation it should be (-4.0/3.0), but
// in Mooney-Rivlin it is (-2.0/3.0)

    pt.m_F = F_prev ;

	return c;
}



//-----------------------------------------------------------------------------
mat3ds CartilageGnR::GnRDevFiberStress(const mat3d& Fe, const vec3d& n0)
{

    double J = Fe.det();

    double Jm23 = pow(J, -2.0/3.0);

    // calculate deviatoric right Cauchy-Green tensor
    mat3ds C;
    C.xx() = Jm23*(Fe[0][0]*Fe[0][0]+Fe[1][0]*Fe[1][0]+Fe[2][0]*Fe[2][0]); // = C[0][0]
    C.yy() = Jm23*(Fe[0][1]*Fe[0][1]+Fe[1][1]*Fe[1][1]+Fe[2][1]*Fe[2][1]); // = C[1][1]
    C.zz() = Jm23*(Fe[0][2]*Fe[0][2]+Fe[1][2]*Fe[1][2]+Fe[2][2]*Fe[2][2]); // = C[2][2]
    C.xy() = Jm23*(Fe[0][0]*Fe[0][1]+Fe[1][0]*Fe[1][1]+Fe[2][0]*Fe[2][1]); // = C[0][1]
    C.yz() = Jm23*(Fe[0][1]*Fe[0][2]+Fe[1][1]*Fe[1][2]+Fe[2][1]*Fe[2][2]); // = C[1][2]
    C.xz() = Jm23*(Fe[0][0]*Fe[0][2]+Fe[1][0]*Fe[1][2]+Fe[2][0]*Fe[2][2]); // = C[0][2]

    mat3ds s;

    // Calculate In = n0*C*n0
    double In_1 = n0*(C*n0) - 1.0;

    // only take fibers in tension into consideration
    const double eps = m_epsf* std::numeric_limits<double>::epsilon();
    if (In_1 >= eps)
    {
        // get the global spatial fiber direction in current configuration
        vec3d nt = Fe*n0;

        // calculate the outer product of nt
        mat3ds N = dyad(nt);

        // calculate strain energy derivative
        double Wl = m_k1*pow(In_1, m_beta - 1.0)*exp(m_k2*pow(In_1, m_beta));

        // calculate the fiber stress
        s = N*(2.0*Wl / J);

    }
    else
    {
        s.zero();
    }

    return s.dev();
}

//-----------------------------------------------------------------------------
tens4ds CartilageGnR::GnRDevFiberTangent(const mat3d& Fe, const vec3d& n0)
{
    double J = Fe.det();

    double Jm23 = pow(J, -2.0/3.0);

    // calculate deviatoric right Cauchy-Green tensor
    mat3ds C;
    C.xx() = Jm23*(Fe[0][0]*Fe[0][0]+Fe[1][0]*Fe[1][0]+Fe[2][0]*Fe[2][0]); // = C[0][0]
    C.yy() = Jm23*(Fe[0][1]*Fe[0][1]+Fe[1][1]*Fe[1][1]+Fe[2][1]*Fe[2][1]); // = C[1][1]
    C.zz() = Jm23*(Fe[0][2]*Fe[0][2]+Fe[1][2]*Fe[1][2]+Fe[2][2]*Fe[2][2]); // = C[2][2]
    C.xy() = Jm23*(Fe[0][0]*Fe[0][1]+Fe[1][0]*Fe[1][1]+Fe[2][0]*Fe[2][1]); // = C[0][1]
    C.yz() = Jm23*(Fe[0][1]*Fe[0][2]+Fe[1][1]*Fe[1][2]+Fe[2][1]*Fe[2][2]); // = C[1][2]
    C.xz() = Jm23*(Fe[0][0]*Fe[0][2]+Fe[1][0]*Fe[1][2]+Fe[2][0]*Fe[2][2]); // = C[0][2]


    const double eps = m_epsf*std::numeric_limits<double>::epsilon();

    mat3ds s;
    tens4ds c;

    // Calculate In = n0*C*n0
    double In_1 = n0*(C*n0) - 1.0;

    // only take fibers in tension into consideration
    if (In_1 >= eps)
    {
        // get the global spatial fiber direction in current configuration
        vec3d nt = Fe*n0;

        // calculate the outer product of nt
        mat3ds N = dyad(nt);
        tens4ds NxN = dyad1s(N);

        // calculate strain energy derivative
        double Wl = m_k1*pow(In_1, m_beta - 1.0)*exp(m_k2*pow(In_1, m_beta));

        // calculate the fiber stress
        s = N*(2.0*Wl / J);

        // calculate strain energy 2nd derivative
        double tmp = m_k2*pow(In_1, m_beta);
        double Wll = m_k1*pow(In_1, m_beta - 2.0)*((tmp + 1)*m_beta - 1.0)*exp(tmp);

        // calculate the fiber tangent
        c = NxN*(4.0*Wll / J);

    }
    else
    {
        c.zero();
    }

    // This is the final value of the elasticity tensor
    mat3dd I(1);
    tens4ds IxI = dyad1s(I);
    tens4ds I4 = dyad4s(I);
    c += ((I4 + IxI / 3.0)*s.tr() - dyad1s(I, s))*(2. / 3.) - (ddots(IxI, c) - IxI*(c.tr() / 3.)) / 3.;

    return c;
}

//-----------------------------------------------------------------------------
double CartilageGnR::IntegratedFiberDensity(FEMaterialPoint &mp)
{
//    m_IFD = 0;
    // get the local coordinate systems
    mat3d QT = GetLocalCS(mp).transpose();
    double IFD = 0;

    FEFiberIntegrationSchemeIterator* it = m_pFint->GetIterator(0);
    if (it->IsValid())
    {
        do
        {
            // get fiber direction
            // BUG: why is this fiber not rotated to the local coordinates like in the other functions
            vec3d& n0e = it->m_fiber;

            // rotate to local configuration to evaluate ellipsoidally distributed material coefficients
            vec3d n0a = QT * n0e;
            // evaluate local fiber distribution
            double R = m_pFDD->FiberDensity(mp, n0a);

            // integrate the fiber distribution
            IFD += R*it->m_weight;
        }
        while (it->Next());
    }

    // don't forget to delete the iterator
    delete it;

    // just in case
    if (IFD == 0.0) IFD = 1.0;

    return IFD;
}

