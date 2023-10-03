#pragma once
//#include "FEUncoupledMaterial.h"
#include <FECore/FEModelParam.h>
#include <FECore/FEMaterialPoint.h>
#include "FECore/FEMaterial.h"
#include "FEElasticMaterial.h"
#include "FEElasticFiberMaterial.h"
#include "FEFiberDensityDistribution.h"
#include "FEFiberIntegrationScheme.h"
#include "FEFiberMaterialPoint.h"
#include "CartilageGnRMaterialPoint.h"

//-----------------------------------------------------------------------------
//! Cartilage Growth and Remodeling material

class CartilageGnR : public FEElasticMaterial
{
public:
	CartilageGnR(FEModel* pfem);
    ~CartilageGnR();

    // Initialization
    bool Init() override;

    //! get the elastic base material
//    FEElasticMaterial* GetBaseMaterial() { return m_pBase; }
    FESolidMaterial* GetBaseMaterial() { return m_pBase; }
public:
    double	m_mu;	//!< Neo-Hookean shear modulus
    double	m_K;	//!< Bulk modulus

    double m_k1;    //!< Fiber parameter k_1
    double m_k2;    //!< Fiber parameter k_2

    double m_epsf;
    double m_beta;


    double UJ;
    double UJJ;
    double nvolume;
    double n_FIB;
    double norm_col_fun;
    double norm_pg;

    double nuhat_col;
    double nuhat_pg;

    double m_v;
    double J_cp;

    int m_p;
    int vgnr;

public:
    // returns a pointer to a new material point object
    FEMaterialPointData* CreateMaterialPointData() override;

	//! calculate stress at material point
	mat3ds Stress(FEMaterialPoint& mp) override;

	//! calculate deviatoric tangent stiffness at material point
	tens4ds Tangent(FEMaterialPoint& mp) override;

    //! calculate growth deviatoric stress for a single fiber
    mat3ds GnRDevFiberStress(const mat3d& Fe, const vec3d& n0);

    //! calculate deviatoric tangent stiffness for a single fiber
    tens4ds GnRDevFiberTangent(const mat3d& Fe, const vec3d& n0);



private:
    double IntegratedFiberDensity(FEMaterialPoint& mp);

//protected:
//    FEElasticMaterial*          m_pMat;     // pointer to elastic material
//    FEElasticMaterial*  m_pBase;        //!< pointer to elastic solid material
    FESolidMaterial*  m_pBase;        //!< pointer to elastic solid material
    FEElasticFiberMaterial*     m_pFmat;    // pointer to fiber material
    FEFiberDensityDistribution* m_pFDD;     // pointer to fiber density distribution
    FEFiberIntegrationScheme*   m_pFint;    // pointer to fiber integration scheme


	// declare the parameter list
	DECLARE_FECORE_CLASS();
};
