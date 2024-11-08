#ifndef COMBINED_ANY_LOCAL_MLIP
#define COMBINED_ANY_LOCAL_MLIP

#include "../src/pair_potentials.h"
#include "../src/basic_mlip.h"
#include "mtpr.h"
#include "eam.h"

class CombinedAnyLocalMLIP : public AnyLocalMLIP
{

protected:
    std::ostream* p_logstrm;

public:
    void CalcSiteEnergyDers(const Neighborhood&) override;
    void CalcE(Configuration&) override; 
    void CalcEFS(Configuration&) override;	
    
    void Load(const std::string& filename) override;
    void Save(const std::string& filename) override;										

    AnyLocalMLIP* p_pot1 = nullptr;
    AnyLocalMLIP* p_pot2 = nullptr;

    CombinedAnyLocalMLIP(AnyLocalMLIP* _p_pot1, AnyLocalMLIP* _p_pot2, std::ostream* _p_logstrm = nullptr)
    { 
        if (_p_pot1 == nullptr) 
            ERROR("CombinedAnyLocalMLIP: attempting to initialize with p_pot1 = nullptr!");
        if (_p_pot2 == nullptr) 
            ERROR("CombinedAnyLocalMLIP: attempting to initialize with p_pot2 = nullptr!");

        // Potential with higher cutoff is p_pot1, another one is p_pot2
        if (_p_pot2->CutOff() > _p_pot1->CutOff())
        {
            p_pot1 = _p_pot2;
            p_pot2 = _p_pot1;
        }
        else
        {    
            p_pot1 = _p_pot1;
            p_pot2 = _p_pot2;
        }
        p_RadialBasis = p_pot1->p_RadialBasis;
        p_logstrm = _p_logstrm;
    };
    

    int CoeffCount() //!< number of coefficients
    {
        return p_pot1->CoeffCount();
    }
    
    double* Coeff() //!< coefficients themselves
    {
        return &p_pot1->Coeff()[0];
    }
    
    void AccumulateCombinationGrad(const Neighborhood&, std::vector<double>&,
				    const double se_weight = 0.0,
				    const Vector3* se_ders_weights = nullptr) override;

    ~CombinedAnyLocalMLIP() {};
										
};

#endif // COMBINED_ANY_LOCAL_MLIP
