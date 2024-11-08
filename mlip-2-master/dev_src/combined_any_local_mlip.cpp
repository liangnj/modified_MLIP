#include "combined_any_local_mlip.h"
#ifdef MLIP_MPI
#	include <mpi.h>
#endif

using namespace std;

void CombinedAnyLocalMLIP::Load(const std::string& filename)
{

}

void CombinedAnyLocalMLIP::Save(const std::string& filename)
{

}


void CombinedAnyLocalMLIP::CalcSiteEnergyDers(const Neighborhood& big_nbh)
{
    Neighborhood small_nbh = Neighborhood();
    small_nbh.my_type = big_nbh.my_type;

    for (int j = 0; j < big_nbh.count; j++)
    {
        if (big_nbh.dists[j] < p_pot2->CutOff())
        {
            small_nbh.count += 1;
            small_nbh.inds.push_back(big_nbh.inds[j]);
            small_nbh.vecs.push_back(big_nbh.vecs[j]);
            small_nbh.dists.push_back(big_nbh.dists[j]);
            small_nbh.types.push_back(big_nbh.types[j]);
        }
    }

    p_pot1->CalcSiteEnergyDers(big_nbh);
    p_pot2->CalcSiteEnergyDers(small_nbh);

    buff_site_energy_ = p_pot1->buff_site_energy_ + p_pot2->buff_site_energy_;
    buff_site_energy_ders_.resize(big_nbh.count);
    FillWithZero(buff_site_energy_ders_); 

    int k = 0;
    for (int j = 0; j < big_nbh.count; j++) 
    {
        buff_site_energy_ders_[j] = p_pot1->buff_site_energy_ders_[j];

        if (k < small_nbh.count && big_nbh.vecs[j] == small_nbh.vecs[k])
            {
                buff_site_energy_ders_[j] += p_pot2->buff_site_energy_ders_[k];
                k++;
            }
    }
    // if (k != small_nbh.count - 1)
    //     cout << "Error occured during force summation     " << ind << endl; 
}


void CombinedAnyLocalMLIP::CalcE(Configuration& cfg)
{
    cfg.energy = 0.0;

    Neighborhoods big_neighborhoods(cfg, CutOff());

    for (int ind = 0; ind < cfg.size(); ind++)
    {
        const Neighborhood& big_nbh = big_neighborhoods[ind];
        Neighborhood small_nbh = Neighborhood();
        small_nbh.my_type = big_nbh.my_type;

        for (int j = 0; j < big_nbh.count; j++)
        {
            if (big_nbh.dists[j] < p_pot2->CutOff())
            {
                small_nbh.count += 1;
                small_nbh.inds.push_back(big_nbh.inds[j]);
                small_nbh.vecs.push_back(big_nbh.vecs[j]);
                small_nbh.dists.push_back(big_nbh.dists[j]);
                small_nbh.types.push_back(big_nbh.types[j]);
            }
        }

        // Neighborhoods small_neighborhoods(cfg, p_pot2->CutOff());
        // if (small_nbh.inds == small_neighborhoods[ind].inds && small_nbh.vecs == small_neighborhoods[ind].vecs && small_nbh.dists == small_neighborhoods[ind].dists && small_nbh.types == small_neighborhoods[ind].types)
        // {
        //     // cout << ind << "Neighborhoods are similar";
        // }
        // else
        //     cout << "ERROR in neighborhood constructing";

        p_pot1->CalcSiteEnergyDers(big_nbh);
        p_pot2->CalcSiteEnergyDers(small_nbh);


        cfg.energy += p_pot1->buff_site_energy_ + p_pot2->buff_site_energy_;
    }
}

void CombinedAnyLocalMLIP::CalcEFS(Configuration& cfg)
{
    ResetEFS(cfg);

    cfg.has_energy(true);
	cfg.has_forces(true);
	cfg.has_stresses(true);

    Neighborhoods big_neighborhoods(cfg, CutOff());


    for (int ind = 0; ind < cfg.size(); ind++)
    {
        const Neighborhood& big_nbh = big_neighborhoods[ind];
        Neighborhood small_nbh = Neighborhood();
        small_nbh.my_type = big_nbh.my_type;

        for (int j = 0; j < big_nbh.count; j++)
        {
            if (big_nbh.dists[j] < p_pot2->CutOff())
            {
                small_nbh.count += 1;
                small_nbh.inds.push_back(big_nbh.inds[j]);
                small_nbh.vecs.push_back(big_nbh.vecs[j]);
                small_nbh.dists.push_back(big_nbh.dists[j]);
                small_nbh.types.push_back(big_nbh.types[j]);
            }
        }

        
        // Neighborhoods small_neighborhoods(cfg, p_pot2->CutOff());
        // if (small_nbh.inds == small_neighborhoods[ind].inds && small_nbh.vecs == small_neighborhoods[ind].vecs && small_nbh.dists == small_neighborhoods[ind].dists && small_nbh.types == small_neighborhoods[ind].types)
        // {
        //     cout << ind << "Neighborhoods are similar";
        // }
        // else
        //     cout << "ERROR in neighborhood constructing";

        p_pot1->CalcSiteEnergyDers(big_nbh);
        p_pot2->CalcSiteEnergyDers(small_nbh);

        cfg.energy += p_pot1->buff_site_energy_ + p_pot2->buff_site_energy_;

        int k = 0;

        for (int j = 0; j < big_nbh.count; j++) 
        {

            Vector3 curr_force = p_pot1->buff_site_energy_ders_[j];

            if (k < small_nbh.count && big_nbh.vecs[j] == small_nbh.vecs[k])
            {
                curr_force += p_pot2->buff_site_energy_ders_[k];
                k++;
            }

			cfg.force(ind) += curr_force;
			cfg.force(big_nbh.inds[j]) -= curr_force;

            for (int a = 0; a < 3; a++)
				for (int b = 0; b < 3; b++)
					cfg.stresses[a][b] -= curr_force[a] * big_nbh.vecs[j][b];
		}
        // if (k != small_nbh.count)
        //     cout << "Error occurred during force summation     " << ind << endl; 
    }
}



void CombinedAnyLocalMLIP::AccumulateCombinationGrad(const Neighborhood& nbh,
				            std::vector<double>& out_grad_accumulator,
				            const double se_weight,
				            const Vector3* se_ders_weights)
{
} 
