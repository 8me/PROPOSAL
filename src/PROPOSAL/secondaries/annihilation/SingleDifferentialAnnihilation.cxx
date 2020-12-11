#include "PROPOSAL/secondaries/annihilation/SingleDifferentialAnnihilation.h"
#include "PROPOSAL/Constants.h"

#include <cmath>
#include <stdexcept>
#include <tuple>

using std::fmod;
using std::make_tuple;
using std::sqrt;

using namespace PROPOSAL;

double secondaries::SingleDifferentialAnnihilation::CalculateRho(
    double energy, double rnd, const Component& comp)
{
    for (auto& it : dndx) {
        if (comp.GetName() == it.first->GetName()) {
            auto rate = it.second->Calculate(energy);
            return it.second->GetUpperLimit(energy, rnd * rate);
        }
    }
    std::ostringstream s;
    s << "Component (" << comp.GetName()
      << ") can not be found in the precalculated annihilation tables.";
    throw std::out_of_range(s.str());
}

tuple<Vector3D, Vector3D>
secondaries::SingleDifferentialAnnihilation::CalculateDirections(
    Vector3D primary_dir, double energy, double rho, double rnd)
{
    auto com_energy = energy + ME; // center of mass energy
    auto kin_energy = energy - ME; // kinetic energy
    auto cosphi0 = (com_energy * (1. - rho) - ME)
        / ((1. - rho) * sqrt(com_energy * kin_energy));
    auto cosphi1
        = (com_energy * rho - ME) / (rho * sqrt(com_energy * kin_energy));
    auto rnd_theta = rnd * 2. * PI;
    auto dir_1 = primary_dir;
    auto dir_2 = primary_dir;
    dir_1.deflect(cosphi0, rnd_theta);
    dir_2.deflect(cosphi1, fmod(rnd_theta + PI, 2. * PI));
    return make_tuple(dir_1, dir_2);
}

tuple<double, double>
secondaries::SingleDifferentialAnnihilation::CalculateEnergy(
    double energy, double rho)
{
    assert(rho >= 0);
    assert(rho <= 1);
    auto energy_1 = (energy + ME) * (1 - rho);
    auto energy_2 = (energy + ME) * rho;
    return make_tuple(energy_1, energy_2);
}

vector<ParticleState>
secondaries::SingleDifferentialAnnihilation::CalculateSecondaries(
        StochasticLoss loss, const Component& comp, vector<double> &rnd)
{
    auto rho = CalculateRho(loss.parent_particle_energy, rnd[0], comp);
    auto secondary_energies = CalculateEnergy(loss.parent_particle_energy, rho);
    auto secondary_dir = CalculateDirections(
            loss.direction, loss.parent_particle_energy, rho, rnd[1]);
    auto sec = std::vector<ParticleState>();
    sec.emplace_back(ParticleType::Gamma, loss.position, get<0>(secondary_dir),
            get<0>(secondary_energies), loss.time, 0.);
    sec.emplace_back(ParticleType::Gamma, loss.position, get<1>(secondary_dir),
            get<1>(secondary_energies), loss.time, 0.);
    return sec;
}
