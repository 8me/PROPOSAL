#include <functional>
#include <cmath>
#include <cassert>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/decay/LeptonicDecayChannel.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include "PROPOSAL/math/MathMethods.h"

template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args)
{
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

using namespace PROPOSAL;

/******************************************************************************
*                         LeptonicDecayChannelApprox                         *
******************************************************************************/


const std::string LeptonicDecayChannelApprox::name_ = "LeptonicDecayChannelApprox";

// ------------------------------------------------------------------------- //
LeptonicDecayChannelApprox::LeptonicDecayChannelApprox(const ParticleDef& lepton,
                                           const ParticleDef& neutrino,
                                           const ParticleDef& anti_neutrino)
    : DecayChannel()
    , massive_lepton_(lepton)
    , neutrino_(neutrino)
    , anti_neutrino_(anti_neutrino)
{
}

// ------------------------------------------------------------------------- //
LeptonicDecayChannelApprox::~LeptonicDecayChannelApprox() {}

// ------------------------------------------------------------------------- //
LeptonicDecayChannelApprox::LeptonicDecayChannelApprox(const LeptonicDecayChannelApprox& mode)
    : DecayChannel(mode)
    , massive_lepton_(mode.massive_lepton_)
    , neutrino_(mode.neutrino_)
    , anti_neutrino_(mode.anti_neutrino_)
{
}

// ------------------------------------------------------------------------- //
bool LeptonicDecayChannelApprox::compare(const DecayChannel& channel) const
{
    const LeptonicDecayChannelApprox* leptonic = dynamic_cast<const LeptonicDecayChannelApprox*>(&channel);

    if (!leptonic)
        return false;
    else if (massive_lepton_ != leptonic->massive_lepton_)
        return false;
    else if (neutrino_ != leptonic->neutrino_)
        return false;
    else if (anti_neutrino_ != leptonic->anti_neutrino_)
        return false;
    else
        return true;
}

// ------------------------------------------------------------------------- //
double LeptonicDecayChannelApprox::DecayRate(double x, double parent_mass, double E_max, double right_side)
{
    (void)parent_mass;
    (void)E_max;

    return x * x * x * (1. - 0.5 * x) - right_side;
}

// ------------------------------------------------------------------------- //
double LeptonicDecayChannelApprox::DifferentialDecayRate(double x, double parent_mass, double E_max)
{
    (void)parent_mass;
    (void)E_max;

    return (3 - 2 * x) * x * x;
}

// ------------------------------------------------------------------------- //
double LeptonicDecayChannelApprox::FindRoot(double min, double parent_mass, double E_max, double right_side)
{
    double max        = 1;
    double x_start    = 0.5;

    return NewtonRaphson(std::bind(&LeptonicDecayChannelApprox::DecayRate, this, std::placeholders::_1, parent_mass, E_max, right_side),
                         std::bind(&LeptonicDecayChannelApprox::DifferentialDecayRate, this, std::placeholders::_1, parent_mass, E_max),
                         min, max, x_start, 40, 1e-3);


}

// ------------------------------------------------------------------------- //
double LeptonicDecayChannelApprox::CZSampler(double x, double P, std::string particle_name, double uy)
{
    //P is the magnitude of the polarization vector.
    double a = 0.5 * P * (1. - 2*x);
    double b = (3. - 2.*x);

    if (particle_name == "TauPlus")
        a = -a;

    double c = b - a - 2*(3.-2*x) * uy;

    double result = (-b + sqrt(b*b - 4 * a * c))/( 2* a);

    return result;
}
Cartesian3D LeptonicDecayChannelApprox::RotateX(Cartesian3D vec, double theta){

    double x = vec.GetX();
    double y = vec.GetY();
    double z = vec.GetZ();

    double cth = std::cos(theta);
    double sth = std::sin(theta);
    Cartesian3D rotated_vec = Cartesian3D(x, y*cth - z*sth, y*sth + z*cth);

    return rotated_vec;
}

Cartesian3D LeptonicDecayChannelApprox::RotateY(Cartesian3D vec, double theta){

    double x = vec.GetX();
    double y = vec.GetY();
    double z = vec.GetZ();

    double cth = std::cos(theta);
    double sth = std::sin(theta);
    Cartesian3D rotated_vec = Cartesian3D(x*cth + z*sth, y, -x*sth + z*cth);

    return rotated_vec;
}

Cartesian3D LeptonicDecayChannelApprox::RotateZ(Cartesian3D vec, double theta){

    double x = vec.GetX();
    double y = vec.GetY();
    double z = vec.GetZ();

    double cth = std::cos(theta);
    double sth = std::sin(theta);
    Cartesian3D rotated_vec = Cartesian3D(x*cth - y*sth, x*sth + y*cth, z);

    return rotated_vec;
}

std::vector<ParticleState> LeptonicDecayChannelApprox::Decay(const ParticleDef& p_def, const ParticleState& p_condition)
{
    assert (p_condition.direction.magnitude() > 0);
    // Sample energy from decay rate
    double emax       = (p_def.mass * p_def.mass + massive_lepton_.mass * massive_lepton_.mass) / (2 * p_def.mass);
    double x_min      = massive_lepton_.mass / emax;

    double f_min      = DecayRate(x_min, p_def.mass, emax, 0.0);
    double f_max      = DecayRate(1.0, p_def.mass, emax, 0.0);
    double right_side = f_min + (f_max - f_min) * RandomGenerator::Get().RandomDouble();

    double find_root = FindRoot(x_min, p_def.mass, emax, right_side);

    double lepton_energy   = std::max(find_root * emax, massive_lepton_.mass);
    double lepton_momentum = std::sqrt((lepton_energy - massive_lepton_.mass) * (lepton_energy + massive_lepton_.mass));

    Cartesian3D lepton_direction;
    if (massive_lepton_.name == "TauMinus" or
        massive_lepton_.name == "TauPlus") {

        double lh_probability = lepton_momentum / (lepton_energy + 1776.0);

        Cartesian3D spin;
        if ( (massive_lepton_.name == "TauMinus") &&
             (lh_probability < RandomGenerator::Get().RandomDouble()) ) {
            spin = -p_condition.direction;
        } else {
            spin = p_condition.direction;
        }

        double sampled_cth = CZSampler(find_root, 1, massive_lepton_.name, RandomGenerator::Get().RandomDouble());
        double sampled_phi= 2.0 * PI * RandomGenerator::Get().RandomDouble();
        double cphi = std::cos(sampled_phi);
        double sphi = std::sin(sampled_phi);
        double sampled_sth = sqrt( 1. - sampled_cth*sampled_cth);

        //makes the direction vector in the frame where the spin is aligned with the z-axis.
        Cartesian3D vec_off_spin = Cartesian3D(cphi * sampled_sth, sphi * sampled_sth, sampled_cth);
        // vec_off_spin.CalculateSphericalCoordinates();
        // spin.CalculateSphericalCoordinates();
        std::array<double, 3> spin_spherical = spin.GetSphericalCoordinates();
        lepton_direction = RotateZ(RotateY(vec_off_spin, spin_spherical[2]), spin_spherical[1]);
    } else {
        lepton_direction = GenerateRandomDirection();

    }
    // Sample directions For the massive letpon
    ParticleState massive_lepton((ParticleType)massive_lepton_.particle_type,
                                 p_condition.position,
                                 lepton_direction,
                                 lepton_energy,
                                 p_condition.time,
                                 0.);

    // Sample directions For the massless letpon
    double energy_neutrinos   = p_def.mass - lepton_energy;
    double virtual_mass       = std::sqrt((energy_neutrinos - lepton_momentum) * (energy_neutrinos + lepton_momentum));
    double momentum_neutrinos = 0.5 * virtual_mass;


    auto direction = GenerateRandomDirection();

    ParticleState neutrino((ParticleType)neutrino_.particle_type,
                           p_condition.position,
                           direction,
                           momentum_neutrinos,
                           p_condition.time,
                           0.);

    Cartesian3D opposite_direction = -direction;

    ParticleState anti_neutrino((ParticleType)anti_neutrino_.particle_type,
                                p_condition.position,
                                opposite_direction,
                                momentum_neutrinos,
                                p_condition.time,
                                0.);

    // Boost neutrinos to lepton frame
    // double beta = lepton_momentum / energy_neutrinos;
    double gamma = energy_neutrinos / virtual_mass;
    double betagamma = lepton_momentum / virtual_mass;


    Boost(neutrino, massive_lepton.direction, gamma, betagamma);
    Boost(anti_neutrino, massive_lepton.direction, gamma, betagamma);


    std::vector<ParticleState> secondaries;
    secondaries.push_back(massive_lepton);
    secondaries.push_back(neutrino);
    secondaries.push_back(anti_neutrino);

    // Get Momentum is not defined for pseudo particle decay, so it must be
    // calculated manually
    double primary_momentum = std::sqrt(std::max((p_condition.energy + p_def.mass) * (p_condition.energy - p_def.mass), 0.0));
    // Boost all products in Lab frame (the reason, why the boosting goes in the negative direction of the particle)
    Boost(secondaries, -p_condition.direction, p_condition.energy/p_def.mass, primary_momentum/p_def.mass);

    return secondaries;
}


// ------------------------------------------------------------------------- //
// Print
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
void LeptonicDecayChannelApprox::print(std::ostream& os) const
{
    os << "Massive lepton:\n" << massive_lepton_ << '\n';
    os << "Neutrino:\n" << neutrino_ << '\n';
    os << "Anti neutrino:\n" << anti_neutrino_ << '\n';
}

/******************************************************************************
 *                          LeptonicDecayChannel                              *
 ******************************************************************************/

const std::string LeptonicDecayChannel::name_ = "LeptonicDecayChannel";

// ------------------------------------------------------------------------- //
LeptonicDecayChannel::LeptonicDecayChannel(const ParticleDef& lepton,
                                           const ParticleDef& neutrino,
                                           const ParticleDef& anti_neutrino)
    : LeptonicDecayChannelApprox(lepton, neutrino, anti_neutrino)
{
}

// ------------------------------------------------------------------------- //
LeptonicDecayChannel::~LeptonicDecayChannel() {}

// ------------------------------------------------------------------------- //
LeptonicDecayChannel::LeptonicDecayChannel(const LeptonicDecayChannel& mode)
    : LeptonicDecayChannelApprox(mode)
{
}

// ------------------------------------------------------------------------- //
double LeptonicDecayChannel::DecayRate(double x, double M, double E_max, double right_side)
{
    double M2 = M * M;
    double m  = massive_lepton_.mass;
    double m2 = m * m;

    double E_l     = E_max * x;
    double sqrt_EM = std::sqrt(E_l * E_l - m2);

    return 1.5 * m2 * m2 * M * std::log(sqrt_EM + E_l) +
           sqrt_EM * ((M2 + m2 - M * E_l) * (E_l * E_l - m2) - 1.5 * M * E_l * m2) - right_side;
}

// ------------------------------------------------------------------------- //
double LeptonicDecayChannel::DifferentialDecayRate(double x, double M, double E_max)
{
    double m   = massive_lepton_.mass;
    double E_l = E_max * x;

    return E_max * std::sqrt(E_l * E_l - m * m) * (M * E_l * (3.0 * M - 4.0 * E_l) + m * m * (3.0 * E_l - 2 * M));
}
