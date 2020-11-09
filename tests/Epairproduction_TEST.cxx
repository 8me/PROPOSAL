
#include "gtest/gtest.h"

#include <fstream>
#include "PROPOSAL/crosssection/parametrization/EpairProduction.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/medium/MediumFactory.h"

using namespace PROPOSAL;

ParticleDef getParticleDef(const std::string& name)
{
    if (name == "MuMinus")
    {
        return MuMinusDef();
    } else if (name == "TauMinus")
    {
        return TauMinusDef();
    } else
    {
        return EMinusDef();
    }
}

const std::string testfile_dir = "bin/TestFiles/";

// TEST(Comparison, Comparison_equal)
// {
//     ParticleDef particle_def = MuMinusDef::Get();
//     auto medium = std::make_shared<const Water>();
//     EnergyCutSettings ecuts;
//     double multiplier = 1.;
//     bool lpm          = true;

//     EpairProductionRhoIntegral* EpairInt_A =
//         new EpairKelnerKokoulinPetrukhin(particle_def, medium, ecuts, multiplier, lpm);
//     Parametrization* EpairInt_B = new EpairKelnerKokoulinPetrukhin(particle_def, medium, ecuts, multiplier, lpm);
//     EXPECT_TRUE(*EpairInt_A == *EpairInt_B);

//     EpairKelnerKokoulinPetrukhin param_int(particle_def, medium, ecuts, multiplier, lpm);
//     EXPECT_TRUE(param_int == *EpairInt_A);

//     EpairIntegral* Int_A        = new EpairIntegral(param_int);
//     CrossSectionIntegral* Int_B = new EpairIntegral(param_int);
//     EXPECT_TRUE(*Int_A == *Int_B);

//     InterpolationDef InterpolDef;
//     EpairProductionRhoInterpolant<EpairKelnerKokoulinPetrukhin>* EpairInterpol_A =
//         new EpairProductionRhoInterpolant<EpairKelnerKokoulinPetrukhin>(particle_def, medium, ecuts, multiplier, lpm, InterpolDef);
//     Parametrization* EpairInterpol_B =
//         new EpairProductionRhoInterpolant<EpairKelnerKokoulinPetrukhin>(particle_def, medium, ecuts, multiplier, lpm, InterpolDef);
//     EXPECT_TRUE(*EpairInterpol_A == *EpairInterpol_B);

//     EpairProductionRhoInterpolant<EpairKelnerKokoulinPetrukhin> param_interpol(particle_def, medium, ecuts, multiplier, lpm, InterpolDef);
//     EXPECT_TRUE(param_interpol == *EpairInterpol_A);

//     EpairInterpolant* Interpol_A        = new EpairInterpolant(param_interpol, InterpolDef);
//     CrossSectionInterpolant* Interpol_B = new EpairInterpolant(param_interpol, InterpolDef);
//     EXPECT_TRUE(*Interpol_A == *Interpol_B);

//     delete EpairInt_A;
//     delete EpairInt_B;
//     delete Int_A;
//     delete Int_B;
//     delete EpairInterpol_A;
//     delete EpairInterpol_B;
//     delete Interpol_A;
//     delete Interpol_B;
// }

// TEST(Comparison, Comparison_not_equal)
// {
//     ParticleDef mu_def  = MuMinusDef::Get();
//     ParticleDef tau_def = TauMinusDef::Get();
//     auto medium_1 = std::make_shared<const Water>();
//     auto medium_2 = std::make_shared<const Ice>();
//     EnergyCutSettings ecuts_1(500, -1);
//     EnergyCutSettings ecuts_2(-1, 0.05);
//     double multiplier_1 = 1.;
//     double multiplier_2 = 2.;
//     bool lpm_1          = true;
//     bool lpm_2          = false;

//     EpairKelnerKokoulinPetrukhin EpairInt_A(mu_def, medium_1, ecuts_1, multiplier_1, lpm_1);
//     EpairKelnerKokoulinPetrukhin EpairInt_B(tau_def, medium_1, ecuts_1, multiplier_1, lpm_1);
//     EpairKelnerKokoulinPetrukhin EpairInt_C(mu_def, medium_2, ecuts_1, multiplier_1, lpm_1);
//     EpairKelnerKokoulinPetrukhin EpairInt_D(mu_def, medium_1, ecuts_2, multiplier_1, lpm_1);
//     EpairKelnerKokoulinPetrukhin EpairInt_E(mu_def, medium_1, ecuts_1, multiplier_2, lpm_1);
//     EpairKelnerKokoulinPetrukhin EpairInt_F(mu_def, medium_1, ecuts_1, multiplier_1, lpm_2);
//     EXPECT_TRUE(EpairInt_A != EpairInt_B);
//     EXPECT_TRUE(EpairInt_A != EpairInt_C);
//     EXPECT_TRUE(EpairInt_A != EpairInt_D);
//     EXPECT_TRUE(EpairInt_A != EpairInt_E);
//     EXPECT_TRUE(EpairInt_A != EpairInt_F);

//     EpairIntegral Int_A(EpairInt_A);
//     EpairIntegral Int_B(EpairInt_B);
//     EXPECT_TRUE(Int_A != Int_B);

//     InterpolationDef InterpolDef;
//     EpairProductionRhoInterpolant<EpairKelnerKokoulinPetrukhin> EpairInterpol_A(mu_def, medium_1, ecuts_1, multiplier_1, lpm_1, InterpolDef);
//     EpairProductionRhoInterpolant<EpairKelnerKokoulinPetrukhin> EpairInterpol_B(tau_def, medium_1, ecuts_1, multiplier_1, lpm_1, InterpolDef);
//     EpairProductionRhoInterpolant<EpairKelnerKokoulinPetrukhin> EpairInterpol_C(mu_def, medium_2, ecuts_1, multiplier_1, lpm_1, InterpolDef);
//     EpairProductionRhoInterpolant<EpairKelnerKokoulinPetrukhin> EpairInterpol_D(mu_def, medium_1, ecuts_2, multiplier_1, lpm_1, InterpolDef);
//     EpairProductionRhoInterpolant<EpairKelnerKokoulinPetrukhin> EpairInterpol_E(mu_def, medium_1, ecuts_1, multiplier_2, lpm_1, InterpolDef);
//     EpairProductionRhoInterpolant<EpairKelnerKokoulinPetrukhin> EpairInterpol_F(mu_def, medium_1, ecuts_1, multiplier_1, lpm_2, InterpolDef);
//     EXPECT_TRUE(EpairInterpol_A != EpairInterpol_B);
//     EXPECT_TRUE(EpairInterpol_A != EpairInterpol_C);
//     EXPECT_TRUE(EpairInterpol_A != EpairInterpol_D);
//     EXPECT_TRUE(EpairInterpol_A != EpairInterpol_E);
//     EXPECT_TRUE(EpairInterpol_A != EpairInterpol_F);

//     EpairInterpolant Interpol_A(EpairInterpol_A, InterpolDef);
//     EpairInterpolant Interpol_B(EpairInterpol_B, InterpolDef);
//     EXPECT_TRUE(Interpol_A != Interpol_B);
// }

// TEST(Assignment, Copyconstructor)
// {
//     ParticleDef particle_def = MuMinusDef::Get();
//     auto medium = std::make_shared<const Water>();
//     EnergyCutSettings ecuts;
//     double multiplier = 1.;
//     bool lpm          = true;

//     EpairKelnerKokoulinPetrukhin EpairInt_A(particle_def, medium, ecuts, multiplier, lpm);
//     EpairKelnerKokoulinPetrukhin EpairInt_B = EpairInt_A;
//     EXPECT_TRUE(EpairInt_A == EpairInt_B);

//     EpairIntegral Int_A(EpairInt_A);
//     EpairIntegral Int_B = Int_A;
//     EXPECT_TRUE(Int_A == Int_B);

//     InterpolationDef InterpolDef;
//     EpairProductionRhoInterpolant<EpairKelnerKokoulinPetrukhin> EpairInterpol_A(particle_def, medium, ecuts, multiplier, lpm, InterpolDef);
//     EpairProductionRhoInterpolant<EpairKelnerKokoulinPetrukhin> EpairInterpol_B = EpairInterpol_A;
//     EXPECT_TRUE(EpairInterpol_A == EpairInterpol_B);

//     EpairInterpolant Interpol_A(EpairInterpol_A, InterpolDef);
//     EpairInterpolant Interpol_B = Interpol_A;
//     EXPECT_TRUE(Interpol_A == Interpol_B);
// }

// TEST(Assignment, Copyconstructor2)
// {
//     ParticleDef particle_def = MuMinusDef::Get();
//     auto medium = std::make_shared<const Water>();
//     EnergyCutSettings ecuts;
//     double multiplier = 1.;
//     bool lpm          = true;

//     EpairKelnerKokoulinPetrukhin EpairInt_A(particle_def, medium, ecuts, multiplier, lpm);
//     EpairKelnerKokoulinPetrukhin EpairInt_B(EpairInt_A);
//     EXPECT_TRUE(EpairInt_A == EpairInt_B);

//     EpairIntegral Int_A(EpairInt_A);
//     EpairIntegral Int_B(Int_A);
//     EXPECT_TRUE(Int_A == Int_B);

//     InterpolationDef InterpolDef;
//     EpairProductionRhoInterpolant<EpairKelnerKokoulinPetrukhin> EpairInterpol_A(particle_def, medium, ecuts, multiplier, lpm, InterpolDef);
//     EpairProductionRhoInterpolant<EpairKelnerKokoulinPetrukhin> EpairInterpol_B(EpairInterpol_A);
//     EXPECT_TRUE(EpairInterpol_A == EpairInterpol_B);

//     EpairInterpolant Interpol_A(EpairInterpol_A, InterpolDef);
//     EpairInterpolant Interpol_B(Interpol_A);
//     EXPECT_TRUE(Interpol_A == Interpol_B);
// }

// in polymorphism an assignmant and swap operator doesn't make sense

TEST(Epairproduction, Test_of_dEdx)
{
    std::string filename = testfile_dir + "Epair_dEdx.txt";
    std::ifstream in{filename};
    EXPECT_TRUE(in.good()) << "Test resource file '" << filename << "' could not be opened";

    std::string particleName;
    std::string mediumName;
    double ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    bool lpm;
    double energy;
    std::string parametrization;
    double dEdx_stored;
    double dEdx_new;

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> lpm >> energy >> parametrization >> dEdx_stored)
    {
        parametrization.erase(0,5);
        if (vcut == -1)
            vcut = 1;
        if (ecut == -1)
            ecut = INF;

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        auto cross = crosssection::make_epairproduction(particle_def,
            *medium,
            ecuts,
            false,
            lpm,
            parametrization);

        dEdx_new = cross->CalculatedEdx(energy) * medium->GetMassDensity();

        EXPECT_NEAR(dEdx_new, dEdx_stored, 1e-10 * dEdx_stored);
    }
}

TEST(Epairproduction, Test_of_dNdx)
{
    std::string filename = testfile_dir + "Epair_dNdx.txt";
    std::ifstream in{filename};
    EXPECT_TRUE(in.good()) << "Test resource file '" << filename << "' could not be opened";

    std::string particleName;
    std::string mediumName;
    double ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    bool lpm;
    double energy;
    std::string parametrization;
    double dNdx_stored;
    double dNdx_new;

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> lpm >> energy >> parametrization >> dNdx_stored)
    {
        parametrization.erase(0,5);
        if (vcut == -1)
            vcut = 1;
        if (ecut == -1)
            ecut = INF;

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        auto cross = crosssection::make_epairproduction(particle_def,
            *medium,
            ecuts,
            false,
            lpm,
            parametrization);

        dNdx_new = cross->CalculatedNdx(energy) * medium->GetMassDensity();

        ASSERT_NEAR(dNdx_new, dNdx_stored, 1e-10 * dNdx_stored);
    }
}

TEST(Epairproduction, Test_Stochastic_Loss)
{
    std::string filename = testfile_dir + "Epair_e.txt";
    std::ifstream in{filename};
    EXPECT_TRUE(in.good()) << "Test resource file '" << filename << "' could not be opened";

    std::string particleName;
    std::string mediumName;
    double ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    bool lpm;
    double energy;
    std::string parametrization;
    double rnd1;
    double rnd2;
    double stochastic_loss_stored;
    double stochastic_loss_new;

    std::cout.precision(16);
    RandomGenerator::Get().SetSeed(0);

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> lpm >> energy >> parametrization >> rnd1 >> rnd2 >> stochastic_loss_stored)
    {
        parametrization.erase(0,5);
        if (vcut == -1)
            vcut = 1;
        if (ecut == -1)
            ecut = INF;

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        auto cross = crosssection::make_epairproduction(particle_def,
            *medium,
            ecuts,
            false,
            lpm,
            parametrization);

        auto dNdx_full = cross->CalculatedNdx(energy);
        auto components = cross->GetTargets();
        double sum = 0;

        for (auto comp : components)
        {
            double dNdx_for_comp = cross->CalculatedNdx(energy, comp);
            sum += dNdx_for_comp;
            if (sum > dNdx_full * (1. - rnd2)) {
                double rate_new = dNdx_for_comp * rnd1;
                if (ecut == INF and vcut == 1 ) {
                    EXPECT_DEATH(cross->CalculateStochasticLoss(comp, energy, rate_new), "");
                } else {
                    stochastic_loss_new = energy * cross->CalculateStochasticLoss(comp, energy, rate_new);
                    EXPECT_NEAR(stochastic_loss_new, stochastic_loss_stored, 1E-6 * stochastic_loss_stored);
                    break;
                }
            }
        }
    }
}

TEST(Epairproduction, Test_of_dEdx_Interpolant)
{
    std::string filename = testfile_dir + "Epair_dEdx_interpol.txt";
    std::ifstream in{filename};
    EXPECT_TRUE(in.good()) << "Test resource file '" << filename << "' could not be opened";

    std::string particleName;
    std::string mediumName;
    double ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    bool lpm;
    double energy;
    std::string parametrization;
    double dEdx_stored;
    double dEdx_new;

    InterpolationDef InterpolDef;

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> lpm >> energy >> parametrization >> dEdx_stored)
    {
        parametrization.erase(0,5);
        if (vcut == -1)
            vcut = 1;
        if (ecut == -1)
            ecut = INF;

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        auto cross = crosssection::make_epairproduction(particle_def,
            *medium,
            ecuts,
            true,
            lpm,
            parametrization);

        dEdx_new = cross->CalculatedEdx(energy) * medium->GetMassDensity();

        EXPECT_NEAR(dEdx_new, dEdx_stored, 1e-10 * dEdx_stored);
    }
}

TEST(Epairproduction, Test_of_dNdx_Interpolant)
{
    std::string filename = testfile_dir + "Epair_dNdx_interpol.txt";
    std::ifstream in{filename};
    EXPECT_TRUE(in.good()) << "Test resource file '" << filename << "' could not be opened";

    std::string particleName;
    std::string mediumName;
    double ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    bool lpm;
    double energy;
    std::string parametrization;
    double dNdx_stored;
    double dNdx_new;

    InterpolationDef InterpolDef;

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> lpm >> energy >> parametrization >> dNdx_stored)
    {
        parametrization.erase(0,5);
        if (vcut == -1)
            vcut = 1;
        if (ecut == -1)
            ecut = INF;

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        auto cross = crosssection::make_epairproduction(particle_def,
            *medium,
            ecuts,
            true,
            lpm,
            parametrization);

        dNdx_new = cross->CalculatedNdx(energy) * medium->GetMassDensity();

        EXPECT_NEAR(dNdx_new, dNdx_stored, 1e-10 * dNdx_stored);
    }
}

TEST(Epairproduction, Test_of_e_interpol)
{
    std::string filename = testfile_dir + "Epair_e_interpol.txt";
    std::ifstream in{filename};
    EXPECT_TRUE(in.good()) << "Test resource file '" << filename << "' could not be opened";

    std::string particleName;
    std::string mediumName;
    double ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    bool lpm;
    double energy;
    std::string parametrization;
    double rnd1;
    double rnd2;
    double stochastic_loss_stored;
    double stochastic_loss_new;

    InterpolationDef InterpolDef;

    RandomGenerator::Get().SetSeed(0);

    while (in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> lpm >> energy >> parametrization >> rnd1 >> rnd2 >> stochastic_loss_stored)
    {        parametrization.erase(0,5);
        if (vcut == -1)
            vcut = 1;
        if (ecut == -1)
            ecut = INF;

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        auto cross = crosssection::make_epairproduction(particle_def,
            *medium,
            ecuts,
            true,
            lpm,
            parametrization);

        auto dNdx_full = cross->CalculatedNdx(energy);
        auto components = cross->GetTargets();
        double sum = 0;

        for (auto comp : components)
        {
            double dNdx_for_comp = cross->CalculatedNdx(energy, comp);
            sum += dNdx_for_comp;
            if (sum > dNdx_full * (1. - rnd2)) {
                double rate_new = dNdx_for_comp * rnd1;
                if (ecut == INF and vcut == 1 ) {
                    EXPECT_DEATH(cross->CalculateStochasticLoss(comp, energy, rate_new), "");
                } else {
                    stochastic_loss_new = energy * cross->CalculateStochasticLoss(comp, energy, rate_new);
                    EXPECT_NEAR(stochastic_loss_new, stochastic_loss_stored, 1E-6 * stochastic_loss_stored);
                    break;
                }
            }
        }
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
