
// #include <iostream>
// #include <string>

#include "gtest/gtest.h"

#include <fstream>
#include "PROPOSAL/crosssection/parametrization/Compton.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/medium/MediumFactory.h"

using namespace PROPOSAL;

const std::string testfile_dir = "bin/TestFiles/";

// TEST(Comparison, Comparison_equal)
// {
// ParticleDef particle_def = GammaDef::Get();
// auto medium = std::make_shared<const Water>();
// EnergyCutSettings ecuts;
// double multiplier = 1.;

// ComptonKleinNishina* Compton_A =
//         new ComptonKleinNishina(particle_def, medium, ecuts, multiplier);
// Parametrization* Compton_B = new ComptonKleinNishina(particle_def, medium, ecuts, multiplier);
// EXPECT_TRUE(*Compton_A == *Compton_B);

// ComptonKleinNishina param(particle_def, medium, ecuts, multiplier);
// EXPECT_TRUE(param == *Compton_A);

// ComptonIntegral* Int_A        = new ComptonIntegral(param);
// CrossSectionIntegral* Int_B = new ComptonIntegral(param);
// EXPECT_TRUE(*Int_A == *Int_B);

// InterpolationDef InterpolDef;
// ComptonInterpolant* Interpol_A        = new ComptonInterpolant(param, InterpolDef);
// CrossSectionInterpolant* Interpol_B = new ComptonInterpolant(param, InterpolDef);
// EXPECT_TRUE(*Interpol_A == *Interpol_B);

// delete Compton_A;
// delete Compton_B;
// delete Int_A;
// delete Int_B;
// delete Interpol_A;
// delete Interpol_B;
// }

// TEST(Comparison, Comparison_not_equal)
// {
// ParticleDef gamma_def  = GammaDef::Get();
// auto medium_1 = std::make_shared<const Water>();
// auto medium_2 = std::make_shared<const Ice>();
// EnergyCutSettings ecuts_1(500, -1);
// EnergyCutSettings ecuts_2(-1, 0.05);
// double multiplier_1 = 1.;
// double multiplier_2 = 2.;

// ComptonKleinNishina Compton_A(gamma_def, medium_1, ecuts_1, multiplier_1);
// ComptonKleinNishina Compton_B(gamma_def, medium_2, ecuts_1, multiplier_1);
// ComptonKleinNishina Compton_C(gamma_def, medium_1, ecuts_2, multiplier_1);
// ComptonKleinNishina Compton_D(gamma_def, medium_1, ecuts_1, multiplier_2);
// EXPECT_TRUE(Compton_A != Compton_B);
// EXPECT_TRUE(Compton_A != Compton_C);
// EXPECT_TRUE(Compton_A != Compton_C);
// }

// TEST(Assignment, Copyconstructor)
// {
// ParticleDef particle_def = GammaDef::Get();
// auto medium = std::make_shared<const Water>();
// EnergyCutSettings ecuts;
// double multiplier = 1.;

// ComptonKleinNishina Compton_A(particle_def, medium, ecuts, multiplier);
// ComptonKleinNishina Compton_B = Compton_A;
// EXPECT_TRUE(Compton_A == Compton_B);

// ComptonIntegral Int_A(Compton_A);
// ComptonIntegral Int_B = Int_A;
// EXPECT_TRUE(Int_A == Int_B);

// InterpolationDef InterpolDef;
// ComptonInterpolant Interpol_A(Compton_A, InterpolDef);
// ComptonInterpolant Interpol_B = Interpol_A;
// EXPECT_TRUE(Interpol_A == Interpol_B);
// }

// TEST(Assignment, Copyconstructor2)
// {
// ParticleDef particle_def = GammaDef::Get();
// auto medium = std::make_shared<const Water>();
// EnergyCutSettings ecuts;
// double multiplier = 1.;

// ComptonKleinNishina Compton_A(particle_def, medium, ecuts, multiplier);
// ComptonKleinNishina Compton_B(Compton_A);
// EXPECT_TRUE(Compton_A == Compton_B);

// ComptonIntegral Int_A(Compton_A);
// ComptonIntegral Int_B(Int_A);
// EXPECT_TRUE(Int_A == Int_B);

// InterpolationDef InterpolDef;
// ComptonInterpolant Interpol_A(Compton_A, InterpolDef);
// ComptonInterpolant Interpol_B(Interpol_A);
// EXPECT_TRUE(Interpol_A == Interpol_B);
// }

// in polymorphism an assignment and swap operator doesn't make sense

TEST(Compton, Test_of_dEdx)
{
    std::string filename = testfile_dir + "Compton_dEdx.txt";
    std::ifstream in{filename};
    EXPECT_TRUE(in.good()) << "Test resource file '" << filename << "' could not be opened";

    std::string mediumName;
    double ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    std::string parametrization;
    double energy;
    double dEdx_stored;
    double dEdx_new;

    std::cout.precision(16);

    while (in >> mediumName >> ecut >> vcut >> multiplier >> energy >> dEdx_stored >> parametrization)
    {
        parametrization.erase(0,7);

        if (vcut == -1)
            vcut = 1;
        if (ecut == -1)
            ecut = INF;

        ParticleDef particle_def = GammaDef();
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        auto cross = crosssection::make_compton(particle_def,
            *medium,
            ecuts,
            false,
            parametrization);

        dEdx_new = cross->CalculatedEdx(energy) * medium->GetMassDensity();

        ASSERT_NEAR(dEdx_new, dEdx_stored, 1e-3 * dEdx_stored);
    }
}

TEST(Compton, Test_of_dNdx)
{
    std::string filename = testfile_dir + "Compton_dNdx.txt";
    std::ifstream in{filename};
    EXPECT_TRUE(in.good()) << "Test resource file '" << filename << "' could not be opened";

    std::string mediumName;
    double ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    std::string parametrization;
    double energy;
    double dNdx_stored;
    double dNdx_new;

    std::cout.precision(16);

    while (in >> mediumName >> ecut >> vcut >> multiplier >> energy >> dNdx_stored >> parametrization)
    {
        parametrization.erase(0,7);

        if (vcut == -1)
            vcut = 1;
        if (ecut == -1)
            ecut = INF;

        ParticleDef particle_def = GammaDef();
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        auto cross = crosssection::make_compton(particle_def,
            *medium,
            ecuts,
            false,
            parametrization);

        dNdx_new = cross->CalculatedNdx(energy) * medium->GetMassDensity();

        ASSERT_NEAR(dNdx_new, dNdx_stored, 1e-3 * dNdx_stored);
    }
}

TEST(Compton, Test_of_e)
{
    std::string filename = testfile_dir + "Compton_e.txt";
    std::ifstream in{filename};
    EXPECT_TRUE(in.good()) << "Test resource file '" << filename << "' could not be opened";

    std::string mediumName;
    double ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    std::string parametrization;
    double energy;
    double rnd1;
    double rnd2;
    double stochastic_loss_stored;
    double stochastic_loss_new;

    std::cout.precision(16);

    RandomGenerator::Get().SetSeed(0);

    while (in >> mediumName >> ecut >> vcut >> multiplier >> energy >> rnd1 >> rnd2 >> stochastic_loss_stored >> parametrization)
    {
        parametrization.erase(0,7);

        if (vcut == -1)
            vcut = 1;
        if (ecut == -1)
            ecut = INF;

        ParticleDef particle_def = GammaDef();
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        auto cross = crosssection::make_compton(particle_def,
            *medium,
            ecuts,
            false,
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
                    #ifndef NDEBUG
                    EXPECT_DEATH(cross->CalculateStochasticLoss(comp, energy, rate_new), "");
                    #endif
                } else {
                    stochastic_loss_new = energy * cross->CalculateStochasticLoss(comp, energy, rate_new);
                    EXPECT_NEAR(stochastic_loss_new, stochastic_loss_stored, 1E-6 * stochastic_loss_stored);
                    break;
                }
            }
        }
    }
}

TEST(Compton, Test_of_dEdx_Interpolant)
{
    std::string filename = testfile_dir + "Compton_dEdx_interpol.txt";
    std::ifstream in{filename};
    EXPECT_TRUE(in.good()) << "Test resource file '" << filename << "' could not be opened";

    std::string mediumName;
    double ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    std::string parametrization;
    double energy;
    double dEdx_stored;
    double dEdx_new;

    std::cout.precision(16);

    while (in >> mediumName >> ecut >> vcut >> multiplier >> energy >> dEdx_stored >> parametrization)
    {
        parametrization.erase(0,7);

        if (vcut == -1)
            vcut = 1;
        if (ecut == -1)
            ecut = INF;

        ParticleDef particle_def = GammaDef();
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        auto cross = crosssection::make_compton(particle_def,
            *medium,
            ecuts,
            false,
            parametrization);

        dEdx_new = cross->CalculatedEdx(energy) * medium->GetMassDensity();

        ASSERT_NEAR(dEdx_new, dEdx_stored, 1e-3 * dEdx_stored);
    }
}

TEST(Compton, Test_of_dNdx_Interpolant)
{
    std::string filename = testfile_dir + "Compton_dNdx_interpol.txt";
    std::ifstream in{filename};
    EXPECT_TRUE(in.good()) << "Test resource file '" << filename << "' could not be opened";

    std::string mediumName;
    double ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    std::string parametrization;
    double energy;
    double dNdx_stored;
    double dNdx_new;

    std::cout.precision(16);

    while (in >> mediumName >> ecut >> vcut >> multiplier >> energy >> dNdx_stored >> parametrization)
    {
        parametrization.erase(0,7);

        if (vcut == -1)
            vcut = 1;
        if (ecut == -1)
            ecut = INF;

        ParticleDef particle_def = GammaDef();
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        auto cross = crosssection::make_compton(particle_def,
            *medium,
            ecuts,
            false,
            parametrization);

        dNdx_new = cross->CalculatedNdx(energy) * medium->GetMassDensity();

        EXPECT_NEAR(dNdx_new, dNdx_stored, 1e-3 * dNdx_stored);
    }
}

TEST(Compton, Test_of_e_Interpolant)
{
    std::string filename = testfile_dir + "Compton_e_interpol.txt";
    std::ifstream in{filename};
    EXPECT_TRUE(in.good()) << "Test resource file '" << filename << "' could not be opened";

    std::string mediumName;
    double ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    std::string parametrization;
    double energy;
    double rnd1;
    double rnd2;
    double stochastic_loss_stored;
    double stochastic_loss_new;

    std::cout.precision(16);

    RandomGenerator::Get().SetSeed(0);

    while (in >> mediumName >> ecut >> vcut >> multiplier >> energy >> rnd1 >> rnd2 >> stochastic_loss_stored >> parametrization)
    {
        parametrization.erase(0,7);

        if (vcut == -1)
            vcut = 1;
        if (ecut == -1)
            ecut = INF;

        ParticleDef particle_def = GammaDef();
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        auto cross = crosssection::make_compton(particle_def,
            *medium,
            ecuts,
            false,
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
                    #ifndef NDEBUG
                    EXPECT_DEATH(cross->CalculateStochasticLoss(comp, energy, rate_new), "");
                    #endif
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
