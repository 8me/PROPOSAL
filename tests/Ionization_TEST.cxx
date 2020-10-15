
#include "gtest/gtest.h"

#include <fstream>
// #include "PROPOSAL/Constants.h"
#include "PROPOSAL/crosssection/parametrization/Ionization.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/medium/MediumFactory.h"
// #include "PROPOSAL/methods.h"

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
//     ParticleDef particle_def = MuMinusDef();
//     auto medium = std::make_shared<const Water>();
//     EnergyCutSettings ecuts;
//     double multiplier = 1.;

//     IonizBetheBlochRossi* Ioniz_A = new IonizBetheBlochRossi(particle_def, medium, ecuts, multiplier);
//     Parametrization* Ioniz_B = new IonizBetheBlochRossi(particle_def, medium, ecuts, multiplier);
//     EXPECT_TRUE(*Ioniz_A == *Ioniz_B);

//     IonizBetheBlochRossi param(particle_def, medium, ecuts, multiplier);
//     EXPECT_TRUE(param == *Ioniz_A);

//     IonizIntegral* Int_A        = new IonizIntegral(param);
//     CrossSectionIntegral* Int_B = new IonizIntegral(param);
//     EXPECT_TRUE(*Int_A == *Int_B);

//     InterpolationDef InterpolDef;
//     IonizInterpolant* Interpol_A        = new IonizInterpolant(param, InterpolDef);
//     CrossSectionInterpolant* Interpol_B = new IonizInterpolant(param, InterpolDef);
//     EXPECT_TRUE(*Interpol_A == *Interpol_B);

//     delete Ioniz_A;
//     delete Ioniz_B;
//     delete Int_A;
//     delete Int_B;
//     delete Interpol_A;
//     delete Interpol_B;
// }

// TEST(Comparison, Comparison_not_equal)
// {
//     ParticleDef mu_def  = MuMinusDef();
//     ParticleDef tau_def = TauMinusDef();
//     auto medium_1 = std::make_shared<const Water>();
//     auto medium_2 = std::make_shared<const Ice>();
//     EnergyCutSettings ecuts_1(500, -1);
//     EnergyCutSettings ecuts_2(-1, 0.05);
//     double multiplier_1 = 1.;
//     double multiplier_2 = 2.;

//     IonizBetheBlochRossi Ioniz_A(mu_def, medium_1, ecuts_1, multiplier_1);
//     IonizBetheBlochRossi Ioniz_B(tau_def, medium_1, ecuts_1, multiplier_1);
//     IonizBetheBlochRossi Ioniz_C(mu_def, medium_2, ecuts_1, multiplier_1);
//     IonizBetheBlochRossi Ioniz_D(mu_def, medium_1, ecuts_2, multiplier_1);
//     IonizBetheBlochRossi Ioniz_E(mu_def, medium_1, ecuts_1, multiplier_2);
//     EXPECT_TRUE(Ioniz_A != Ioniz_B);
//     EXPECT_TRUE(Ioniz_A != Ioniz_C);
//     EXPECT_TRUE(Ioniz_A != Ioniz_D);
//     EXPECT_TRUE(Ioniz_A != Ioniz_E);

//     IonizIntegral Int_A(Ioniz_A);
//     IonizIntegral Int_B(Ioniz_B);
//     EXPECT_TRUE(Int_A != Int_B);

//     InterpolationDef InterpolDef;
//     IonizInterpolant Interpol_A(Ioniz_A, InterpolDef);
//     IonizInterpolant Interpol_B(Ioniz_B, InterpolDef);
//     EXPECT_TRUE(Interpol_A != Interpol_B);
// }

// TEST(Assignment, Copyconstructor)
// {
//     ParticleDef mu_def = MuMinusDef();
//     auto medium = std::make_shared<const Water>();
//     EnergyCutSettings ecuts(500, -1);
//     double multiplier = 1.;

//     IonizBetheBlochRossi Ioniz_A(mu_def, medium, ecuts, multiplier);
//     IonizBetheBlochRossi Ioniz_B = Ioniz_A;

//     IonizIntegral Int_A(Ioniz_A);
//     IonizIntegral Int_B = Int_A;
//     EXPECT_TRUE(Int_A == Int_B);

//     InterpolationDef InterpolDef;
//     IonizInterpolant Interpol_A(Ioniz_A, InterpolDef);
//     IonizInterpolant Interpol_B = Interpol_A;
//     EXPECT_TRUE(Interpol_A == Interpol_B);
// }

// TEST(Assignment, Copyconstructor2)
// {

//     ParticleDef mu_def = MuMinusDef();
//     auto medium = std::make_shared<const Water>();
//     EnergyCutSettings ecuts(500, -1);
//     double multiplier = 1.;

//     IonizBetheBlochRossi Ioniz_A(mu_def, medium, ecuts, multiplier);
//     IonizBetheBlochRossi Ioniz_B(Ioniz_A);

//     IonizIntegral Int_A(Ioniz_A);
//     IonizIntegral Int_B(Int_A);
//     EXPECT_TRUE(Int_A == Int_B);

//     InterpolationDef InterpolDef;
//     IonizInterpolant Interpol_A(Ioniz_A, InterpolDef);
//     IonizInterpolant Interpol_B(Interpol_A);
//     EXPECT_TRUE(Interpol_A == Interpol_B);
// }

// in polymorphism an assignmant and swap operator doesn't make sense

TEST(Ionization, Test_of_dEdx)
{
    std::string filename = testfile_dir + "Ioniz_dEdx.txt";
	std::ifstream in{filename};
	EXPECT_TRUE(in.good()) << "Test resource file '" << filename << "' could not be opened";

    char firstLine[256];
    in.getline(firstLine, 256);

    std::string particleName;
    std::string mediumName;
    double ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    std::string parametrization;
    double energy;
    double dEdx_stored;
    double dEdx_new;

    while (in.good())
    {
        in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> energy >> dEdx_stored >> parametrization;

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        auto cross = crosssection::make_ionization(particle_def,
            *medium,
            ecuts,
            false,
            parametrization);

        dEdx_new = cross->CalculatedEdx(energy);

        ASSERT_NEAR(dEdx_new, dEdx_stored, 1e-10 * dEdx_stored);

    }
}

TEST(Ionization, Test_of_dNdx)
{
    std::string filename = testfile_dir + "Ioniz_dNdx.txt";
    std::ifstream in{filename};
    EXPECT_TRUE(in.good()) << "Test resource file '" << filename << "' could not be opened";

    char firstLine[256];
    in.getline(firstLine, 256);

    std::string particleName;
    std::string mediumName;
    double ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    std::string parametrization;
    double energy;
    double dNdx_stored;
    double dNdx_new;

    while (in.good())
    {
        in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> energy >> dNdx_stored >> parametrization;

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        auto cross = crosssection::make_ionization(particle_def,
            *medium,
            ecuts,
            false,
            parametrization);

        dNdx_new = cross->CalculatedNdx(energy);

        ASSERT_NEAR(dNdx_new, dNdx_stored, 1e-10 * dNdx_stored);
    }
}

TEST(Ionization, Test_Stochastic_Loss)
{
    std::string filename = testfile_dir + "Ioniz_e.txt";
    std::ifstream in{filename};
    EXPECT_TRUE(in.good()) << "Test resource file '" << filename << "' could not be opened";

    char firstLine[256];
    in.getline(firstLine, 256);

    std::string particleName;
    std::string mediumName;
    double ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    std::string parametrization;
    double energy;
    double rate;
    double rnd1, rnd2;
    double stochastic_loss_stored;
    double stochastic_loss_new;

    RandomGenerator::Get().SetSeed(0);

    while (in.good())
    {
        in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> energy >> rnd1 >> rnd2 >>
            stochastic_loss_stored >> parametrization;

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        auto cross = crosssection::make_ionization(particle_def,
            *medium,
            ecuts,
            false,
            parametrization);

        auto components = medium->GetComponents();
        for (auto comp : components)
        {
            auto tmp = std::make_shared<const Component>(comp);
            rate = cross->CalculatedNdx(energy, tmp);
            stochastic_loss_new = cross->CalculateStochasticLoss(
                tmp, energy, rate);

            ASSERT_NEAR(stochastic_loss_new, stochastic_loss_stored, 1E-6 * stochastic_loss_stored);
        }

    }
}

TEST(Ionization, Test_of_dEdx_Interpolant)
{
    std::string filename = testfile_dir + "Ioniz_dEdx_interpol.txt";
	std::ifstream in{filename};
	EXPECT_TRUE(in.good()) << "Test resource file '" << filename << "' could not be opened";

    char firstLine[256];
    in.getline(firstLine, 256);

    std::string particleName;
    std::string mediumName;
    double ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    std::string parametrization;
    double energy;
    double dEdx_stored;
    double dEdx_new;

    while (in.good())
    {
        in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> energy >> dEdx_stored >> parametrization;

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        auto cross = crosssection::make_ionization(particle_def,
            *medium,
            ecuts,
            true,
            parametrization);

        dEdx_new = cross->CalculatedEdx(energy);

        ASSERT_NEAR(dEdx_new, dEdx_stored, 1e-10 * dEdx_stored);

    }
}

TEST(Ionization, Test_of_dNdx_Interpolant)
{
    std::string filename = testfile_dir + "Ioniz_dNdx_interpol.txt";
	std::ifstream in{filename};
	EXPECT_TRUE(in.good()) << "Test resource file '" << filename << "' could not be opened";

    char firstLine[256];
    in.getline(firstLine, 256);

    std::string particleName;
    std::string mediumName;
    double ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    std::string parametrization;
    double energy;
    double dNdx_stored;
    double dNdx_new;

    // InterpolationDef InterpolDef;

    while (in.good())
    {
        in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> energy >> dNdx_stored >> parametrization;

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        auto cross = crosssection::make_ionization(particle_def,
            *medium,
            ecuts,
            true,
            parametrization);

        dNdx_new = cross->CalculatedNdx(energy);

        ASSERT_NEAR(dNdx_new, dNdx_stored, 1e-10 * dNdx_stored);
    }
}

TEST(Ionization, Test_of_e_interpol)
{
    std::string filename = testfile_dir + "Ioniz_e_interpol.txt";
	std::ifstream in{filename};
	EXPECT_TRUE(in.good()) << "Test resource file '" << filename << "' could not be opened";

    char firstLine[256];
    in.getline(firstLine, 256);

    std::string particleName;
    std::string mediumName;
    double ecut;
    double vcut;
    bool cont_rand = false;
    double multiplier;
    std::string parametrization;
    double energy;
    double rate;
    double rnd1, rnd2;
    double stochastic_loss_stored;
    double stochastic_loss_new;

    // InterpolationDef InterpolDef;

    RandomGenerator::Get().SetSeed(0);

    while (in.good())
    {
        in >> particleName >> mediumName >> ecut >> vcut >> multiplier >> energy >> rnd1 >> rnd2 >>
            stochastic_loss_stored >> parametrization;

        ParticleDef particle_def = getParticleDef(particleName);
        auto medium = CreateMedium(mediumName);
        auto ecuts = std::make_shared<EnergyCutSettings>(ecut, vcut, cont_rand);

        auto cross = crosssection::make_ionization(particle_def,
            *medium,
            ecuts,
            true,
            parametrization);

        auto components = medium->GetComponents();
        for (auto comp : components)
        {
            auto tmp = std::make_shared<const Component>(comp);
            rate = cross->CalculatedNdx(energy, tmp);
            stochastic_loss_new = cross->CalculateStochasticLoss(
                tmp, energy, rate);

            ASSERT_NEAR(stochastic_loss_new, stochastic_loss_stored, 1E-6 * stochastic_loss_stored);
        }

    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
