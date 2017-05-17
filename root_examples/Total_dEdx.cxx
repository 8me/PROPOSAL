/** \file
 * This is an example to calculate dE/dx for every cross section.\n
 * If you have ROOT installed this file will be compilied when running make\n
 * and you can find an executable in your build directory in root_examples:\n
 * -->   Total_dEdx   <--\n
 * A ROOT-file named Total_dEdx.root will be generated,
 * which contains a TTree with dE/dx for every parametrization and particle type
 * and three different media: ice, hydrogen and uranium.
 * @brief Example to calculate dE/dx for every parametrization
 * @author Jan-Hendrik Köhne
 */

#include "PROPOSAL/Bremsstrahlung.h"
#include "PROPOSAL/Photonuclear.h"
#include "PROPOSAL/Ionization.h"
#include "PROPOSAL/Epairproduction.h"
#include "PROPOSAL/PROPOSALParticle.h"
#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/Medium.h"
#include "PROPOSAL/CrossSections.h"
#include "TFile.h"
#include "TTree.h"
#include <cmath>
#include <sstream>
#include "TCanvas.h"
#include "TGraph.h"
#include "TMultiGraph.h"

using namespace std;
using namespace PROPOSAL;

int main()
{
    cout<<"\n------------------------------------------------------------------\n"
        <<"This is an example to calculate dE/dx for every cross section.\n"
        <<"A ROOT-file named Total_dEdx.root will be generated,\n"
        <<"which contains a TTree with dE/dx for every parametrization and particle type\n"
        <<"and three different media: ice, hydrogen and uranium.\n"
        <<"------------------------------------------------------------------\n"
        <<endl;


    TFile *file     =   new TFile("Total_dEdx.root","RECREATE");

    PROPOSALParticle *mu    =   new PROPOSALParticle(ParticleType::MuMinus);
    PROPOSALParticle *tau   =   new PROPOSALParticle(ParticleType::TauMinus);
    PROPOSALParticle *e     =   new PROPOSALParticle(ParticleType::EMinus);

    Medium  *med1   =   new Medium("hydrogen",1.);
    Medium  *med2   =   new Medium("ice",1.);
    Medium  *med3   =   new Medium("uranium",1.);

    EnergyCutSettings*  cuts    =   new EnergyCutSettings(-1,-1);

    vector<CrossSections*> cros;

    cros.push_back(new Bremsstrahlung(mu,med1,cuts));
    cros.push_back(new Bremsstrahlung(mu,med2,cuts));
    cros.push_back(new Bremsstrahlung(mu,med3,cuts));
    cros.push_back(new Bremsstrahlung(tau,med1,cuts));
    cros.push_back(new Bremsstrahlung(tau,med2,cuts));
    cros.push_back(new Bremsstrahlung(tau,med3,cuts));
    cros.push_back(new Bremsstrahlung(e,med1,cuts));
    cros.push_back(new Bremsstrahlung(e,med2,cuts));
    cros.push_back(new Bremsstrahlung(e,med3,cuts));

    cros.push_back(new Bremsstrahlung(mu,med1,cuts));
    cros.push_back(new Bremsstrahlung(mu,med2,cuts));
    cros.push_back(new Bremsstrahlung(mu,med3,cuts));
    cros.push_back(new Bremsstrahlung(tau,med1,cuts));
    cros.push_back(new Bremsstrahlung(tau,med2,cuts));
    cros.push_back(new Bremsstrahlung(tau,med3,cuts));
    cros.push_back(new Bremsstrahlung(e,med1,cuts));
    cros.push_back(new Bremsstrahlung(e,med2,cuts));
    cros.push_back(new Bremsstrahlung(e,med3,cuts));

    cros.push_back(new Bremsstrahlung(mu,med1,cuts));
    cros.push_back(new Bremsstrahlung(mu,med2,cuts));
    cros.push_back(new Bremsstrahlung(mu,med3,cuts));
    cros.push_back(new Bremsstrahlung(tau,med1,cuts));
    cros.push_back(new Bremsstrahlung(tau,med2,cuts));
    cros.push_back(new Bremsstrahlung(tau,med3,cuts));
    cros.push_back(new Bremsstrahlung(e,med1,cuts));
    cros.push_back(new Bremsstrahlung(e,med2,cuts));
    cros.push_back(new Bremsstrahlung(e,med3,cuts));

    cros.push_back(new Bremsstrahlung(mu,med1,cuts));
    cros.push_back(new Bremsstrahlung(mu,med2,cuts));
    cros.push_back(new Bremsstrahlung(mu,med3,cuts));
    cros.push_back(new Bremsstrahlung(tau,med1,cuts));
    cros.push_back(new Bremsstrahlung(tau,med2,cuts));
    cros.push_back(new Bremsstrahlung(tau,med3,cuts));
    cros.push_back(new Bremsstrahlung(e,med1,cuts));
    cros.push_back(new Bremsstrahlung(e,med2,cuts));
    cros.push_back(new Bremsstrahlung(e,med3,cuts));

    for(unsigned int i = 0; i<cros.size() ; i++)
    {
        if(i<9) cros.at(i)->SetParametrization(ParametrizationType::BremsKelnerKokoulinPetrukhin);
        else if(i<18) cros.at(i)->SetParametrization(ParametrizationType::BremsAndreevBezrukovBugaev);
        else if(i<27) cros.at(i)->SetParametrization(ParametrizationType::BremsPetrukhinShestakov);
        else if(i<36) cros.at(i)->SetParametrization(ParametrizationType::BremsCompleteScreeningCase);
    }

    cros.push_back(new Epairproduction(mu,med1,cuts));
    cros.push_back(new Epairproduction(mu,med2,cuts));
    cros.push_back(new Epairproduction(mu,med3,cuts));
    cros.push_back(new Epairproduction(tau,med1,cuts));
    cros.push_back(new Epairproduction(tau,med2,cuts));
    cros.push_back(new Epairproduction(tau,med3,cuts));
    cros.push_back(new Epairproduction(e,med1,cuts));
    cros.push_back(new Epairproduction(e,med2,cuts));
    cros.push_back(new Epairproduction(e,med3,cuts));

    cros.push_back(new Ionization(mu,med1,cuts));
    cros.push_back(new Ionization(mu,med2,cuts));
    cros.push_back(new Ionization(mu,med3,cuts));
    cros.push_back(new Ionization(tau,med1,cuts));
    cros.push_back(new Ionization(tau,med2,cuts));
    cros.push_back(new Ionization(tau,med3,cuts));
    cros.push_back(new Ionization(e,med1,cuts));
    cros.push_back(new Ionization(e,med2,cuts));
    cros.push_back(new Ionization(e,med3,cuts));

    vector<Photonuclear*> photo;

    for(int i =0 ; i<14; i++)
    {
        photo.push_back(new Photonuclear(mu,med1,cuts));
        photo.push_back(new Photonuclear(mu,med2,cuts));
        photo.push_back(new Photonuclear(mu,med3,cuts));
        photo.push_back(new Photonuclear(tau,med1,cuts));
        photo.push_back(new Photonuclear(tau,med2,cuts));
        photo.push_back(new Photonuclear(tau,med3,cuts));
        photo.push_back(new Photonuclear(e,med1,cuts));
        photo.push_back(new Photonuclear(e,med2,cuts));
        photo.push_back(new Photonuclear(e,med3,cuts));
    }

    for(unsigned int i = 0; i<photo.size() ; i++)
    {
        if(i<9) photo.at(i)->SetParametrization(ParametrizationType::PhotoKokoulinShadowBezrukovSoft);
        else if(i<18) photo.at(i)->SetParametrization(ParametrizationType::PhotoKokoulinShadowBezrukovHard);
        else if(i<27) photo.at(i)->SetParametrization(ParametrizationType::PhotoRhodeShadowBezrukovSoft);
        else if(i<36) photo.at(i)->SetParametrization(ParametrizationType::PhotoRhodeShadowBezrukovHard);
        else if(i<45) photo.at(i)->SetParametrization(ParametrizationType::PhotoBezrukovBugaevShadowBezrukovSoft);
        else if(i<54) photo.at(i)->SetParametrization(ParametrizationType::PhotoBezrukovBugaevShadowBezrukovHard);
        else if(i<63) photo.at(i)->SetParametrization(ParametrizationType::PhotoZeusShadowBezrukovSoft);
        else if(i<72) photo.at(i)->SetParametrization(ParametrizationType::PhotoZeusShadowBezrukovHard);
        else if(i<81) photo.at(i)->SetParametrization(ParametrizationType::PhotoAbramowiczLevinLevyMaor91ShadowDutta);
        else if(i<90) photo.at(i)->SetParametrization(ParametrizationType::PhotoAbramowiczLevinLevyMaor91ShadowButkevich);
        else if(i<99) photo.at(i)->SetParametrization(ParametrizationType::PhotoAbramowiczLevinLevyMaor97ShadowDutta);
        else if(i<108) photo.at(i)->SetParametrization(ParametrizationType::PhotoAbramowiczLevinLevyMaor97ShadowButkevich);
        else if(i<117) photo.at(i)->SetParametrization(ParametrizationType::PhotoButkevichMikhailovShadowDutta);
        else if(i<126) photo.at(i)->SetParametrization(ParametrizationType::PhotoButkevichMikhailovShadowButkevich);
    }

    for(unsigned int i =0; i < photo.size();i++)
    {
        cros.push_back(photo.at(i));
    }


    double energy;
    vector<double> dEdx;

    TTree *tree = new TTree("dEdx","dEdx");

    tree->Branch("energy",&energy,"energy/D");

    stringstream ss;
    stringstream ss1;

    dEdx.resize(cros.size());

    for(unsigned int i = 0; i < cros.size(); i++)
    {
        cout<<i<<endl;
        cros.at(i)->EnableDEdxInterpolation("resources/tables",true);

        ss<<cros.at(i)->GetName();

        if(cros.at(i)->GetName().compare("Bremsstrahlung")==0 || cros.at(i)->GetName().compare("Photonuclear")==0)
        {
            ss<<"_"<<cros.at(i)->GetParametrization();
        }

        ss<<"_"<<cros.at(i)->GetParticle()->GetName()<<"_"<<cros.at(i)->GetMedium()->GetName();
        ss1<<ss.str()<<"/D";

        tree->Branch(ss.str().c_str(),&dEdx.at(i),ss1.str().c_str());

        ss.str("");
        ss.clear();
        ss1.str("");
        ss1.clear();
    }

    for(double log_energy = 0 ;  log_energy < 11; log_energy +=0.01)
    {
        for(unsigned int i = 0 ; i < cros.size() ; i++)
        {
            energy  =   cros.at(i)->GetParticle()->GetMass()   +   pow(10,log_energy);

            cros.at(i)->GetParticle()->SetEnergy(energy);

            dEdx.at(i)  =   cros.at(i)->CalculatedEdx();

            tree->Fill();

        }

    }

    vector<TGraph*> mu_ice_graphs;
    vector<TGraph*> mu_uranium_graphs;
    vector<TGraph*> mu_hydrogen_graphs;
    vector<TGraph*> e_ice_graphs;
    vector<TGraph*> e_uranium_graphs;
    vector<TGraph*> e_hydrogen_graphs;
    vector<TGraph*> tau_ice_graphs;
    vector<TGraph*> tau_uranium_graphs;
    vector<TGraph*> tau_hydrogen_graphs;

    TMultiGraph* mu_ice_graph;
    TMultiGraph* mu_uranium_graph;
    TMultiGraph* mu_hydrogen_graph;
    TMultiGraph* e_ice_graph;
    TMultiGraph* e_uranium_graph;
    TMultiGraph* e_hydrogen_graph;
    TMultiGraph* tau_ice_graph;
    TMultiGraph* tau_uranium_graph;
    TMultiGraph* tau_hydrogen_graph;

    double Ioniz;
    double Epair;
    double Photo;
    double Brems;
    double Total_dEdx;

    int point   =   0;

    stringstream ss2;
    stringstream graph_name;

    tree->SetBranchAddress("energy",&energy);


    vector<ParametrizationType::Enum> brems_parametrization_vector;
    brems_parametrization_vector.push_back(ParametrizationType::BremsKelnerKokoulinPetrukhin);
    brems_parametrization_vector.push_back(ParametrizationType::BremsAndreevBezrukovBugaev);
    brems_parametrization_vector.push_back(ParametrizationType::BremsPetrukhinShestakov);
    brems_parametrization_vector.push_back(ParametrizationType::BremsCompleteScreeningCase);

    vector<ParametrizationType::Enum> photo_parametrization_vector;
    photo_parametrization_vector.push_back(ParametrizationType::PhotoKokoulinShadowBezrukovSoft);
    photo_parametrization_vector.push_back(ParametrizationType::PhotoKokoulinShadowBezrukovHard);
    photo_parametrization_vector.push_back(ParametrizationType::PhotoRhodeShadowBezrukovSoft);
    photo_parametrization_vector.push_back(ParametrizationType::PhotoRhodeShadowBezrukovHard);
    photo_parametrization_vector.push_back(ParametrizationType::PhotoBezrukovBugaevShadowBezrukovSoft);
    photo_parametrization_vector.push_back(ParametrizationType::PhotoBezrukovBugaevShadowBezrukovHard);
    photo_parametrization_vector.push_back(ParametrizationType::PhotoZeusShadowBezrukovSoft);
    photo_parametrization_vector.push_back(ParametrizationType::PhotoZeusShadowBezrukovHard);
    photo_parametrization_vector.push_back(ParametrizationType::PhotoAbramowiczLevinLevyMaor91ShadowDutta);
    photo_parametrization_vector.push_back(ParametrizationType::PhotoAbramowiczLevinLevyMaor91ShadowButkevich);
    photo_parametrization_vector.push_back(ParametrizationType::PhotoAbramowiczLevinLevyMaor97ShadowDutta);
    photo_parametrization_vector.push_back(ParametrizationType::PhotoAbramowiczLevinLevyMaor97ShadowButkevich);
    photo_parametrization_vector.push_back(ParametrizationType::PhotoButkevichMikhailovShadowDutta);
    photo_parametrization_vector.push_back(ParametrizationType::PhotoButkevichMikhailovShadowButkevich);

    //Muons ice
    for(unsigned int j = 0; j < photo_parametrization_vector.size(); j++)
    {
        tree->SetBranchStatus("*",0);
        ss<<"Photonuclear_"<<photo_parametrization_vector.at(j)<<"_mu_ice";

        tree->SetBranchAddress(ss.str().c_str(),&Photo);
        tree->SetBranchStatus(ss.str().c_str(),1);
        tree->SetBranchStatus("energy",1);

        ss.str("");
        ss.clear();

        for(unsigned int k = 0; k < brems_parametrization_vector.size(); k++)
        {
            ss<<"Bremsstrahlung_"<<brems_parametrization_vector.at(k)<<"_mu_ice";
            tree->SetBranchAddress(ss.str().c_str(),&Brems);
            tree->SetBranchAddress("Ionization_mu_ice",&Ioniz);
            tree->SetBranchAddress("Epairproduction_mu_ice",&Epair);

            tree->SetBranchStatus(ss.str().c_str(),1);
            tree->SetBranchStatus("Ionization_mu_ice",1);
            tree->SetBranchStatus("Epairproduction_mu_ice",1);

            graph_name<<"mu_ice_photo_"<<photo_parametrization_vector.at(j)<<"_brems_"<<brems_parametrization_vector.at(k);

            mu_ice_graphs.push_back(new TGraph(graph_name.str().c_str(),graph_name.str().c_str()));

            if(photo_parametrization_vector.at(j) == ParametrizationType::PhotoAbramowiczLevinLevyMaor97ShadowButkevich
                && brems_parametrization_vector.at(k) == ParametrizationType::BremsKelnerKokoulinPetrukhin)
            {
                mu_ice_graphs.back()->SetLineColor(kBlack);
            }
            else
            {
                mu_ice_graphs.back()->SetLineColor(kRed);
            }

            point=0;

            for(unsigned int l =0 ; l < tree->GetEntries(); l++)
            {
                tree->GetEntry(l);
                Total_dEdx = Ioniz + Brems + Epair + Photo;

                mu_ice_graphs.back()->SetPoint(point,energy,Total_dEdx);
                point++;
            }
            mu_ice_graph->Add(mu_ice_graphs.back(),"l");
            ss.str("");
            ss.clear();
        }
    }

    TCanvas *mu_ice = new TCanvas("mu_ice","Muons in ice",1400,1050);
    mu_ice->cd();
    mu_ice_graph->Draw("A");


    mu_ice->Write();
    tree->Write();
    file->Close();

}
