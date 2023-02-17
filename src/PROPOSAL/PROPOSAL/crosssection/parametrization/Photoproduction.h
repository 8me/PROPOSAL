#pragma once

#include "PROPOSAL/crosssection/parametrization/ParametrizationDirect.h"
#include "PROPOSAL/crosssection/parametrization/Parametrization.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/methods.h"

namespace PROPOSAL {
    class Component;

    namespace crosssection {
        class Photoproduction : public ParametrizationDirect {
            virtual double PhotonAtomCrossSection(double, const Component&);
            double ShadowingFactor(double, const Component&);
        protected:
            double GetCutOff(const Component& comp) const;
        public:
            Photoproduction() = default;

            virtual double PhotonNucleonCrossSection(double, const Component&) = 0;

            double CalculatedNdx(double, const ParticleDef&, const Medium&, cut_ptr) override;
            double CalculatedNdx(double, size_t, const ParticleDef&, const Medium&, cut_ptr) override;
            std::vector<std::pair<size_t, double>> CalculatedNdx_PerTarget(
                    double, const ParticleDef&, const Medium&, cut_ptr) override;

            // no continuous losses
            double CalculatedEdx(double, const ParticleDef&, const Medium&, cut_ptr) override { return 0.; };
            double CalculatedE2dx(double, const ParticleDef&, const Medium&, cut_ptr) override { return 0.; };

            double CalculateCumulativeCrosssection(double, size_t, double, const ParticleDef&, const Medium&, cut_ptr) override {
                throw std::logic_error("No cumulative crosssection defined for Photoproduction.");
            };

            double CalculateStochasticLoss(size_t, double, double, const ParticleDef&, const Medium&, cut_ptr) override {
                return 1.; // all energy is always lost, i.e. v=1
            };

            double GetLowerEnergyLim(const ParticleDef&, const Medium&, cut_ptr) const override;

            size_t GetHash(const ParticleDef&, const Medium& m, cut_ptr) const noexcept override;

            InteractionType GetInteractionType() const noexcept override;
        };

        class PhotoproductionZeus : public Photoproduction {
        public:
            PhotoproductionZeus();
            std::unique_ptr<ParametrizationDirect> clone() const final;
            double PhotonNucleonCrossSection(double, const Component&) override;
        };

        template <> struct ParametrizationName<PhotoproductionZeus> {
            static constexpr auto value = "Zeus";
        };

        class PhotoproductionBezrukovBugaev : virtual public Photoproduction {
        public:
            PhotoproductionBezrukovBugaev();
            std::unique_ptr<ParametrizationDirect> clone() const override;
            double PhotonNucleonCrossSection(double, const Component&) override;
        };

        template <> struct ParametrizationName<PhotoproductionBezrukovBugaev> {
            static constexpr auto value = "BezrukovBugaev";
        };

        class PhotoproductionCaldwell : virtual public Photoproduction {
        public:
            PhotoproductionCaldwell();
            std::unique_ptr<ParametrizationDirect> clone() const override;
            double PhotonNucleonCrossSection(double, const Component&) override;
        };

        template <> struct ParametrizationName<PhotoproductionCaldwell> {
            static constexpr auto value = "Caldwell";
        };


        class PhotoproductionKokoulin : public PhotoproductionBezrukovBugaev, public PhotoproductionCaldwell {
        public:
            PhotoproductionKokoulin();
            std::unique_ptr<ParametrizationDirect> clone() const final;
            double PhotonNucleonCrossSection(double, const Component&) override;
        };

        template <> struct ParametrizationName<PhotoproductionKokoulin> {
            static constexpr auto value = "Kokoulin";
        };

        class PhotoproductionRhode : public PhotoproductionCaldwell {
            std::shared_ptr<Interpolant> interpolant_;
        public:
            PhotoproductionRhode();
            std::unique_ptr<ParametrizationDirect> clone() const override;
            double PhotonNucleonCrossSection(double, const Component&) override;
        };

        template <> struct ParametrizationName<PhotoproductionRhode> {
            static constexpr auto value = "Rhode";
        };

        class PhotoproductionHeck : virtual public Photoproduction {
            double sigma_bw(double s, double E, double M, double Gamma, double Sigma0);
            double Q(double E, double E_th, double w);
            const std::array<double, 3> resonances_M = {1.231, 1.515, 1.680};
            const std::array<double, 3> resonances_Gamma = {0.11, 0.11, 0.125};
            const std::array<double, 3> resonances_Sigma0 = {31.125, 25.567, 17.508};
            const std::array<double, 3> resonances_w = {0.17, 0.38, 0.38};

        public:
            PhotoproductionHeck();
            std::unique_ptr<ParametrizationDirect> clone() const override;
            double PhotonNucleonCrossSection(double, const Component&) override;
        };

        template <> struct ParametrizationName<PhotoproductionHeck> {
            static constexpr auto value = "Heck";
        };

        class PhotoproductionHeckC7Shadowing : virtual public PhotoproductionHeck {
        public:
            PhotoproductionHeckC7Shadowing();
            double PhotonAtomCrossSection(double energy, const Component& comp) override;
            std::unique_ptr<ParametrizationDirect> clone() const override;
        };

        template <> struct ParametrizationName<PhotoproductionHeckC7Shadowing> {
            static constexpr auto value = "HeckC7Shadowing";
        };

    }
}

