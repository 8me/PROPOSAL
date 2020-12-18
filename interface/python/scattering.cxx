
#include "PROPOSAL/scattering/Scattering.h"
#include "PROPOSAL/propagation_utility/PropagationUtility.h"
#include "PROPOSAL/scattering/multiple_scattering/Highland.h"
#include "PROPOSAL/scattering/multiple_scattering/HighlandIntegral.h"
#include "PROPOSAL/scattering/multiple_scattering/Moliere.h"
#include "PROPOSAL/scattering/multiple_scattering/ScatteringFactory.h"

#include "PROPOSAL/scattering/stochastic_deflection/bremsstrahlung/NaivBremsstrahlung.h"

#include "pyBindings.h"

#include <pybind11/stl.h>

namespace py = pybind11;
using namespace PROPOSAL;

void init_scattering(py::module& m)
{
    py::module m_sub = m.def_submodule("scattering");

    py::class_<multiple_scattering::Parametrization,
        std::shared_ptr<multiple_scattering::Parametrization>>(
        m_sub, "MultipleScattering")
        .def("scatter", &multiple_scattering::Parametrization::Scatter,
            py::arg("grammage"), py::arg("e_i"), py::arg("e_f"),
            py::arg("position"), py::arg("random_numbers"),
            R"pbdoc(
                Calculate a random averaged scatterangle `u` alonge the given
                displacement`dr` and the particle directions after distance
                `n_i`.

                Args:
                    dr(double): displacement of particle
                    ei(double): inital energy
                    ef(double): final energy
            )pbdoc");

    auto scattering_doc
        = R"pbdoc(Factory method for a MultipleScattering object.

                Args:
                    name(string): name of the scattering method
                    particle(ParticleDef): particle related constants
                    medium(Medium): medium related constants

                Returns: MultipleScattering
    )pbdoc";

    m.def(
        "make_multiple_scattering",
        [](std::string const& n, ParticleDef const& p, Medium const& m,
            crosssection_list_t<ParticleDef, Medium> c, bool i) {
            return std::shared_ptr<multiple_scattering::Parametrization>(
                make_multiple_scattering(n, p, m, c, i));
        },
        py::arg("name"), py::arg("particle"), py::arg("medium"),
        py::arg("cross"), py::arg("interpolate"), scattering_doc);

    m.def(
        "make_multiple_scattering",
        [](std::string const& n, ParticleDef const& p, Medium const& m) {
            return std::shared_ptr<multiple_scattering::Parametrization>(
                make_multiple_scattering(n, p, m));
        },
        py::arg("name"), py::arg("particle"), py::arg("medium"),
        scattering_doc);

    py::class_<stochastic_deflection::Parametrization,
        std::shared_ptr<stochastic_deflection::Parametrization>>(
        m_sub, "StochasticDeflection")
        .def("n_rnd",
            &stochastic_deflection::Parametrization::RequiredRandomNumbers,
            R"pbdoc(Required random numbers for a stochastic deflection calculation)pbdoc")
        .def("type",
            &stochastic_deflection::Parametrization::GetInteractionType,
            R"pbdoc(Interaction type which causes the stochastic deflection calculation)pbdoc")
        .def("stochastic_deflection",
            &stochastic_deflection::Parametrization::
                CalculateStochasticDeflection,
            py::arg("initial_energy"), py::arg("loss_energy"),
            py::arg("random_numbers"),
            R"pbdoc(TODO: Doc is missing because it's not clear if second argument
            should be the lost energy or the final energy. Please contact the
            maintainers if required.)pbdoc");

    using deflect_ptr = std::shared_ptr<stochastic_deflection::Parametrization>;
    using deflect_list_t = std::vector<deflect_ptr>;

    m.def(
        "make_stochastic_deflection",
        [](std::vector<InteractionType> t, ParticleDef const& p,
            Medium const& m) {
            auto v_shared = deflect_list_t();
            for (auto& v_i : make_stochastic_deflection(t, p, m))
                v_shared.emplace_back(v_i->clone());
            return v_shared;
        },
        py::arg("type"), py::arg("particle"), py::arg("medium"));

    using multiple_scattering_t = multiple_scattering::Parametrization;

    py::class_<Scattering, std::shared_ptr<Scattering>>(m_sub, "Scattering")
        .def(
            py::init([](multiple_scattering_t const& s,
                         deflect_list_t const& d) { return Scattering(s, d); }))
        .def("n_rnd_mulitple_scatter",
            &Scattering::MultipleScatteringRandomNumbers)
        .def("n_rnd_stochastic_deflect",
            &Scattering::StochasticDeflectionRandomNumbers)
        .def("stochastic_deflection",
            &Scattering::CalculateStochasticDeflection<double, double,
                std::vector<double> const&>)
        .def("multiple_scattering",
            &Scattering::CalculateMultipleScattering<double, double, double,
                const Vector3D&, const std::array<double, 4>&>);
}
