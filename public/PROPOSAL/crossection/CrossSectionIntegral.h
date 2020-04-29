
/******************************************************************************
 *                                                                            *
 * This file is part of the simulation tool PROPOSAL.                         *
 *                                                                            *
 * Copyright (C) 2017 TU Dortmund University, Department of Physics,          *
 *                    Chair Experimental Physics 5b                           *
 *                                                                            *
 * This software may be modified and distributed under the terms of a         *
 * modified GNU Lesser General Public Licence version 3 (LGPL),               *
 * copied verbatim in the file "LICENSE".                                     *
 *                                                                            *
 * Modifcations to the LGPL License:                                          *
 *                                                                            *
 *      1. The user shall acknowledge the use of PROPOSAL by citing the       *
 *         following reference:                                               *
 *                                                                            *
 *         J.H. Koehne et al.  Comput.Phys.Commun. 184 (2013) 2070-2090 DOI:  *
 *         10.1016/j.cpc.2013.04.001                                          *
 *                                                                            *
 *      2. The user should report any bugs/errors or improvments to the       *
 *         current maintainer of PROPOSAL or open an issue on the             *
 *         GitHub webpage                                                     *
 *                                                                            *
 *         "https://github.com/tudo-astroparticlephysics/PROPOSAL"            *
 *                                                                            *
 ******************************************************************************/

#pragma once

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crossection/CrossSection.h"
#include "PROPOSAL/math/Integral.h"
#include <array>
#include <functional>
#include <utility>
#include <vector>

using std::function;
using std::vector;

namespace PROPOSAL {
class CrossSectionIntegral : public virtual CrossSection {
protected:
    Integral integral_;

    vector<function<double(double, double)>> dndx_integral_;
    vector<function<double(double)>> dedx_integral_;
    vector<function<double(double)>> de2dx_integral_;

public:
    template <typename T>
    CrossSectionIntegral(T&& param, shared_ptr<const EnergyCutSettings> cut);
    virtual ~CrossSectionIntegral() = default;

    enum {TOTAL_RATE, SAMPLED_RATE};
    virtual double dndx_integral(double energy, double rnd);
    virtual double dedx_integral(double energy);
    virtual double de2dx_integral(double energy);

    virtual double CalculatedEdx(double energy);
    virtual double CalculatedE2dx(double energy);
    virtual double CalculatedNdx(double energy, double rnd = 0);
};
} // namespace PROPOSAL



namespace PROPOSAL {
template <typename T>
CrossSectionIntegral::CrossSectionIntegral(
    T&& param, shared_ptr<const EnergyCutSettings> cut)
    : CrossSection(param, cut)
{
    for (const auto& comp : parametrization_->GetComponents()) {
        parametrization_->current_component_(comp);
        dndx_integral_.emplace_back(&CrossSectionIntegral::dndx_integral);
        dedx_integral_.emplace_back(&CrossSectionIntegral::dedx_integral);
        de2dx_integral_.emplace_back(&CrossSectionIntegral::de2dx_integral);
    }
}
} // namespace PROPOSAL
