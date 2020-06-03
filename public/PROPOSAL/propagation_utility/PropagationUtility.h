
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

#include <vector>

#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include "PROPOSAL/propagation_utility/ContRand.h"
#include "PROPOSAL/propagation_utility/Decay.h"
#include "PROPOSAL/propagation_utility/Displacement.h"
#include "PROPOSAL/propagation_utility/Interaction.h"
#include "PROPOSAL/propagation_utility/Time.h"
#include "PROPOSAL/scattering/Scattering.h"

using std::tuple;
namespace PROPOSAL {

class PropagationUtility {
public:
    struct Collection {

        bool operator==(const Collection& lhs);

        // obligatory pointers
        std::shared_ptr<Interaction> interaction_calc;
        std::shared_ptr<Displacement> displacement_calc;
        std::shared_ptr<Time> time_calc;

        // optional pointers
        std::shared_ptr<Scattering> scattering;
        std::shared_ptr<Decay> decay_calc;
        std::shared_ptr<ContRand> cont_rand;
    };

    PropagationUtility(Collection collection);
    // PropagationUtility(const PropagationUtility&);

    tuple<InteractionType, double> EnergyStochasticloss(double, double);
    double EnergyDecay(double, std::function<double()>);
    double EnergyInteraction(double, std::function<double()>);
    double EnergyRandomize(double, double, std::function<double()>);
    double EnergyDistance(double, double);
    double LengthContinuous(double, double);
    double TimeElapsed(double, double, double);

    // TODO: return value doesn't tell what it include. Maybe it would be better
    // to give a tuple of two directions back. One is the mean over the
    // displacement and the other is the actual direction. With a get method
    // there could be a possible access with the position of the object stored
    // in an enum.

    tuple<Vector3D, Vector3D> DirectionsScatter(
        double, double, double, const Vector3D&, const std::array<double, 4>&);
    // std::pair<double, double> DirectionDeflect(CrossSection&, double, double);

    Collection collection;
};
} // namespace PROPOSAL
