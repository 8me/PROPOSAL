#pragma once

#include "PROPOSAL/scattering/Scattering.h"

namespace PROPOSAL {
class ScatteringMultiplier : public Scattering {
    double multiple_scatt = 1;
    std::vector<std::pair<InteractionType, double>> stochastic_deflect;

    DirectionChangeAngular _scale_deflect(DirectionChangeAngular& angles, InteractionType t) override
    {
        for (auto m : stochastic_deflect) {
            if (m.first == t) {
                angles.zenith *= m.second;
                assert(angles.zenith <= PI);
                return angles;
            }
        }
        return angles;
    }

    std::array<double, 4> _scale_scatter(std::array<double, 4>& angles) override
    {
        for (auto& a : angles)
            a *= multiple_scatt;
        return angles;
    }

public:
    /**
     * @brief A wrapper class for handling with scattering multiplier.
     * Multiplier are linear factors for scattering angles.
     *
     * @tparam T1 multiple_scattering::Parametrization or nullptr_t
     * @tparam T2 container of stochastic_deflection::Parametrization or
     * @param _m Multiple scattering calculator to take deflections caused by
     * continuous losses into account
     * @param _s list of deflection calculator to take stochastic deflections
     * @param _mm multiple scattering factor
     * @param _sm interaction type dependent stochastic deflection factor
     */
    template <typename T1, typename T2>
    ScatteringMultiplier(T1&& _m, T2&& _s, double _mm,
        std::vector<std::pair<InteractionType, double>> _sm)
        : Scattering(std::forward<T1>(_m), std::forward<T2>(_s))
        , multiple_scatt(_mm)
        , stochastic_deflect(_sm)
    {
    }
};
} // namespace PROPOSAL
