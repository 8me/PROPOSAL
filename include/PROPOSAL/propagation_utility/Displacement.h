#pragma once

#include "PROPOSAL/crosssection/CrossSectionVector.h"
#include "PROPOSAL/math/InterpolantBuilder.h"
#include <vector>

namespace PROPOSAL {
    class CrossSectionBase;
}

namespace PROPOSAL {

class Displacement {
protected:
    using crossbase_list_t = std::vector<std::shared_ptr<CrossSectionBase>>;
    crossbase_list_t cross_list;
    double lower_lim;
    size_t hash;

public:
    Displacement() = default;

    template <typename Cross>
    Displacement(Cross const& cross)
        : cross_list(std::begin(cross), std::end(cross))
        , lower_lim(CrossSectionVector::GetLowerLim(cross))
        , hash(CrossSectionVector::GetHash(cross))
    {
        if (cross.size() < 1)
            throw std::invalid_argument(
                "At least one crosssection is required.");
    }
    virtual ~Displacement() = default;

    static Interpolant1DBuilder::Definition interpol_def;

    double FunctionToIntegral(double);
    virtual double SolveTrackIntegral(double, double) = 0;
    virtual double UpperLimitTrackIntegral(double, double) = 0;

    auto GetHash() const noexcept { return hash; }
    auto GetLowerLim() const noexcept { return lower_lim; }
};

} // namespace PROPOSAL
