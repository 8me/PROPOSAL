#pragma once
#include <memory>
#include "PROPOSAL/json_fwd.hpp"

namespace PROPOSAL {
    class CrossSectionBase;
    struct ParticleDef;
    class Medium;

    using cross_ptr = std::unique_ptr<CrossSectionBase>;
}

namespace PROPOSAL {
    cross_ptr make_weakinteraction(const ParticleDef&, const Medium&, bool,
                                   const nlohmann::json&);
}
