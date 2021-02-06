#pragma once
#include <memory>
#include "PROPOSAL/json_fwd.hpp"

namespace PROPOSAL {
    class CrossSectionBase;
    struct ParticleDef;
    class Medium;
    class EnergyCutSettings;

    using cross_ptr = std::unique_ptr<CrossSectionBase>;
}

namespace PROPOSAL {
    cross_ptr make_epairproduction(const ParticleDef&, const Medium&,
                                  std::shared_ptr<const EnergyCutSettings>,
                                  bool, const nlohmann::json&,
                                  double density_correction = 1.0);
}
