#pragma once

namespace ga {

    struct DefaultPolicies {
        using Scalar = float;

        static constexpr Scalar epsilon() noexcept { return static_cast<Scalar>(1e-6); }
    };

    using Policies = DefaultPolicies;

} // namespace ga
