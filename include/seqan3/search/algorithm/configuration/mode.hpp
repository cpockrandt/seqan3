// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ============================================================================

/*!\file
 * \brief Provides the search strategy configuration "all".
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/algorithm/all.hpp>
#include <seqan3/core/metafunction/basic.hpp>
#include <seqan3/core/metafunction/template_inspection.hpp>
#include <seqan3/search/algorithm/configuration/utility.hpp>

/*!\addtogroup search
 * \{
 */

namespace seqan3::detail
{

struct search_mode_all {};

struct search_mode_all_best {};

struct search_mode_best {};

} // namespace seqan3::detail

namespace seqan3::search_cfg
{

inline detail::search_mode_all constexpr all;
inline detail::search_mode_all_best constexpr all_best;
inline detail::search_mode_best constexpr best;

/*!\brief A strong type of underlying type `uint8_t` or `double` that represents the number or rate of total errors.
 * \tparam value_t The underlying type
 * \ingroup search_configuration
 */
struct strata : detail::strong_type<uint8_t, strata, detail::strong_type_skill::convert>
{
    using detail::strong_type<uint8_t, strata, detail::strong_type_skill::convert>::strong_type;
};

} // namespace seqan3::search_cfg

namespace seqan3::detail
{
/*!\brief Configuration element to receive all possible hits.
 * \ingroup search_configuration
 */
template <typename mode_t>
struct search_config_mode
{
    //!\cond
    mode_t value;
    //!\endcond
};

/*!\brief The seqan3::search_cfg::mode adaptor enabling pipe notation.
 * \ingroup search_configuration
 */
template <template <typename ...> typename search_config_mode_type>
struct search_config_mode_adaptor : public configuration_fn_base<search_config_mode_adaptor<search_config_mode>>
{

    /*!\brief Adds to the configuration the seqan3::search_cfg::mode configuration element.
     * \param[in] cfg The configuration to be extended.
     * \returns A new configuration containing the seqan3::search_cfg::mode configuration element.
     */
    template <typename configuration_t, typename mode_t>
    //!\cond
        requires is_algorithm_configuration_v<remove_cvref_t<configuration_t>>
    //!\endcond
    // TODO: require strong type
    constexpr auto invoke(configuration_t && cfg, mode_t mode) const
    {
        static_assert(is_valid_search_configuration_v<search_cfg::id::mode, remove_cvref_t<configuration_t>>,
                      SEQAN3_INVALID_CONFIG(search_cfg::id::mode));

        return std::forward<configuration_t>(cfg).push_front(search_config_mode<mode_t>{std::move(mode)});
    }
};

//!\brief Helper template meta-function associated with detail::search_config_mode.
//!\ingroup search_configuration
template <>
struct on_search_config<search_cfg::id::mode>
{
    //!\brief Type alias used by meta::find_if
    template <config_element_concept t>
    using invoke = typename is_type_specialisation_of<t, search_config_mode>::type;
};

//!\brief Mapping from the detail::search_config_mode type to it's corresponding seqan3::search_cfg::id.
//!\ingroup search_configuration
template <typename mode_t>
struct search_config_type_to_id<search_config_mode<mode_t>>
{
    //!\brief The associated seqan3::search_cfg::id.
    static constexpr search_cfg::id value = search_cfg::id::mode;
};
} // namespace seqan3::detail

namespace seqan3::search_cfg
{
/*!\brief Configuration element to receive all possible hits.
 * \ingroup search_configuration
 */
inline detail::search_config_mode_adaptor<seqan3::detail::search_config_mode> constexpr mode;

} // namespace seqan3::search_cfg

//!\}
