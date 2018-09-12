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
 * \brief Provides the configuration for maximum number of insertion errors.
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 */

#pragma once

#include <seqan3/search/configuration/utility.hpp>
#include <seqan3/core/algorithm/all.hpp>
#include <seqan3/core/metafunction/basic.hpp>
#include <seqan3/core/metafunction/template_inspection.hpp>

namespace seqan3::detail
{
/*!\brief A configuration element for the maximum number of insertion errors.
 * \ingroup search_configuration
 */
struct search_config_max_insertion_error
{
    //!\brief The actual value.
    uint8_t value;
};

/*!\brief The max_insertion_error adaptor enabling pipe notation.
 * \ingroup search_configuration
 */
struct search_config_max_insertion_error_adaptor : public configuration_fn_base<search_config_max_insertion_error_adaptor>
{

    /*!\brief Adds to the configuration a max_insertion_error configuration element.
     * \relates seqan3::search_config_max_insertion_error
     * \param[in] cfg  The configuration to be extended.
     * \param[in] nbr The number of maximum insertion errors used for the algorithm.
     * \returns A new configuration containing the max_insertion_error configuration element.
     */
    template <typename configuration_t>
    //!\cond
        requires is_algorithm_configuration_v<remove_cvref_t<configuration_t>>
    //!\endcond
    constexpr auto invoke(configuration_t && cfg, uint8_t const nbr) const
    {
        static_assert(is_valid_search_configuration_v<search_cfg::id::max_insertion_error, remove_cvref_t<configuration_t>>,
                      SEQAN3_INVALID_CONFIG(search_cfg::id::max_insertion_error));

        search_config_max_insertion_error tmp{nbr};
        return std::forward<configuration_t>(cfg).push_front(std::move(tmp));
    }
};

//!\brief Helper template meta-function associated with detail::search_config_max_insertion_error.
//!\ingroup search_configuration
template <>
struct on_search_config<search_cfg::id::max_insertion_error>
{
    //!\brief Type alias used by meta::find_if
    template <config_element_concept t>
    using invoke = typename std::is_same<t, search_config_max_insertion_error>::type;
};

//!\brief Mapping from the detail::search_config_max_insertion_error type to it's corresponding seqan3::search_cfg::id.
//!\ingroup search_configuration
template <>
struct search_config_type_to_id<search_config_max_insertion_error>
{
    //!\brief The associated seqan3::search_cfg::id.
    static constexpr search_cfg::id value = search_cfg::id::max_insertion_error;
};
} // namespace seqan3::detail

namespace seqan3::search_cfg
{
/*!\brief A configuration element for the maximum number of insertion errors.
 * \ingroup search_configuration
 */
inline constexpr detail::search_config_max_insertion_error_adaptor max_insertion_error;

} // namespace seqan3::search_cfg
