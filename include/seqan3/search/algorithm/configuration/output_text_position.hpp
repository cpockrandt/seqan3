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
 * \brief Provides the configuration for returning positions in the text.
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
/*!\brief Configuration element to specify the return type of the hits as positions in the indexed text.
 * \ingroup search_configuration
 */
struct search_config_output_text_position
{
    //!\cond
    static constexpr bool value{true};
    //!\endcond
};

/*!\brief The seqan3::search_cfg::output_text_position adaptor enabling pipe notation.
 * \ingroup search_configuration
 */
struct search_config_output_text_position_adaptor : public configuration_fn_base<search_config_output_text_position_adaptor>
{

    /*!\brief Adds to the configuration the seqan3::search_cfg::output_text_position configuration element.
     * \param[in] cfg The configuration to be extended.
     * \returns A new configuration containing the seqan3::search_cfg::output_text_position configuration element.
     */
    template <typename configuration_t>
    //!\cond
        requires is_algorithm_configuration_v<remove_cvref_t<configuration_t>>
    //!\endcond
    constexpr auto invoke(configuration_t && cfg) const
    {
        static_assert(is_valid_search_configuration_v<search_cfg::id::output_text_position, remove_cvref_t<configuration_t>>,
                      SEQAN3_INVALID_CONFIG(search_cfg::id::output_text_position));

        return std::forward<configuration_t>(cfg).push_front(search_config_output_text_position{});
    }
};

//!\brief Helper template meta-function associated with detail::search_config_output_text_position.
//!\ingroup search_configuration
template <>
struct on_search_config<search_cfg::id::output_text_position>
{
    //!\brief Type alias used by meta::find_if
    template <config_element_concept t>
    using invoke = typename std::is_same<t, search_config_output_text_position>::type;
};

//!\brief Mapping from the detail::search_config_output_text_position type to it's corresponding seqan3::search_cfg::id.
//!\ingroup search_configuration
template <>
struct search_config_type_to_id<search_config_output_text_position>
{
    //!\brief The associated seqan3::search_cfg::id.
    static constexpr search_cfg::id value = search_cfg::id::output_text_position;
};
} // namespace seqan3::detail

namespace seqan3::search_cfg
{
/*!\brief Configuration element to specify the return type of the hits as positions in the indexed text.
 * \ingroup search_configuration
 */
inline constexpr detail::search_config_output_text_position_adaptor output_text_position;

} // namespace seqan3::search_cfg

//!\}