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
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 * \brief
 */

#pragma once

#include <seqan3/range/view/persist.hpp>
#include <seqan3/search/algorithm/detail/search.hpp>
#include <seqan3/search/fm_index/all.hpp>

using namespace seqan3::search_cfg;

namespace seqan3
{

template <class a>
struct nothing;

//!\brief \todo Document!
template <typename index_t, typename queries_t, typename config_t>
//!\cond
    requires
        (std::ranges::RandomAccessRange<queries_t> ||
            (std::ranges::ForwardRange<queries_t> && std::ranges::RandomAccessRange<value_type_t<queries_t>>)) &&
        detail::is_algorithm_configuration_v<remove_cvref_t<config_t>>
//!\endcond
inline auto search(index_t const & index, queries_t const & queries, config_t const & cfg)
{
    // set default error types when total is set but no specific error type
    // bool
    //     substitution = contains<id::max_substitution_error_rate>(cfg) || contains<id::max_substitution_error>(cfg),
    //     insertion    = contains<id::max_insertion_error_rate>(cfg) || contains<id::max_insertion_error>(cfg),
    //     deletion     = contains<id::max_deletion_error_rate>(cfg) || contains<id::max_deletion_error>(cfg);

    // TODO!
    // if constexpr (!substitution && !insertion && !deletion)
    // {
    //     if constexpr (contains<id::max_total_error>(cfg))
    //     {
    //         uint8_t const total_error = get<id::max_total_error>(cfg);
    //         detail::configuration const cfg2 = cfg | max_substitution_error(total_error)
    //                                                | max_insertion_error(total_error)
    //                                                | max_deletion_error(total_error);
    //         return detail::search<true, true, true>(index, queries, cfg2);
    //     }
    //     else if constexpr (contains<id::max_total_error_rate>(cfg))
    //     {
    //         double const total_error_rate = get<id::max_total_error_rate>(cfg);
    //         detail::configuration const cfg2 = cfg | max_substitution_error_rate(total_error_rate)
    //                                                | max_insertion_error_rate(total_error_rate)
    //                                                | max_deletion_error_rate(total_error_rate);
    //         return detail::search<true, true, true>(index, queries, cfg2);
    //     }
    // }
    if constexpr (contains<id::mode>(cfg))
    {
        return detail::_search(index, queries, cfg);
    }
    else
    {
        // TODO: overload pipe operator for empty config object
        if constexpr (std::Same<remove_cvref_t<decltype(cfg)>, detail::configuration<>>)
        {
            detail::configuration const cfg2 = mode(all);
            return detail::_search(index, queries, cfg2);
        }
        else
        {
            detail::configuration const cfg2 = cfg | mode(all);
            return detail::_search(index, queries, cfg2);
        }
    }
}

// TODO: const &, &&, etc.: use forwarding references and add a view::persist
// DOC: insertion/deletion are with resp. to the query. i.e. an insertion is the insertion of a base into the query
// that does not occur in the text at the position
//!\brief \todo Document!
template <typename index_t, typename queries_t>
//!\cond
    requires std::ranges::RandomAccessRange<queries_t> ||
             (std::ranges::ForwardRange<queries_t> && std::ranges::RandomAccessRange<value_type_t<queries_t>>)
//!\endcond
inline auto search(index_t const & index, queries_t const & queries)
{
    // TODO: auto const queries_lvalue = queries | view::persist;

    // set default configuration
    detail::configuration const cfg; // max_total_error(0) | return_text_position() | strategy_all();
    return search(index, queries/*_lvalue*/, cfg);
}

} // namespace seqan3
