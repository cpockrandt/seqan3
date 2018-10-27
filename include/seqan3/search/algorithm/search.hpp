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
 * \brief Provides the public interface for search algorithms.
 */

#pragma once

#include <seqan3/range/view/persist.hpp>
#include <seqan3/search/algorithm/detail/search.hpp>
#include <seqan3/search/fm_index/all.hpp>

namespace seqan3
{

/*!\brief Search a query or a range of queries in an index.
 * \param[in] index String index to be searched.
 * \param[in] queries A single query or a (forward) range of queries. A query must model
                      std::ranges::RandomAccessRange.
 * \param[in] cfg A configuration object specifying the search parameters (e.g. number of errors, error types,
 *                output format, etc.).
 *
 * ### Complexity
 *
 * Exponential in the numbers of errors.
 *
 * ### Exceptions
 *
 * No-throw guarantee.
 */
template <typename index_t, typename queries_t, typename config_t>
//!\cond
    requires
        (std::ranges::RandomAccessRange<queries_t> ||
            (std::ranges::ForwardRange<queries_t> &&
             std::ranges::RandomAccessRange<value_type_t<queries_t>>)) &&
        detail::is_algorithm_configuration_v<remove_cvref_t<config_t>>
//!\endcond
inline auto search(index_t const & index, queries_t && queries, config_t const & cfg)
{
    if constexpr (contains<search_cfg::id::max_error>(cfg))
    {
        auto & [total, subs, ins, del] = get<search_cfg::id::max_error>(cfg);
        if (subs > total)
            throw std::invalid_argument("Substitution error threshold is higher than the total error "
                                        "threshold.");
        if (ins > total)
            throw std::invalid_argument("Insertion error threshold is higher than the total error threshold.");
        if (del > total)
            throw std::invalid_argument("Deletion error threshold is higher than the total error threshold.");
    }
    else if constexpr (contains<search_cfg::id::max_error_rate>(cfg))
    {
        auto & [total, subs, ins, del] = get<search_cfg::id::max_error_rate>(cfg);
        if (subs > total)
            throw std::invalid_argument("Substitution error threshold is higher than the total error "
                                        "threshold.");
        if (ins > total)
            throw std::invalid_argument("Insertion error threshold is higher than the total error threshold.");
        if (del > total)
            throw std::invalid_argument("Deletion error threshold is higher than the total error threshold.");
    }

    if constexpr (contains<search_cfg::id::mode>(cfg))
    {
        if constexpr (contains<search_cfg::id::output>(cfg))
            return detail::search_all(index, queries, cfg);
        else
            return detail::search_all(index, queries, cfg | search_cfg::output(search_cfg::text_position));
    }
    else
    {
        detail::configuration const cfg2 = cfg | search_cfg::mode(search_cfg::all);
        if constexpr (contains<search_cfg::id::output>(cfg))
            return detail::search_all(index, queries, cfg2);
        else
            return detail::search_all(index, queries, cfg2 | search_cfg::output(search_cfg::text_position));
    }
}

/*!\brief Search a query or a range of queries in an index.
 * \param[in] index String index to be searched.
 * \param[in] queries A single query or a (forward) range of queries. A query must model
                      std::ranges::RandomAccessRange.
 * \param[in] cfg A configuration object specifying the search parameters (e.g. number of errors, error types,
 *                output format, etc.).
 *
 * ### Complexity
 *
 * Exponential in the numbers of errors.
 *
 * ### Exceptions
 *
 * No-throw guarantee.
 */
template <typename index_t, typename queries_t>
//!\cond
    requires std::ranges::RandomAccessRange<queries_t> ||
             (std::ranges::ForwardRange<queries_t> && std::ranges::RandomAccessRange<value_type_t<queries_t>>)
//!\endcond
inline auto search(index_t const & index, queries_t && queries)
{
    detail::configuration const default_cfg = search_cfg::max_error(search_cfg::total{0},
                                                                    search_cfg::substitution{0},
                                                                    search_cfg::insertion{0},
                                                                    search_cfg::deletion{0})
                                            | search_cfg::output(search_cfg::text_position)
                                            | search_cfg::mode(search_cfg::all);
    return search(index, queries, default_cfg);
}

} // namespace seqan3
