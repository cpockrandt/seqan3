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

#include <iostream>

#include <seqan3/core/metafunction/pre.hpp>
#include <seqan3/search/algorithm/detail/search_trivial.hpp>

using namespace seqan3::search_cfg;

//!\cond

namespace seqan3::detail
{

// template <fm_index_iterator_concept iterator_t>
// inline void _filter_hits(std::vector<iterator_t> & hits)
// {
//     std::vector<iterator_t> hits_filtered;
//     hits_filtered.reserve(hits.size());
//
//     auto comp = [](iterator_t const & i1, iterator_t const & i2)
//     {
//         auto const & sa1 = i1.get_sa_range();
//         auto const & sa2 = i2.get_sa_range();
//         return (sa1.i1 < sa2.i1) || (sa1.i1 == sa2.i1 && sa1.i2 < sa2.i2);
//     };
//     std::sort(hits.begin(), hits.end(), comp);
//
//     for (auto it = hits.begin(); it != hits.end(); ++it)
//     {
//
//     }
//
//     return hits_filtered;
// }

// template <seqan3::search_cfg::id id_error, seqan3::search_cfg::id id_error_rate>
// inline uint8_t _compute_errors(auto const & cfg, auto const query_size)
// {
//     // NOTE: Casting doubles rounds towards zero (i.e. floor for positive numbers). Thus given a rate of 10% and a read
//     // length of 101 the maxiumum number of errors is correctly casted from 10.1 errors to 10
//
//     if constexpr (contains<id_error>(cfg))
//     {
//         return get<id_error>(cfg);
//     }
//     else if constexpr (contains<id_error_rate>(cfg))
//     {
//         return static_cast<uint8_t>(get<id_error_rate>(cfg) * query_size);
//     }
//     else
//     {
//         return 0;
//     }
// }

template <typename index_t, typename query_t, typename config_t>
inline auto _search_single(index_t const & index, query_t const & query, config_t const & cfg)
{
    // retrieve error numbers / rates
    detail::search_params max_error{0, 0, 0, 0};
    if constexpr (contains<id::max_error>(cfg))
    {
        auto const t = get<id::max_error>(cfg);
        max_error.total = std::get<0>(t);
        max_error.substitution = std::get<1>(t);
        max_error.insertion = std::get<2>(t);
        max_error.deletion = std::get<3>(t);
    }
    else if constexpr (contains<id::max_error_rate>(cfg))
    {
        auto const t = get<id::max_error_rate>(cfg);
        max_error.total = std::get<0>(t) * query.size();
        max_error.substitution = std::get<1>(t) * query.size();
        max_error.insertion = std::get<2>(t) * query.size();
        max_error.deletion = std::get<3>(t) * query.size();
    }

    // TODO: this would be a lot nicer if all errors would either be ints or doubles and nothing mixed
    // total not set but other error types
    // if constexpr (!contains<id::max_total_error_rate>(cfg) && !contains<id::max_total_error>(cfg) &&
    //               (substitution || insertion || deletion))
    // {
    //     max_error.total = max_error.deletion + max_error.substitution + max_error.insertion;
    // }
    // else if constexpr (contains<id::max_total_error_rate>(cfg) || contains<id::max_total_error>(cfg))
    // {
    //     if (max_error.total == 0 && (max_error.deletion > 0 || max_error.substitution > 0 || max_error.insertion > 0))
    //         throw std::invalid_argument("The total number of errors is set to zero while there is a positive number "
    //                                     "of errors for a specific error type.");
    //     // else if (max_error.total > 0) // not a good idea
    //     // {
    //     //     if constexpr (!deletion)
    //     //         max_error.deletion = max_error.total;
    //     //     if constexpr (!insertion)
    //     //         max_error.insertion = max_error.total;
    //     //     if constexpr (!substitution)
    //     //         max_error.substitution = max_error.total;
    //     // }
    // }

    // TODO: throw exception when any error number or rate is higher than the total error number/rate

    // construct internal delegate for collecting hits for later filtering (if necessary)
    // TODO: pass "it" by reference or value?
    // std::vector<std::tuple<typename index_t::iterator_type, detail::search_params>> internal_hits;
    std::vector<typename index_t::iterator_type> internal_hits;
    auto internal_delegate = [&internal_hits, &max_error](auto const & it/*, detail::search_params const error_left*/)
    {
        internal_hits.push_back(it);
        // internal_hits.push_back({it, max_error - error_left});
    };

    // choose strategy
    if constexpr (contains<id::strategy_best>(cfg))
    {
        detail::search_params max_error2{max_error};
        max_error2.total = 0;
        while (internal_hits.empty() && max_error2.total <= max_error.total)
        {
            detail::search_trivial<true>(index, query, max_error2, internal_delegate);
            max_error2.total++;
        }
    }
    else if constexpr (contains<id::strategy_all_best>(cfg))
    {
        detail::search_params max_error2{max_error};
        max_error2.total = 0;
        while (internal_hits.empty() && max_error2.total <= max_error.total)
        {
            detail::search_trivial<false>(index, query, max_error2, internal_delegate);
            max_error2.total++;
        }
    }
    else if constexpr (contains<id::strategy_strata>(cfg))
    {
        detail::search_params max_error2{max_error};
        max_error2.total = 0;
        while (internal_hits.empty() && max_error2.total <= max_error.total)
        {
            detail::search_trivial<true>(index, query, max_error2, internal_delegate);
            max_error2.total++;
        }
        if (!internal_hits.empty())
        {
            internal_hits.clear(); // don't clear when using Optimum Search Schemes with lower error bounds
            uint8_t const s = seqan3::get<id::strategy_strata>(cfg);
            max_error2.total += s - 1;
            detail::search_trivial<false>(index, query, max_error2, internal_delegate);
        }
    }
    else // "strategy_all" or not specified
    {
        detail::search_trivial<false>(index, query, max_error, internal_delegate);
    }

    // TODO: filter hits and only do it when necessary (depending on error types)
    // _filter_hits(hits);

    // output iterators or text_positions
    if constexpr (contains<id::output_index_iterator>(cfg))
    {
        return internal_hits;
    }
    else
    {
        std::vector<typename index_t::size_type> hits;
        // std::vector<std::tuple<typename index_t::size_type, detail::search_params>> hits;
        // for (auto const & [it, error] : internal_hits)
        if constexpr (contains<id::strategy_best>(cfg))
        {
            // only one iterator is reported but it might contain more than one text position
            if (!internal_hits.empty())
            {
                auto const & text_pos = internal_hits[0].lazy_locate();
                hits.push_back(text_pos[0]);
            }
        }
        else
        {
            for (auto const & it : internal_hits)
                for (auto const & text_pos : it.locate())
                    hits.push_back(text_pos);
                    // hits.push_back({text_pos, error});
        }
        return hits;
    }
}

template <typename index_t, typename queries_t, typename config_t>
    requires
        (std::ranges::RandomAccessRange<queries_t> ||
            (std::ranges::ForwardRange<queries_t> && std::ranges::RandomAccessRange<value_type_t<queries_t>>)) &&
        detail::is_algorithm_configuration_v<remove_cvref_t<config_t>>
inline auto _search(index_t const & index, queries_t const & queries, config_t const & cfg)
{
    // return type: for each query: a vector of text_position (or iterators) and number of errors spent
    // delegate params: text_position (or iterator), number of errors spent and query id. (TODO: or return vector)
    //                  we will withhold all hits of one query anyway to filter duplicates. more efficient to call delegate once with one vector instead of calling delegate for each hit separately at once.
    using hit_t = std::conditional_t<contains<id::output_index_iterator>(cfg),
                                     typename index_t::iterator_type,
                                     typename index_t::size_type>;

    if constexpr (std::ranges::ForwardRange<queries_t> && std::ranges::RandomAccessRange<value_type_t<queries_t>>)
    {
        // TODO: if constexpr (contains<id::on_hit>(cfg))
        // std::vector<std::vector<std::tuple<hit_t, detail::search_params>>> hits;
        std::vector<std::vector<hit_t>> hits;
        hits.reserve(queries.size());
        for (auto query_iter = queries.begin(); query_iter != queries.end(); query_iter++)
        {
            hits.push_back(_search_single(index, *query_iter, cfg));
        }
        return hits;
    }
    else // std::ranges::RandomAccessRange<queries_t>
    {
        // TODO: if constexpr (contains<id::on_hit>(cfg))
        return _search_single(index, queries, cfg);
    }

    // as we allow different error rates for different error types, we cannot really disallow this error configuration
    // if (!substitution && insertion && deletion)
    //     throw std::invalid_argument("Searching with insertions and deletions without substitutions is not allowed. "
    //                                 "Please turn on "); // or static_assert

    // case error_type_enum::insertion | error_type_enum::deletion: // Illegal error_type configuration in the search module. Allowing insertions and deletions without mismatches is illogical and not supported (An insertion followed by a deletion and vice versa corresponds to a mismatch).
    //     if constexpr (std::ranges::ForwardRange<queries_t> &&
    //                   std::ranges::RandomAccessRange<value_type_t<queries_t>>)
    //     {
    //         return std::vector<std::vector<typename index_t::iterator_type>>{};
    //     }
    //     else // std::ranges::RandomAccessRange<queries_t>
    //     {
    //         return std::vector<typename index_t::iterator_type>{};
    //     }
}

} // namespace seqan3::detail

//!\endcond
