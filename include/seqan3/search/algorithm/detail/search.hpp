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

#include <seqan3/core/metafunction/pre.hpp>
#include <seqan3/search/algorithm/detail/search_trivial.hpp>
#include <seqan3/search/algorithm/detail/search_scheme_algorithm.hpp>
#include <seqan3/search/fm_index/concept.hpp>

namespace seqan3::detail
{

//!\brief \todo Docs missing
template <typename index_t, typename query_t, typename config_t>
inline auto _search_single(index_t const & index, query_t const & query, config_t const & cfg)
{
    // retrieve error numbers / rates
    detail::search_params max_error{0, 0, 0, 0};
    auto & [total, subs, ins, del] = max_error;
    if constexpr (contains<search_cfg::id::max_error>(cfg))
    {
        std::tie(total, subs, ins, del) = get<search_cfg::id::max_error>(cfg);
    }
    else if constexpr (contains<search_cfg::id::max_error_rate>(cfg))
    {
        // NOTE: Casting doubles rounds towards zero (i.e. floor for positive numbers). Thus given a rate of 10% and a
        // read length of 101 the maxiumum number of errors is correctly casted from 10.1 errors to 10
        std::tie(total, subs, ins, del) = std::apply([& query](auto && ... args)
            {
                return std::tuple{(args * query.size())...};
            }, get<search_cfg::id::max_error_rate>(cfg));
    }

    // TODO: if total not set: max_error.total = max_error.deletion + max_error.substitution + max_error.insertion;
    // TODO: throw exception when any error number or rate is higher than the total error number/rate
    //         throw std::invalid_argument("The total number of errors is set to zero while there is a positive number "
    //                                     "of errors for a specific error type.");

    // construct internal delegate for collecting hits for later filtering (if necessary)
    // TODO: pass "it" by reference or value?
    std::vector<typename index_t::iterator_type> internal_hits;
    auto internal_delegate = [&internal_hits, &max_error](auto const & it)
    {
        internal_hits.push_back(it); // TODO: report max_error-error_left
    };

    // choose mode
    auto const & selected_mode = seqan3::get<search_cfg::id::mode>(cfg);
    if constexpr (std::Same<remove_cvref_t<decltype(selected_mode)>, detail::search_mode_best>)
    {
        detail::search_params max_error2{max_error};
        max_error2.total = 0;
        while (internal_hits.empty() && max_error2.total <= max_error.total)
        {
            detail::search_trivial<true>(index, query, max_error2, internal_delegate);
            max_error2.total++;
        }
    }
    else if constexpr (std::Same<remove_cvref_t<decltype(selected_mode)>, detail::search_mode_all_best>)
    {
        detail::search_params max_error2{max_error};
        max_error2.total = 0;
        while (internal_hits.empty() && max_error2.total <= max_error.total)
        {
            detail::search_trivial<false>(index, query, max_error2, internal_delegate);
            max_error2.total++;
        }
    }
    else if constexpr (std::Same<remove_cvref_t<decltype(selected_mode)>, search_cfg::strata>)
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
            uint8_t const s = selected_mode;
            max_error2.total += s - 1;
            detail::search_trivial<false>(index, query, max_error2, internal_delegate);
        }
    }
    else // detail::search_mode_all
    {
        detail::search_trivial<false>(index, query, max_error, internal_delegate);
    }

    // TODO: filter hits and only do it when necessary (depending on error types)

    // output iterators or text_positions
    auto const & output = seqan3::get<search_cfg::id::output>(cfg);
    if constexpr (std::Same<remove_cvref_t<decltype(output)>, detail::search_output_index_iterator>)
    {
        return internal_hits;
    }
    else
    {
        std::vector<typename index_t::size_type> hits;
        auto const & selected_mode = seqan3::get<search_cfg::id::mode>(cfg);
        if constexpr (std::Same<remove_cvref_t<decltype(selected_mode)>, detail::search_mode_best>)
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
            {
                for (auto const & text_pos : it.locate())
                    hits.push_back(text_pos);
                std::sort(hits.begin(), hits.end());
                hits.erase(std::unique(hits.begin(), hits.end()), hits.end());
            }
        }
        return hits;
    }
}

//!\brief \todo Docs missing
template <typename index_t, typename queries_t, typename config_t>
    requires
        (std::ranges::RandomAccessRange<queries_t> ||
            (std::ranges::ForwardRange<queries_t> && std::ranges::RandomAccessRange<value_type_t<queries_t>>)) &&
        detail::is_algorithm_configuration_v<remove_cvref_t<config_t>>
inline auto _search(index_t const & index, queries_t const & queries, config_t const & cfg)
{
    // return type: for each query: a vector of text_position (or iterators) and number of errors spent
    // delegate params: text_position (or iterator), number of errors spent and query id. (TODO: or return vector)
    //                  we will withhold all hits of one query anyway to filter duplicates. more efficient to call
    //                  delegate once with one vector instead of calling delegate for each hit separately at once.
    auto const & output = seqan3::get<search_cfg::id::output>(cfg);
    using hit_t = std::conditional_t<std::Same<remove_cvref_t<decltype(output)>, detail::search_output_index_iterator>,
                                     typename index_t::iterator_type,
                                     typename index_t::size_type>;

    if constexpr (std::ranges::ForwardRange<queries_t> && std::ranges::RandomAccessRange<value_type_t<queries_t>>)
    {
        // TODO: if constexpr (contains<search_cfg::id::on_hit>(cfg))
        std::vector<std::vector<hit_t>> hits;
        hits.reserve(queries.size());
        for (auto const query : queries)
        {
            hits.push_back(_search_single(index, query, cfg));
        }
        return hits;
    }
    else // std::ranges::RandomAccessRange<queries_t>
    {
        // TODO: if constexpr (contains<search_cfg::id::on_hit>(cfg))
        return _search_single(index, queries, cfg);
    }
}

} // namespace seqan3::detail
