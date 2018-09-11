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

#include <seqan3/search/detail/search_trivial.hpp>

using namespace seqan3::search_cfg;

namespace seqan3
{

inline void _filter_hits()
{

}

template <bool substitution, bool insertion, bool deletion, typename index_t, typename query_t, typename config_t>
inline auto _search_single(index_t const & index, query_t const & query, config_t const & cfg)
{
    constexpr bool to_text_position = false;

    // Retrieve max_error:
    // If error_type is used, the error number must be positive and vice versa.
    uint8_t max_error = 0;
    if constexpr (seqan3::contains<id::max_error>(cfg))
    {
        max_error = seqan3::get<id::max_error>(cfg);
        assert((substitution || insertion || deletion) == (max_error > 0));
    }
    else if constexpr (seqan3::contains<id::max_error_rate>(cfg))
    {
        auto const max_error_rate = seqan3::get<id::max_error_rate>(cfg);
        max_error = static_cast<uint8_t>(max_error_rate * query.size()); // cast rounds towards zero (i.e. floor for positive numbers)
        // Do not compare max_error_rate since the range of queries could contain some queries that are too short
        // for the given error rate to allow for any errors. This should not trigger any errors.
        assert((substitution || insertion || deletion) == (max_error_rate > 0));
    }
    else
    {
        assert(!(substitution || insertion || deletion));
    }

    constexpr bool const search_strategy_best = seqan3::contains<id::strategy_best>(cfg);
    constexpr bool const search_strategy_strata = seqan3::contains<id::strategy_strata>(cfg);
    constexpr bool const search_strategy_all_best = seqan3::contains<id::strategy_all_best>(cfg);

    using hit_t = std::conditional_t<to_text_position, uint64_t, typename index_t::iterator_type>;

    std::vector<hit_t> hits;
    auto internal_delegate = [&hits](auto const & it, uint8_t const errors) // TODO: pass it by reference or value?
    {
        if constexpr (to_text_position)
        {
            for (auto const & text_pos : it.locate())
                hits.push_back(text_pos);
        }
        else
        {
            hits.push_back(it);
        }
    };

    // choose strategy
    if constexpr (search_strategy_best)
    {
        uint8_t errors = 0;
        while (hits.empty() && errors <= max_error)
        {
            detail::search_trivial<substitution, insertion, deletion, true>(index, query, errors, internal_delegate);
            ++errors;
        }
    }
    else if constexpr (search_strategy_all_best)
    {
        uint8_t errors = 0;
        while (hits.empty() && errors <= max_error)
        {
            detail::search_trivial<substitution, insertion, deletion, false>(index, query, errors, internal_delegate);
            ++errors;
        }
    }
    else if constexpr (search_strategy_strata)
    {
        uint8_t errors = 0;
        while (hits.empty() && errors <= max_error)
        {
            detail::search_trivial<substitution, insertion, deletion, true>(index, query, errors, internal_delegate);
            ++errors;
        }
        if (!hits.empty())
        {
            hits.clear();
            uint8_t const s = seqan3::get<id::strategy_strata>(cfg);
            detail::search_trivial<substitution, insertion, deletion, false>(index, query, errors - 1 + s, internal_delegate);
        }
    }
    else // "strategy_all" or not specified
    {
        detail::search_trivial<substitution, insertion, deletion, false>(index, query, max_error, internal_delegate);
    }

    // TODO: filtering: all_best, all (and strata): filtering when non-disjoint error_type

    // TODO: output iterators or text_positions

    return hits;
}

template <bool substitution, bool insertion, bool deletion, typename index_t, typename queries_t, typename config_t>
inline auto _search(index_t const & index, queries_t const & queries, config_t const & cfg)
{
    // return type: for each query: a vector of text_position (or iterators) and number of errors spent
    // delegate params: text_position (or iterator), number of errors spent and query id. (TODO: or return vector)
    //                  we will withhold all hits of one query anyway to filter duplicates. more efficient to call delegate once with one vector instead of calling delegate for each hit separately at once.
    constexpr bool to_text_position = false;
    using hit_t = std::conditional_t<to_text_position, uint64_t, typename index_t::iterator_type>;

    if constexpr (std::ranges::ForwardRange<queries_t> && std::ranges::RandomAccessRange<typename queries_t::value_type>)
    {
        // TODO: call delegate instead of return if on_hit cfg is being passed
        std::vector<std::vector<hit_t>> hits;
        hits.reserve(queries.size());
        for (auto query_iter = queries.begin(); query_iter != queries.end(); query_iter++)
        {
            hits.push_back(_search_single<substitution, insertion, deletion>(index, *query_iter, cfg));
        }
        return hits;
    }
    else // std::ranges::RandomAccessRange<queries_t>
    {
        // TODO: call delegate instead of return if on_hit cfg is being passed
        return _search_single<substitution, insertion, deletion>(index, queries, cfg);
    }
}

template <typename index_t, typename queries_t, typename config_t>
//!\cond
    // TODO: ForwardRange does not guarantee member type "value_type"
    requires
        (std::ranges::RandomAccessRange<queries_t> ||
            (std::ranges::ForwardRange<queries_t> && std::ranges::RandomAccessRange<typename queries_t::value_type>)) //&&
        // TODO: doesn't work detail::is_algorithm_configuration_v<remove_cvref_t<config_t>>
//!\endcond
inline auto search(index_t const & index, queries_t const & queries, config_t const & cfg)
{
    error_type_enum error_type = error_type_enum::none;
    if constexpr (seqan3::contains<id::error_type>(cfg))
    {
        error_type = seqan3::get<id::error_type>(cfg);
    }

    // <substitution, insertion, deletion>
    switch (error_type)
    {
        case error_type_enum::none:
            return _search<false, false, false>(index, queries, cfg);
        case error_type_enum::deletion:
            return _search<false, false, true>(index, queries, cfg);
        case error_type_enum::insertion:
            return _search<false, true, false>(index, queries, cfg);
        case error_type_enum::substitution:
            return _search<true, false, false>(index, queries, cfg);
        case error_type_enum::substitution | error_type_enum::deletion:
            return _search<true, false, true>(index, queries, cfg);
        case error_type_enum::substitution | error_type_enum::insertion:
            return _search<true, true, false>(index, queries, cfg);
        case error_type_enum::substitution | error_type_enum::insertion | error_type_enum::deletion:
            return _search<true, true, true>(index, queries, cfg);
        case error_type_enum::insertion | error_type_enum::deletion: // Illegal error_type configuration in the search module. Allowing insertions and deletions without mismatches is illogical and not supported (An insertion followed by a deletion and vice versa corresponds to a mismatch).
        default: // Illegal error_type configuration in the search module.
            assert(false);
            if constexpr (std::ranges::ForwardRange<queries_t> &&
                          std::ranges::RandomAccessRange<typename queries_t::value_type>)
            {
                return std::vector<std::vector<typename index_t::iterator_type>>{};
            }
            else // std::ranges::RandomAccessRange<queries_t>
            {
                return std::vector<typename index_t::iterator_type>{};
            }
    }
}

template <typename index_t, typename queries_t>
//!\cond
    // TODO: ForwardRange does not guarantee member type "value_type"
    requires std::ranges::RandomAccessRange<queries_t> || (std::ranges::ForwardRange<queries_t> && std::ranges::RandomAccessRange<typename queries_t::value_type>)
//!\endcond
inline auto search(index_t const & index, queries_t const & queries)
{
    detail::configuration cfg = max_error(0) | strategy_all();
    return search(index, queries, cfg);
}

} // namespace seqan3
