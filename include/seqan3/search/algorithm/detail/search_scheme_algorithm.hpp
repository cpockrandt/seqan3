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
 * \brief Provides the algorithm to search in an index using optimum search schemes.
 */

#pragma once

#include <type_traits>

#include <seqan3/search/algorithm/detail/search_scheme_precomputed.hpp>

/*!\addtogroup search
 * \{
 */

namespace seqan3::detail
{

//!\brief Computes a (not optimal) search scheme.
inline auto compute_search_scheme(uint8_t const min_error, uint8_t const max_error)
{
    // TODO: Simple backtracking. Replace this at least by the pigeonhole principle or even better by 01*0 schemes.
    std::vector<search_dyn> scheme{{{1}, {min_error}, {max_error}}};
    return scheme;
}

//!\brief Returns for each search the cumulative length of blocks in the order of blocks in each search.
inline auto search_scheme_block_info(auto const & search_scheme, uint64_t const query_length)
{
    // TODO: use std::array for optimum_search_scheme, std::vector otherwise
    std::vector<std::tuple<std::vector<uint64_t>, uint64_t>> result(search_scheme.size());

    uint8_t  const blocks     {search_scheme[0].blocks()};
    uint64_t const blocklength{query_length / blocks};
    uint8_t  const rest       {query_length - blocks * blocklength};

    std::vector<uint64_t> blocklengths(blocks, blocklength);
    for (uint8_t block_id = 0; block_id < rest; ++block_id)
        ++blocklengths[block_id];

    for (uint8_t search_id = 0; search_id < search_scheme.size(); ++search_id)
    {
        auto const & search = search_scheme[search_id];

        auto & [search_blocklengths, start_pos] = result[search_id];

        start_pos = 0;
        search_blocklengths[0] = blocklengths[search.pi[0] - 1];
        for (uint8_t i = 1; i < blocks; ++i)
        {
            search_blocklengths[i] = blocklengths[search.pi[i] - 1] + search_blocklengths[i - 1];
            if (search.pi[i] < search.pi[0])
                start_pos += search_blocklengths[i] - search_blocklengths[i - 1];
        }
    }

    return result;
}

template <bool abort_on_hit>
inline bool search_scheme_single_search(auto const & index, auto & query, uint64_t const lb, uint64_t const rb,
                                        uint8_t const block_id, auto const & search, auto const & blocklength,
                                        search_params const error_left, auto && delegate)
{
    // TODO
    return false;
}

template <bool abort_on_hit>
inline void search_search_scheme(auto const & index, auto & query, search_params error_left, auto const & search_scheme,
                                 auto && delegate)
{
    auto const block_info = search_scheme_block_info(search_scheme, query.size());

    for (uint8_t search_id = 0; search_id < search_scheme.size(); ++search_id)
    {
        auto const & search = search_scheme[search_id];
        auto const & [blocklength, start_pos] = block_info[search_id];

        bool const hit = search_scheme_single_search<abort_on_hit>(
                            index.begin(), query,     // data
                            start_pos, start_pos + 1, // infix range
                            0,                        // current block id
                            search, blocklength,      // search scheme information
                            error_left,
                            delegate);

        if (abort_on_hit && hit)
            return;
    }
}

}

//!\}
