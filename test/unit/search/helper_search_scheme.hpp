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
 * \brief Provides helper functions for testing search schemes.
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 */

#pragma once

#include <algorithm>
#include <vector>

#include <seqan3/io/stream/debug_stream.hpp>
#include <seqan3/search/algorithm/detail/search_scheme_precomputed.hpp>

namespace seqan3
{

//!\brief Order the elements in a vector in the same order as `search.pi`.
template <typename T>
inline void order_search_vector(std::vector<T> & v, auto const & search)
{
    std::vector<T> v_tmp{v};
    for (uint8_t i = 0; i < search.blocks(); ++i)
    {
        uint8_t const index = search.pi[i] - 1;
        v[index] = v_tmp[i];
    }
}

//!\brief Reorder the search and blocks_length s.t. all fields are in the order from left to right.
inline void get_ordered_search(auto const & search, auto const & blocks_length,
                               auto & ordered_search, auto & ordered_blocks_length)
{
    ordered_blocks_length.resize(search.blocks());
    for (uint8_t i = 0; i < search.blocks(); ++i)
    {
        uint8_t const index = search.pi[i] - 1;
        ordered_search.pi[index] = search.pi[i]; // i.e. index + 1
        ordered_search.l[index] = search.l[i];
        ordered_search.u[index] = search.u[i];
        ordered_blocks_length[index] = blocks_length[i] - ((i > 0) ? blocks_length[i - 1] : 0);
    }
}

//!\brief Helper function for search_error_distribution(res, search).
inline void search_error_distribution(std::vector<std::vector<uint8_t> > & res, auto l, auto u, uint8_t const e)
{
    if (l.size() == 0)
    {
        res.emplace_back(std::forward<std::vector<uint8_t>>(std::vector<uint8_t>{}));
        return;
    }

    uint8_t l0 = l[0];
    uint8_t u0 = u[0];
    l.erase(l.begin());
    u.erase(u.begin());

    for (uint8_t i = std::max(e, l0); i <= u0; ++i)
    {
        std::vector<std::vector<uint8_t> > res_recursive;
        search_error_distribution(res_recursive, l, u, i);
        for (auto & res_elem : res_recursive)
        {
            res_elem.insert(res_elem.begin(), i - e);
            res.push_back(res_elem);
        }
    }
}

//!\brief Compute all possible error distributions given a search.
//        The result is in the same order as the search (i.e. search.pi).
inline void search_error_distribution(std::vector<std::vector<uint8_t> > & res, auto const & search)
{
    res.clear();
    std::vector<uint8_t> l(search.l.begin(), search.l.end());
    std::vector<uint8_t> u(search.u.begin(), search.u.end());
    search_error_distribution(res, l, u, 0u);
}

//!\brief Compute all possible error distributions for each search given a search scheme.
inline void search_scheme_error_distribution(std::vector<std::vector<uint8_t> > & res, auto const & search_scheme)
{
    res.clear();
    for (auto && search : search_scheme)
    {
        std::vector<std::vector<uint8_t> > res_search;
        search_error_distribution(res_search, search);
        for (auto & single_error_distribution : res_search)
            order_search_vector(single_error_distribution, search);
        res.insert(res.end(), res_search.begin(), res_search.end());
    }
}

//!\brief Construct a trivial search scheme with one search and a specified number of blocks.
inline detail::search_scheme_dyn_type trivial_search_scheme(uint8_t const min_error, uint8_t const max_error,
                                                            uint8_t const blocks)
{
    detail::search_scheme_dyn_type search_scheme(1);
    auto & search = search_scheme.front();
    search.pi.resize(blocks);
    search.l.resize(blocks);
    search.u.resize(blocks);
    for (uint8_t i = 0; i < blocks; ++i)
    {
        search.pi[i] = i + 1;
        search.l[i] = 0;
        search.u[i] = max_error;
    }
    search.l[blocks - 1] = min_error;
    return search_scheme;
}

} // namespace std
