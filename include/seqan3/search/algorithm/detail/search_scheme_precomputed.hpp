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
 * \brief Provides the data structures and precomputed instances for (optimum) search schemes.
 */

#pragma once

/*!\addtogroup search
 * \{
 */

namespace seqan3::detail
{

//!\brief Object storing information for a search (of a search scheme).
//        Number of blocks have to be known at compile time.
template <uint8_t nbr_blocks>
struct search
{
    //!\brief Type for storing the length of blocks
    typedef std::array<uint64_t, nbr_blocks> blocks_length_type;

    //!\brief Order of blocks
    std::array<uint8_t, nbr_blocks> pi;
    //!\brief Lower error bound for each block (accumulated values)
    std::array<uint8_t, nbr_blocks> l;
    //!\brief Upper error bound for each block (accumulated values)
    std::array<uint8_t, nbr_blocks> u;

    //!\brief Returns the number of blocks
    constexpr uint8_t blocks() const noexcept
    {
        return nbr_blocks;
    }
};

//!\brief Object storing information for a search (of a search scheme).
//        Number of blocks do not have to be known at compile time.
struct search_dyn
{
    //!\brief Type for storing the length of blocks
    typedef std::vector<uint64_t> blocks_length_type;

    //!\brief Order of blocks
    std::vector<uint8_t> pi;
    //!\brief Lower error bound for each block (accumulated values)
    std::vector<uint8_t> l;
    //!\brief Upper error bound for each block (accumulated values)
    std::vector<uint8_t> u;

    //!\brief Returns the number of blocks
    uint8_t blocks() const noexcept
    {
        return pi.size();
    }
};

//!\brief Type for storing search schemes. Number of blocks have to be known at compile time.
template <uint8_t nbr_searches, uint8_t nbr_blocks>
using search_scheme_type = std::array<search<nbr_blocks>, nbr_searches>;

//!\brief Type for storing search schemes. Number of blocks do not have to be known at compile time.
using search_scheme_dyn_type = std::vector<search_dyn>;

/*!\brief Search scheme that is optimal in the running time for the specified lower and upper error bound.
 * \tparam min_error Lower bound of errors.
 * \tparam max_error Upper bound of errors.
 * \detail Please note that the searches within each search scheme are sorted by their asymptotical running time
 *         (i.e. upper error bound string). s.t. easy to compute searches come first. This improves the running time of
 *         algorithms that abort after the first hit (e.g. search mode: best). Even though it is not guaranteed, this
 *         seems to be a good greedy approach.
 */
template <uint8_t min_error, uint8_t max_error>
inline int constexpr optimum_search_scheme;

//!\cond

template <>
inline search_scheme_type<1, 1> constexpr optimum_search_scheme<0, 0>
{{
    {{1}, {0}, {0}}
}};

template <>
inline search_scheme_type<2, 2> constexpr optimum_search_scheme<0, 1>
{{
    {{1, 2}, {0, 0}, {0, 1}},
    {{2, 1}, {0, 1}, {0, 1}}
}};

template <>
inline search_scheme_type<2, 2> constexpr optimum_search_scheme<1, 1>
{{
    {{1, 2}, {0, 1}, {0, 1}},
    {{2, 1}, {0, 1}, {0, 1}}
}};

template <>
inline search_scheme_type<3, 4> constexpr optimum_search_scheme<0, 2>
{{
    {{1, 2, 3, 4}, {0, 0, 1, 1}, {0, 0, 2, 2}},
    {{3, 2, 1, 4}, {0, 0, 0, 0}, {0, 1, 1, 2}},
    {{4, 3, 2, 1}, {0, 0, 0, 2}, {0, 1, 2, 2}}
}};

template <>
inline search_scheme_type<3, 4> constexpr optimum_search_scheme<1, 2>
{{
    {{1, 2, 3, 4}, {0, 0, 0, 1}, {0, 0, 2, 2}},
    {{3, 2, 1, 4}, {0, 0, 1, 1}, {0, 1, 1, 2}},
    {{4, 3, 2, 1}, {0, 0, 0, 2}, {0, 1, 2, 2}}
}};

template <>
inline search_scheme_type<3, 4> constexpr optimum_search_scheme<2, 2>
{{
    {{4, 3, 2, 1}, {0, 0, 1, 2}, {0, 0, 2, 2}},
    {{2, 3, 4, 1}, {0, 0, 0, 2}, {0, 1, 1, 2}},
    {{1, 2, 3, 4}, {0, 0, 0, 2}, {0, 1, 2, 2}}
}};

template <>
inline search_scheme_type<4, 5> constexpr optimum_search_scheme<0, 3>
{{
    {{5, 4, 3, 2, 1}, {0, 0, 0, 0, 0}, {0, 0, 3, 3, 3}},
    {{3, 4, 5, 2, 1}, {0, 0, 1, 1, 1}, {0, 1, 1, 2, 3}},
    {{2, 3, 4, 5, 1}, {0, 0, 0, 2, 2}, {0, 1, 2, 2, 3}},
    {{1, 2, 3, 4, 5}, {0, 0, 0, 0, 3}, {0, 2, 2, 3, 3}}
}};

// TODO: add the following missing optimum search schemes (computation has not finished yet)
// optimum_search_scheme<i, 3>, 1 < i <= 3
// optimum_search_scheme<i, 4>, 0 < i <= 4

//!\endcond

}

//!\}
