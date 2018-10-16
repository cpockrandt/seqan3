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
 * \brief Provides the data structures and precomputed instances for optimum search schemes.
 */

#pragma once

/*!\addtogroup search
 * \{
 */

namespace seqan3::detail
{

template <uint8_t nbr_blocks>
struct search
{
    std::array<uint8_t, nbr_blocks> pi;
    std::array<uint8_t, nbr_blocks> l;
    std::array<uint8_t, nbr_blocks> u;

    // TODO?: search() : pi(), l(), u() {};

    uint8_t blocks() const noexcept
    {
        return pi.size();
    }
};

template <uint8_t nbr_searches, uint8_t nbr_blocks>
using search_scheme = std::array<search<nbr_blocks>, nbr_searches>;

struct search_dyn
{
    std::vector<uint8_t> pi;
    std::vector<uint8_t> l;
    std::vector<uint8_t> u;

    uint8_t blocks() const noexcept
    {
        return pi.size();
    }
};

using search_scheme_dyn = std::vector<search_dyn>;

template <uint8_t min_errors, uint8_t max_errors>
struct optimum_search_scheme;

template <>
struct optimum_search_scheme<0, 0>
{
    static constexpr search_scheme<1, 3> value
    {{
        {{2, 1, 3}, {0, 0, 0}, {0, 0, 0}}
    }};
};

template <>
struct optimum_search_scheme<0, 1>
{
    static constexpr search_scheme<2, 2> value
    {{
        {{1, 2}, {0, 0}, {0, 1}},
        {{2, 1}, {0, 1}, {0, 1}}
    }};
};

template <>
struct optimum_search_scheme<1, 1>
{
    static constexpr search_scheme<2, 2> value
    {{
        {{1, 2}, {0, 1}, {0, 1}},
        {{2, 1}, {0, 1}, {0, 1}}
    }};
};

template <>
struct optimum_search_scheme<0, 2>
{
    static constexpr search_scheme<3, 4> value
    {{
        {{1, 2, 3, 4}, {0, 0, 1, 1}, {0, 0, 2, 2}},
        {{3, 2, 1, 4}, {0, 0, 0, 0}, {0, 1, 1, 2}},
        {{4, 3, 2, 1}, {0, 0, 0, 2}, {0, 1, 2, 2}}
    }};
};

template <>
struct optimum_search_scheme<0, 3>
{
    static constexpr search_scheme<4, 5> value
    {{
        {{1, 2, 3, 4, 5}, {0, 0, 0, 0, 3}, {0, 2, 2, 3, 3}},
        {{2, 3, 4, 5, 1}, {0, 0, 0, 2, 2}, {0, 1, 2, 2, 3}},
        {{3, 4, 5, 2, 1}, {0, 0, 1, 1, 1}, {0, 1, 1, 2, 3}},
        {{5, 4, 3, 2, 1}, {0, 0, 0, 0, 0}, {0, 0, 3, 3, 3}}
    }};
};

// TODO: OSS with min error != 0
// TODO: compute OSS for 4 errors

}

//!\}
