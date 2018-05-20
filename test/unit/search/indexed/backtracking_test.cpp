// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2017, Knut Reinert, FU Berlin
// Copyright (c) 2016-2017, Knut Reinert & MPI Molekulare Genetik
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
// ==========================================================================

#include <sstream>

#include <gtest/gtest.h>

#include <seqan3/search/indexed/search_schemes.hpp>

using namespace seqan3;

struct my_search_parameters
{
    // static constexpr uint8_t min_errors = 0u;
    static constexpr uint8_t max_errors = 1u;
    static constexpr search_parameters_metric metric {search_parameters_metric::hamming};
    static constexpr bool output_alignments = false;
};

void random_string(sdsl::int_vector<8> & iv, uint64_t length)
{
    iv.resize(length);
    sdsl::util::set_random_bits(iv);
}

TEST(backtracking, hamming)
{
    sdsl::int_vector<8> text, pattern;
    random_string(text, 10000);
    random_string(pattern, 5);

    typedef sdsl::csa_wt<sdsl::wt_blcd<sdsl::bit_vector, sdsl::rank_support_v<>, sdsl::select_support_scan<>, sdsl::select_support_scan<0> >, 10, 10000000> index_t;

    index_t index;
    construct_im(index, text, 0);

    // TODO: test nicht zusammenhängende alphabete. alles mit und ohne sentinels! alles mit allen möglichen indices!

    std::vector<std::pair<uint64_t, uint8_t> > hits; // (position, errors)
    auto callback = [&hits](uint64_t l, uint64_t r, uint8_t errors)
    {
        for (unsigned i = l; i < r; ++i)
            hits.push_back({i, errors}); // TODO: i ist SA-pos, nicht text-pos!!!
    };

    const my_search_parameters params;
    detail::search_backtracking(index, pattern.begin(), pattern.end(), params, callback);

    // remove duplicates
    std::sort(hits.begin(), hits.end());
    hits.erase(std::unique(hits.begin(), hits.end()), hits.end());

    uint64_t hit_pos = 0;
    for (uint64_t i = 0; (i < text.size() - pattern.size() + 1) && hit_pos < hits.size(); ++i)
    {
        uint8_t errors = 0;
        for (uint64_t j = 0; j < pattern.size() && errors <= params.max_errors; ++j)
        {
            errors += text[i + j] != pattern[i];
        }

        if (errors <= params.max_errors) // expecting a hit!
        {
            EXPECT_EQ(i, hits[hit_pos].first);
            EXPECT_EQ(errors, hits[hit_pos].second);
            ++hit_pos;
        }
    }
    EXPECT_EQ(hit_pos, hits.size());
}
