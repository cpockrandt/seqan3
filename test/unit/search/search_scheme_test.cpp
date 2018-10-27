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

#include <algorithm>
#include <type_traits>

#include "helper_search_scheme.hpp"

#include <seqan3/search/algorithm/detail/search_scheme_algorithm.hpp>

#include <gtest/gtest.h>

using namespace seqan3;

template <uint8_t min_error, uint8_t max_error, bool precomputed_scheme>
void error_distributions(auto & expected, auto & actual)
{
    if constexpr (precomputed_scheme)
    {
        auto const & oss{detail::optimum_search_scheme<min_error, max_error>};
        search_scheme_error_distribution(actual, oss);
        search_scheme_error_distribution(expected, trivial_search_scheme(min_error, max_error, oss.front().blocks()));
    }
    else
    {
        auto const & ss{detail::compute_ss(min_error, max_error)};
        search_scheme_error_distribution(actual, ss);
        search_scheme_error_distribution(expected, trivial_search_scheme(min_error, max_error, ss.front().blocks()));
    }
    std::sort(expected.begin(), expected.end());
    std::sort(actual.begin(), actual.end());
}

TEST(search_scheme_test, error_distribution_coverage_optimum_search_schemes)
{
    std::vector<std::vector<uint8_t> > expected, actual;

    error_distributions<0, 0, true>(expected, actual);
    EXPECT_EQ(actual, expected);

    error_distributions<0, 1, true>(expected, actual);
    EXPECT_EQ(actual, expected);

    error_distributions<1, 1, true>(expected, actual);
    EXPECT_EQ(actual, expected);

    error_distributions<0, 2, true>(expected, actual);
    EXPECT_EQ(actual, expected);

    error_distributions<1, 2, true>(expected, actual);
    EXPECT_EQ(actual, expected);

    error_distributions<2, 2, true>(expected, actual);
    EXPECT_EQ(actual, expected);

    error_distributions<0, 3, true>(expected, actual);
    EXPECT_EQ(actual, expected);
}

TEST(search_scheme_test, error_distribution_coverage_computed_search_schemes)
{
    std::vector<std::vector<uint8_t> > expected, actual;

    error_distributions<0, 0, false>(expected, actual);
    EXPECT_EQ(actual, expected);

    error_distributions<0, 1, false>(expected, actual);
    EXPECT_EQ(actual, expected);

    error_distributions<1, 1, false>(expected, actual);
    EXPECT_EQ(actual, expected);

    error_distributions<0, 2, false>(expected, actual);
    EXPECT_EQ(actual, expected);

    error_distributions<1, 2, false>(expected, actual);
    EXPECT_EQ(actual, expected);

    error_distributions<2, 2, false>(expected, actual);
    EXPECT_EQ(actual, expected);

    error_distributions<0, 3, false>(expected, actual);
    EXPECT_EQ(actual, expected);

    error_distributions<1, 3, false>(expected, actual);
    EXPECT_EQ(actual, expected);

    error_distributions<2, 3, false>(expected, actual);
    EXPECT_EQ(actual, expected);

    error_distributions<3, 3, false>(expected, actual);
    EXPECT_EQ(actual, expected);

    error_distributions<3, 5, false>(expected, actual);
    EXPECT_EQ(actual, expected);
    error_distributions<0, 6, false>(expected, actual);
    EXPECT_EQ(actual, expected);
    error_distributions<7, 7, false>(expected, actual);
    EXPECT_EQ(actual, expected);
}

template <uint8_t min_error, uint8_t max_error, bool precomputed_scheme>
bool check_disjoint_search_scheme()
{
    std::vector<std::vector<uint8_t> > error_distributions;

    auto const & oss{detail::optimum_search_scheme<min_error, max_error>};
    search_scheme_error_distribution(error_distributions, oss);
    uint64_t size = error_distributions.size();
    std::sort(error_distributions.begin(), error_distributions.end());
    error_distributions.erase(std::unique(error_distributions.begin(), error_distributions.end()), error_distributions.end());
    return size == error_distributions.size();
}

TEST(search_scheme_test, error_distribution_disjoint_optimum_search_schemes)
{
    bool ret;

    ret = check_disjoint_search_scheme<0, 0, false>();
    EXPECT_TRUE(ret);
    ret = check_disjoint_search_scheme<0, 1, false>();
    EXPECT_TRUE(ret);
    ret = check_disjoint_search_scheme<0, 2, false>();
    EXPECT_TRUE(ret);
    ret = check_disjoint_search_scheme<0, 3, false>();
    EXPECT_TRUE(ret);
}
