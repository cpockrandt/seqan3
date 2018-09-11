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

#include <type_traits>

#include <seqan3/search/all.hpp>
#include <seqan3/test/comparison.hpp>

#include <gtest/gtest.h>

using namespace seqan3;
using namespace seqan3::test;
using namespace seqan3::literal;
using namespace seqan3::search_cfg;

template <typename T>
class search_test : public ::testing::Test
{
public:
    std::vector<dna4> text{"ACGTACGT"_dna4};
    T index{text};
};

using fm_index_types = ::testing::Types<fm_index<std::vector<dna4>>, bi_fm_index<std::vector<dna4>>>;

TYPED_TEST_CASE(search_test, fm_index_types);

template <typename text_position_t>
inline void compare_hits(std::vector<uint64_t> const & actual, std::vector<text_position_t> const & expected)
{
    EXPECT_TRUE(is_set_equal(actual, expected));
}

template <fm_index_iterator_concept iterator_t, typename text_position_t>
inline void compare_hits(std::vector<iterator_t> const & actual, std::vector<text_position_t> const & expected)
{
    // TODO: replace this with a nice view
    std::vector<text_position_t> actual_located;
    actual_located.reserve(expected.size());
    for (iterator_t const & iter : actual)
    {
        auto const & pos = iter.locate();
        actual_located.insert(actual_located.end(), pos.begin(), pos.end());
    }
    EXPECT_TRUE(is_set_equal(actual_located, expected));
}

template <typename actual_t, typename text_position_t>
inline void compare_hits(std::vector<std::vector<actual_t>> const & actual,
                         std::vector<std::vector<text_position_t>> const & expected)
{
    EXPECT_EQ(actual.size(), expected.size());
    for (uint64_t i = 0; i < actual.size(); ++i)
    {
        compare_hits(actual[i], expected[i]);
    }
}

TYPED_TEST(search_test, error_free)
{
    std::vector<dna4> query{"ACGT"_dna4};
    std::vector<uint64_t> hits{0, 4};

    {
        compare_hits(search(this->index, query), hits);
    }

    {
        detail::configuration const cfg;
        compare_hits(search(this->index, query), hits);
    }

    {
        detail::configuration const cfg = max_error(0);
        compare_hits(search(this->index, query, cfg), hits);
    }

    {
        detail::configuration const cfg = max_error_rate(.0);
        compare_hits(search(this->index, query, cfg), hits);
    }
}

TYPED_TEST(search_test, multiple_queries)
{
    std::vector<std::vector<dna4>> queries{{"ACGT"_dna4, "GG"_dna4, "CGTA"_dna4}};
    std::vector<std::vector<uint64_t>> hits{{0, 4}, {}, {1}};

    detail::configuration const cfg = max_error_rate(.0);
    compare_hits(search(this->index, queries, cfg), hits);
}

TYPED_TEST(search_test, error_substitution)
{
    std::vector<std::vector<dna4>> queries{{"ACGT"_dna4, "ACGGACG"_dna4, "CGTC"_dna4, "CGG"_dna4}};
    std::vector<std::vector<uint64_t>> hits{{0, 4}, {0}, {1}, {}};

    detail::configuration const cfg = max_error_rate(.25) | error_type(error_type_enum::substitution);
    compare_hits(search(this->index, queries, cfg), hits);
}

TYPED_TEST(search_test, error_insertion)
{
}

TYPED_TEST(search_test, error_deletion)
{
}

TYPED_TEST(search_test, error_levenshtein)
{
}

TYPED_TEST(search_test, search_strategy_all)
{
}

TYPED_TEST(search_test, search_strategy_best)
{
}

TYPED_TEST(search_test, search_strategy_all_best)
{
}

TYPED_TEST(search_test, search_strategy_strata)
{
}

TYPED_TEST(search_test, return_iterator)
{
}

TYPED_TEST(search_test, on_hit)
{
}
