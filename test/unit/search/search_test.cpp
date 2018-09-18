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

#include "helper.hpp"

#include <seqan3/search/algorithm/all.hpp>

#include <gtest/gtest.h>

using namespace seqan3;
using namespace seqan3::literal;
using namespace seqan3::search_cfg;

template <typename T>
class search_test : public ::testing::Test
{
public:
    std::vector<dna4> text{"ACGTACGTACGT"_dna4};
    T index{text};
};

using fm_index_types = ::testing::Types<fm_index<std::vector<dna4>>/*, bi_fm_index<std::vector<dna4>>*/>;

TYPED_TEST_CASE(search_test, fm_index_types);

TYPED_TEST(search_test, error_free)
{
    using hits_result_t = std::vector<typename TypeParam::size_type>;

    {
        // successful and unsuccesful exact search without cfg
        EXPECT_EQ(sort(search(this->index, "ACGT"_dna4)), (hits_result_t{0, 4, 8}));
        EXPECT_EQ(sort(search(this->index, "ACGG"_dna4)), (hits_result_t{}));
    }

    {
        // successful and unsuccesful exact search with empty cfg
        detail::configuration const cfg;
        EXPECT_EQ(sort(search(this->index, "ACGT"_dna4, cfg)), (hits_result_t{0, 4, 8}));
        EXPECT_EQ(sort(search(this->index, "ACGG"_dna4, cfg)), (hits_result_t{}));
    }

    {
        // successful and unsuccesful exact search using max_total_error
        detail::configuration const cfg = max_error(total{0}, substitution{0}, insertion{0}, deletion{0});
        EXPECT_EQ(sort(search(this->index, "ACGT"_dna4, cfg)), (hits_result_t{0, 4, 8}));
        EXPECT_EQ(sort(search(this->index, "ACGG"_dna4, cfg)), (hits_result_t{}));
    }

    {
        // successful and unsuccesful exact search using max_total_error_rate
        detail::configuration const cfg = max_error_rate(total{.0}, substitution{.0}, insertion{.0}, deletion{.0});
        EXPECT_EQ(sort(search(this->index, "ACGT"_dna4, cfg)), (hits_result_t{0, 4, 8}));
        EXPECT_EQ(sort(search(this->index, "ACGG"_dna4, cfg)), (hits_result_t{}));
    }
}

TYPED_TEST(search_test, multiple_queries)
{
    using hits_result_t = std::vector<std::vector<typename TypeParam::size_type>>;
    std::vector<std::vector<dna4>> const queries{{"GG"_dna4, "ACGTACGTACGT"_dna4, "ACGTA"_dna4}};

    detail::configuration const cfg = max_error_rate(total{.0}, substitution{.0}, insertion{.0}, deletion{.0});
    EXPECT_EQ(sort(search(this->index, queries, cfg)), (hits_result_t{{}, {0}, {0, 4}})); // 0, 1 and 2 hits
}

TYPED_TEST(search_test, error_substitution)
{
    using hits_result_t = std::vector<typename TypeParam::size_type>;

    {
        detail::configuration const cfg = max_error_rate(total{.25}, substitution{.25}, insertion{.0}, deletion{.0});

        EXPECT_EQ(sort(search(this->index, "ACGT"_dna4    , cfg)), (hits_result_t{0, 4, 8})); // exact match
        EXPECT_EQ(sort(search(this->index, "CGG"_dna4     , cfg)), (hits_result_t{}));        // not enough mismatches
        EXPECT_EQ(sort(search(this->index, "CGTC"_dna4    , cfg)), (hits_result_t{1, 5}));    // 1 mismatch
        EXPECT_EQ(sort(search(this->index, "ACGGACG"_dna4 , cfg)), (hits_result_t{0, 4}));    // 1 mismatch
        EXPECT_EQ(sort(search(this->index, "ACGGACGG"_dna4, cfg)), (hits_result_t{0, 4}));    // 2 mismatches
    }

    {
        // TODO: remove insertion and deletion
        detail::configuration const cfg = max_error(total{1}, substitution{1}, insertion{0}, deletion{0});

        EXPECT_EQ(sort(search(this->index, "ACGT"_dna4    , cfg)), (hits_result_t{0, 4, 8})); // exact match
        EXPECT_EQ(sort(search(this->index, "CGTTT"_dna4   , cfg)), (hits_result_t{}));        // not enough mismatches
        EXPECT_EQ(sort(search(this->index, "CGG"_dna4     , cfg)), (hits_result_t{1, 5, 9})); // 1 mismatch
        EXPECT_EQ(sort(search(this->index, "ACGGACG"_dna4 , cfg)), (hits_result_t{0, 4}));    // 1 mismatch
        EXPECT_EQ(sort(search(this->index, "CGTCCGTA"_dna4, cfg)), (hits_result_t{1}));       // 1 mismatch
    }

    {
        detail::configuration const cfg = max_error(total{1}, substitution{1}, insertion{0}, deletion{0});

        EXPECT_EQ(sort(search(this->index, "ACGT"_dna4    , cfg)), (hits_result_t{0, 4, 8})); // exact match
        EXPECT_EQ(sort(search(this->index, "CGTTT"_dna4   , cfg)), (hits_result_t{}));        // not enough mismatches
        EXPECT_EQ(sort(search(this->index, "CGG"_dna4     , cfg)), (hits_result_t{1, 5, 9})); // 1 mismatch
        EXPECT_EQ(sort(search(this->index, "ACGGACG"_dna4 , cfg)), (hits_result_t{0, 4}));    // 1 mismatch
        EXPECT_EQ(sort(search(this->index, "CGTCCGTA"_dna4, cfg)), (hits_result_t{1}));       // 1 mismatch
    }
}

TYPED_TEST(search_test, error_insertion)
{
    using hits_result_t = std::vector<typename TypeParam::size_type>;

    {
        detail::configuration const cfg = max_error_rate(total{.25}, substitution{.0}, insertion{.25}, deletion{.0});

        // exact match and insertion at the beginning of the query
        EXPECT_EQ(sort(search(this->index, "ACGT"_dna4    , cfg)), (hits_result_t{0, 1, 4, 5, 8, 9}));
        // 1 insertion
        EXPECT_EQ(sort(search(this->index, "CCGT"_dna4    , cfg)), (hits_result_t{1, 5, 9}));
        // 2 insertions
        EXPECT_EQ(sort(search(this->index, "ACCGGTAC"_dna4, cfg)), (hits_result_t{0, 4}));
        // 2 insertions necessary, only 1 allowed
        EXPECT_EQ(sort(search(this->index, "ACCGG"_dna4   , cfg)), (hits_result_t{}));
        // deletion necessary, not allowed
        EXPECT_EQ(sort(search(this->index, "ACTACGT"_dna4 , cfg)), (hits_result_t{}));
    }

    {
        detail::configuration const cfg = max_error(total{1}, substitution{0}, insertion{1}, deletion{0});

        // exact match and insertion at the beginning of the query
        EXPECT_EQ(sort(search(this->index, "ACGT"_dna4    , cfg)), (hits_result_t{0, 1, 4, 5, 8, 9}));
        // 1 insertion
        EXPECT_EQ(sort(search(this->index, "CCGT"_dna4    , cfg)), (hits_result_t{1, 5, 9}));
        // 2 insertions necessary, only 1 allowed
        EXPECT_EQ(sort(search(this->index, "ACCGGTAC"_dna4, cfg)), (hits_result_t{}));
        // deletion necessary, not allowed
        EXPECT_EQ(sort(search(this->index, "ACTACGT"_dna4 , cfg)), (hits_result_t{}));
    }
}

TYPED_TEST(search_test, error_deletion)
{
    using hits_result_t = std::vector<typename TypeParam::size_type>;

    {
        detail::configuration const cfg = max_error_rate(total{.25}, substitution{.0}, insertion{.0}, deletion{.25});

        // exact match, no deletion
        EXPECT_EQ(sort(search(this->index, "ACGT"_dna4    , cfg)), (hits_result_t{0, 4, 8}));
        // not enough max errors
        EXPECT_EQ(sort(search(this->index, "AGT"_dna4     , cfg)), (hits_result_t{}));
        // one deletion (C)
        EXPECT_EQ(sort(search(this->index, "AGTA"_dna4    , cfg)), (hits_result_t{0, 4}));
        // two deletion (C)
        EXPECT_EQ(sort(search(this->index, "AGTAGTAC"_dna4, cfg)), (hits_result_t{0}));
        // no deletion at beginning. 0 and 4 cannot be reported
        EXPECT_EQ(sort(search(this->index, "CGTACGT"_dna4 , cfg)), (hits_result_t{1, 5}));
    }

    {
        detail::configuration const cfg = max_error(total{1}, substitution{0}, insertion{0}, deletion{1});

        // exact match, no deletion
        EXPECT_EQ(sort(search(this->index, "ACGT"_dna4    , cfg)), (hits_result_t{0, 4, 8}));
        // one deletion (C)
        EXPECT_EQ(sort(search(this->index, "AGTA"_dna4    , cfg)), (hits_result_t{0, 4}));
        // 2 deletions necessary, only 1 allowed
        EXPECT_EQ(sort(search(this->index, "AGTAGTAC"_dna4, cfg)), (hits_result_t{}));
        // no deletion at beginning. 0 and 4 cannot be reported
        EXPECT_EQ(sort(search(this->index, "CGTACGT"_dna4 , cfg)), (hits_result_t{1, 5}));
    }
}

TYPED_TEST(search_test, error_levenshtein)
{
    using hits_result_t = std::vector<typename TypeParam::size_type>;

    {
        // TODO: remove everything except total{2}
        detail::configuration const cfg = max_error(total{2}, substitution{2}, insertion{2}, deletion{2});
        EXPECT_EQ(sort(search(this->index, "CCGT"_dna4, cfg)), (hits_result_t{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}));
    }

    {
        detail::configuration const cfg = max_error(total{2}, substitution{2}, insertion{2}, deletion{2});
        EXPECT_EQ(sort(search(this->index, "CCGT"_dna4, cfg)), (hits_result_t{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}));
    }
}

TYPED_TEST(search_test, search_strategy_all)
{
    using hits_result_t = std::vector<typename TypeParam::size_type>;

    {
        detail::configuration const cfg = max_error(total{1}, substitution{1}, insertion{1}, deletion{1});
        EXPECT_EQ(sort(search(this->index, "ACGT"_dna4, cfg)), (hits_result_t{0, 1, 4, 5, 8, 9}));
    }

    {
        detail::configuration const cfg = max_error(total{1}, substitution{1}, insertion{1}, deletion{1})
                                        | strategy_all();
        EXPECT_EQ(sort(search(this->index, "ACGT"_dna4, cfg)), (hits_result_t{0, 1, 4, 5, 8, 9}));
    }
}

TYPED_TEST(search_test, search_strategy_best)
{
    using hits_result_t = std::vector<typename TypeParam::size_type>;

    {
        detail::configuration const cfg = max_error(total{1}, substitution{1}, insertion{1}, deletion{1})
                                        | strategy_best();

        hits_result_t possible_hits{0, 4, 8}; // any of 0, 4, 8 ... 1, 5, 9 are not best hits
        hits_result_t result = search(this->index, "ACGT"_dna4, cfg);
        ASSERT_EQ(result.size(), 1);
        EXPECT_TRUE(std::find(possible_hits.begin(), possible_hits.end(), result[0]) != possible_hits.end());

        EXPECT_EQ(search(this->index, "AAAA"_dna4, cfg), (hits_result_t{})); // no hit
    }
}

TYPED_TEST(search_test, search_strategy_all_best)
{
    using hits_result_t = std::vector<typename TypeParam::size_type>;

    {
        detail::configuration const cfg = max_error(total{1}, substitution{1}, insertion{1}, deletion{1})
                                        | strategy_all_best();

        EXPECT_EQ(sort(search(this->index, "ACGT"_dna4, cfg)), (hits_result_t{0, 4, 8})); // 1, 5, 9 are not best hits

        EXPECT_EQ(search(this->index, "AAAA"_dna4, cfg), (hits_result_t{})); // no hit
    }
}

TYPED_TEST(search_test, search_strategy_strata)
{
    using hits_result_t = std::vector<typename TypeParam::size_type>;

    {
        detail::configuration const cfg = max_error(total{1}, substitution{1}, insertion{1}, deletion{1})
                                        | strategy_strata(0);
        EXPECT_EQ(sort(search(this->index, "ACGT"_dna4, cfg)), (hits_result_t{0, 4, 8}));
    }

    {
        detail::configuration const cfg = max_error(total{1}, substitution{1}, insertion{1}, deletion{1})
                                        | strategy_strata(1);
        EXPECT_EQ(sort(search(this->index, "ACGT"_dna4, cfg)), (hits_result_t{0, 1, 4, 5, 8, 9}));
    }

    {
        detail::configuration const cfg = max_error(total{1}, substitution{1}, insertion{1}, deletion{1})
                                        | strategy_strata(1);
        EXPECT_EQ(search(this->index, "AAAA"_dna4, cfg), (hits_result_t{})); // no hit
    }

    // {
    //     // best hit ACGT with 1 error, i.e. 1+1
    //     detail::configuration const cfg = max_total_error(1) | strategy_strata(1);
    //     EXPECT_EQ(sort(search(this->index, "CCGT"_dna4, cfg)), (hits_result_t{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}));
    // }

    // {
    //     // best hit ACGT with 1 error, i.e. 1+1
    //     detail::configuration const cfg = max_total_error(1) | strategy_strata(1);
    //     EXPECT_EQ(sort(search(this->index, "CCGT"_dna4, cfg)), (hits_result_t{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}));
    // }
}

// TYPED_TEST(search_test, return_iterator_index)
// {
// }
//
// TYPED_TEST(search_test, on_hit)
// {
// }
