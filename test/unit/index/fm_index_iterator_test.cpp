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

#include <range/v3/algorithm/equal.hpp>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/search/fm_index/bi_fm_index.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/test/comparison.hpp>

#include <gtest/gtest.h>

using namespace seqan3;
using namespace seqan3::literal;
using namespace seqan3::test;

struct fm_index_byte_alphabet_traits
{
    using sdsl_index_type = sdsl::csa_wt<
        sdsl::wt_blcd<
            sdsl::bit_vector,
            sdsl::rank_support_v<>,
            sdsl::select_support_scan<>,
            sdsl::select_support_scan<0>
        >,
        16,
        10000000,
        sdsl::sa_order_sa_sampling<>,
        sdsl::isa_sampling<>,
        sdsl::byte_alphabet
    >;
};

struct bi_fm_index_byte_alphabet_traits
{
   using fm_index_traits = fm_index_byte_alphabet_traits;
   using rev_fm_index_traits = fm_index_byte_alphabet_traits;
};

template <typename T>
class fm_index_iterator_test : public ::testing::Test
{};

using fm_index_iterator_types = ::testing::Types<
        fm_index_iterator<fm_index<std::vector<dna4>, fm_index_default_traits>>,
        fm_index_iterator<fm_index<std::vector<dna4>, fm_index_byte_alphabet_traits>>,
        bi_fm_index_iterator<bi_fm_index<std::vector<dna4>, bi_fm_index_default_traits>>,
        bi_fm_index_iterator<bi_fm_index<std::vector<dna4>, bi_fm_index_byte_alphabet_traits>>>;

TYPED_TEST_CASE(fm_index_iterator_test, fm_index_iterator_types);

TYPED_TEST(fm_index_iterator_test, ctr)
{
    typename TypeParam::index_type::text_type text{"ACGACG"_dna4};
    typename TypeParam::index_type fm{text};

    // custom constructor
    TypeParam it0{fm};
    EXPECT_EQ(it0.query_length(), 0);
    EXPECT_EQ(it0.locate().size(), fm.size());

    // default construction (does not initialize the iterator)
    TypeParam it1;

    // copy construction
    TypeParam it2{it0};
    EXPECT_EQ(it0, it2);

    // copy assignment
    TypeParam it3 = it0;
    EXPECT_EQ(it0, it3);

    // move construction
    TypeParam it4{std::move(it0)};
    EXPECT_EQ(it0, it4);

    // move assigment
    TypeParam it5 = std::move(it0);
    EXPECT_EQ(it0, it5);
}

TYPED_TEST(fm_index_iterator_test, begin)
{
    typename TypeParam::index_type::text_type text{"ACGACG"_dna4};
    typename TypeParam::index_type fm{text};

    // begin
    TypeParam it(fm);
    EXPECT_TRUE(is_set_equal(it.locate(), (std::vector<uint64_t>{0, 1, 2, 3, 4, 5, 6}))); // sentinel position included
    EXPECT_EQ(it.query_length(), 0);
    EXPECT_EQ(it.count(), 7);
}

TYPED_TEST(fm_index_iterator_test, extend_right_range)
{
    typename TypeParam::index_type::text_type text{"ACGACG"_dna4};
    typename TypeParam::index_type fm{text};

    // successful extend_right(range)
    TypeParam it(fm);
    EXPECT_TRUE(it.extend_right("CG"_dna4));
    EXPECT_TRUE(is_set_equal(it.locate(), (std::vector<uint64_t>{1, 4})));
    EXPECT_EQ(it.query_length(), 2);
    EXPECT_EQ(it.count(), 2);

    EXPECT_TRUE(it.extend_right("A"_dna4));
    EXPECT_EQ(it.locate(), (std::vector<uint64_t>{1}));
    EXPECT_EQ(it.query_length(), 3);
    EXPECT_EQ(it.count(), 1);

    // unsuccessful extend_right(range), it remains untouched
    TypeParam it_cpy = it;
    EXPECT_FALSE(it.extend_right("A"_dna4));
    EXPECT_EQ(it, it_cpy);

    // extend_right(empty range)
    it_cpy = it;
    EXPECT_TRUE(it.extend_right(""_dna4));
    EXPECT_EQ(it, it_cpy);
}

// TODO: doesn't work with the current structure of typed tests
// TYPED_TEST(fm_index_iterator_test, extend_right_convertible_range)
// {
//     typename TypeParam::index_type::text_type text{"ANGACGNN"_dna5};
//     typename TypeParam::index_type fm{text};
//
//     // successful extend_right(range) using a different alphabet
//     TypeParam it(fm);
//     EXPECT_TRUE(it.extend_right("GA"_dna4));
//     EXPECT_EQ(it.locate(), (std::vector<uint64_t>{2}));
//     EXPECT_EQ(it.query_length(), 2);
// }

TYPED_TEST(fm_index_iterator_test, extend_right_char)
{
    typename TypeParam::index_type::text_type text{"ACGACG"_dna4};
    typename TypeParam::index_type fm{text};

    // successful extend_right(char)
    TypeParam it(fm);
    EXPECT_TRUE(it.extend_right(dna4::A));
    EXPECT_TRUE(is_set_equal(it.locate(), (std::vector<uint64_t>{0, 3})));
    EXPECT_EQ(it.query_length(), 1);

    EXPECT_TRUE(it.extend_right(dna4::C));
    EXPECT_TRUE(is_set_equal(it.locate(), (std::vector<uint64_t>{0, 3})));
    EXPECT_EQ(it.query_length(), 2);

    // unsuccessful extend_right(char), it remains untouched
    TypeParam it_cpy = it;
    EXPECT_FALSE(it.extend_right(dna4::C));
    EXPECT_EQ(it, it_cpy);
}

// TODO: doesn't work with the current structure of typed tests
// TYPED_TEST(fm_index_iterator_test, extend_right_convertible_char)
// {
//     typename TypeParam::index_type::text_type text{"ANGACGNN"_dna5};
//     typename TypeParam::index_type fm{text};
//
//     // successful extend_right(char) using a different alphabet
//     TypeParam it(fm);
//     EXPECT_TRUE(it.extend_right(dna4::A));
//     EXPECT_TRUE(is_set_equal(it.locate(), (std::vector<uint64_t>{0, 3})));
//     EXPECT_EQ(it.query_length(), 1);
// }

TYPED_TEST(fm_index_iterator_test, extend_right_range_and_cycle)
{
    typename TypeParam::index_type::text_type text{"ACGAACGC"_dna4};
    typename TypeParam::index_type fm{text};

    // successful extend_right() and cycle_back()
    TypeParam it(fm);
    EXPECT_TRUE(it.extend_right("ACGA"_dna4));
    EXPECT_EQ(it.locate(), (std::vector<uint64_t>{0}));
    EXPECT_EQ(it.query_length(), 4);

    EXPECT_TRUE(it.cycle_back());
    EXPECT_EQ(it.locate(), (std::vector<uint64_t>{4}));
    EXPECT_EQ(it.query_length(), 4);
}

TYPED_TEST(fm_index_iterator_test, extend_right_char_and_cycle)
{
    typename TypeParam::index_type::text_type text{"ACGAACGC"_dna4};
    typename TypeParam::index_type fm{text};

    // successful extend_right() and cycle_back()
    TypeParam it(fm);
    EXPECT_TRUE(it.extend_right(dna4::A));
    EXPECT_TRUE(is_set_equal(it.locate(), (std::vector<uint64_t>{0, 3, 4})));
    EXPECT_EQ(it.query_length(), 1);

    EXPECT_TRUE(it.cycle_back());
    EXPECT_TRUE(is_set_equal(it.locate(), (std::vector<uint64_t>{1, 5, 7})));
    EXPECT_EQ(it.query_length(), 1);
}

TYPED_TEST(fm_index_iterator_test, extend_right_and_cycle)
{
    typename TypeParam::index_type::text_type text{"ACGACG"_dna4};
    typename TypeParam::index_type fm{text};

    // successful extend_right() and cycle_back()
    TypeParam it(fm);
    EXPECT_TRUE(it.extend_right());
    EXPECT_TRUE(is_set_equal(it.locate(), (std::vector<uint64_t>{0, 3})));
    EXPECT_EQ(it.query_length(), 1);

    EXPECT_TRUE(it.cycle_back());
    EXPECT_TRUE(is_set_equal(it.locate(), (std::vector<uint64_t>{1, 4})));
    EXPECT_EQ(it.query_length(), 1);

    EXPECT_TRUE(it.extend_right());
    EXPECT_TRUE(is_set_equal(it.locate(), (std::vector<uint64_t>{1, 4})));
    EXPECT_EQ(it.query_length(), 2);

    // unsuccessful cycle_back(), it remains untouched
    TypeParam it_cpy = it;
    EXPECT_FALSE(it.cycle_back());
    EXPECT_EQ(it, it_cpy);

    // unsuccessful extend_right(), it remains untouched
    it = TypeParam(fm);
    EXPECT_TRUE(it.extend_right("GACG"_dna4));
    it_cpy = it;
    EXPECT_FALSE(it.extend_right());
    EXPECT_EQ(it, it_cpy);

    // cycle_back() cannot be called on the root node
    it = TypeParam(fm);
#ifndef NDEBUG
    EXPECT_DEATH(it.cycle_back(), "");
#endif
    EXPECT_EQ(it, TypeParam(fm));
}

TYPED_TEST(fm_index_iterator_test, query)
{
    typename TypeParam::index_type::text_type text{"ACGACG"_dna4};
    typename TypeParam::index_type fm{text};

    // query()
    TypeParam it(fm);
    EXPECT_TRUE(it.extend_right("ACG"_dna4));
    EXPECT_TRUE(ranges::equal(*it, "ACG"_dna4));
    EXPECT_TRUE(ranges::equal(it.query(), "ACG"_dna4));
}

// TODO: test last_char()
TYPED_TEST(fm_index_iterator_test, incomplete_alphabet)
{
    // search a char that does not occur in the text (higher rank than largest char occurring in text)
    {
        typename TypeParam::index_type::text_type text{"ACGACG"_dna4};
        typename TypeParam::index_type fm{text};
        TypeParam it = TypeParam(fm);
        EXPECT_FALSE(it.extend_right(dna4::T));
        EXPECT_EQ(it, TypeParam(fm));
    }

    // search a char that does not occur in the text (smaller rank than smallest char occurring in text)
    {
        typename TypeParam::index_type::text_type text{"CGTCGT"_dna4};
        typename TypeParam::index_type fm{text};
        TypeParam it = TypeParam(fm);
        EXPECT_FALSE(it.extend_right(dna4::A));
        EXPECT_EQ(it, TypeParam(fm));
    }

    // search a char that does not occur in the text
    // (some rank that is neither the smallest nor the highest occurring in text)
    {
        typename TypeParam::index_type::text_type text{"ATATAT"_dna4};
        typename TypeParam::index_type fm{text};
        TypeParam it = TypeParam(fm);
        EXPECT_FALSE(it.extend_right(dna4::C));
        EXPECT_FALSE(it.extend_right(dna4::G));
        EXPECT_FALSE(it.extend_right("ACGT"_dna4));
        EXPECT_FALSE(it.extend_right("G"_dna4));
        EXPECT_EQ(it, TypeParam(fm));

        EXPECT_TRUE(it.extend_right(dna4::A));
        EXPECT_TRUE(it.cycle_back());
        EXPECT_TRUE(ranges::equal(it.query(), "T"_dna4));
    }
}

TYPED_TEST(fm_index_iterator_test, lazy_locate)
{
    typename TypeParam::index_type::text_type text{"ACGTACGT"_dna4};
    typename TypeParam::index_type fm{text};

    TypeParam it = TypeParam(fm);
    it.extend_right("ACG"_dna4);

    EXPECT_TRUE(ranges::equal(it.locate(), it.lazy_locate()));
}

TEST(fm_index, concepts)
{
    EXPECT_TRUE(fm_index_iterator_concept<fm_index_iterator<fm_index<std::vector<dna4>>>>);
    EXPECT_TRUE(fm_index_iterator_concept<bi_fm_index_iterator<bi_fm_index<std::vector<dna4>>>>);
    EXPECT_TRUE(bi_fm_index_iterator_concept<bi_fm_index_iterator<bi_fm_index<std::vector<dna4>>>>);
}
