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

#include <seqan3/index/fm_index.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
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

template <typename T>
class fm_index_iterator_test : public ::testing::Test
{};

using fm_index_iterator_types = ::testing::Types<
        fm_index_iterator<fm_index<std::vector<dna4>, fm_index_default_traits>>,
        fm_index_iterator<fm_index<std::vector<dna4>, fm_index_byte_alphabet_traits>>>;

TYPED_TEST_CASE(fm_index_iterator_test, fm_index_iterator_types);

TYPED_TEST(fm_index_iterator_test, ctr)
{
    typename TypeParam::index_type::text_type text{"ACGACG"_dna4};
    typename TypeParam::index_type fm{text};

    // custom constructor
    TypeParam it0{fm};
    EXPECT_EQ(it0.depth(), 0);
    EXPECT_EQ(it0.locate().size(), fm.size());

    // default construction (does not set the iterator to the root node)
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

TYPED_TEST(fm_index_iterator_test, root_node)
{
    typename TypeParam::index_type::text_type text{"ACGACG"_dna4};
    typename TypeParam::index_type fm{text};

    // root
    TypeParam it(fm);
    EXPECT_TRUE(is_set_equal(it.locate(), (std::vector<uint64_t>{0, 1, 2, 3, 4, 5, 6}))); // sentinel position included
    EXPECT_EQ(it.depth(), 0);
    EXPECT_EQ(it.count(), 7);
}

TYPED_TEST(fm_index_iterator_test, down_range)
{
    typename TypeParam::index_type::text_type text{"ACGACG"_dna4};
    typename TypeParam::index_type fm{text};

    // successful down(range)
    TypeParam it(fm);
    EXPECT_TRUE(it.down("CG"_dna4));
    EXPECT_TRUE(is_set_equal(it.locate(), (std::vector<uint64_t>{1, 4})));
    EXPECT_EQ(it.depth(), 2);
    EXPECT_EQ(it.count(), 2);

    EXPECT_TRUE(it.down("A"_dna4));
    EXPECT_EQ(it.locate(), (std::vector<uint64_t>{1}));
    EXPECT_EQ(it.depth(), 3);
    EXPECT_EQ(it.count(), 1);

    // unsuccessful down(range), it remains untouched
    TypeParam it_cpy = it;
    EXPECT_FALSE(it.down("A"_dna4));
    EXPECT_EQ(it, it_cpy);

    // down(range) does not take an empty range
    it_cpy = it;
    ASSERT_DEATH(it.down(""_dna4), "");
    EXPECT_EQ(it, it_cpy);
}

// TODO: doesn't work with the current structure of typed tests
// TYPED_TEST(fm_index_iterator_test, down_convertible_range)
// {
//     typename TypeParam::index_type::text_type text{"ANGACGNN"_dna5};
//     typename TypeParam::index_type fm{text};
//
//     // successful down(range) using a different alphabet
//     TypeParam it(fm);
//     EXPECT_TRUE(it.down("GA"_dna4));
//     EXPECT_EQ(it.locate(), (std::vector<uint64_t>{2}));
//     EXPECT_EQ(it.depth(), 2);
// }

TYPED_TEST(fm_index_iterator_test, down_char)
{
    typename TypeParam::index_type::text_type text{"ACGACG"_dna4};
    typename TypeParam::index_type fm{text};

    // successful down(char)
    TypeParam it(fm);
    EXPECT_TRUE(it.down(dna4::A));
    EXPECT_TRUE(is_set_equal(it.locate(), (std::vector<uint64_t>{0, 3})));
    EXPECT_EQ(it.depth(), 1);

    EXPECT_TRUE(it.down(dna4::C));
    EXPECT_TRUE(is_set_equal(it.locate(), (std::vector<uint64_t>{0, 3})));
    EXPECT_EQ(it.depth(), 2);

    // unsuccessful down(char), it remains untouched
    TypeParam it_cpy = it;
    EXPECT_FALSE(it.down(dna4::C));
    EXPECT_EQ(it, it_cpy);
}

// TODO: doesn't work with the current structure of typed tests
// TYPED_TEST(fm_index_iterator_test, down_convertible_char)
// {
//     typename TypeParam::index_type::text_type text{"ANGACGNN"_dna5};
//     typename TypeParam::index_type fm{text};
//
//     // successful down(char) using a different alphabet
//     TypeParam it(fm);
//     EXPECT_TRUE(it.down(dna4::A));
//     EXPECT_TRUE(is_set_equal(it.locate(), (std::vector<uint64_t>{0, 3})));
//     EXPECT_EQ(it.depth(), 1);
// }

TYPED_TEST(fm_index_iterator_test, down_range_and_right)
{
    typename TypeParam::index_type::text_type text{"ACGAACGC"_dna4};
    typename TypeParam::index_type fm{text};

    // successful down() and right()
    TypeParam it(fm);
    EXPECT_TRUE(it.down("ACGA"_dna4));
    EXPECT_EQ(it.locate(), (std::vector<uint64_t>{0}));
    EXPECT_EQ(it.depth(), 4);

    EXPECT_TRUE(it.right());
    EXPECT_EQ(it.locate(), (std::vector<uint64_t>{4}));
    EXPECT_EQ(it.depth(), 4);
}

TYPED_TEST(fm_index_iterator_test, down_char_and_right)
{
    typename TypeParam::index_type::text_type text{"ACGAACGC"_dna4};
    typename TypeParam::index_type fm{text};

    // successful down() and right()
    TypeParam it(fm);
    EXPECT_TRUE(it.down(dna4::A));
    EXPECT_TRUE(is_set_equal(it.locate(), (std::vector<uint64_t>{0, 3, 4})));
    EXPECT_EQ(it.depth(), 1);

    EXPECT_TRUE(it.right());
    EXPECT_TRUE(is_set_equal(it.locate(), (std::vector<uint64_t>{1, 5, 7})));
    EXPECT_EQ(it.depth(), 1);
}

TYPED_TEST(fm_index_iterator_test, down_and_right)
{
    typename TypeParam::index_type::text_type text{"ACGACG"_dna4};
    typename TypeParam::index_type fm{text};

    // successful down() and right()
    TypeParam it(fm);
    EXPECT_TRUE(it.down());
    EXPECT_TRUE(is_set_equal(it.locate(), (std::vector<uint64_t>{0, 3})));
    EXPECT_EQ(it.depth(), 1);

    EXPECT_TRUE(it.right());
    EXPECT_TRUE(is_set_equal(it.locate(), (std::vector<uint64_t>{1, 4})));
    EXPECT_EQ(it.depth(), 1);

    EXPECT_TRUE(it.down());
    EXPECT_TRUE(is_set_equal(it.locate(), (std::vector<uint64_t>{1, 4})));
    EXPECT_EQ(it.depth(), 2);

    // unsuccessful right(), it remains untouched
    TypeParam it_cpy = it;
    EXPECT_FALSE(it.right());
    EXPECT_EQ(it, it_cpy);

    // unsuccessful down(), it remains untouched
    it = TypeParam(fm);
    EXPECT_TRUE(it.down("GACG"_dna4));
    it_cpy = it;
    EXPECT_FALSE(it.down());
    EXPECT_EQ(it, it_cpy);

    // right() cannot be called on the root node
    it = TypeParam(fm);
    ASSERT_DEATH(it.right(), "");
    EXPECT_EQ(it, TypeParam(fm));
}

template <typename iterator_type, typename index_type>
auto get_all_child_iterators(iterator_type it, index_type fm)
{
    std::array<iterator_type, alphabet_size_v<typename iterator_type::index_type::char_type>> result;

    uint8_t i = 0;
    if (it.down())
    {
        do
        {
            result[i++] = it;
        } while (it.right());
    }

    // fill the rest with iterators pointing to root
    while (i < iterator_type::index_type::char_type::value_size)
    {
        result[i++] = iterator_type(fm);
    }

    return result;
}

TYPED_TEST(fm_index_iterator_test, children)
{
    typename TypeParam::index_type::text_type text{"ACGTAGGT"_dna4};
    typename TypeParam::index_type fm{text};

    // all children
    TypeParam it(fm);
    EXPECT_EQ(it.children(), get_all_child_iterators(it, fm));

    // some children
    it.down(dna4::A);
    EXPECT_EQ(it.children(), get_all_child_iterators(it, fm));

    // one child
    it.down(dna4::G);
    EXPECT_EQ(it.children(), get_all_child_iterators(it, fm));

    // no children
    it.down("GT"_dna4);
    EXPECT_EQ(it.children(), get_all_child_iterators(it, fm));
}

TYPED_TEST(fm_index_iterator_test, path_label)
{
    typename TypeParam::index_type::text_type text{"ACGACG"_dna4};
    typename TypeParam::index_type fm{text};

    // path_label()
    TypeParam it(fm);
    EXPECT_TRUE(it.down("ACG"_dna4));
    EXPECT_TRUE(ranges::equal(*it, "ACG"_dna4));
    EXPECT_TRUE(ranges::equal(it.path_label(), "ACG"_dna4));
}

TYPED_TEST(fm_index_iterator_test, incomplete_alphabet)
{
    // search a char that does not occur in the text (higher rank than largest char occurring in text)
    {
        typename TypeParam::index_type::text_type text{"ACGACG"_dna4};
        typename TypeParam::index_type fm{text};
        TypeParam it = TypeParam(fm);
        EXPECT_FALSE(it.down(dna4::T));
        EXPECT_EQ(it, TypeParam(fm));
    }

    // search a char that does not occur in the text (smaller rank than smallest char occurring in text)
    {
        typename TypeParam::index_type::text_type text{"CGTCGT"_dna4};
        typename TypeParam::index_type fm{text};
        TypeParam it = TypeParam(fm);
        EXPECT_FALSE(it.down(dna4::A));
        EXPECT_EQ(it, TypeParam(fm));
    }

    // search a char that does not occur in the text
    // (some rank that is neither the smallest nor the highest occurring in text)
    {
        typename TypeParam::index_type::text_type text{"ATATAT"_dna4};
        typename TypeParam::index_type fm{text};
        TypeParam it = TypeParam(fm);
        EXPECT_FALSE(it.down(dna4::C));
        EXPECT_FALSE(it.down(dna4::G));
        EXPECT_FALSE(it.down("ACGT"_dna4));
        EXPECT_FALSE(it.down("G"_dna4));
        EXPECT_EQ(it, TypeParam(fm));

        EXPECT_TRUE(it.down(dna4::A));
        EXPECT_TRUE(it.right());
        EXPECT_TRUE(ranges::equal(it.path_label(), "T"_dna4));
    }
}

TYPED_TEST(fm_index_iterator_test, lazy_locate)
{
    typename TypeParam::index_type::text_type text{"ACGTACGT"_dna4};
    typename TypeParam::index_type fm{text};

    TypeParam it = TypeParam(fm);
    it.down("ACG"_dna4);

    EXPECT_TRUE(ranges::equal(it.locate(), it.lazy_locate()));
}

// TODO: tests for other alphabets

TEST(fm_index, concepts)
{
    EXPECT_TRUE(fm_index_iterator_concept<fm_index_iterator<fm_index<std::vector<dna4>>>>);
}
