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

#include <seqan3/test/comparison.hpp>

#include <range/v3/algorithm/equal.hpp>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/search/fm_index/bi_fm_index.hpp>

#include <gtest/gtest.h>

using namespace seqan3;
using namespace seqan3::literal;
using namespace seqan3::test;

template <typename T>
class bi_fm_index_iterator_test : public ::testing::Test
{};

using bi_fm_index_iterator_types = ::testing::Types<bi_fm_index_iterator<bi_fm_index<std::vector<dna4>>>>;

TYPED_TEST_CASE(bi_fm_index_iterator_test, bi_fm_index_iterator_types);

TYPED_TEST(bi_fm_index_iterator_test, begin)
{
    using text_t = typename TypeParam::index_type::text_type;
    text_t text{"AACGATCGGA"_dna4};
    auto rev_text = ranges::view::reverse(text);
    // using rev_text_t = decltype(ranges::view::reverse(text_t{}));

    using fm_fwd_t = fm_index<text_t>;
    using fm_rev_t = fm_index<decltype(rev_text)>;

    typename TypeParam::index_type bi_fm{text};
    fm_fwd_t fm_fwd{text};
    fm_rev_t fm_rev{rev_text};

    TypeParam bi_it = bi_fm.begin();
    EXPECT_TRUE(is_set_equal(bi_it.locate(), bi_fm.fwd_begin().locate()));
    EXPECT_TRUE(is_set_equal(bi_it.locate(), bi_fm.rev_begin().locate()));
}

TYPED_TEST(bi_fm_index_iterator_test, extend)
{
    typename TypeParam::index_type::text_type text{"ACGGTAGGACG"_dna4};
    typename TypeParam::index_type bi_fm{text};

    auto it = bi_fm.begin();
    EXPECT_TRUE(it.extend_right()); // "A"
    EXPECT_TRUE(is_set_equal(it.locate(), std::vector<uint64_t>{0, 5, 8}));
    EXPECT_TRUE(it.extend_left()); // "GA"
    EXPECT_TRUE(is_set_equal(it.locate(), std::vector<uint64_t>{7}));
    EXPECT_TRUE(it.extend_right()); // "GAC"
    EXPECT_TRUE(is_set_equal(it.locate(), std::vector<uint64_t>{7}));
    EXPECT_TRUE(it.extend_right()); // "GACG"
    EXPECT_TRUE(is_set_equal(it.locate(), std::vector<uint64_t>{7}));
    EXPECT_FALSE(it.extend_right()); // "GACG"
    EXPECT_TRUE(is_set_equal(it.locate(), std::vector<uint64_t>{7}));
    EXPECT_TRUE(it.extend_left()); // "GGACG"
    EXPECT_TRUE(is_set_equal(it.locate(), std::vector<uint64_t>{6}));
}

TYPED_TEST(bi_fm_index_iterator_test, extend_char)
{
    typename TypeParam::index_type::text_type text{"ACGGTAGGACG"_dna4};
    typename TypeParam::index_type bi_fm{text};

    auto it = bi_fm.begin();
    EXPECT_TRUE(it.extend_left(dna4::G)); // "G"
    EXPECT_TRUE(is_set_equal(it.locate(), std::vector<uint64_t>{2, 3, 6, 7, 10}));
    EXPECT_TRUE(it.extend_left(dna4::C)); // "CG"
    EXPECT_TRUE(is_set_equal(it.locate(), std::vector<uint64_t>{1, 9}));
    EXPECT_FALSE(it.extend_left(dna4::C)); // "CG"
    EXPECT_TRUE(is_set_equal(it.locate(), std::vector<uint64_t>{1, 9}));
    EXPECT_FALSE(it.extend_left(dna4::G)); // "CG"
    EXPECT_TRUE(is_set_equal(it.locate(), std::vector<uint64_t>{1, 9}));
    EXPECT_FALSE(it.extend_right(dna4::T)); // "CG"
    EXPECT_TRUE(is_set_equal(it.locate(), std::vector<uint64_t>{1, 9}));
    EXPECT_TRUE(it.extend_right(dna4::G)); // "CGG"
    EXPECT_TRUE(is_set_equal(it.locate(), std::vector<uint64_t>{1}));
    EXPECT_TRUE(it.extend_right(dna4::T)); // "CGGT"
    EXPECT_TRUE(is_set_equal(it.locate(), std::vector<uint64_t>{1}));
    EXPECT_TRUE(it.extend_right(dna4::A)); // "CGGTA"
    EXPECT_TRUE(is_set_equal(it.locate(), std::vector<uint64_t>{1}));
    EXPECT_TRUE(it.extend_left(dna4::A)); // "ACGGTA"
    EXPECT_TRUE(is_set_equal(it.locate(), std::vector<uint64_t>{0}));
    EXPECT_FALSE(it.extend_left(dna4::A)); // "ACGGTA"
    EXPECT_TRUE(is_set_equal(it.locate(), std::vector<uint64_t>{0}));
}

TYPED_TEST(bi_fm_index_iterator_test, extend_range)
{
    typename TypeParam::index_type::text_type text{"ACGGTAGGACG"_dna4};
    typename TypeParam::index_type bi_fm{text};

    auto it = bi_fm.begin();
    EXPECT_FALSE(it.extend_left("CAG"_dna4)); // ""
    // sentinel position included
    EXPECT_TRUE(is_set_equal(it.locate(), std::vector<uint64_t>{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}));
    EXPECT_TRUE(it.extend_left("CG"_dna4)); // "CG"
    EXPECT_TRUE(is_set_equal(it.locate(), std::vector<uint64_t>{1, 9}));
    EXPECT_TRUE(it.extend_right("GTA"_dna4)); // "CGGTA"
    EXPECT_TRUE(is_set_equal(it.locate(), std::vector<uint64_t>{1}));
    EXPECT_FALSE(it.extend_left("TA"_dna4)); // "CGGTA"
    EXPECT_TRUE(is_set_equal(it.locate(), std::vector<uint64_t>{1}));
    EXPECT_TRUE(it.extend_left("A"_dna4)); // "ACGGTA"
    EXPECT_TRUE(is_set_equal(it.locate(), std::vector<uint64_t>{0}));
}

TYPED_TEST(bi_fm_index_iterator_test, extend_and_cycle)
{
    typename TypeParam::index_type::text_type text{"ACGGTAGGACG"_dna4};
    typename TypeParam::index_type bi_fm{text};

    auto it = bi_fm.begin();
    EXPECT_TRUE(it.extend_right()); // "A"
#ifndef NDEBUG
    EXPECT_DEATH(it.cycle_front(), "");
#endif
    EXPECT_TRUE(it.extend_left()); // "GA"
    EXPECT_TRUE(is_set_equal(it.locate(), std::vector<uint64_t>{7}));
#ifndef NDEBUG
    EXPECT_DEATH(it.cycle_back(), "");
#endif
    EXPECT_TRUE(it.cycle_front()); // "TA"
    EXPECT_TRUE(is_set_equal(it.locate(), std::vector<uint64_t>{4}));
    EXPECT_FALSE(it.cycle_front()); // "TA"
    EXPECT_TRUE(is_set_equal(it.locate(), std::vector<uint64_t>{4}));
}

TYPED_TEST(bi_fm_index_iterator_test, extend_range_and_cycle)
{
    typename TypeParam::index_type::text_type text{"ACGGTAGGACGTAG"_dna4};
    typename TypeParam::index_type bi_fm{text};

    auto it = bi_fm.begin();
    EXPECT_TRUE(it.extend_right("AC"_dna4)); // "AC"
    EXPECT_TRUE(is_set_equal(it.locate(), std::vector<uint64_t>{0, 8}));
#ifndef NDEBUG
    EXPECT_DEATH(it.cycle_front(), "");
#endif
    EXPECT_TRUE(it.cycle_back()); // "AG"
    EXPECT_TRUE(is_set_equal(it.locate(), std::vector<uint64_t>{5, 12}));
#ifndef NDEBUG
    EXPECT_DEATH(it.cycle_front(), "");
#endif
    EXPECT_FALSE(it.extend_left("TT"_dna4)); // "AG"
    EXPECT_TRUE(it.extend_left("CGT"_dna4)); // "CGTAG"
    EXPECT_TRUE(is_set_equal(it.locate(), std::vector<uint64_t>{9}));
#ifndef NDEBUG
    EXPECT_DEATH(it.cycle_back(), "");
#endif
    EXPECT_TRUE(it.cycle_front()); // "GGTAG"
    EXPECT_TRUE(is_set_equal(it.locate(), std::vector<uint64_t>{2}));
}

TYPED_TEST(bi_fm_index_iterator_test, to_fwd_iterator)
{
    typename TypeParam::index_type::text_type text{"ACGGTAGGACGTAGC"_dna4};
    typename TypeParam::index_type bi_fm{text};

    {
        auto it = bi_fm.begin();
        EXPECT_TRUE(it.extend_right("GTAGC"_dna4)); // "GTAGC"
        EXPECT_TRUE(is_set_equal(it.locate(), std::vector<uint64_t>{10}));

        auto fwd_it = it.to_fwd_iterator();
        EXPECT_TRUE(fwd_it.cycle_back()); // "GTAGG"
        EXPECT_TRUE(is_set_equal(fwd_it.locate(), std::vector<uint64_t>{3}));
        EXPECT_TRUE(ranges::equal(*fwd_it, "GTAGG"_dna4));
        EXPECT_FALSE(fwd_it.cycle_back());
    }

    {
        auto it = bi_fm.begin();
        EXPECT_TRUE(it.extend_left("GTAG"_dna4)); // "GTAG"
        EXPECT_TRUE(is_set_equal(it.locate(), std::vector<uint64_t>{3, 10}));

        auto fwd_it = it.to_fwd_iterator();
    #ifndef NDEBUG
        EXPECT_DEATH(fwd_it.cycle_back(), "");
    #endif
        EXPECT_TRUE(fwd_it.extend_right());
        EXPECT_TRUE(is_set_equal(fwd_it.locate(), std::vector<uint64_t>{10}));
        EXPECT_TRUE(ranges::equal(*fwd_it, "GTAGC"_dna4));
        EXPECT_TRUE(fwd_it.cycle_back());
        EXPECT_TRUE(is_set_equal(fwd_it.locate(), std::vector<uint64_t>{3}));
        EXPECT_TRUE(ranges::equal(*fwd_it, "GTAGG"_dna4));
    }
}

TYPED_TEST(bi_fm_index_iterator_test, to_rev_iterator)
{
    typename TypeParam::index_type::text_type text{"ACGGTAGGACGTAGC"_dna4};
    typename TypeParam::index_type bi_fm{text};

    {
        auto it = bi_fm.begin();
        EXPECT_TRUE(it.extend_left("CGTAG"_dna4)); // "CGTAG"
        EXPECT_TRUE(is_set_equal(it.locate(), std::vector<uint64_t>{9}));

        auto rev_it = it.to_rev_iterator(); // text "CGATGCAGGATGGCA"
        EXPECT_TRUE(is_set_equal(rev_it.locate(), std::vector<uint64_t>{1}));
        EXPECT_TRUE(ranges::equal(*rev_it, "GATGC"_dna4));
        EXPECT_TRUE(rev_it.cycle_back()); // "GATGG"
        EXPECT_TRUE(is_set_equal(rev_it.locate(), std::vector<uint64_t>{8}));
        EXPECT_TRUE(ranges::equal(*rev_it, "GATGG"_dna4));
        EXPECT_FALSE(rev_it.cycle_back());
    }

    {
        auto it = bi_fm.begin();
        EXPECT_TRUE(it.extend_right("GTAG"_dna4)); // "GTAG"
        EXPECT_TRUE(is_set_equal(it.locate(), std::vector<uint64_t>{3, 10}));

        auto rev_it = it.to_rev_iterator(); // text "CGATGCAGGATGGCA"
    #ifndef NDEBUG
        EXPECT_DEATH(rev_it.cycle_back(), "");
    #endif
        EXPECT_TRUE(rev_it.extend_right()); // "CGTAG" resp. "GATGC"
        EXPECT_TRUE(is_set_equal(rev_it.locate(), std::vector<uint64_t>{1}));
        EXPECT_TRUE(ranges::equal(*rev_it, "GATGC"_dna4));
        EXPECT_TRUE(rev_it.cycle_back()); // "GGTAG" resp. "GATGG"
        EXPECT_TRUE(is_set_equal(rev_it.locate(), std::vector<uint64_t>{8}));
        EXPECT_TRUE(ranges::equal(*rev_it, "GATGG"_dna4));
    }
}
