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

#include <seqan3/search/fm_index/all.hpp>
#include <seqan3/test/tmp_filename.hpp>

#include <gtest/gtest.h>

using namespace seqan3;
using namespace seqan3::literal;

// TODO: EXPECT_EQ is not supported by sdsl

template <typename T>
class fm_index_test : public ::testing::Test
{};

using fm_index_types = ::testing::Types<fm_index<std::vector<dna4>>,
                                        bi_fm_index<std::vector<dna4>>,
                                        bi_fm_index<std::vector<aa27>>,
                                        bi_fm_index<std::vector<char>>>;

TYPED_TEST_CASE(fm_index_test, fm_index_types);

TYPED_TEST(fm_index_test, ctr)
{
    typename TypeParam::text_type text(10); // initialized with smallest char

    // default/zero construction
    TypeParam fm0;
    fm0.construct(text);

    // copy construction
    TypeParam fm1{fm0};
    EXPECT_EQ(fm0.size(), fm1.size());

    // copy assignment
    TypeParam fm2 = fm0;
    EXPECT_EQ(fm0.size(), fm2.size());

    // move construction
    TypeParam fm3{std::move(fm0)};
    EXPECT_EQ(fm0.size(), fm3.size());

    // move assigment
    TypeParam fm4 = std::move(fm0);
    EXPECT_EQ(fm0.size(), fm4.size());

    // container contructor
    TypeParam fm5{text};
    EXPECT_EQ(fm0.size(), fm5.size());
}

TYPED_TEST(fm_index_test, swap)
{
    typename TypeParam::text_type textA(10);
    typename TypeParam::text_type textB(20);

    TypeParam fm0{textA};
    TypeParam fm1{textB};
    TypeParam fm2{fm0};
    TypeParam fm3{fm1};

    EXPECT_EQ(fm0.size(), fm2.size());
    EXPECT_EQ(fm1.size(), fm3.size());
    EXPECT_NE(fm0.size(), fm1.size());

    std::swap(fm1, fm2);

    EXPECT_EQ(fm0.size(), fm1.size());
    EXPECT_EQ(fm2.size(), fm3.size());
    EXPECT_NE(fm0.size(), fm2.size());
}

TYPED_TEST(fm_index_test, size)
{
    TypeParam fm;
    EXPECT_TRUE(fm.empty());

    typename TypeParam::text_type test(8);
    fm.construct(test);
    EXPECT_EQ(fm.size(), 9); // including a sentinel character
}

TYPED_TEST(fm_index_test, serialization)
{
    typename TypeParam::text_type text(8);
    TypeParam fm0{text};

    test::tmp_filename filename{"fm_index"};
    auto const & path = filename.get_path();

    EXPECT_TRUE(fm0.store(path));

    TypeParam fm1{};
    EXPECT_TRUE(fm1.load(path));

    EXPECT_EQ(fm1.size(), 9);
}

TEST(fm_index_test, concepts)
{
    EXPECT_TRUE(fm_index_concept<fm_index<std::vector<dna4>>>);
    EXPECT_TRUE(fm_index_traits_concept<fm_index_default_traits>);

    EXPECT_TRUE(bi_fm_index_concept<bi_fm_index<std::vector<dna4>>>);
    EXPECT_TRUE(bi_fm_index_traits_concept<bi_fm_index_default_traits>);
}
