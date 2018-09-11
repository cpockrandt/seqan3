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

#include <gtest/gtest.h>

using namespace seqan3;
using namespace seqan3::literal;

template <typename T>
class search_configuration_test : public ::testing::Test
{};

// TODO: this should go to a typed configuration test that also checks the alignment configuration
TEST(search_configuration_test, symmetric_configuration)
{
    for (uint8_t i = 0; i < static_cast<uint8_t>(search_cfg::id::SIZE); ++i)
    {
        // no element can occur twice in a configuration
        EXPECT_FALSE(detail::search_config_validation_matrix[i][i])
            << "There is a TRUE value on the diagonal of the search configuration matrix.";
        for (uint8_t j = 0; j < i; ++j)
        {
            // symmetric matrix
            EXPECT_EQ(detail::search_config_validation_matrix[i][j], detail::search_config_validation_matrix[j][i])
                << "Search configuration matrix is not symmetric.";
        }
    }
}

TEST(search_configuration_test, illegal_runtime_configurations)
{
    std::vector<dna4> text{"ACGT"_dna4}, query{"ACG"_dna4};
    fm_index<std::vector<dna4>> fm{text};

    // max_error* without error_type
    search(fm, query, max_error(0) | error_type(error_type_enum::none));
    EXPECT_DEATH(search(fm, query, detail::configuration {max_error(1)}), "");
    EXPECT_DEATH(search(fm, query, detail::configuration {max_error_rate(.1)}), "");
    EXPECT_DEATH(search(fm, query, max_error(1) | error_type(error_type_enum::none)), "");

    // error_type without max_error*
    search(fm, query, detail::configuration {error_type(error_type_enum::none)});
    search(fm, query, error_type(error_type_enum::none) | max_error_rate(.0));
    EXPECT_DEATH(search(fm, query, detail::configuration {error_type(error_type_enum::substitution)}), "");
    EXPECT_DEATH(search(fm, query, error_type(error_type_enum::substitution) | max_error(0)), "");
    EXPECT_DEATH(search(fm, query, error_type(error_type_enum::substitution) | max_error_rate(.0)), "");

    // error_type with max_error*
    search(fm, query, error_type(error_type_enum::substitution) | max_error(1));
    search(fm, query, error_type(error_type_enum::substitution) | max_error_rate(.1));
}
