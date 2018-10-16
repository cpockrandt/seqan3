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

#include <cstdint>

#include <algorithm>
#include <type_traits>

#include "helper.hpp"
#include "helper_search_scheme.hpp"

#include <seqan3/io/stream/debug_stream.hpp>
#include <seqan3/search/algorithm/detail/search_scheme_algorithm.hpp>
#include <seqan3/search/algorithm/detail/search_trivial.hpp>
#include <seqan3/search/fm_index/all.hpp>
#include <seqan3/range/view/convert.hpp>

#include <range/v3/view/slice.hpp>

#include <gtest/gtest.h>

using namespace seqan3;

auto trivial_search_scheme(uint8_t const min_error, uint8_t const max_error, uint8_t const blocks)
{
    detail::search_dyn simple_search;
    simple_search.pi.resize(blocks);
    simple_search.l.resize(blocks);
    simple_search.u.resize(blocks);
    for (unsigned i = 0; i < blocks; ++i)
    {
        simple_search.pi[i] = i + 1;
        simple_search.l[i] = min_error;
        simple_search.u[i] = max_error;
    }
    return simple_search;
}

template <uint8_t min_error, uint8_t max_error, bool precomputed_scheme>
void error_distributions(auto & expected, auto & actual)
{
    if constexpr (precomputed_scheme)
    {
        auto const & oss{detail::optimum_search_scheme<min_error, max_error>::value};
        search_scheme_error_distribution(actual, oss);

        auto const simple_search = trivial_search_scheme(min_error, max_error, oss[0].blocks());
        search_error_distribution(expected, simple_search);
    }
    else
    {
        auto const & ss{detail::compute_search_scheme(min_error, max_error)};
        search_scheme_error_distribution(actual, ss);

        auto const simple_search = trivial_search_scheme(min_error, max_error, ss[0].blocks());
        search_error_distribution(expected, simple_search);
    }

    std::sort(expected.begin(), expected.end());
    std::sort(actual.begin(), actual.end());

    // debug_stream << "------------------\n";
    // debug_stream << expected << '\n';
    // debug_stream << actual << '\n';
}

TEST(search_scheme_test, error_distribution_coverage_optimum_search_schemes)
{
    std::vector<std::vector<uint8_t> > expected, actual;

    error_distributions<0, 0, true>(expected, actual);
    EXPECT_EQ(actual, expected);

    error_distributions<0, 1, true>(expected, actual);
    EXPECT_EQ(actual, expected);

    error_distributions<0, 2, true>(expected, actual);
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

    auto const & oss{detail::optimum_search_scheme<min_error, max_error>::value};
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

TEST(search_scheme_test, error_distribution_disjoint_computed_search_schemes)
{

}














inline void test_search(auto it,
                        auto const & text,
                        auto const & search,
                        uint64_t const needle_length,
                        std::vector<uint8_t> const & error_distribution,
                        // time_t const seed,
                        auto const & blocklength, auto const & ordered_blocklength, uint64_t const start_pos)
{
    using char_t = dna4;

    uint64_t const pos = std::rand() % (text.size() - needle_length + 1);
    std::vector<dna4> const orig_needle = text | ranges::view::slice(pos, pos + needle_length);

    // Modify needle s.t. it has errors matching error_distribution.
    auto needle = orig_needle;
    uint64_t cumulativeBlocklength = 0;
    for (uint8_t block = 0; block < search.blocks(); ++block)
    {
        uint64_t single_block_length = ordered_blocklength[block];
        if (error_distribution[block] > single_block_length)
        {
            debug_stream << "Error in block " << block << "(+ 1): " << error_distribution[block]
                         << " errors cannot fit into a block of length " << single_block_length << "." << '\n'
                         << "Error Distribution: " << error_distribution << '\n';
            print_search(search, blocklength);
            // print_search(orderd_search);
            exit(1);
        }

        // choose random positions in needle that will be a mismatch/indel
        // repeat until all error positions are unique
        std::vector<uint8_t> error_positions(error_distribution[block]);
        do
        {
            error_positions.clear();
            for (uint8_t error = 0; error < error_distribution[block]; ++error)
                error_positions.push_back(std::rand() % single_block_length);
            sort(error_positions.begin(), error_positions.end());
        } while (std::adjacent_find(error_positions.begin(), error_positions.end()) != error_positions.end());

        // construct needle with chosen error positions
        for (unsigned error = 0; error < error_positions.size(); ++error)
        {
            unsigned pos = error_positions[error] + cumulativeBlocklength;
            char_t newChar;
            do
            {
                assign_rank(newChar, std::rand() % 4);
            } while(needle[pos] == newChar);
            needle[pos] = newChar;
        }
        cumulativeBlocklength += single_block_length;
    }

    uint8_t max_errors = search.u.back();

    std::vector<uint64_t> hits, expectedHitsSS, expectedHitsTrivial;
    auto delegate = [&hits](auto const & it)
    {
        auto const & hits_tmp = it.locate();
        hits.insert(hits.end(), hits_tmp.begin(), hits_tmp.end());
    };

    using namespace seqan3::literal;

    // Find all hits using search schemes.
    /*if (text == "GCCCATTAAG"_dna4 && needle == "GCCCA"_dna4)
    {
        std::cout << "yoooo\n";
    }*/

    detail::search_params error_left{max_errors, max_errors, max_errors, max_errors};

    detail::search_scheme_single_search<false>(
        it, needle,     // data
        start_pos, start_pos + 1, // infix range
        0,                        // errors spent
        0,                        // current block id
        true,                     // go to the right
        search, blocklength,      // search scheme information
        error_left,
        delegate);

    for (uint64_t hit : hits)
    {
        dna4_vector matched_seq = text | ranges::view::slice(hit, hit + needle_length);
        if (matched_seq == orig_needle)
            expectedHitsSS.push_back(hit);
    }

    // Find all hits using trivial backtracking.
    hits.clear();
    detail::_search_trivial<false>(it, needle, 0, error_left, delegate);

    for (uint64_t hit : hits)
    {
        // filter only correct error distributions
        dna4_vector matched_seq = text | ranges::view::slice(hit, hit + needle_length);
        if (orig_needle == matched_seq)
        {
            bool distributionOkay = true;
            uint64_t leftRange = 0, rightRange = 0;
            for (uint8_t block = 0; block < search.blocks(); ++block)
            {
                rightRange += ordered_blocklength[block];

                uint8_t errors = 0;
                for (unsigned i = leftRange; i < rightRange; ++i)
                    if (hit + i >= text.size())
                        ++errors;
                    else
                        errors += needle[i] != text[hit + i];
                if (errors != error_distribution[block])
                    distributionOkay = false;
                leftRange += ordered_blocklength[block];
            }
            if (distributionOkay || (error_left.insertion > 0 || error_left.deletion > 0))
                expectedHitsTrivial.push_back(hit);
        }
    }

    // eliminate duplicates
    std::sort(expectedHitsSS.begin(), expectedHitsSS.end());
    std::sort(expectedHitsTrivial.begin(), expectedHitsTrivial.end());
    expectedHitsSS.erase(std::unique(expectedHitsSS.begin(), expectedHitsSS.end()), expectedHitsSS.end());
    expectedHitsTrivial.erase(std::unique(expectedHitsTrivial.begin(), expectedHitsTrivial.end()), expectedHitsTrivial.end());

    if (expectedHitsSS != expectedHitsTrivial)
    {
        debug_stream //<< "Seed: " << seed << '\n'
                     << "Text: " << text << '\n'
                     << "Error_distribution: " << error_distribution << '\n'
                     << "Original: " << orig_needle << '\n'
                     << "Modified: " << needle << '\n'
                     << "ExpectedHitsSS: " << expectedHitsSS << '\n'
                     << "ExpectedHitsTrivial: " << expectedHitsTrivial << '\n';
        print_search(search, blocklength);
        exit(1);
    }
}

template <typename search_scheme_t>
inline void test_search_scheme(search_scheme_t const & search_scheme)
{
    using index_t = bi_fm_index<dna4_vector>;

    search_scheme_t ordered_search_scheme;
    std::vector<std::vector<std::vector<uint8_t> > > error_distributions(search_scheme.size());

    // Calculate all error distributions and sort each of them (from left to right).
    uint8_t errors = 0;
    for (uint8_t search_id = 0; search_id < search_scheme.size(); ++search_id)
    {
        ordered_search_scheme[search_id] = search_scheme[search_id];
        search_error_distribution(error_distributions[search_id], search_scheme[search_id]);
        for (std::vector<uint8_t> & resElem : error_distributions[search_id])
            order_search_vector(resElem, search_scheme[search_id]);
        errors = std::max(errors, search_scheme[search_id].u.back());
    }

    for (uint64_t text_length = 10; text_length < 10000; text_length *= 10)
    {
        std::vector<dna4> text;
        random_text(text, text_length);
        //debug_stream << text << '\n';
        index_t index(text);
        for (uint64_t needle_length = std::max<uint64_t>(5, search_scheme[0].blocks() * errors); needle_length < std::min<uint64_t>(16, text_length); ++needle_length)
        {
            auto const block_info = search_scheme_block_info(search_scheme, needle_length);
            for (uint8_t search_id = 0; search_id < search_scheme.size(); ++search_id)
            {
                /*using namespace seqan3::literal;
                if (needle_length == 5 && search_id == 1 && text == "TCAGCCAAGCGGATCG"_dna4)
                {
                    debug_stream << needle_length << '\n';
                    debug_stream << search_id << '\n';
                }*/
                auto const & [blocklength, start_pos] = block_info[search_id];

                std::vector<uint64_t> ordered_blocklength;
                get_ordered_search(search_scheme[search_id], blocklength, ordered_search_scheme[search_id], ordered_blocklength);

                for (auto const & error_distribution : error_distributions[search_id])
                {
                    //std::cout << "-----------------------------------\n";
                    test_search(index.begin(), text, search_scheme[search_id], needle_length, error_distribution, /*seed,*/ blocklength, ordered_blocklength, start_pos);
                }
            }
        }
    }
}

TEST(search_scheme_test, search_scheme_mismatches)
{
    time_t seed = std::time(nullptr);
    std::srand(seed);
    std::cout << "seed = " << seed << std::endl;

    for (uint64_t i = 0; i < 1000; ++i)
    {
        test_search_scheme(detail::optimum_search_scheme<0, 0>::value);
        test_search_scheme(detail::optimum_search_scheme<0, 1>::value);
        test_search_scheme(detail::optimum_search_scheme<1, 1>::value);
        test_search_scheme(detail::optimum_search_scheme<0, 2>::value);
        test_search_scheme(detail::optimum_search_scheme<0, 3>::value);
        if (i > 0 && i % 20 == 0)
            std::cout << std::endl;
        std::cout << ".";
        std::cout.flush();
    }
}
