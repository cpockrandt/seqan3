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
#include "helper_search_scheme.hpp"

#include <seqan3/io/stream/debug_stream.hpp>
#include <seqan3/search/algorithm/detail/search_scheme_algorithm.hpp>
#include <seqan3/search/algorithm/detail/search_trivial.hpp>
#include <seqan3/search/fm_index/all.hpp>

#include <range/v3/view/slice.hpp>

#include <gtest/gtest.h>

#ifdef NDEBUG
#define ITERATIONS 1000
#else
#define ITERATIONS 10
#endif

using namespace seqan3;

template <typename text_t>
inline void test_search_hamming(auto it, text_t const & text, auto const & search, uint64_t const query_length,
                                std::vector<uint8_t> const & error_distribution, time_t const seed,
                                auto const & blocks_length, auto const & ordered_blocks_length,
                                uint64_t const start_pos)
{
    using char_t = typename text_t::value_type;

    uint64_t const pos = std::rand() % (text.size() - query_length + 1);
    text_t const orig_query = text | ranges::view::slice(pos, pos + query_length);

    // Modify query s.t. it has errors matching error_distribution.
    auto query = orig_query;
    uint64_t current_blocks_length = 0;
    for (uint8_t block = 0; block < search.blocks(); ++block)
    {
        uint64_t const single_block_length = ordered_blocks_length[block];
        EXPECT_LE(error_distribution[block], single_block_length);
        if (error_distribution[block] > single_block_length)
        {
            debug_stream << "Error in block " << block << "(+ 1): " << error_distribution[block]
                         << " errors cannot fit into a block of length " << single_block_length << "." << '\n'
                         << "Error Distribution: " << error_distribution << '\n';
            exit(1);
        }

        // Choose random positions in the query sequence for substitutions. Repeat until all error positions are unique.
        std::vector<uint8_t> error_positions(error_distribution[block]);
        do
        {
            error_positions.clear();
            for (uint8_t error = 0; error < error_distribution[block]; ++error)
                error_positions.push_back(std::rand() % single_block_length);
            std::sort(error_positions.begin(), error_positions.end());
        } while (std::adjacent_find(error_positions.begin(), error_positions.end()) != error_positions.end());

        // Construct query sequence with chosen error positions.
        for (uint8_t error = 0; error < error_positions.size(); ++error)
        {
            uint64_t const pos = error_positions[error] + current_blocks_length;
            // Decrease alphabet size by one because we don't want to replace query[pos], with the same character.
            uint8_t new_rank = std::rand() % (alphabet_size_v<char_t> - 1);
            // If it is a match now, it can't be the highest rank of the alphabet. Thus we set it to the highest rank.
            if (new_rank == to_rank(query[pos]))
                new_rank = alphabet_size_v<char_t> - 1;
            assign_rank(query[pos], new_rank);
        }
        current_blocks_length += single_block_length;
    }

    std::vector<uint64_t> hits_trivial, hits_ss;

    auto delegate_trivial = [&hits_trivial](auto const & it)
    {
        auto const & hits_tmp = it.locate();
        hits_trivial.insert(hits_trivial.end(), hits_tmp.begin(), hits_tmp.end());
    };

    auto delegate_ss = [&hits_ss](auto const & it)
    {
        auto const & hits_tmp = it.locate();
        hits_ss.insert(hits_ss.end(), hits_tmp.begin(), hits_tmp.end());
    };

    uint8_t const total        = search.u.back();
    uint8_t const substitution = std::rand() % (total + 1);

    detail::search_param error_left{total, substitution, 0, 0};

    // Find all hits using search schemes.
    detail::search_ss<false>(it, query, start_pos, start_pos + 1, 0, 0, true, search, blocks_length, error_left,
                             delegate_ss);
    // Find all hits using trivial backtracking.
    detail::search_trivial<false>(it, query, 0, error_left, delegate_trivial);

    // Eliminate hits that we are not interested in (based on the search and chosen error distribution)
    hits_ss.erase(std::remove_if(hits_ss.begin(), hits_ss.end(),
                                 [&text, &orig_query, query_length](uint64_t const hit)
                                 {
                                     dna4_vector matched_seq = text | ranges::view::slice(hit, hit + query_length);
                                     return (matched_seq != orig_query);
                                 }), hits_ss.end());
    hits_trivial.erase(std::remove_if(hits_trivial.begin(), hits_trivial.end(),
                                      [&](uint64_t const hit)
                                      {
                                          // filter only correct error distributions
                                          dna4_vector matched_seq = text | ranges::view::slice(hit, hit + query_length);
                                          if (orig_query != matched_seq)
                                              return true;

                                          uint64_t lb = 0, rb = 0;
                                          uint8_t total_errors = 0;
                                          for (uint8_t block = 0; block < search.blocks(); ++block)
                                          {
                                              rb += ordered_blocks_length[block];

                                              uint8_t errors = 0;
                                              for (uint64_t i = lb; i < rb; ++i)
                                                  if (hit + i >= text.size())
                                                      ++errors;
                                                  else
                                                      errors += query[i] != text[hit + i];
                                              total_errors += errors;
                                              if (errors != error_distribution[block])
                                                  return true;
                                              lb += ordered_blocks_length[block];
                                          }
                                          return false;
                                      }), hits_trivial.end());

    // Eliminate duplicates
    hits_ss = sort(hits_ss);
    hits_trivial = sort(hits_trivial);

    EXPECT_EQ(hits_ss, hits_trivial);
    if (hits_ss != hits_trivial)
    {
        debug_stream << "Seed: " << seed << '\n'
                     << "Text: " << text << '\n'
                     << "Query: " << query << '\n'
                     << "Errors: " << total << ", " << substitution << '\n'
                     << "SS hits: " << hits_ss << '\n'
                     << "Trivial hits: " << hits_trivial << '\n';
    }
}

template <typename search_scheme_t>
inline void test_search_scheme_hamming(search_scheme_t const & search_scheme, time_t const seed,
                                       uint64_t const iterations)
{
    using text_t = dna4_vector;

    text_t text;

    search_scheme_t ordered_search_scheme;
    std::vector<std::vector<std::vector<uint8_t> > > error_distributions(search_scheme.size());

    // Calculate all error distributions and sort each of them (from left to right).
    uint8_t max_error = 0;
    for (uint8_t search_id = 0; search_id < search_scheme.size(); ++search_id)
    {
        ordered_search_scheme[search_id] = search_scheme[search_id];
        search_error_distribution(error_distributions[search_id], search_scheme[search_id]);
        for (std::vector<uint8_t> & resElem : error_distributions[search_id])
            order_search_vector(resElem, search_scheme[search_id]);
        max_error = std::max(max_error, search_scheme[search_id].u.back());
    }

    for (uint64_t text_length = 10; text_length < 10000; text_length *= 10)
    {
        uint64_t const query_length_min = std::max<uint64_t>(3, search_scheme.front().blocks() * max_error);
        uint64_t const query_length_max = std::min<uint64_t>(16, text_length);

        random_text(text, text_length);
        bi_fm_index index(text);

        for (uint64_t i = 0; i < iterations; ++i)
        {
            for (uint64_t query_length = query_length_min; query_length < query_length_max; ++query_length)
            {
                auto const block_info = search_scheme_block_info(search_scheme, query_length);
                for (uint8_t search_id = 0; search_id < search_scheme.size(); ++search_id)
                {
                    auto const & [blocks_length, start_pos] = block_info[search_id];

                    std::vector<uint64_t> ordered_blocks_length;
                    get_ordered_search(search_scheme[search_id], blocks_length,
                                       ordered_search_scheme[search_id], ordered_blocks_length);

                    for (auto && error_distribution : error_distributions[search_id])
                    {
                        test_search_hamming(index.begin(), text, search_scheme[search_id], query_length,
                                            error_distribution, seed, blocks_length, ordered_blocks_length, start_pos);
                    }
                }
            }
        }
    }
}

template <typename search_scheme_t>
inline void test_search_scheme_edit(search_scheme_t const & search_scheme, time_t const seed, uint64_t const iterations)
{
    using text_t = dna4_vector;

    text_t text, query;

    // retrieve maximum number of errors from search_scheme
    uint8_t max_error = 0;
    for (auto const & search : search_scheme)
        max_error = std::max(max_error, search.u.back());

    for (uint64_t text_length = 10; text_length < 10000; text_length *= 10)
    {
        uint64_t const query_length_min = std::max<uint64_t>(3, search_scheme.front().blocks() * max_error);
        uint64_t const query_length_max = std::min<uint64_t>(16, text_length);

        random_text(text, text_length);
        bi_fm_index index(text);

        uint8_t const substitution = std::rand() % (max_error + 1);
        uint8_t const insertion    = std::rand() % (max_error + 1);
        uint8_t const deletion     = std::rand() % (max_error + 1);
        detail::search_param error_left{max_error, substitution, insertion, deletion};

        for (uint64_t i = 0; i < iterations; ++i)
        {
            for (uint64_t query_length = query_length_min; query_length < query_length_max; ++query_length)
            {
                random_text(query, query_length);

                std::vector<uint64_t> hits_trivial, hits_ss;

                auto delegate_trivial = [&hits_trivial](auto const & it)
                {
                    auto const & hits_tmp = it.locate();
                    hits_trivial.insert(hits_trivial.end(), hits_tmp.begin(), hits_tmp.end());
                };

                auto delegate_ss = [&hits_ss](auto const & it)
                {
                    auto const & hits_tmp = it.locate();
                    hits_ss.insert(hits_ss.end(), hits_tmp.begin(), hits_tmp.end());
                };

                // Find all hits using search schemes.
                detail::search_ss<false>(index, query, error_left, search_scheme, delegate_ss);
                // Find all hits using trivial backtracking.
                detail::search_trivial<false>(index, query, error_left, delegate_trivial);

                // Eliminate duplicates
                hits_ss = sort(hits_ss);
                hits_trivial = sort(hits_trivial);

                EXPECT_EQ(hits_ss, hits_trivial);
                if (hits_ss != hits_trivial)
                {
                    debug_stream << "Seed: " << seed << '\n'
                                 << "Text: " << text << '\n'
                                 << "Query: " << query << '\n'
                                 << "Errors: " << max_error << ", " << substitution << ", "
                                               << insertion << ", " << deletion << '\n'
                                 << "SS hits: " << hits_ss << '\n'
                                 << "Trivial hits: " << hits_trivial << '\n';
                }
            }
        }
    }
}

TEST(search_scheme_test, search_scheme_hamming)
{
    time_t seed = std::time(nullptr);
    std::srand(seed);

    test_search_scheme_hamming(detail::optimum_search_scheme<0, 0>, seed, ITERATIONS);
    test_search_scheme_hamming(detail::optimum_search_scheme<0, 1>, seed, ITERATIONS);
    test_search_scheme_hamming(detail::optimum_search_scheme<1, 1>, seed, ITERATIONS);
    test_search_scheme_hamming(detail::optimum_search_scheme<0, 2>, seed, ITERATIONS);
    test_search_scheme_hamming(detail::optimum_search_scheme<0, 3>, seed, ITERATIONS);
}

TEST(search_scheme_test, search_scheme_edit)
{
    time_t seed = std::time(nullptr);
    std::srand(seed);

    test_search_scheme_edit(detail::optimum_search_scheme<0, 0>, seed, ITERATIONS);
    test_search_scheme_edit(detail::optimum_search_scheme<0, 1>, seed, ITERATIONS);
    test_search_scheme_edit(detail::optimum_search_scheme<0, 2>, seed, ITERATIONS);
    test_search_scheme_edit(detail::optimum_search_scheme<0, 3>, seed, ITERATIONS);
}
