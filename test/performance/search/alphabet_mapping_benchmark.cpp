// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
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
// ==========================================================================

#include <algorithm>
#include <cmath>
#include <cstring>

#include <benchmark/benchmark.h>

#include <seqan3/index/fm_index.hpp>

using namespace seqan3;
using namespace seqan3::literal;

template <typename csa_alphabet_stragy_t>
struct fm_index_mapping_traits
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
        csa_alphabet_stragy_t
    >;
};

// ============================================================================
//  alphabet_mapping_bench
// ============================================================================

template <bool complete_alphabet, typename char_type>
inline void create_random_string(std::vector<char_type> & str, const uint64_t length, std::mt19937_64 & rng)
{
    if constexpr(complete_alphabet)
    {
        for (uint64_t i = 0; i < length; ++i)
            str[i].assign_rank(rng() % alphabet_size_v<char_type>);
    }
    else
    {
        // create a subset of characters that is the same for each call of this function
        std::mt19937_64 rng2(42);
        uint64_t subset_chars_size = alphabet_size_v<char_type> / 2;
        std::vector<uint8_t> subset_chars{};
        while (subset_chars.size() != subset_chars_size)
        {
            auto const rank = rng2();
            if (std::find(subset_chars.begin(), subset_chars.end(), rank) == subset_chars.end())
                subset_chars.push_back(rank);
        }

        for (uint64_t i = 0; i < length; ++i)
        {
            auto const rank = rng() % subset_chars_size;
            str[i].assign_rank(subset_chars[rank]);
        }
    }
}

template <typename char_type, typename alphabet_mapping_t, bool complete_alphabet>
static void alphabet_mapping_bench(benchmark::State& state)
{
    using traits = fm_index_mapping_traits<alphabet_mapping_t>;

    uint64_t const query_length = 10;
    uint64_t const text_length = 1000000;

	std::mt19937_64 rng;
	rng.seed(42);

    std::vector<char_type> text(text_length);
    create_random_string<complete_alphabet>(text, text_length, rng);

    fm_index<std::vector<char_type>, traits> sa{text};

    std::vector<std::vector<char_type>> queries{};
    for (unsigned i = 0; i < 10000; ++i)
    {
        std::vector<char_type> query(query_length);
        create_random_string<complete_alphabet>(query, query_length, rng);
        queries.push_back(query);
    }

    for (auto _ : state)
    {
        [[maybe_unused]] volatile uint64_t count_total = 0;
        for (auto it = queries.cbegin(); it != queries.cend(); ++it)
        {
            auto sa_it = sa.root();
            sa_it.down(*it);
            count_total += sa_it.count();
        }
    }
}

// BENCHMARK_TEMPLATE(alphabet_mapping_bench, dna4, sdsl::byte_alphabet, false);
// // BENCHMARK_TEMPLATE(alphabet_mapping_bench, dna4, sdsl::succinct_byte_alphabet<sdsl::bit_vector, sdsl::bit_vector::rank_1_type, sdsl::bit_vector::select_1_type, sdsl::int_vector<> >, false);
// BENCHMARK_TEMPLATE(alphabet_mapping_bench, dna4, sdsl::plain_byte_alphabet, false);

BENCHMARK_TEMPLATE(alphabet_mapping_bench, dna4, sdsl::byte_alphabet, true);
// BENCHMARK_TEMPLATE(alphabet_mapping_bench, dna4, sdsl::succinct_byte_alphabet<sdsl::bit_vector, sdsl::bit_vector::rank_1_type, sdsl::bit_vector::select_1_type, sdsl::int_vector<> >, true);
BENCHMARK_TEMPLATE(alphabet_mapping_bench, dna4, sdsl::plain_byte_alphabet, true);

// BENCHMARK_TEMPLATE(alphabet_mapping_bench, dna15, sdsl::byte_alphabet, false);
// // BENCHMARK_TEMPLATE(alphabet_mapping_bench, dna15, sdsl::succinct_byte_alphabet<sdsl::bit_vector, sdsl::bit_vector::rank_1_type, sdsl::bit_vector::select_1_type, sdsl::int_vector<> >, false);
// BENCHMARK_TEMPLATE(alphabet_mapping_bench, dna15, sdsl::plain_byte_alphabet, false);

BENCHMARK_TEMPLATE(alphabet_mapping_bench, dna15, sdsl::byte_alphabet, true);
// BENCHMARK_TEMPLATE(alphabet_mapping_bench, dna15, sdsl::succinct_byte_alphabet<sdsl::bit_vector, sdsl::bit_vector::rank_1_type, sdsl::bit_vector::select_1_type, sdsl::int_vector<> >, true);
BENCHMARK_TEMPLATE(alphabet_mapping_bench, dna15, sdsl::plain_byte_alphabet, true);

// ============================================================================
//  run
// ============================================================================

BENCHMARK_MAIN();
