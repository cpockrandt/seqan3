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

#include <seqan3/search/all.hpp>
#include <seqan3/io/stream/debug_stream.hpp>

using namespace seqan3;
using namespace seqan3::literal;

// ============================================================================
//  alphabet_mapping_bench
// ============================================================================

template <typename char_t>
inline void create_random_string(std::vector<char_t> & str, const uint64_t length, std::mt19937_64 & rng)
{
    for (uint64_t i = 0; i < length; ++i)
        str[i].assign_rank(rng() % alphabet_size_v<char_t>);
}

// uint64_t const query_length = 10;
// uint64_t const text_length = 10000000;
//
// std::mt19937_64 rng;
// rng.seed(42);
//
// std::vector<dna4> text(text_length);
// create_random_string(text, text_length, rng);
//
// fm_index sa{text};
//
// std::vector<std::vector<dna4>> queries{};
// for (unsigned i = 0; i < 1000; ++i)
// {
//     std::vector<dna4> query(query_length);
//     create_random_string(query, query_length, rng);
//     queries.push_back(query);
// }

uint64_t const query_length = 10;
uint64_t const text_length = 10000000;

fm_index<std::vector<dna4>> sa;
std::vector<dna4> text;
std::vector<std::vector<dna4>> queries;

static void init(benchmark::State& /*state*/)
{
	std::mt19937_64 rng;
	rng.seed(42);

    text.resize(text_length);
    create_random_string(text, text_length, rng);
    sa.construct(text);

    for (unsigned i = 0; i < 10; ++i)
    {
        std::vector<dna4> query(query_length);
        create_random_string(query, query_length, rng);
        queries.push_back(query);
    }
}

static void bench1(benchmark::State& state)
{
    detail::configuration cfg = max_total_error(2)
                              | max_substitution_error(2);

    for (auto _ : state)
    {
        auto res = search(sa, queries, cfg);
        [[maybe_unused]] volatile size_t fin = res.size();
    }
}

static void bench2(benchmark::State& state)
{
    detail::configuration cfg = max_total_error(2)
                              | max_substitution_error(2)
                              | max_insertion_error(0)
                              | max_deletion_error(0)
                              ;

    for (auto _ : state)
    {
        auto res = search(sa, queries, cfg);
        [[maybe_unused]] volatile size_t fin = res.size();
    }
}

BENCHMARK(init);
BENCHMARK(bench1);
BENCHMARK(bench2);

// ============================================================================
//  run
// ============================================================================

BENCHMARK_MAIN();
