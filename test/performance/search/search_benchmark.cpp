// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <ctime>

#include <benchmark/benchmark.h>

#include <seqan3/range/view/slice.hpp>

#include <seqan3/test/performance/units.hpp>

#include <seqan3/alphabet/aminoacid/aa10murphy.hpp>
#include <seqan3/search/fm_index/all.hpp>
#include <seqan3/io/sequence_file/all.hpp>

using namespace seqan3;
using namespace seqan3::test;

struct sequence_file_input_default_traits_aa10 : sequence_file_input_default_traits_dna
{
     using sequence_alphabet = aa10murphy;
     using sequence_legal_alphabet = aa10murphy;
};

constexpr uint64_t genome_len{1'000'000'000};
constexpr uint64_t read_len{100};
constexpr uint64_t read_count{1'000'000};

std::filesystem::path root{"/home/pocki/epr_bench/"};

// sequence_file_input<sequence_file_input_default_traits_aa10> fin{"/home/pocki/Downloads/murphy10.fa"};
// std::vector<aa10murphy> genome;
// // for (auto & rec : fin)
// // {
//     genome = get<field::SEQ>(*fin.begin());
// //     break;
// // }
//
// std::cout << genome.size() << " genome bases." << std::endl;

static void fm_index_build_data(benchmark::State& /*state*/) {

    std::vector<aa10murphy> genome(genome_len);

    {
        // build and store file
        std::srand(std::time(nullptr));
        for (uint64_t i = 0; i < genome.size(); ++i)
            genome[i].assign_rank(std::rand() % 10);

        std::filesystem::path p1{root / std::filesystem::path{"murphy10.fa"}};
        sequence_file_output fout{p1};
        fout.options.fasta_letters_per_line = 120;
        fout.emplace_back(genome, "seq0");
        std::cout << "Genome built." << std::endl;

        // sample reads
        std::filesystem::path p2{root / "reads.fa"};
        sequence_file_output freads{p2};
        freads.options.fasta_letters_per_line = 0;
        for (uint64_t i = 0; i < read_count; ++i)
        {
            // max startPos is: genome_length - read_len
            uint64_t const startPos{std::rand() % (genome_len - read_len + 1)};
            std::vector<aa10murphy> read{genome | ranges::view::slice(startPos, startPos + read_len - 1)};
            std::string id{"seq" + std::to_string(i)};
            freads.push_back(std::tie(read, id));
        }
        std::cout << "Reads written." << std::endl;
    }

    {
        fm_index wt(genome);

        std::filesystem::path p3{root / "index.wt"};
        std::ofstream os{p3, std::ios::binary};
        cereal::BinaryOutputArchive oarchive{os};
        oarchive(wt);

        std::cout << "WT Index built." << std::endl;
    }

    {
        fm_index<text_layout::single, sdsl_epr_index_type> epr(genome);

        std::filesystem::path p4{root / "index.epr"};
        std::ofstream os{p4, std::ios::binary};
        cereal::BinaryOutputArchive oarchive{os};
        oarchive(epr);

        std::cout << "EPR Index built." << std::endl;
    }

    exit(0);
}

// static void fm_index_wt_benchmark(benchmark::State& state) {
//     sequence_file_input<sequence_file_input_default_traits_aa10> fin{root / "reads.fa"};
//     std::vector<std::vector<aa10murphy>> reads;
//     reads.reserve(read_count);
//
//     for (auto & rec : fin)
//         reads.emplace_back(get<field::SEQ>(rec));
//
//     std::cout << reads.size() << " reads read." << std::endl;
//
//     fm_index wt;
//
//     std::filesystem::path p3{root / "index.wt"};
//     std::ifstream os{p3, std::ios::binary};
//     cereal::BinaryInputArchive iarchive{os};
//     iarchive(wt);
//
//     uint64_t count = 0;
//     for (auto _ : state)
//     {
//         for (uint64_t i = 0; i < reads.size(); ++i)
//         {
//             auto it = wt.begin();
//             it.extend_right(reads[i]);
//             count += it.count();
//         }
//     }
//     std::cout << "Hits: " << count << std::endl;
// }
//
// static void fm_index_epr_benchmark(benchmark::State& state) {
//     sequence_file_input<sequence_file_input_default_traits_aa10> fin{root / "reads.fa"};
//     std::vector<std::vector<aa10murphy>> reads;
//     reads.reserve(read_count);
//
//     for (auto & rec : fin)
//         reads.emplace_back(get<field::SEQ>(rec));
//
//     std::cout << reads.size() << " reads read." << std::endl;
//
//     fm_index<text_layout::single, sdsl_epr_index_type> epr;
//
//     std::filesystem::path p3{root / "index.epr"};
//     std::ifstream os{p3, std::ios::binary};
//     cereal::BinaryInputArchive iarchive{os};
//     iarchive(epr);
//
//     uint64_t count = 0;
//     for (auto _ : state)
//     {
//         for (uint64_t i = 0; i < reads.size(); ++i)
//         {
//             auto it = epr.begin();
//             it.extend_right(reads[i]);
//             count += it.count();
//         }
//     }
//     std::cout << "Hits: " << count << std::endl;
// }

// Register the function as a benchmark
BENCHMARK(fm_index_build_data);
// BENCHMARK(fm_index_wt_benchmark)->Unit(benchmark::kMillisecond);
// BENCHMARK(fm_index_epr_benchmark)->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN();
