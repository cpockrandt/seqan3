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

 /*!\file
  * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
  * \brief Meta-header for the FM index module.
  *
  * \defgroup submodule_fm_index FM Index
  * \ingroup search
  *
  * ## FM Indices
  *
  * FM indices are text indices similar to suffix trees or suffix arrays which are based on the Burrow Wheeler
  * transform and a sampled suffix array. FM indices are significantly smaller in space without sacrificing speed. They
  * also allow for adjusting the speed respectively space by choosing the underlying data structures accordingly or
  * modyfing the sampling rate of the sampled suffix array.
  *
  * The FM indices are based on the <a href="https://github.com/xxsds/sdsl-lite">SDSL 3</a> (succinct data structure
  * library). You are able to specify the underlying implementation of the SDSL to adjust it to your needs as well as
  * choose one of the preconfigured indices that are suitable for common applications in sequence analysis.
  *
  * For technical reasons you can currently only build indices over a seqan3::alphabet_concept if its
  * seqan3::alphabet_size is smaller or equal 256.
  *
  * You can choose between unidirectional and bidirectional FM indices (which can be thought of suffix trees
  * and affix trees, i.e. a combination of suffix and prefix trees being able to search a pattern from left to
  * right, right to left and character by character in any arbitrary order). Roughly speaking bidirectional
  * FM indices are more powerful for approximate string matching at the cost of a higher space consumption
  * \todo (between a factor of X and Y depending on the configuration).
  *
  * ## FM Index Iterators
  *
  * Index Iterators are lightweight objects, i.e. they are cheap to copy.
  *
  * Note that although the SeqAn3 index iterators are called "iterators", they don't model any of the standard library
  * iterator concepts, not even std::Iterator.
  *
  */

#pragma once

#include <seqan3/search/fm_index/bi_fm_index.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
