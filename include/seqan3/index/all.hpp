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
  * \brief Meta-header for the index module.
  *
  * \defgroup index Index
  *
  * ## Introduction
  *
  * Indices are a core component for searching large amounts of data and are used for tools such as read mappers,
  * assemblers or protein search tools. There are currently two major kind of indices: FM indices and k-mer indices
  * (also known as q-gram indices).
  *
  * Besides searching the index by yourself using the iterator interfaces, SeqAn3 also provides a very powerful
  * search module that makes using iterators and implementing your own index-based search algorithms
  * superfluous.
  *
  * ## FM Indices
  *
  * You can choose between unidirectional and bidirectional FM indices (which can be thought of suffix trees
  * and affix trees, i.e. a combination of suffix and prefix trees being able to search a pattern from left to
  * right, right to left and character by character in any arbitrary order). Roughly speaking bidirectional
  * FM indices are more powerful for approximate string matching at the cost of a higher space consumption
  * \todo (between a factor of X and Y depending on the configuration).
  *
  * The FM indices are based on the <a href="https://github.com/xxsds/sdsl-lite">SDSL 3</a> (succinct data structure
  * library). You are able to specify the underlying implementation of the SDSL to adjust it to your needs as well as
  * choose one of the preconfigured indices that are suitable for common applications in sequence analysis.
  *
  * Even though the SDSL supports both byte and integer alphabets, SeqAn3 is optimised for byte alphabets. For
  * integer alphabets you currently cannot use any of the index interfaces of SeqAn3.
  *
  * All FM indices have a suffix-tree-like interface. Even though FM indices are actually prefix trees, they
  * can be searched like a suffix tree for convenience, i.e. there is no need to reverse the text before
  * indexing, the pattern before searching or recomputing text positions afterwards.
  *
  * ## k-mer Indices
  *
  * Coming soon. Stay tuned!
  *
  */

#pragma once

#include <seqan3/index/fm_index.hpp>
