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
 * \brief Meta-header for the search module.
 *
 * \defgroup search Search
 *
 * ## Introduction
 *
 * Searching is a key component in many sequence analysis tools. The search module is a powerful and easy way to search
 * sequences in a large text or an arbitrary nested collection of texts. When it comes to searching, indices are a core
 * component for searching large amounts of data and are used for tools such as read mappers, assemblers or protein
 * search tools. There are currently two major kind of indices: FM indices and k-mer indices (also known as q-gram
 * indices).
 *
 * \todo Elaborate on that (space consumption for growing k, maybe a rule of thumb).
 *
 * Generally speaking k-mer indices support very fast searching of exact k-mers (strings of length k) or k-mers with
 * predefined wildcard positions that do not have to match. FM indices on the other hand are more versatile and work
 * with arbitrary pattern lengths and error numbers / positions.
 *
 * SeqAn3 currently supports very fast FM indices. For more details visit the \ref submodule_fm_index
 * "FM index submodule".
 *
 * \todo k-mer indices are coming soon. Stay tuned!
 */

#pragma once

#include <seqan3/search/fm_index/all.hpp>
