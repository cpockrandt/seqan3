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
 * \brief Provides the internal representation of a node of the seqan3::fm_index_iterator.
 */

#pragma once

#include <type_traits>

#include <seqan3/index/concept.hpp>
#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{

/*!\addtogroup index
 * \{
 */

//!\privatesection

/*!\brief Internal representation of the node of an FM index iterator.
 * \ingroup fm_index
 * \tparam index_t The type of the underlying index; must satisfy seqan3::fm_index_concept.
 */
template <typename index_t>
struct fm_index_iterator_node
{
    //!\brief Type for representing positions in the indexed text.
    using size_type = typename index_t::size_type;
    /*!\brief The type of the reduced alphabet type. (The reduced alphabet might be smaller than the original alphabet
     *        in case not all possible characters occur in the indexed text.)
     */
    using sdsl_char_type = typename index_t::sdsl_char_type;

    //!\brief Left suffix array bound.
    size_type lb;
    //!\brief Right suffix array bound.
    size_type rb;
    //!\brief Depth of the node in the suffix tree, i.e. length of the searched sequence.
    size_type depth;
    //!\brief Label of the last edge moved down. Needed for right().
    sdsl_char_type last_char;

    //!\brief Comparison of two iterator nodes.
    bool operator==(fm_index_iterator_node const & rhs) const
    {
        // NOTE: last_char is implementation specific for right().
        // lb, rb and depth already determine the node in the suffix tree.
        // Thus there is no need to compare last_char.
        return lb == rhs.lb && rb == rhs.rb && depth == rhs.depth;
    }

    //!\brief Comparison of two iterator nodes.
    bool operator!=(fm_index_iterator_node const & rhs) const
    {
        return !(*this == rhs);
    }
};

//!\publicsection

//!\}

}
