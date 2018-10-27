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
 * \brief Provides an approximate string matching algorithm based on simple backtracking.
 *        This should only be used as a reference for unit testing.
 */

#pragma once

#include <type_traits>

#include <range/v3/view/drop_exactly.hpp>

#include <seqan3/range/concept.hpp>
#include <seqan3/search/algorithm/detail/search_common.hpp>

namespace seqan3::detail
{

/*!\brief Searches a query sequence in an index using trivial backtracking.
 * \tparam abort_on_hit If the flag is set, the search algorithm aborts on the first hit.
 * \param[in] it Iterator of atring index built on the text that will be searched.
 * \param[in] query Query sequence to be searched with the iterator.
 * \param[in] query_pos Position in the query sequence indicating the prefix that has already been searched.
 * \param[in] error_left Number of errors left for matching the remaining suffix of the query sequence.
 * \param[in] delegate Function that is called on every hit. Takes `index::iterator_type` as argument.
 */
template <bool abort_on_hit>
inline bool search_trivial(auto it, auto & query, uint64_t const query_pos, search_param const error_left,
                           auto && delegate)
{
    // Exact case (end of query sequence or no errors left)
    if (query_pos == query.size() || error_left.total == 0)
    {
        // If not at end of query sequence, try searching the remaining suffix without any errors.
        if (query_pos == query.size() || it.extend_right(ranges::view::drop_exactly(query, query_pos)))
        {
            // TODO: pre-merge: only call delegate if the search was not finished with a deletion
            delegate(it);
            return true;
        }
    }
    // Approximate case
    else
    {
        // Insertion
        if (error_left.insertion > 0)
        {
            search_param error_left2{error_left};
            error_left2.insertion--;
            error_left2.total--;

            // always perform a recursive call. Abort recursion if and only if recursive call found a hit and
            // abort_on_hit is set to true.
            if (search_trivial<abort_on_hit>(it, query, query_pos + 1, error_left2, delegate) && abort_on_hit)
                return true;
        }

        // do not allow deletions at the beginning of the query sequence
        if (((query_pos > 0 && error_left.deletion > 0) || error_left.substitution > 0) && it.extend_right())
        {
            do
            {
                // Match (when error_left.substitution > 0) and Mismatch
                if (error_left.substitution > 0)
                {
                    bool delta = it.last_char() != query[query_pos];
                    search_param error_left2{error_left};
                    error_left2.total -= delta;
                    error_left2.substitution -= delta;

                    if (search_trivial<abort_on_hit>(it, query, query_pos + 1, error_left2, delegate) && abort_on_hit)
                        return true;
                }

                // Deletion
                // do not allow deletions at the beginning of the query sequence
                if (query_pos > 0)
                {
                    // Match (when error_left.substitution == 0)
                    if (error_left.substitution == 0 && it.last_char() == query[query_pos])
                    {
                        if (search_trivial<abort_on_hit>(it, query, query_pos + 1, error_left, delegate) &&
                            abort_on_hit)
                        {
                            return true;
                        }
                    }

                    if (error_left.deletion > 0)
                    {
                        search_param error_left2{error_left};
                        error_left2.total--;
                        error_left2.deletion--;

                        if (search_trivial<abort_on_hit>(it, query, query_pos, error_left2, delegate) && abort_on_hit)
                            return true;
                    }
                }
            } while (it.cycle_back());
        }
        else
        {
            // Match (when error_left.substitution == 0)
            if (it.extend_right(query[query_pos]))
            {
                if (search_trivial<abort_on_hit>(it, query, query_pos + 1, error_left, delegate) && abort_on_hit)
                    return true;
            }
        }
    }

    return false;
}

/*!\brief Searches a query sequence in an index using trivial backtracking.
* \tparam abort_on_hit If the flag is set, the search algorithm aborts on the first hit.
* \param[in] index String index built on the text that will be searched.
* \param[in] query Query sequence to be searched in the index.
* \param[in] error_left Number of errors left for matching the remaining suffix of the query sequence.
* \param[in] delegate Function that is called on every hit. Takes `index::iterator_type` as argument.
*
* ### Complexity
*
* Exponential in the numbers of errors.
*
* ### Exceptions
*
* No-throw guarantee.
*/
template <bool abort_on_hit>
inline void search_trivial(auto const & index, auto & query, search_param const error_left, auto && delegate)
{
    search_trivial<abort_on_hit>(index.begin(), query, 0, error_left, delegate);
}

}
