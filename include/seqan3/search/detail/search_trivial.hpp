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

// NOTE: this file will not be documented since it is only a prototype and will be replaced by search schemes later.

/*!\file
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 * \brief
 */

#pragma once

#include <type_traits>

#include <range/v3/view/drop_exactly.hpp>

#include <seqan3/range/concept.hpp>

//!\cond

namespace seqan3::detail
{
    template <bool substitution, bool insertion, bool deletion, bool allow_insertion, bool allow_deletion,
              bool abort_on_hit>
    auto _search_trivial(auto it, auto const & query, uint64_t const query_pos,
                         uint8_t const errors, uint8_t const max_errors, auto && delegate)
    {
        // Exact case
        if (query_pos == query.size() || errors == max_errors)
        {
            if (query_pos == query.size() || it.extend_right(ranges::view::drop_exactly(query, query_pos)))
            {
                delegate(it, errors);
                if constexpr (abort_on_hit)
                    return true;
                else
                    return; // specify return type before recursive calls
            }
        }
        // Approximate case
        else
        {
            // // Base case
            // if (query_pos == query.size())
            // {
            //     delegate(it, query, errors);
            //     if constexpr (abort_on_hit)
            //         return true;
            // }
            // else
            // // Recursive case

            // Insertion
            if constexpr (insertion && allow_insertion)
            {
                // do not allow deletion in the next step
                if constexpr (abort_on_hit)
                {
                    bool ret = _search_trivial<substitution, insertion, deletion, true, false, abort_on_hit>(it,
                        query, query_pos + 1, errors + 1, max_errors, delegate);
                    if (ret)
                        return true;
                }
                else
                {
                    _search_trivial<substitution, insertion, deletion, true, false, abort_on_hit>(it,
                        query, query_pos + 1, errors + 1, max_errors, delegate);
                }
            }

            if constexpr (deletion || substitution)
            {
                if (it.extend_right())
                {
                    do
                    {
                        // Match / Mismatch
                        if constexpr (substitution)
                        {
                            bool delta = it.last_char() != query[query_pos];
                            if constexpr (abort_on_hit)
                            {
                                bool ret = _search_trivial<substitution, insertion, deletion, true, true, abort_on_hit>(it,
                                    query, query_pos + 1, errors + delta, max_errors, delegate);
                                if (ret)
                                    return true;
                            }
                            else
                            {
                                _search_trivial<substitution, insertion, deletion, true, true, abort_on_hit>(it,
                                    query, query_pos + 1, errors + delta, max_errors, delegate);
                            }
                        }

                        // Deletion
                        if constexpr (deletion && allow_deletion)
                        {
                            // do not allow deletion in the next step
                            if constexpr (abort_on_hit)
                            {
                                bool ret = _search_trivial<substitution, insertion, deletion, false, true, abort_on_hit>(it,
                                    query, query_pos, errors + 1, max_errors, delegate);
                                if (ret)
                                    return true;
                            }
                            else
                            {
                                _search_trivial<substitution, insertion, deletion, false, true, abort_on_hit>(it,
                                    query, query_pos, errors + 1, max_errors, delegate);
                            }
                        }
                    }
                    while (it.cycle_back());
                }
            }
            else
            {
                // Match
                if (it.extend_right(query[query_pos]))
                {
                    if constexpr (abort_on_hit)
                    {
                        bool ret = _search_trivial<substitution, insertion, deletion, true, true, abort_on_hit>(it,
                                query, query_pos + 1, errors, max_errors, delegate);
                        if (ret)
                            return true;
                    }
                    else
                    {
                        _search_trivial<substitution, insertion, deletion, true, true, abort_on_hit>(it,
                            query, query_pos + 1, errors, max_errors, delegate);
                    }
                }
            }
        }
    }

    template <bool substitution, bool insertion, bool deletion, bool abort_on_hit, typename TIndex, typename TQuery, typename TDelegate>
    inline void search_trivial(TIndex const & index, TQuery const & query, uint8_t const max_errors, TDelegate && delegate)
    {
        // <substitution, insertion, deletion, allow_insertion, allow_deletion>
        // do not allow deletions at the beginning of the query
        _search_trivial<substitution, insertion, deletion, true, false, abort_on_hit>(index.begin(),
            query, 0, 0, max_errors, delegate);
    }

}

//!\endcond
