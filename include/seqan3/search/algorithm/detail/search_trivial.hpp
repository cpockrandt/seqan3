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

struct search_params
{
    uint8_t total, substitution, insertion, deletion;

    // search_params & operator-=(search_params const & rhs)
    // {
    //     total -= rhs.total;
    //     substitution -= rhs.substitution;
    //     insertion -= rhs.insertion;
    //     deletion -= rhs.deletion;
    //
    //     return *this;
    // }
    //
    // search_params operator-(search_params const & rhs)
    // {
    //     search_params result{};
    //     result.total = total - rhs.total;
    //     result.substitution = total - rhs.substitution;
    //     result.insertion = total - rhs.insertion;
    //     result.deletion = total - rhs.deletion;
    //
    //     return result;
    // }
};

} // namespace seqan3::detail

namespace seqan3::detail
{

template <bool abort_on_hit>
auto _search_trivial(auto it, auto const & query, uint64_t const query_pos, search_params error_left,
                     auto && delegate)
{
    // Exact case
    if (query_pos == query.size() || error_left.total == 0)
    {
        if (query_pos == query.size() || it.extend_right(ranges::view::drop_exactly(query, query_pos)))
        {
            delegate(it/*, error_left.total*/);
            if constexpr (abort_on_hit)
                return true;
            else
                return; // specify return type "void" before recursive calls
        }
    }
    // Approximate case
    else
    {
        // Insertion
        if (error_left.insertion > 0)
        {
            search_params error_left2{error_left};
            error_left2.insertion--;
            error_left2.total--;

            // do not allow deletion in the next step
            if constexpr (abort_on_hit)
            {
                if (_search_trivial<abort_on_hit>(it, query, query_pos + 1, error_left2, delegate))
                    return true;
            }
            else
            {
                _search_trivial<abort_on_hit>(it, query, query_pos + 1, error_left2, delegate);
            }
        }

        if (((query_pos > 0 && error_left.deletion > 0) || error_left.substitution > 0) && it.extend_right())
        {
            do
            {
                // Match / Mismatch
                if (error_left.substitution > 0)
                {
                    bool delta = it.last_char() != query[query_pos];
                    search_params error_left2{error_left};
                    error_left2.total -= delta;
                    error_left2.substitution -= delta;

                    if constexpr (abort_on_hit)
                    {
                        if (_search_trivial<abort_on_hit>(it, query, query_pos + 1, error_left2, delegate))
                            return true;
                    }
                    else
                    {
                        _search_trivial<abort_on_hit>(it, query, query_pos + 1, error_left2, delegate);
                    }
                }

                // Deletion
                if (query_pos > 0) // do not allow deletions at the beginning of the query
                {
                    // Match
                    if (error_left.substitution == 0 && it.last_char() == query[query_pos])
                    {
                        if constexpr (abort_on_hit)
                        {
                            if (_search_trivial<abort_on_hit>(it, query, query_pos + 1, error_left, delegate))
                                return true;
                        }
                        else
                        {
                            _search_trivial<abort_on_hit>(it, query, query_pos + 1, error_left, delegate);
                        }
                    }

                    if (error_left.deletion > 0)
                    {
                        search_params error_left2{error_left};
                        error_left2.total--;
                        error_left2.deletion--;

                        // do not allow deletion in the next step
                        if constexpr (abort_on_hit)
                        {
                            if (_search_trivial<abort_on_hit>(it, query, query_pos, error_left2, delegate))
                                return true;
                        }
                        else
                        {
                            _search_trivial<abort_on_hit>(it, query, query_pos, error_left2, delegate);
                        }
                    }
                }
            }
            while (it.cycle_back());
        }
        else
        {
            // Match
            if (it.extend_right(query[query_pos]))
            {
                if constexpr (abort_on_hit)
                {
                    if (_search_trivial<abort_on_hit>(it, query, query_pos + 1, error_left, delegate))
                        return true;
                }
                else
                {
                    _search_trivial<abort_on_hit>(it, query, query_pos + 1, error_left, delegate);
                }
            }
        }
    }
}

template <bool abort_on_hit, typename TIndex, typename TQuery, typename TDelegate>
inline void search_trivial(TIndex const & index, TQuery const & query, search_params error_left, TDelegate && delegate)
{
    // <abort_on_hit>
    _search_trivial<abort_on_hit>(index.begin(), query, 0, error_left, delegate);
}

}

//!\endcond
