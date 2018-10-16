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
 * \brief Provides the algorithm to search in an index using optimum search schemes.
 */

#pragma once

#include <type_traits>

#include <range/v3/view/slice.hpp>

#include <seqan3/search/algorithm/detail/search_common.hpp>
#include <seqan3/search/algorithm/detail/search_scheme_precomputed.hpp>

/*!\addtogroup search
 * \{
 */

namespace seqan3::detail
{

//!\brief Computes a (not optimal) search scheme.
inline auto compute_search_scheme(uint8_t const min_error, uint8_t const max_error)
{
    // TODO: Simple backtracking. Replace this at least by the pigeonhole principle or even better by 01*0 schemes.
    std::vector<search_dyn> scheme{{{1}, {min_error}, {max_error}}};
    return scheme;
}

//!\brief Returns for each search the cumulative length of blocks in the order of blocks in each search.
inline auto search_scheme_block_info(auto const & search_scheme, uint64_t const query_length)
{
    // TODO: use std::array for optimum_search_scheme, std::vector otherwise
    std::vector<std::tuple<std::vector<uint64_t>, uint64_t>> result(search_scheme.size());

    uint8_t  const blocks     {search_scheme[0].blocks()};
    uint64_t const blocklength{query_length / blocks};
    uint8_t  const rest       {query_length - blocks * blocklength};

    std::vector<uint64_t> blocklengths(blocks, blocklength);
    for (uint8_t block_id = 0; block_id < rest; ++block_id)
        ++blocklengths[block_id];

    for (uint8_t search_id = 0; search_id < search_scheme.size(); ++search_id)
    {
        auto const & search = search_scheme[search_id];

        auto & [search_blocklengths, start_pos] = result[search_id];

        start_pos = 0;
        search_blocklengths.resize(blocks);
        search_blocklengths[0] = blocklengths[search.pi[0] - 1];
        for (uint8_t i = 1; i < blocks; ++i)
        {
            search_blocklengths[i] = blocklengths[search.pi[i] - 1] + search_blocklengths[i - 1];
            if (search.pi[i] < search.pi[0])
                start_pos += search_blocklengths[i] - search_blocklengths[i - 1];
        }
    }

    return result;
}

// forward declaration
template <bool abort_on_hit>
inline bool search_scheme_single_search(auto it, auto & query, uint64_t const lb, uint64_t const rb,
                                        uint8_t const errors_spent, uint8_t const block_id, bool const go_right,
                                        auto const & search, auto const & blocklength,
                                        search_params const error_left, auto && delegate);

template <bool abort_on_hit>
inline bool search_scheme_single_search_exact(auto it, auto & query, uint64_t const lb, uint64_t const rb,
                                              uint8_t const errors_spent, uint8_t const block_id, bool const go_right,
                                              auto const & search, auto const & blocklength,
                                              search_params const error_left, auto && delegate)
{
    uint8_t const block_id2 = std::min(block_id + 1, static_cast<uint8_t>(search.blocks()) - 1);
    bool const go_right2 = (block_id < search.blocks() - 1) && search.pi[block_id + 1] > search.pi[block_id];
    // NOTE: from sven:
    // bool const goToRight2 = (blockIndex < s.pi.size() - 1) ? s.pi[blockIndex + 1] > s.pi[blockIndex] : s.pi[blockIndex] > s.pi[blockIndex - 1];

    if (go_right)
    {
        uint64_t const infix_lb = rb - 1; // inclusive
        uint64_t const infix_rb = lb + blocklength[block_id] - 1; // exclusive

        if (!it.extend_right(query | ranges::view::slice(infix_lb, infix_rb + 1)))
            return false;

        if (search_scheme_single_search<abort_on_hit>(it, query, lb, infix_rb + 2, errors_spent, block_id2, go_right2, search, blocklength, error_left, delegate) && abort_on_hit)
            return true;
    }
    else
    {
        uint64_t const infix_lb = rb - blocklength[block_id] - 1; // inclusive
        uint64_t const infix_rb = lb - 1; // inclusive

        if (!it.extend_left(query | ranges::view::slice(infix_lb, infix_rb + 1)))
            return false;

        if (search_scheme_single_search<abort_on_hit>(it, query, infix_lb, rb, errors_spent, block_id2, go_right2, search, blocklength, error_left, delegate) && abort_on_hit)
            return true;
    }
    return false;
}

template <bool abort_on_hit>
inline bool search_scheme_single_search_deletion(auto it, auto & query, uint64_t const lb, uint64_t const rb,
                                                 uint8_t const errors_spent, uint8_t const block_id, bool const go_right,
                                                 auto const & search, auto const & blocklength,
                                                 search_params const error_left, auto && delegate)
{
    uint8_t const max_errors_left_in_block = search.u[block_id] - errors_spent;
    uint8_t const min_errors_left_in_block = std::max(search.l[block_id] - errors_spent, 0); // NOTE: changed

    if (min_errors_left_in_block == 0)
    {
        uint8_t const block_id2 = std::min(block_id + 1, static_cast<uint8_t>(search.blocks()) - 1);
        bool const go_right2 = search.pi[block_id2] > search.pi[block_id2 - 1];

        if (search_scheme_single_search<abort_on_hit>(it, query, lb, rb, errors_spent, block_id2, go_right2, search, blocklength, error_left, delegate) && abort_on_hit)
            return true;
    }

    if (!(search.pi[block_id] == 1 && !go_right) && max_errors_left_in_block > 0 && error_left.total > 0 && error_left.deletion > 0 && ((go_right && it.extend_right()) || (!go_right && it.extend_left())))
    {
        search_params error_left2{error_left};
        error_left2.total--;
        error_left2.deletion--;
        do
        {
            if (search_scheme_single_search_deletion<abort_on_hit>(it, query, lb, rb, errors_spent + 1, block_id, go_right, search, blocklength, error_left2, delegate) && abort_on_hit)
                return true;
        } while ((go_right && it.cycle_back()) || (!go_right && it.cycle_front()));
    }
    return false;
}

template <bool abort_on_hit>
inline bool search_scheme_single_search_children(auto it, auto & query, uint64_t const lb, uint64_t const rb,
                                                 uint8_t const errors_spent, uint8_t const block_id, bool const go_right,
                                                 uint8_t const min_errors_left_in_block,
                                                 auto const & search, auto const & blocklength,
                                                 search_params const error_left, auto && delegate)
{
    if ((go_right && it.extend_right()) || (!go_right && it.extend_left()))
    {
        uint64_t const chars_left = blocklength[block_id] - (rb - lb - 1);

        uint64_t lb2 = lb - !go_right;
        uint64_t rb2 = rb + go_right;

        do
        {
            bool const delta = it.last_char() != query[(go_right ? rb : lb) - 1];

            // NOTE: move that outside the if / do-while struct
            // NOTE: check that as well before doing non-edit-distance steps inside the loop
            if (!(/*error_left.insertion > 0 || */error_left.deletion > 0) && // TODO: remove insertion case
                min_errors_left_in_block > 0 && chars_left + delta < min_errors_left_in_block + 1u) // charsLeft - 1 < minErrorsLeftInBlock - delta
            {
                continue;
            }

            if (!delta || error_left.substitution > 0)
            {
                search_params error_left2{error_left};
                error_left2.total -= delta;
                error_left2.substitution -= delta;

                if (rb - lb == blocklength[block_id])
                {
                    // leave the possibility for one or multiple deletions! therefore, don't change direction, etc!
                    if (/*error_left.insertion > 0 || */error_left.deletion > 0) // TODO: remove insertion case
                    {
                        if (search_scheme_single_search_deletion<abort_on_hit>(it, query, lb2, rb2, errors_spent + delta, block_id, go_right, search, blocklength, error_left2, delegate) && abort_on_hit)
                            return true;
                    }
                    else
                    {
                        // NOTE: call std::min<uint8_t>
                        uint8_t const block_id2 = std::min(block_id + 1, static_cast<uint8_t>(search.blocks()) - 1);
                        bool const go_right2 = search.pi[block_id2] > search.pi[block_id2 - 1];

                        if (search_scheme_single_search<abort_on_hit>(it, query, lb2, rb2, errors_spent + delta, block_id2, go_right2, search, blocklength, error_left2, delegate) && abort_on_hit)
                            return true;
                    }
                }
                else
                {
                    if (search_scheme_single_search<abort_on_hit>(it, query, lb2, rb2, errors_spent + delta, block_id, go_right, search, blocklength, error_left2, delegate) && abort_on_hit)
                        return true;
                }
            }

            // Deletion
            // TODO: check in search_scheme_single_search when calling delegate that the last error was not a deletion
            // (no matter whether its the first or last block???)
            if (error_left.deletion > 0) // TODO: is this correct?
            {
                search_params error_left3{error_left};
                error_left3.total--;
                error_left3.deletion--;
                search_scheme_single_search<abort_on_hit>(it, query, lb, rb, errors_spent + 1, block_id, go_right, search, blocklength, error_left3, delegate); // NOTE: previous call to _deletion. wrong transfer from SeqAn2
            }

        } while ((go_right && it.cycle_back()) || (!go_right && it.cycle_front()));
    }
    return false;
}

template <bool abort_on_hit>
inline bool search_scheme_single_search(auto it, auto & query, uint64_t const lb, uint64_t const rb,
                                        uint8_t const errors_spent, uint8_t const block_id, bool const go_right,
                                        auto const & search, auto const & blocklength,
                                        search_params const error_left, auto && delegate)
{
    /*using namespace seqan3::literal;
    if (query == "AGC"_dna4)
    {
        std::cout << "OOO\n";
    }*/
    uint8_t const max_errors_left_in_block = search.u[block_id] - errors_spent;
    uint8_t const min_errors_left_in_block = std::max(search.l[block_id] - errors_spent, 0); // NOTE: changed

    // Done.
    if (min_errors_left_in_block == 0 && lb == 0 && rb == query.size() + 1)
    {
        delegate(it);
        return true;
    }
    // Exact search in current block.
    else if (max_errors_left_in_block == 0 && rb - lb - 1 != blocklength[block_id] || (error_left.total == 0 && min_errors_left_in_block == 0))
    {
        if (search_scheme_single_search_exact<abort_on_hit>(it, query, lb, rb, errors_spent, block_id, go_right, search, blocklength, error_left, delegate) && abort_on_hit);
            return true;
    }
    // Approximate search in current block.
    else if (error_left.total > 0)// (blocklength[block_id] - (rb - lb - (lb != rb)) >= min_errors_left_in_block)
    {
        // Insertion
        if (error_left.insertion > 0)
        {
            uint64_t const lb2 = lb - !go_right;
            uint64_t const rb2 = rb + go_right;

            search_params error_left2{error_left};
            error_left2.total--;
            error_left2.insertion--;
            if (rb - lb == blocklength[block_id])
            {
                // leave the possibility for one or multiple deletions! therefore, don't change direction, etc!
                // TODO: can the first "if statement" lead to a duplicate path in backtracking? should we enforce a
                // going down edges before calling this function?
                if (search_scheme_single_search_deletion<abort_on_hit>(it, query, lb2, rb2, errors_spent + 1, block_id, go_right, search, blocklength, error_left2, delegate) && abort_on_hit)
                    return true;
            }
            else
            {
                if (search_scheme_single_search<abort_on_hit>(it, query, lb2, rb2, errors_spent + 1, block_id, go_right, search, blocklength, error_left2, delegate) && abort_on_hit)
                    return true;
            }
        }
        if (search_scheme_single_search_children<abort_on_hit>(it, query, lb, rb, errors_spent, block_id, go_right, min_errors_left_in_block, search, blocklength, error_left, delegate) && abort_on_hit)
            return true;
    }
    return false;
}

template <bool abort_on_hit>
inline void search_search_scheme(auto const & index, auto & query, search_params const errors_left, auto const & search_scheme,
                                 auto && delegate)
{
    auto const block_info = search_scheme_block_info(search_scheme, query.size());

    for (uint8_t search_id = 0; search_id < search_scheme.size(); ++search_id)
    {
        auto const & search = search_scheme[search_id];
        auto const & [blocklength, start_pos] = block_info[search_id];

        bool const hit = search_scheme_single_search<abort_on_hit>(
                            index.begin(), query,     // data
                            start_pos, start_pos + 1, // infix range
                            0,                        // errors spent
                            0,                        // current block id
                            true,                     // go to the right
                            search, blocklength,      // search scheme information
                            errors_left,
                            delegate);

        if (abort_on_hit && hit)
            return;
    }
}

}

//!\}
