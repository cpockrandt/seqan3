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
 * \brief Provides the algorithm to search in an index using search schemes.
 */

#pragma once

#include <type_traits>

#include <range/v3/view/slice.hpp>

#include <seqan3/core/metafunction/transformation_trait_or.hpp>
#include <seqan3/search/algorithm/detail/search_common.hpp>
#include <seqan3/search/algorithm/detail/search_scheme_precomputed.hpp>
#include <seqan3/search/algorithm/detail/search_trivial.hpp>
#include <seqan3/search/fm_index/concept.hpp>

/*!\addtogroup search
 * \{
 */

namespace seqan3::detail
{

/*!\brief Computes a (non-optimal) search scheme. Currently the generated search scheme represents trivial
 *        backtracking.
 * \param[in] min_error Minimum number of errors allowed.
 * \param[in] max_error Maximum number of errors allowed.
 */
inline std::vector<search_dyn> compute_ss(uint8_t const min_error, uint8_t const max_error)
{
    // TODO: Replace this at least by the pigeonhole principle or even better by 01*0 schemes.
    std::vector<search_dyn> scheme{{{1}, {min_error}, {max_error}}};
    return scheme;
}

/*!\brief Returns for each search the cumulative length of blocks in the order of blocks in each search and the
 *        starting position of the first block in the query sequence.
 * \param[in] search_scheme Search scheme that will be used for searching.
 * \param[in] query_length Length of the query that will be searched in an index.
 */
template <typename search_scheme_t>
inline auto search_scheme_block_info(search_scheme_t const & search_scheme, uint64_t const query_length)
{
    using blocks_length_type = typename search_scheme_t::value_type::blocks_length_type;

    std::vector<std::tuple<blocks_length_type, uint64_t>> result(search_scheme.size());

    uint8_t  const blocks      {search_scheme[0].blocks()};
    uint64_t const block_length{query_length / blocks};
    uint8_t  const rest        {query_length - blocks * block_length};

    blocks_length_type blocks_length;
    // set all blocks_length values to block_length
    // resp. block_length + 1 for the first `rest = block_length % blocks` values
    if constexpr (std::is_same_v<blocks_length_type, std::vector<uint64_t>>)
        blocks_length.resize(blocks, block_length);
    else
        blocks_length.fill(block_length);

    for (uint8_t block_id = 0; block_id < rest; ++block_id)
        ++blocks_length[block_id];

    for (uint8_t search_id = 0; search_id < search_scheme.size(); ++search_id)
    {
        auto const & search = search_scheme[search_id];

        auto & [search_blocks_length, start_pos] = result[search_id];

        // compute cumulative blocks_length and starting position
        start_pos = 0;
        if constexpr (std::is_same_v<blocks_length_type, std::vector<uint64_t>>)
            search_blocks_length.resize(blocks);
        search_blocks_length[0] = blocks_length[search.pi[0] - 1];
        for (uint8_t i = 1; i < blocks; ++i)
        {
            search_blocks_length[i] = blocks_length[search.pi[i] - 1] + search_blocks_length[i - 1];
            if (search.pi[i] < search.pi[0])
                start_pos += search_blocks_length[i] - search_blocks_length[i - 1];
        }
    }

    return result;
}

//!\cond
// forward declaration
template <bool abort_on_hit>
inline bool search_ss(auto it, auto & query, uint64_t const lb, uint64_t const rb, uint8_t const errors_spent,
                      uint8_t const block_id, bool const go_right, auto const & search,
                      auto const & blocks_length, search_param const error_left, auto && delegate);
//!\endcond

/*!\brief Searches a query sequence in a bidirectional index using a single search of a search schemes.
 *        Sub-function for searching the remaining part of the current block without any errors.
 * \tparam abort_on_hit If the flag is set, the search aborts on the first hit.
 * \param[in] it Iterator of atring index built on the text that will be searched.
 * \param[in] query Query sequence to be searched with the iterator.
 * \param[in] lb Left bound of the infix of `query` already searched (exlusive).
 * \param[in] rb Right bound of the infix of `query` already searched (exlusive).
 * \param[in] errors_spent Number of errors spent while searching the infix of `query`.
 * \param[in] block_id Id of the block that infix is extended in next.
 * \param[in] go_right The infix will be extended to the right if the flag is set to true.
 * \param[in] search Search of a search scheme to be used for searching.
 * \param[in] blocks_length Cumulative block lengths of the search.
 * \param[in] error_left Number of errors left for matching the remaining suffix of the query sequence.
 * \param[in] delegate Function that is called on every hit. Takes `index::iterator_type` as argument.
 */
template <bool abort_on_hit>
inline bool search_ss_exact(auto it, auto & query, uint64_t const lb, uint64_t const rb,
                            uint8_t const errors_spent, uint8_t const block_id, bool const go_right,
                            auto const & search, auto const & blocks_length, search_param const error_left,
                            auto && delegate)
{
    uint8_t const block_id2 = std::min<uint8_t>(block_id + 1, search.blocks() - 1);
    bool const go_right2 = (block_id < search.blocks() - 1) && search.pi[block_id + 1] > search.pi[block_id];

    if (go_right)
    {
        uint64_t const infix_lb = rb - 1; // inclusive
        uint64_t const infix_rb = lb + blocks_length[block_id] - 1; // exclusive

        if (!it.extend_right(query | ranges::view::slice(infix_lb, infix_rb + 1)))
            return false;

        if (search_ss<abort_on_hit>(it, query, lb, infix_rb + 2, errors_spent, block_id2, go_right2, search,
                                    blocks_length, error_left, delegate) && abort_on_hit)
        {
            return true;
        }
    }
    else
    {
        uint64_t const infix_lb = rb - blocks_length[block_id] - 1; // inclusive
        uint64_t const infix_rb = lb - 1; // inclusive

        if (!it.extend_left(query | ranges::view::slice(infix_lb, infix_rb + 1)))
            return false;

        if (search_ss<abort_on_hit>(it, query, infix_lb, rb, errors_spent, block_id2, go_right2, search,
                                    blocks_length, error_left, delegate) && abort_on_hit)
        {
            return true;
        }
    }
    return false;
}

/*!\brief Searches a query sequence in a bidirectional index using a single search of a search schemes.
 *        Sub-function for deletions at the end of a block.
 * \tparam abort_on_hit If the flag is set, the search aborts on the first hit.
 * \param[in] it Iterator of atring index built on the text that will be searched.
 * \param[in] query Query sequence to be searched with the iterator.
 * \param[in] lb Left bound of the infix of `query` already searched (exlusive).
 * \param[in] rb Right bound of the infix of `query` already searched (exlusive).
 * \param[in] errors_spent Number of errors spent while searching the infix of `query`.
 * \param[in] block_id Id of the block that infix is extended in next.
 * \param[in] go_right The infix will be extended to the right if the flag is set to true.
 * \param[in] search Search of a search scheme to be used for searching.
 * \param[in] blocks_length Cumulative block lengths of the search.
 * \param[in] error_left Number of errors left for matching the remaining suffix of the query sequence.
 * \param[in] delegate Function that is called on every hit. Takes `index::iterator_type` as argument.
 */
template <bool abort_on_hit>
inline bool search_ss_deletion(auto it, auto & query, uint64_t const lb, uint64_t const rb,
                               uint8_t const errors_spent, uint8_t const block_id, bool const go_right,
                               auto const & search, auto const & blocks_length, search_param const error_left,
                               auto && delegate)
{
    uint8_t const max_error_left_in_block = search.u[block_id] - errors_spent;
    uint8_t const min_error_left_in_block = std::max(search.l[block_id] - errors_spent, 0);

    // Switch to the next block when the min number of errors is reached
    if (min_error_left_in_block == 0)
    {
        uint8_t const block_id2 = std::min<uint8_t>(block_id + 1, search.blocks() - 1);
        bool const go_right2 = search.pi[block_id2] > search.pi[block_id2 - 1];

        if (search_ss<abort_on_hit>(it, query, lb, rb, errors_spent, block_id2, go_right2, search,
                                    blocks_length, error_left, delegate) && abort_on_hit)
        {
            return true;
        }
    }

    // Insert deletions into the current block as long as possible
    if (!(search.pi[block_id] == 1 && !go_right) && max_error_left_in_block > 0 &&
        error_left.total > 0 && error_left.deletion > 0 &&
        ((go_right && it.extend_right()) || (!go_right && it.extend_left())))
    {
        search_param error_left2{error_left};
        error_left2.total--;
        error_left2.deletion--;
        do
        {
            if (search_ss_deletion<abort_on_hit>(it, query, lb, rb, errors_spent + 1, block_id, go_right,
                                                 search, blocks_length, error_left2, delegate) && abort_on_hit)
            {
                return true;
            }
        } while ((go_right && it.cycle_back()) || (!go_right && it.cycle_front()));
    }
    return false;
}

/*!\brief Searches a query sequence in a bidirectional index using a single search of a search schemes.
 *        Sub-function for approximate search step (iterating over all children in a conceptual suffix tree).
 * \tparam abort_on_hit If the flag is set, the search aborts on the first hit.
 * \param[in] it Iterator of atring index built on the text that will be searched.
 * \param[in] query Query sequence to be searched with the iterator.
 * \param[in] lb Left bound of the infix of `query` already searched (exlusive).
 * \param[in] rb Right bound of the infix of `query` already searched (exlusive).
 * \param[in] errors_spent Number of errors spent while searching the infix of `query`.
 * \param[in] block_id Id of the block that infix is extended in next.
 * \param[in] go_right The infix will be extended to the right if the flag is set to true.
 * \param[in] min_error_left_in_block Number of remaining errors that need to be spent in the current block.
 * \param[in] search Search of a search scheme to be used for searching.
 * \param[in] blocks_length Cumulative block lengths of the search.
 * \param[in] error_left Number of errors left for matching the remaining suffix of the query sequence.
 * \param[in] delegate Function that is called on every hit. Takes `index::iterator_type` as argument.
 */
template <bool abort_on_hit>
inline bool search_ss_children(auto it, auto & query, uint64_t const lb, uint64_t const rb,
                               uint8_t const errors_spent, uint8_t const block_id, bool const go_right,
                               uint8_t const min_error_left_in_block, auto const & search,
                               auto const & blocks_length, search_param const error_left, auto && delegate)
{
    if ((go_right && it.extend_right()) || (!go_right && it.extend_left()))
    {
        uint64_t const chars_left = blocks_length[block_id] - (rb - lb - 1);

        uint64_t lb2 = lb - !go_right;
        uint64_t rb2 = rb + go_right;

        do
        {
            bool const delta = it.last_char() != query[(go_right ? rb : lb) - 1];

            // skip if there are more min errors left in the current block than characters in the block
            // i.e. charsLeft - 1 < minErrorsLeftInBlock - delta
            // TODO: move that outside the if / do-while struct
            // TODO: incorporate error_left.deletion into formula and simplify a bit
            if (error_left.deletion == 0 && min_error_left_in_block > 0 &&
                chars_left + delta < min_error_left_in_block + 1u)
            {
                continue;
            }

            if (!delta || error_left.substitution > 0)
            {
                search_param error_left2{error_left};
                error_left2.total -= delta;
                error_left2.substitution -= delta;

                // At the end of the current block
                if (rb - lb == blocks_length[block_id])
                {
                    // Leave the possibility for one or multiple deletions at the end of a block.
                    // Thus do not change the direction (go_right) yet.
                    if (error_left.deletion > 0)
                    {
                        if (search_ss_deletion<abort_on_hit>(it, query, lb2, rb2, errors_spent + delta,
                                                             block_id, go_right, search, blocks_length,
                                                             error_left2, delegate) && abort_on_hit)
                        {
                            return true;
                        }
                    }
                    else
                    {
                        uint8_t const block_id2 = std::min<uint8_t>(block_id + 1, search.blocks() - 1);
                        bool const go_right2 = search.pi[block_id2] > search.pi[block_id2 - 1];

                        if (search_ss<abort_on_hit>(it, query, lb2, rb2, errors_spent + delta, block_id2,
                                                    go_right2, search, blocks_length, error_left2, delegate) &&
                            abort_on_hit)
                        {
                            return true;
                        }
                    }
                }
                else
                {
                    if (search_ss<abort_on_hit>(it, query, lb2, rb2, errors_spent + delta, block_id, go_right,
                                                search, blocks_length, error_left2, delegate) && abort_on_hit)
                    {
                        return true;
                    }
                }
            }

            // Deletion
            // TODO: check in search_ss when calling delegate that the last error was not a deletion
            // (no matter whether its the first or last block???)
            if (error_left.deletion > 0) // TODO: is this correct?
            {
                search_param error_left3{error_left};
                error_left3.total--;
                error_left3.deletion--;
                search_ss<abort_on_hit>(it, query, lb, rb, errors_spent + 1, block_id, go_right, search,
                                        blocks_length, error_left3, delegate);
            }
        } while ((go_right && it.cycle_back()) || (!go_right && it.cycle_front()));
    }
    return false;
}

/*!\brief Searches a query sequence in a bidirectional index using a single search of a search schemes.
 *        Kianfar, K., Pockrandt, C., Torkamandi, B., Luo, H., & Reinert, K. (2018).
 *        Optimum Search Schemes for Approximate String Matching Using Bidirectional FM-Index. bioRxiv, 301085.
 * \tparam abort_on_hit If the flag is set, the search aborts on the first hit.
 * \param[in] it Iterator of atring index built on the text that will be searched.
 * \param[in] query Query sequence to be searched with the iterator.
 * \param[in] lb Left bound of the infix of `query` already searched (exlusive).
 * \param[in] rb Right bound of the infix of `query` already searched (exlusive).
 * \param[in] errors_spent Number of errors spent while searching the infix of `query`.
 * \param[in] block_id Id of the block that infix is extended in next.
 * \param[in] go_right The infix will be extended to the right if the flag is set to true.
 * \param[in] search Search of a search scheme to be used for searching.
 * \param[in] blocks_length Cumulative block lengths of the search.
 * \param[in] error_left Number of errors left for matching the remaining suffix of the query sequence.
 * \param[in] delegate Function that is called on every hit. Takes `index::iterator_type` as argument.
 */
template <bool abort_on_hit>
inline bool search_ss(auto it, auto & query, uint64_t const lb, uint64_t const rb, uint8_t const errors_spent,
                      uint8_t const block_id, bool const go_right, auto const & search,
                      auto const & blocks_length, search_param const error_left, auto && delegate)
{
    uint8_t const max_error_left_in_block = search.u[block_id] - errors_spent;
    uint8_t const min_error_left_in_block = std::max(search.l[block_id] - errors_spent, 0); // NOTE: changed

    // Done.
    if (min_error_left_in_block == 0 && lb == 0 && rb == query.size() + 1)
    {
        delegate(it);
        return true;
    }
    // Exact search in current block.
    else if (max_error_left_in_block == 0 && rb - lb - 1 != blocks_length[block_id] ||
             (error_left.total == 0 && min_error_left_in_block == 0))
    {
        if (search_ss_exact<abort_on_hit>(it, query, lb, rb, errors_spent, block_id, go_right, search,
                                          blocks_length, error_left, delegate) && abort_on_hit)
        {
            return true;
        }
    }
    // Approximate search in current block.
    // i.e. blocks_length[block_id] - (rb - lb - (lb != rb)) >= min_error_left_in_block
    else if (error_left.total > 0)
    {
        // Insertion
        if (error_left.insertion > 0)
        {
            uint64_t const lb2 = lb - !go_right;
            uint64_t const rb2 = rb + go_right;

            search_param error_left2{error_left};
            error_left2.total--;
            error_left2.insertion--;
            // At the end of the current block
            if (rb - lb == blocks_length[block_id])
            {
                // Leave the possibility for one or multiple deletions at the end of a block.
                // Thus do not change the direction (go_right) yet.
                if (search_ss_deletion<abort_on_hit>(it, query, lb2, rb2, errors_spent + 1, block_id, go_right,
                                                     search, blocks_length, error_left2, delegate) &&
                    abort_on_hit)
                {
                    return true;
                }
            }
            else
            {
                if (search_ss<abort_on_hit>(it, query, lb2, rb2, errors_spent + 1, block_id, go_right, search,
                                            blocks_length, error_left2, delegate) && abort_on_hit)
                {
                    return true;
                }
            }
        }
        if (search_ss_children<abort_on_hit>(it, query, lb, rb, errors_spent, block_id, go_right,
                                             min_error_left_in_block, search, blocks_length, error_left,
                                             delegate) && abort_on_hit)
        {
            return true;
        }
    }
    return false;
}

/*!\brief Searches a query sequence in a bidirectional index using search schemes.
 * \tparam abort_on_hit If the flag is set, the search aborts on the first hit.
 * \param[in] index String index built on the text that will be searched.
 * \param[in] query Query sequence to be searched in the index.
 * \param[in] error_left Number of errors left for matching the remaining suffix of the query sequence.
 * \param[in] search_scheme Search scheme to be used for searching.
 * \param[in] delegate Function that is called on every hit. Takes `index::iterator_type` as argument.
 */
template <bool abort_on_hit>
inline void search_ss(auto const & index, auto & query, search_param const error_left,
                      auto const & search_scheme, auto && delegate)
{
    // retrieve cumulative block lengths and starting position
    auto const block_info = search_scheme_block_info(search_scheme, query.size());

    for (uint8_t search_id = 0; search_id < search_scheme.size(); ++search_id)
    {
        auto const & search = search_scheme[search_id];
        auto const & [blocks_length, start_pos] = block_info[search_id];

        bool const hit = search_ss<abort_on_hit>(
                             index.begin(),            // iterator on the index
                             query,                    // query to be searched
                             start_pos, start_pos + 1, // infix range already searched (open interval)
                                                       // the first character of `query` has the index 1 (not 0)
                             0,                        // errors spent
                             0,                        // current block id in search scheme
                             true,                     // search the first block from left to right
                             search, blocks_length,     // search scheme information
                             error_left,               // errors left (broken down by error types)
                             delegate                  // delegate function called on hit
                         );

        if (abort_on_hit && hit)
            return;
    }
}

/*!\brief Searches a query sequence in a bidirectional index.
 * \tparam abort_on_hit If the flag is set, the search aborts on the first hit.
 * \param[in] index String index built on the text that will be searched.
 * \param[in] query Query sequence to be searched in the index.
 * \param[in] error_left Number of errors left for matching the remaining suffix of the query sequence.
 * \param[in] delegate Function that is called on every hit. Takes `index::iterator_type` as argument.
 */
template <bool abort_on_hit>
inline void search_algo_bi(auto const & index, auto & query, search_param const error_left, auto && delegate)
{
    switch (error_left.total)
    {
        case 0:
            search_ss<abort_on_hit>(index, query, error_left, optimum_search_scheme<0, 0>, delegate);
            break;
        case 1:
            search_ss<abort_on_hit>(index, query, error_left, optimum_search_scheme<0, 1>, delegate);
            break;
        case 2:
            search_ss<abort_on_hit>(index, query, error_left, optimum_search_scheme<0, 2>, delegate);
            break;
        case 3:
            search_ss<abort_on_hit>(index, query, error_left, optimum_search_scheme<0, 3>, delegate);
            break;
        default:
            auto const & search_scheme{compute_ss(0, error_left.total)};
            search_ss<abort_on_hit>(index, query, error_left, search_scheme, delegate);
            break;
    }
}

/*!\brief Searches a query sequence in a unidirectional index.
 * \tparam abort_on_hit If the flag is set, the search aborts on the first hit.
 * \param[in] index String index built on the text that will be searched.
 * \param[in] query Query sequence to be searched in the index.
 * \param[in] error_left Number of errors left for matching the remaining suffix of the query sequence.
 * \param[in] delegate Function that is called on every hit. Takes `index::iterator_type` as argument.
 */
template <bool abort_on_hit>
inline void search_algo_uni(auto const & index, auto & query, search_param const error_left, auto && delegate)
{
    search_trivial<abort_on_hit>(index, query, error_left, delegate);
}

/*!\brief Searches a query sequence in an index.
 * \tparam abort_on_hit If the flag is set, the search aborts on the first hit.
 * \param[in] index String index built on the text that will be searched.
 * \param[in] query Query sequence to be searched in the index.
 * \param[in] error_left Number of errors left for matching the remaining suffix of the query sequence.
 * \param[in] delegate Function that is called on every hit. Takes `index::iterator_type` as argument.
 */
template <bool abort_on_hit, typename index_t>
inline void search_algo(index_t const & index, auto & query, search_param const error_left, auto && delegate)
{
    if constexpr (bi_fm_index_concept<index_t>)
        search_algo_bi<abort_on_hit>(index, query, error_left, delegate);
    else
        search_algo_uni<abort_on_hit>(index, query, error_left, delegate);
}

}

//!\}
