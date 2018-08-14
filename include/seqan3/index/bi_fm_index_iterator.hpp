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
 * \brief Provides the seqan3::bi_fm_index_iterator for searching in the bidirectional seqan3::bi_fm_index.
 */

#pragma once

#include <array>

#include <sdsl/suffix_trees.hpp>

#include <range/v3/view/iota.hpp>
#include <range/v3/view/slice.hpp>
#include <range/v3/view/transform.hpp>

#include <seqan3/alphabet/all.hpp>
#include <seqan3/index/detail/fm_index_iterator.hpp>
#include <seqan3/index/concept.hpp>
#include <seqan3/core/metafunction/range.hpp>

namespace seqan3
{

/*!\addtogroup index
 * \{
 */

// TODO: naming. better bwd than rev?

template <typename index_t>
class bi_fm_index_iterator
{

public:
    //!\brief Type of the index.
    using index_type = index_t;
    //!\brief Type for representing positions in the indexed text.
    using size_type = typename index_type::size_type;
    // using text_pos_t = typename index_t::size_type; // depends on dimensions of the text

    using fwd_iterator = fm_index_iterator<fm_index<typename index_type::text_type, typename index_type::index_traits::fm_index_traits>>;
    using rev_iterator = fm_index_iterator<fm_index<typename index_type::rev_text_type, typename index_type::index_traits::fm_index_traits>>;

protected:
    //!\privatesection

    //!\brief Type of the representation of characters in the underlying SDSL index.
    using sdsl_char_type = typename index_type::sdsl_char_type;

    //!\brief Type of the underlying FM index.
    index_type const * index;

    /*!\name Suffix array ranges of forward and reverse iterators.
     * \{
     */
     //!\brief Left suffix array range of the forward iterator.
    size_type fwd_lb;
    //!\brief Right suffix array range of the forward iterator.
    size_type fwd_rb;
    //!\brief Left suffix array range of the reverse iterator.
    size_type rev_lb;
    //!\brief Right suffix array range of the reverse iterator.
    size_type rev_rb;
    //\}

    // parent range and last_char only have to be stored for the iterator that has been used last for down() or right()
    // either fwd or rev. i.e. no need to store it twice. once the iterator is switched, the information become invalid
    // anyway
    /*!\name Information for right()
     *       Only stored for the iterator that has been used last to go down an edge because once one iterator is
     *       touched, the others parent information become invalid and cannot be used for right() anymore.
     * \{
     */
    //!\brief Left suffix array range of the parent node.
    size_type parent_lb;
    //!\brief Left suffix array range of the parent node.
    size_type parent_rb;
    //!\brief Label of the last edge moved down. Needed for right().
    sdsl_char_type last_char;
    //\}

    //!\brief Depth of the node in the suffix tree, i.e. length of the searched sequence.
    size_type m_depth; // equal for both iterators. only stored once

    // supports assertions to check whether right() resp. right_rev() is called on the same direction as the last
    // down([...]) resp. down_rev([...])
    #ifndef NDEBUG
        //!\brief Turns on debug mode for assert() in right() and right_rev()
        static constexpr bool is_debug = true;
        //!\brief Stores the information which iterator has been used last for down([...]) to allow for assert()
        //        in right() and right_rev()
        bool fwd_iter_last_used = false;
    #elif
        //!\brief Turns off debug mode for assert() in right() and right_rev()
        static constexpr bool is_debug = false;
    #endif

    //!\brief Helper function to recompute text positions since the indexed text is reversed.
    size_type offset() const noexcept
    {
        return index->size() - depth() - 1; // since the string is reversed during construction
    }

    //!\brief Optimized bidirectional search without alphabet mapping
    template <typename csa_t>
    bool bidirectional_search(csa_t const & csa,
                              size_type const l_fwd, size_type const r_fwd,
                              size_type const l_bwd, size_type const r_bwd,
                              sdsl_char_type const c,
                              size_type & l_fwd_res, size_type & r_fwd_res,
                              size_type & l_bwd_res, size_type & r_bwd_res) const noexcept
    {
        assert(l_fwd <= r_fwd && r_fwd < csa.size());
        assert(r_fwd + 1 >= l_fwd && r_bwd + 1 - l_bwd == r_fwd + 1 - l_fwd);

        if constexpr(std::is_same_v<typename csa_t::alphabet_type, sdsl::plain_byte_alphabet>)
        {
            size_type const c_begin = csa.C[c];

    		if (r_fwd + 1 - l_fwd == csa.size()) // TODO [[unlikely]]
            {
    			l_fwd_res = c_begin;
    			l_bwd_res = c_begin;
    			r_fwd_res = csa.C[c + 1] - 1;
    			r_bwd_res = r_fwd_res;
            }
            else
            {
                auto      const r_s_b   = csa.wavelet_tree.lex_count(l_fwd, r_fwd + 1, c);
                size_type const rank_l  = std::get<0>(r_s_b);
                size_type const s       = std::get<1>(r_s_b), b = std::get<2>(r_s_b);
                size_type const rank_r  = r_fwd - l_fwd - s - b + rank_l;
                l_fwd_res = c_begin + rank_l;
                r_fwd_res = c_begin + rank_r;
                l_bwd_res = l_bwd + s;
                r_bwd_res = r_bwd - b;
            }
        }
        else
        {
            size_type const cc = csa.char2comp[c];
            if (cc == 0 && c > 0) // [[unlikely]]
                return false;

        	size_type const c_begin = csa.C[cc];
    		if (r_fwd + 1 - l_fwd == csa.size()) // [[unlikely]]
            {
    			l_fwd_res = c_begin;
    			l_bwd_res = c_begin;
    			r_fwd_res = csa.C[cc + 1] - 1;
    			r_bwd_res = r_fwd_res;
            	return true;
            }
            else
            {
            	auto      const r_s_b  = csa.wavelet_tree.lex_count(l_fwd, r_fwd + 1, c);
            	size_type const rank_l = std::get<0>(r_s_b);
            	size_type const s      = std::get<1>(r_s_b), b = std::get<2>(r_s_b);
            	size_type const rank_r = r_fwd - l_fwd - s - b + rank_l;
            	l_fwd_res = c_begin + rank_l;
            	r_fwd_res = c_begin + rank_r;
            	l_bwd_res = l_bwd + s;
            	r_bwd_res = r_bwd - b;
            }
        }
        assert(r_fwd_res + 1 >= l_fwd_res && r_bwd_res + 1 - l_bwd_res == r_fwd_res + 1 - l_fwd_res);
        return r_fwd_res >= l_fwd_res;
    }

    //!\brief Optimized bidirectional search for right() and right_rev() without alphabet mapping
    template <typename csa_t>
    bool bidirectional_search_right(csa_t const & csa,
                                    size_type const l_fwd, size_type const r_fwd,
                                    size_type const l_bwd, size_type const r_bwd,
                                    sdsl_char_type const c,
                                    size_type & l_fwd_res, size_type & r_fwd_res,
                                    size_type & l_bwd_res, size_type & r_bwd_res) const noexcept
    {
        assert(l_fwd <= r_fwd && r_fwd < csa.size());

        if constexpr(std::is_same_v<typename csa_t::alphabet_type, sdsl::plain_byte_alphabet>)
        {
            size_type const c_begin = csa.C[c];
            auto const r_s_b = csa.wavelet_tree.lex_count(l_fwd, r_fwd + 1, c);
            size_type const rank_l  = std::get<0>(r_s_b);
            size_type const s = std::get<1>(r_s_b), b = std::get<2>(r_s_b);
            size_type const rank_r = r_fwd - l_fwd - s - b + rank_l;
            l_fwd_res = c_begin + rank_l;
            r_fwd_res = c_begin + rank_r;
            l_bwd_res = r_bwd + 1;
            r_bwd_res = r_bwd + 1 + rank_r - rank_l; // TODO: check maybe +1/-1 missing

            assert(r_fwd_res + 1 >= l_fwd_res && r_bwd_res + 1 - l_bwd_res == r_fwd_res + 1 - l_fwd_res);
            return r_fwd_res + 1 - l_fwd_res; // TODO: r_fwd_res >= l_fwd_res
        }
        else
        {
            size_type const c_begin = csa.C[csa.char2comp[c]];
            auto const r_s_b = csa.wavelet_tree.lex_count(l_fwd, r_fwd + 1, c);
            size_type const rank_l  = std::get<0>(r_s_b);
            size_type const s = std::get<1>(r_s_b), b = std::get<2>(r_s_b);
            size_type const rank_r = r_fwd - l_fwd - s - b + rank_l;
            l_fwd_res = c_begin + rank_l;
            r_fwd_res = c_begin + rank_r;
            l_bwd_res = r_bwd + 1;
            r_bwd_res = r_bwd + 1 + rank_r - rank_l; // TODO: check maybe +1/-1 missing

            assert(r_fwd_res + 1 >= l_fwd_res && r_bwd_res + 1 - l_bwd_res == r_fwd_res + 1 - l_fwd_res);
            return r_fwd_res + 1 - l_fwd_res; // TODO: r_fwd_res >= l_fwd_res
        }
    }

    //!\publicsection

public:

    /*!\name Constructors and destructor
     * \{
     */
    //!\brief Default constructor. Accessing member functions on a default constructed object is undefined behavior.
    //        Default construction is necessary to make this class semi-regular and e.g., to allow construction of
    //        std::array of iterators.
    bi_fm_index_iterator() noexcept = default;
    bi_fm_index_iterator(bi_fm_index_iterator const &) noexcept = default;
    bi_fm_index_iterator & operator=(bi_fm_index_iterator const &) noexcept = default;
    bi_fm_index_iterator(bi_fm_index_iterator &&) noexcept = default;
    bi_fm_index_iterator & operator=(bi_fm_index_iterator &&) noexcept = default;

    bi_fm_index_iterator(index_t const & _index) noexcept : index(&_index), m_depth(0),
                                                            fwd_lb(0), fwd_rb(_index.size() - 1),
                                                            rev_lb(0), rev_rb(_index.size() - 1)
    {}
    //\}

    /*!\brief Compares two iterators.
     * \param[in] rhs Other iterator to compare it to.
     * \returns `true` if both iterators are equal, `false` otherwise.
     *
     * ### Complexity
     *
     * Constant. (\todo: will change when comparing indices, not yet supported by SDSL)
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    bool operator==(bi_fm_index_iterator const & rhs) const noexcept
    {
        assert(index != nullptr);
        assert(!(fwd_lb == rhs.fwd_lb && fwd_rb == rhs.fwd_rb && m_depth == rhs.m_depth)
            || m_depth == 0 || parent_lb == rhs.parent_lb && parent_rb == rhs.parent_rb && last_char == rhs.last_char);

        return fwd_lb == rhs.fwd_lb && fwd_rb == rhs.fwd_rb && m_depth == rhs.m_depth;
    }

    /*!\brief Compares two iterators.
     * \param[in] rhs Other iterator to compare it to.
     * \returns `true` if the iterators are not equal, `false` otherwise.
     *
     * ### Complexity
     *
     * Constant. (\todo: will change when comparing indices, not yet supported by SDSL)
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    bool operator!=(bi_fm_index_iterator const & rhs) const noexcept
    {
        assert(index != nullptr);

        return !(*this == rhs);
    }

    /*!\brief Goes down the leftmost (lexicographically smallest) edge.
     * \returns `true` if the iterator moved down the leftmost edge successfully
     *          (i.e. it was not in a leaf node before).
     *
     * ### Complexity
     *
     * O(sigma) * O(T_BACKWARD_SEARCH). It scans linearly over the alphabet until it finds the smallest character that
     * is represented by an edge.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    bool down() noexcept
    {
        if constexpr(is_debug)
            fwd_iter_last_used = true;

        assert(index != nullptr);

        sdsl_char_type c = 1; // NOTE: start with 0 or 1 depending on implicit_sentintel
        size_type _fwd_lb, _fwd_rb, _rev_lb, _rev_rb;
        while (c < index->fwd_fm.m_index.sigma &&
               !bidirectional_search(index->fwd_fm.m_index, fwd_lb, fwd_rb, rev_lb, rev_rb,
                                     index->fwd_fm.m_index.comp2char[c], _fwd_lb, _fwd_rb, _rev_lb, _rev_rb))
        {
            ++c;
        }

        if (c != index->fwd_fm.m_index.sigma)
        {
            parent_lb = fwd_lb;
            parent_rb = fwd_rb;

            fwd_lb = _fwd_lb;
            fwd_rb = _fwd_rb;
            rev_lb = _rev_lb;
            rev_rb = _rev_rb;

            last_char = c;
            ++m_depth;

            return true;
        }
        return false;
    }

    bool down_rev() noexcept
    {
        if constexpr(is_debug)
            fwd_iter_last_used = false;

        assert(index != nullptr);

        sdsl_char_type c = 1; // NOTE: start with 0 or 1 depending on implicit_sentintel
        size_type _fwd_lb, _fwd_rb, _rev_lb, _rev_rb;
        while (c < index->rev_fm.m_index.sigma &&
               !bidirectional_search(index->rev_fm.m_index, rev_lb, rev_rb, fwd_lb, fwd_rb,
                                     index->rev_fm.m_index.comp2char[c], _rev_lb, _rev_rb, _fwd_lb, _fwd_rb))
        {
            ++c;
        }

        if (c != index->rev_fm.m_index.sigma)
        {
            parent_lb = rev_lb;
            parent_rb = rev_rb;

            fwd_lb = _fwd_lb;
            fwd_rb = _fwd_rb;
            rev_lb = _rev_lb;
            rev_rb = _rev_rb;

            last_char = c;
            ++m_depth;

            return true;
        }
        return false;
    }

    /*!\brief Goes down the edge labelled with `c`.
     * \tparam char_t Type of the character needs to be convertible to the character type `char_type` of the indexed
     *                text.
     * \param[in] c Character to extend the searched sequence with.
     * \returns `true` if the iterator moved down the corresponding edge, `false` otherwise (i.e. the searched sequence
     *          does not occur in the indexed text).
     *
     * ### Complexity
     *
     * O(T_BACKWARD_SEARCH)
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    template <alphabet_concept char_t>
    //!\cond
        requires implicitly_convertible_to_concept<char_t, typename index_t::char_type>
    //!\endcond
    bool down(char_t const c) noexcept
    {
        if constexpr(is_debug)
            fwd_iter_last_used = true;

        assert(index != nullptr);

        size_type _fwd_lb, _fwd_rb, _rev_lb, _rev_rb;

        auto c_char = to_rank(c) + 1;

        if (bidirectional_search(index->fwd_fm.m_index, fwd_lb, fwd_rb, rev_lb, rev_rb,
                                 c_char, _fwd_lb, _fwd_rb, _rev_lb, _rev_rb))
        {
            parent_lb = fwd_lb;
            parent_rb = fwd_rb;

            fwd_lb = _fwd_lb;
            fwd_rb = _fwd_rb;
            rev_lb = _rev_lb;
            rev_rb = _rev_rb;

            last_char = c_char;
            ++m_depth;

            return true;
        }
        return false;
    }

    template <alphabet_concept char_t>
    //!\cond
       requires implicitly_convertible_to_concept<char_t, typename index_t::char_type>
    //!\endcond
    bool down_rev(char_t const c) noexcept
    {
        if constexpr(is_debug)
            fwd_iter_last_used = false;

       assert(index != nullptr);

       size_type _fwd_lb, _fwd_rb, _rev_lb, _rev_rb;

       auto c_char = to_rank(c) + 1;

       if (bidirectional_search(index->rev_fm.m_index, rev_lb, rev_rb, fwd_lb, fwd_rb,
                                c_char, _rev_lb, _rev_rb, _fwd_lb, _fwd_rb))
       {
            parent_lb = rev_lb;
            parent_rb = rev_rb;

            fwd_lb = _fwd_lb;
            fwd_rb = _fwd_rb;
            rev_lb = _rev_lb;
            rev_rb = _rev_rb;

            last_char = c_char;
            ++m_depth;

            return true;
       }
       return false;
    }

    /*!\brief Goes down multiple edges labeled with `pattern`.
     * \tparam char_t Type of the character needs to be convertible to the character type `char_type` of the indexed
     *                text.
     * \param[in] pattern Sequence to extend the searched sequence with.
     * \returns `true` if the iterator moved down *all* corresponding edges, `false` otherwise (i.e. the searched
     *          sequence does not occur in the indexed text).
     *
     * If going down multiple edges fails in the middle of the sequence, all previous computations are rewound to
     * restore the iterator's state before calling this method.
     *
     * ### Complexity
     *
     * |pattern| * O(T_BACKWARD_SEARCH).
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    template <std::ranges::ForwardRange pattern_t>
    //!\cond
        requires implicitly_convertible_to_concept<innermost_value_type_t<pattern_t>, typename index_t::char_type>
    //!\endcond
    bool down(pattern_t && pattern) noexcept
    {
        if constexpr(is_debug)
            fwd_iter_last_used = true;

        auto first = pattern.cbegin();
        auto last = pattern.cend();

        assert(index != nullptr && first != last); // range must not be empty!

        size_type _fwd_lb = fwd_lb, _fwd_rb = fwd_rb, _rev_lb = rev_lb, _rev_rb = rev_rb;
        size_type _parent_lb, _parent_rb;

        sdsl_char_type c;

        for (auto it = first; it != last; ++it)
        {
            c = to_rank(*it) + 1;

            _parent_lb = _fwd_lb;
            _parent_rb = _fwd_rb;
            if (!bidirectional_search(index->fwd_fm.m_index, _fwd_lb, _fwd_rb, _rev_lb, _rev_rb,
                                      c, _fwd_lb, _fwd_rb, _rev_lb, _rev_rb))
                return false;
        }

        fwd_lb = _fwd_lb;
        fwd_rb = _fwd_rb;
        rev_lb = _rev_lb;
        rev_rb = _rev_rb;

        parent_lb = _parent_lb;
        parent_rb = _parent_rb;
        last_char = c;
        m_depth += last - first;

        return true;
    }

    template <std::ranges::ForwardRange pattern_t>
    //!\cond
       requires implicitly_convertible_to_concept<innermost_value_type_t<pattern_t>, typename index_t::char_type>
    //!\endcond
    bool down_rev(pattern_t && pattern) noexcept
    {
        if constexpr(is_debug)
            fwd_iter_last_used = false;

        auto first = pattern.cbegin();
        auto last = pattern.cend();

        assert(index != nullptr && first != last); // range must not be empty!

        size_type _fwd_lb = fwd_lb, _fwd_rb = fwd_rb, _rev_lb = rev_lb, _rev_rb = rev_rb;
        size_type _parent_lb, _parent_rb;

        sdsl_char_type c;

        // TODO: should we iterate from left to right or right to left?
        for (auto it = first; it != last; ++it)
        {
            c = to_rank(*it) + 1;

            _parent_lb = _rev_lb;
            _parent_rb = _rev_rb;
            if (!bidirectional_search(index->rev_fm.m_index, _rev_lb, _rev_rb, _fwd_lb, _fwd_rb,
                                      c, _rev_lb, _rev_rb, _fwd_lb, _fwd_rb))
            return false;
        }

        fwd_lb = _fwd_lb;
        fwd_rb = _fwd_rb;
        rev_lb = _rev_lb;
        rev_rb = _rev_rb;

        parent_lb = _parent_lb;
        parent_rb = _parent_rb;
        last_char = c;
        m_depth += last - first;

        return true;
    }

    /*!\brief Moves the iterator to the right sibling of the current suffix tree node. It would be equivalent to going
     *        up an edge and going down that edge with the smallest character that is smaller than the previous searched
     *        character. Calling right() on an iterator pointing to the root node is undefined behaviour!
     * \returns `true` if the current suffix tree node has a right sibling and the iterator was moved there.
     *
     * Example:
     *
     * ```cpp
     * auto it = index.root(); // create an iterator pointing to the root of the implicit suffix tree.
     * // it.right(); // WRONG! Undefined behaviour!
     * it.down(dna4::A);
     * it.right(); // search the sequence "C", "G" or "T" (the smallest of them that occurs in the text).
     *             // If none occurrs, false is returned.
     * it.down("AAC"_dna4); // search the sequence "AAC".
     * it.right(); // search the sequence "AAG", if not existant "AAT". If not existant either, return false.
     * ```
     *
     * ### Complexity
     *
     * O(sigma) * O(T_BACKWARD_SEARCH). It scans linearly over the alphabet starting from the last searched character
     * (parent edge of the last traversed edge) until it finds the right sibling edge.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    bool right() noexcept
    {
        if constexpr(is_debug)
            assert(fwd_iter_last_used); // right() can only be used into the same direction that has previously been used for down(...)
        assert(index != nullptr && depth() > 0);

        sdsl_char_type c = last_char + 1;
        size_type _fwd_lb, _fwd_rb, _rev_lb, _rev_rb;

        while (c < index->fwd_fm.m_index.sigma && !bidirectional_search_right(index->fwd_fm.m_index,
                   parent_lb, parent_rb, rev_lb, rev_rb, index->fwd_fm.m_index.comp2char[c],
                   _fwd_lb, _fwd_rb, _rev_lb, _rev_rb))
        {
            ++c;
        }

        if (c != index->fwd_fm.m_index.sigma)
        {
            fwd_lb = _fwd_lb;
            fwd_rb = _fwd_rb;
            rev_lb = _rev_lb;
            rev_rb = _rev_rb;

            last_char = c;

            return true;
        }
        return false;
    }

    bool right_rev() noexcept
    {
        if constexpr(is_debug)
            assert(!fwd_iter_last_used); // right() can only be used into the same direction that has previously been used for down(...)

        assert(index != nullptr && depth() > 0);

        sdsl_char_type c = last_char + 1;
        size_type _fwd_lb, _fwd_rb, _rev_lb, _rev_rb;

        while (c < index->rev_fm.m_index.sigma && !bidirectional_search_right(index->rev_fm.m_index,
                   parent_lb, parent_rb, fwd_lb, fwd_rb, index->rev_fm.m_index.comp2char[c],
                   _rev_lb, _rev_rb, _fwd_lb, _fwd_rb))
        {
            ++c;
        }

        if (c != index->rev_fm.m_index.sigma)
        {
            fwd_lb = _fwd_lb;
            fwd_rb = _fwd_rb;
            rev_lb = _rev_lb;
            rev_rb = _rev_rb;

            last_char = c;

            return true;
        }
        return false;
    }

    /*!\brief Returns an array of iterators pointing to the child nodes of the current iterator. Does not modify the
     *        current iterator.
     * \returns Array of iterators of size `alphabet_size(char_type)`, i.e. one iterator for each character. If the
     *          current node does not have an edge for each character, the remaining positions in the array will be
     *          filled with iterators pointing to the root.
     *
     * ### Complexity
     *
     * sigma * O(T_BACKWARD_SEARCH). The asymptotic running time is equal to enumerating all children using `down()` and
     * `right()` but has a better cache performance. (\todo: to be verified)
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    std::array<bi_fm_index_iterator, alphabet_size_v<typename index_type::char_type>> children() const noexcept
    {
        assert(index != nullptr);

        std::array<bi_fm_index_iterator, alphabet_size_v<typename index_type::char_type>> result;

        sdsl_char_type c = 1; // NOTE: start with 0 or 1 depending on implicit_sentintel

        uint8_t i = 0;
        while (c < index->fwd_fm.m_index.sigma)
        {
            size_type _fwd_lb, _fwd_rb, _rev_lb, _rev_rb;
            // TODO: this will be implemented much more efficiently in the future (performed rank queries are almost the
            // same). Rank information for different characters are located in the same cache line.
            if (bidirectional_search_right(index->fwd_fm.m_index, fwd_lb, fwd_rb, rev_lb, rev_rb,
                                           index->fwd_fm.m_index.comp2char[c], _fwd_lb, _fwd_rb, _rev_lb, _rev_rb))
            {
                if constexpr(is_debug)
                    result[i].fwd_iter_last_used = true;

                result[i].fwd_lb = _fwd_lb;
                result[i].fwd_rb = _fwd_rb;
                result[i].rev_lb = _rev_lb;
                result[i].rev_rb = _rev_rb;
                result[i].parent_lb = fwd_lb;
                result[i].parent_rb = fwd_rb;

                result[i].index = index;
                result[i].last_char = c;
                result[i].m_depth = m_depth + 1;
                ++i;
            }
            ++c;
        }

        // fill the rest with iterators pointing to root
        while (i < index_type::char_type::value_size)
        {
            result[i++] = bi_fm_index_iterator(*this->index);
        }

        return result;
    }

    std::array<bi_fm_index_iterator, alphabet_size_v<typename index_type::char_type>> children_rev() const noexcept
    {
        assert(index != nullptr);

        std::array<bi_fm_index_iterator, alphabet_size_v<typename index_type::char_type>> result;

        sdsl_char_type c = 1; // NOTE: start with 0 or 1 depending on implicit_sentintel

        uint8_t i = 0;
        while (c < index->rev_fm.m_index.sigma)
        {
            size_type _fwd_lb, _fwd_rb, _rev_lb, _rev_rb;
            // TODO: this will be implemented much more efficiently in the future (performed rank queries are almost the
            // same). Rank information for different characters are located in the same cache line.
            if (bidirectional_search_right(index->rev_fm.m_index, rev_lb, rev_rb, fwd_lb, fwd_rb,
                                           index->rev_fm.m_index.comp2char[c], _rev_lb, _rev_rb, _fwd_lb, _fwd_rb))
            {
                if constexpr(is_debug)
                    result[i].fwd_iter_last_used = false;

                result[i].fwd_lb = _fwd_lb;
                result[i].fwd_rb = _fwd_rb;
                result[i].rev_lb = _rev_lb;
                result[i].rev_rb = _rev_rb;
                result[i].parent_lb = rev_lb;
                result[i].parent_rb = rev_rb;

                result[i].index = index;
                result[i].last_char = c;
                result[i].m_depth = m_depth + 1;
                ++i;
            }
            ++c;
        }

        // fill the rest with iterators pointing to root
        while (i < index_type::char_type::value_size)
        {
            result[i++] = bi_fm_index_iterator(*this->index);
        }

        return result;
    }

    /*!\brief Returns the depth of the iterator node in the implicit suffix tree, i.e. the length of the sequence
     *        searched.
     * \returns Length of searched sequence.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    size_type depth() const noexcept
    {
        assert(index != nullptr);
        // depth == 0 -> root node
        assert(m_depth != 0 || (fwd_lb == rev_lb && fwd_rb == rev_rb && fwd_lb == 0 && fwd_rb == index->size() - 1));

        return m_depth;
    }

    /*!\brief Checks whether the iterator is at the root node.
     * \returns `true` if the iterator is pointing to the root node, `false` otherwise.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    bool is_root() const noexcept
    {
        assert(index != nullptr);

        return depth() == 0;
    }

    fwd_iterator get_fwd_iterator() const noexcept
    {
        assert(index != nullptr);

        fwd_iterator it(index->fwd_fm);
        it.parent_lb = parent_lb;
        it.parent_rb = parent_rb;
        it.node = {fwd_lb, fwd_rb, m_depth, last_char};

        if constexpr(is_debug)
        {
            if (!fwd_iter_last_used)
            {
                // invalidate parent range
                it.parent_lb = 1;
                it.parent_rb = 0;
            }
        }

        return it;
    }

    rev_iterator get_rev_iterator() const noexcept
    {
        assert(index != nullptr);

        rev_iterator it(index->rev_fm);
        it.parent_lb = parent_lb;
        it.parent_rb = parent_rb;
        it.node = {rev_lb, rev_rb, m_depth, last_char};

        if constexpr(is_debug)
        {
            if (fwd_iter_last_used)
            {
                // invalidate parent range
                it.parent_lb = 1;
                it.parent_rb = 0;
            }
        }

        return it;
    }

    /*!\brief Returns the searched sequence, i.e. a concatenation of all edges from the root node to the iterators
     *        current node.
     * \returns Searched sequence.
     *
     * ### Complexity
     *
     * O(sampling_rate * T_BACKWARD_SEARCH) + |depth()|.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    auto path_label() const noexcept
    {
        assert(index != nullptr && index->text != nullptr);

        size_type const pattern_begin = offset() - index->fwd_fm.m_index[fwd_lb];
        return *index->text | ranges::view::slice(pattern_begin, pattern_begin + depth());
    }

    //!\copydoc path_label()
    auto operator*() const noexcept
    {
       assert(index != nullptr && index->text != nullptr);

       return path_label();
    }

    /*!\brief Counts the number of occurrences of the searched sequence in the text.
     * \returns Number of occurrences of the searched sequence in the text.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    size_type count() const noexcept
    {
        assert(index != nullptr && 1 + fwd_rb - fwd_lb == 1 + rev_rb - rev_lb);

        return 1 + fwd_rb - fwd_lb;
    }

    /*!\brief Locates the occurrences of the searched sequence in the text.
     * \returns Positions in the text.
     *
     * ### Complexity
     *
     * count() * O(T_BACKWARD_SEARCH * SAMPLING_RATE).
     *
     * ### Exceptions
     *
     * Strong exception guarantee (no data is modified in case an exception is thrown).
     */
    std::vector<size_type> locate() const
    {
        assert(index != nullptr);

        std::vector<size_type> occ(count());
        for (size_type i = 0; i < occ.size(); ++i)
        {
            occ[i] = offset() - index->fwd_fm.m_index[fwd_lb + i];
        }
        return occ;
    }


    /*!\brief Locates the occurrences of the searched sequence in the text on demand, i.e. a ranges::view is returned
     *        and every position is located once it is accessed.
     * \returns Positions in the text.
     *
     * ### Complexity
     *
     * count() * O(T_BACKWARD_SEARCH * SAMPLING_RATE).
     *
     * ### Exceptions
     *
     * Strong exception guarantee (no data is modified in case an exception is thrown).
     */
    auto lazy_locate() const
    {
        assert(index != nullptr);

        size_type const _offset = offset();
        return ranges::view::iota(fwd_lb, fwd_lb + count())
               | ranges::view::transform([*this, _offset] (auto sa_pos) {
                   return _offset - index->fwd_fm.m_index[sa_pos];
               });
    }

};

//!\}

} // namespace seqan3
