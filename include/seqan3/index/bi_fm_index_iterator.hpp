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
#include <seqan3/index/bi_fm_index.hpp>
#include <seqan3/core/metafunction/range.hpp>

namespace seqan3
{

/*!\addtogroup index
 * \{
 */

/*!\brief The SeqAn Bidirectional FM Index Iterator.
 * \ingroup bi_fm_index
 * \implements seqan3::bi_fm_index_iterator_concept
 * \tparam index_t The type of the underlying index; must satisfy seqan3::bi_fm_index_concept.
 * \details
 *
 * The iterator's interface provides searching a string both from left to right as well as from right to left in the
 * indexed text. It extends the interface of the unidirectional seqan3::fm_index_iterator.
 * All methods modifying the iterator (e.g. extending by a character with extend_right()) return a `bool` value whether
 * the operation was successful or not. In case of an unsuccessful operation the iterator remains unmodified, i.e. an
 * iterator can never be in an invalid state except default constructed iterators that are always invalid.
 *
 * \cond DEV
 *     The behaviour is equivalent to a prefix and suffix tree with the space and time efficiency of the underlying pure
 *     FM indices. The iterator traverses the implicit prefix and suffix trees beginning at the root node. The implicit
 *     prefix and suffix trees are not compacted, i.e. going down an edge using extend_right(char) will increase the
 *     query by only one character.
 * \endcond
 *
 * The asymptotic running times for using the iterator depend on the SDSL index configuration. To determine the exact
 * running times, you have to additionally look up the running times of the used traits (configuration).
 */
template <typename index_t>
class bi_fm_index_iterator
{

public:
    //!\brief Type of the index.
    using index_type = index_t;
    //!\brief Type for representing positions in the indexed text.
    using size_type = typename index_type::size_type;

    //!\brief Type for the unidirectional iterator on the original text.
    using fwd_iterator = fm_index_iterator<fm_index<typename index_type::text_type, typename index_type::index_traits::fm_index_traits>>;
    //!\brief Type for the unidirectional iterator on the reversed text.
    using rev_iterator = fm_index_iterator<fm_index<typename index_type::rev_text_type, typename index_type::index_traits::fm_index_traits>>;

protected:
    //!\privatesection

    //!\brief Type of the representation of characters in the underlying SDSL index.
    using sdsl_char_type = typename index_type::sdsl_char_type;

    //!\brief Type of the underlying FM index.
    index_type const * m_index;

    /*!\name Suffix array ranges of forward and reverse iterators.
     * \{
     */
     //!\brief Left suffix array range of the forward iterator (for extend_right).
    size_type m_fwd_lb;
    //!\brief Right suffix array range of the forward iterator (for extend_right).
    size_type m_fwd_rb;
    //!\brief Left suffix array range of the reverse iterator (for extend_left).
    size_type m_rev_lb;
    //!\brief Right suffix array range of the reverse iterator (for extend_left).
    size_type m_rev_rb;
    //\}

    /*!\name Information for cycle_back() and cycle_front()
     *       Only stored for the iterator that has been used last to go down an edge because once one iterator is
     *       touched, the others parent information become invalid and cannot be used for cycle_back() anymore.
     * \{
     */

    // m_parent_* and m_last_char only have to be stored for the (unidirectional) iterator that has been used last for
    // extend_right() or cycle_back() resp. extend_left() or cycle_front(), (i.e. either fwd or rev). Thus there is no
    // need to store it twice. Once the iterator is switched, the information becomes invalid anyway.

    //!\brief Left suffix array range of the parent node.
    size_type m_parent_lb;
    //!\brief Left suffix array range of the parent node.
    size_type m_parent_rb;
    //!\brief Label of the last edge moved down. Needed for cycle_back() or cycle_front().
    sdsl_char_type m_last_char;
    //\}

    //!\brief Depth of the node in the suffix tree, i.e. length of the searched query.
    size_type m_depth; // equal for both iterators. only stored once

    // supports assertions to check whether cycle_back() resp. cycle_front() is called on the same direction as the last
    // extend_right([...]) resp. extend_left([...])
    #ifndef NDEBUG
        //!\brief Stores the information which iterator has been used last for extend_*([...]) to allow for assert() in
        //        cycle_back() and cycle_front()
        bool fwd_iter_last_used = false;
    #endif

    //!\brief Helper function to recompute text positions since the indexed text is reversed.
    size_type offset() const noexcept
    {
        return m_index->size() - query_length() - 1; // since the string is reversed during construction
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

    		if (r_fwd + 1 - l_fwd == csa.size()) // [[unlikely]]
            {
    			l_fwd_res = c_begin;
    			l_bwd_res = c_begin;
    			r_fwd_res = csa.C[c + 1] - 1;
    			r_bwd_res = r_fwd_res;
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
        else
        {
            size_type const cc = csa.char2comp[c];
            if (cc == 0 && c > 0) // [[unlikely]]
            {
                return false;
            }

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

    //!\brief Optimized bidirectional search for cycle_back() and cycle_front() without alphabet mapping
    template <typename csa_t>
    bool bidirectional_search_cycle(csa_t const & csa,
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
            return r_fwd_res >= l_fwd_res; // r_fwd_res + 1 - l_fwd_res
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
            return r_fwd_res >= l_fwd_res; // r_fwd_res + 1 - l_fwd_res
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

    bi_fm_index_iterator(index_t const & index) noexcept : m_index(&index), m_depth(0),
                                                           m_fwd_lb(0), m_fwd_rb(index.size() - 1),
                                                           m_rev_lb(0), m_rev_rb(index.size() - 1)
    {}
    //\}

    /*!\brief Compares two iterators.
     * \param[in] rhs Other iterator to compare it to.
     * \returns `true` if both iterators are equal, `false` otherwise.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    bool operator==(bi_fm_index_iterator const & rhs) const noexcept
    {
        assert(m_index != nullptr);
        // equal SA range implies equal parent node information (or both are root nodes)
        assert(!(m_fwd_lb == rhs.m_fwd_lb && m_fwd_rb == rhs.m_fwd_rb && m_depth == rhs.m_depth)
            || m_depth == 0
            || m_parent_lb == rhs.m_parent_lb && m_parent_rb == rhs.m_parent_rb && m_last_char == rhs.m_last_char);

        return m_fwd_lb == rhs.m_fwd_lb && m_fwd_rb == rhs.m_fwd_rb && m_depth == rhs.m_depth;
    }

    /*!\brief Compares two iterators.
     * \param[in] rhs Other iterator to compare it to.
     * \returns `true` if the iterators are not equal, `false` otherwise.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    bool operator!=(bi_fm_index_iterator const & rhs) const noexcept
    {
        assert(m_index != nullptr);

        return !(*this == rhs);
    }

    /*!\brief Tries to extend the query by the smallest possible character to the right such that the query is found in
     *        the text.
     *        \cond DEV
     *            Goes down the leftmost (i.e. lexicographically smallest) edge.
     *        \endcond
     * \returns `true` if the iterator could extend the query successfully.
     *
     * ### Complexity
     *
     * O(SIGMA) * O(T_BACKWARD_SEARCH). It scans linearly over the alphabet until it finds the smallest character that
     * is represented by an edge.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    bool extend_right() noexcept
    {
        #ifndef NDEBUG
            fwd_iter_last_used = true;
        #endif

        assert(m_index != nullptr);

        sdsl_char_type c = 1; // NOTE: start with 0 or 1 depending on implicit_sentintel
        size_type _m_fwd_lb, _m_fwd_rb, _m_rev_lb, _m_rev_rb;
        while (c < m_index->fwd_fm.m_index.sigma &&
               !bidirectional_search(m_index->fwd_fm.m_index, m_fwd_lb, m_fwd_rb, m_rev_lb, m_rev_rb,
                                     m_index->fwd_fm.m_index.comp2char[c], _m_fwd_lb, _m_fwd_rb, _m_rev_lb, _m_rev_rb))
        {
            ++c;
        }

        if (c != m_index->fwd_fm.m_index.sigma)
        {
            m_parent_lb = m_fwd_lb;
            m_parent_rb = m_fwd_rb;

            m_fwd_lb = _m_fwd_lb;
            m_fwd_rb = _m_fwd_rb;
            m_rev_lb = _m_rev_lb;
            m_rev_rb = _m_rev_rb;

            m_last_char = c;
            ++m_depth;

            return true;
        }
        return false;
    }

    /*!\brief Tries to extend the query by the smallest possible character to the left such that the query is found in
     *        the text.
     *        \cond DEV
     *            Goes down the leftmost (i.e. lexicographically smallest) edge in the reverse iterator.
     *        \endcond
     * \returns `true` if the iterator could extend the query successfully.
     *
     * ### Complexity
     *
     * O(SIGMA) * O(T_BACKWARD_SEARCH). It scans linearly over the alphabet until it finds the smallest character that
     * is represented by an edge.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    bool extend_left() noexcept
    {
        #ifndef NDEBUG
            fwd_iter_last_used = false;
        #endif

        assert(m_index != nullptr);

        sdsl_char_type c = 1; // NOTE: start with 0 or 1 depending on implicit_sentintel
        size_type _m_fwd_lb, _m_fwd_rb, _m_rev_lb, _m_rev_rb;
        while (c < m_index->rev_fm.m_index.sigma &&
               !bidirectional_search(m_index->rev_fm.m_index, m_rev_lb, m_rev_rb, m_fwd_lb, m_fwd_rb,
                                     m_index->rev_fm.m_index.comp2char[c], _m_rev_lb, _m_rev_rb, _m_fwd_lb, _m_fwd_rb))
        {
            ++c;
        }

        if (c != m_index->rev_fm.m_index.sigma)
        {
            m_parent_lb = m_rev_lb;
            m_parent_rb = m_rev_rb;

            m_fwd_lb = _m_fwd_lb;
            m_fwd_rb = _m_fwd_rb;
            m_rev_lb = _m_rev_lb;
            m_rev_rb = _m_rev_rb;

            m_last_char = c;
            ++m_depth;

            return true;
        }
        return false;
    }

    /*!\brief Tries to extend the query by the character `c` to the right.
     * \tparam char_t Type of the character needs to be convertible to the character type `char_type` of the indexed
     *                text.
     * \param[in] c Character to extend the query with to the right.
     * \returns `true` if the iterator could extend the query successfully.
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
    bool extend_right(char_t const c) noexcept
    {
        #ifndef NDEBUG
            fwd_iter_last_used = true;
        #endif

        assert(m_index != nullptr);

        size_type _m_fwd_lb, _m_fwd_rb, _m_rev_lb, _m_rev_rb;

        auto c_char = to_rank(c) + 1;

        if (bidirectional_search(m_index->fwd_fm.m_index, m_fwd_lb, m_fwd_rb, m_rev_lb, m_rev_rb,
                                 c_char, _m_fwd_lb, _m_fwd_rb, _m_rev_lb, _m_rev_rb))
        {
            m_parent_lb = m_fwd_lb;
            m_parent_rb = m_fwd_rb;

            m_fwd_lb = _m_fwd_lb;
            m_fwd_rb = _m_fwd_rb;
            m_rev_lb = _m_rev_lb;
            m_rev_rb = _m_rev_rb;

            m_last_char = c_char;
            ++m_depth;

            return true;
        }
        return false;
    }

    /*!\brief Tries to extend the query by the character `c` to the left.
     * \tparam char_t Type of the character needs to be convertible to the character type `char_type` of the indexed
     *                text.
     * \param[in] c Character to extend the query with to the left.
     * \returns `true` if the iterator could extend the query successfully.
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
    bool extend_left(char_t const c) noexcept
    {
        #ifndef NDEBUG
            fwd_iter_last_used = false;
        #endif

       assert(m_index != nullptr);

       size_type _m_fwd_lb, _m_fwd_rb, _m_rev_lb, _m_rev_rb;

       auto c_char = to_rank(c) + 1;

       if (bidirectional_search(m_index->rev_fm.m_index, m_rev_lb, m_rev_rb, m_fwd_lb, m_fwd_rb,
                                c_char, _m_rev_lb, _m_rev_rb, _m_fwd_lb, _m_fwd_rb))
       {
            m_parent_lb = m_rev_lb;
            m_parent_rb = m_rev_rb;

            m_fwd_lb = _m_fwd_lb;
            m_fwd_rb = _m_fwd_rb;
            m_rev_lb = _m_rev_lb;
            m_rev_rb = _m_rev_rb;

            m_last_char = c_char;
            ++m_depth;

            return true;
       }
       return false;
    }

    /*!\brief Tries to extend the query by `seq` to the right.
     * \tparam query_t Type of the character needs to be convertible to the character type `char_type` of the indexed
     *                 text.
     * \param[in] seq Sequence to extend the query with to the right.
     * \returns `true` if the iterator could extend the query successfully.
     *
     * If extending fails in the middle of the sequence, all previous computations are rewound to restore the iterator's
     * state before calling this method.
     *
     * ### Complexity
     *
     * |seq| * O(T_BACKWARD_SEARCH).
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    template <std::ranges::ForwardRange query_t>
    //!\cond
        requires implicitly_convertible_to_concept<innermost_value_type_t<query_t>, typename index_t::char_type>
    //!\endcond
    bool extend_right(query_t && seq) noexcept
    {
        #ifndef NDEBUG
            fwd_iter_last_used = true;
        #endif

        auto first = seq.begin();
        auto last = seq.end();

        assert(m_index != nullptr && first != last); // range must not be empty!

        size_type _m_fwd_lb = m_fwd_lb, _m_fwd_rb = m_fwd_rb, _m_rev_lb = m_rev_lb, _m_rev_rb = m_rev_rb;
        size_type _m_parent_lb, _m_parent_rb;

        sdsl_char_type c;

        for (auto it = first; it != last; ++it)
        {
            c = to_rank(*it) + 1;

            _m_parent_lb = _m_fwd_lb;
            _m_parent_rb = _m_fwd_rb;
            if (!bidirectional_search(m_index->fwd_fm.m_index, _m_fwd_lb, _m_fwd_rb, _m_rev_lb, _m_rev_rb,
                                      c, _m_fwd_lb, _m_fwd_rb, _m_rev_lb, _m_rev_rb))
                return false;
        }

        m_fwd_lb = _m_fwd_lb;
        m_fwd_rb = _m_fwd_rb;
        m_rev_lb = _m_rev_lb;
        m_rev_rb = _m_rev_rb;

        m_parent_lb = _m_parent_lb;
        m_parent_rb = _m_parent_rb;
        m_last_char = c;
        m_depth += last - first;

        return true;
    }

    /*!\brief Tries to extend the query by `seq` to the left.
     * \tparam query_t Type of the character needs to be convertible to the character type `char_type` of the indexed
     *                 text.
     * \param[in] seq Sequence to extend the query with to the left.
     * \returns `true` if the iterator could extend the query successfully.
     *
     * If extending fails in the middle of the sequence, all previous computations are rewound to restore the iterator's
     * state before calling this method.
     *
     * ### Complexity
     *
     * |seq| * O(T_BACKWARD_SEARCH).
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    template <std::ranges::ForwardRange query_t>
    //!\cond
       requires implicitly_convertible_to_concept<innermost_value_type_t<query_t>, typename index_t::char_type>
    //!\endcond
    bool extend_left(query_t && seq) noexcept
    {
        #ifndef NDEBUG
            fwd_iter_last_used = false;
        #endif

        auto first = seq.begin();
        auto last = seq.end();

        assert(m_index != nullptr && first != last); // range must not be empty!

        size_type _m_fwd_lb = m_fwd_lb, _m_fwd_rb = m_fwd_rb, _m_rev_lb = m_rev_lb, _m_rev_rb = m_rev_rb;
        size_type _m_parent_lb, _m_parent_rb;

        sdsl_char_type c;

        // TODO(h-2, eseiler and cpockrandt): should we iterate from left to right or right to left?
        for (auto it = first; it != last; ++it)
        {
            c = to_rank(*it) + 1;

            _m_parent_lb = _m_rev_lb;
            _m_parent_rb = _m_rev_rb;
            if (!bidirectional_search(m_index->rev_fm.m_index, _m_rev_lb, _m_rev_rb, _m_fwd_lb, _m_fwd_rb,
                                      c, _m_rev_lb, _m_rev_rb, _m_fwd_lb, _m_fwd_rb))
            return false;
        }

        m_fwd_lb = _m_fwd_lb;
        m_fwd_rb = _m_fwd_rb;
        m_rev_lb = _m_rev_lb;
        m_rev_rb = _m_rev_rb;

        m_parent_lb = _m_parent_lb;
        m_parent_rb = _m_parent_rb;
        m_last_char = c;
        m_depth += last - first;

        return true;
    }

    /*!\brief Tries to replace the rightmost character of the query by the next lexicographically larger character such
     *        that the query is found in the text.
     *        \cond DEV
     *            Moves the iterator to the right sibling of the current suffix tree node. It would be equivalent to
     *            going up an edge and going down that edge with the smallest character that is larger than the
     *            previous searched character. Calling cycle_*() on an iterator pointing to the root node is undefined
     *            behaviour!
     *        \endcond
     * \returns `true` if there exists a query in the text where the rightmost character of the query is
     *          lexicographically larger than the current rightmost character of the query.
     *
     * Example:
     *
     * \snippet test/snippet/index/bi_fm_index_iterator.cpp cycle
     *
     * ### Complexity
     *
     * O(SIGMA) * O(T_BACKWARD_SEARCH). It scans linearly over the alphabet starting from the rightmost character
     * until it finds the query with a larger rightmost character.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    bool cycle_back() noexcept
    {
        // cycle_back() can only be used into the same direction that has previously been used for down(...)
        #ifndef NDEBUG
            assert(fwd_iter_last_used);
        #endif

        assert(m_index != nullptr && query_length() > 0);

        sdsl_char_type c = m_last_char + 1;
        size_type _m_fwd_lb, _m_fwd_rb, _m_rev_lb, _m_rev_rb;

        while (c < m_index->fwd_fm.m_index.sigma && !bidirectional_search_cycle(m_index->fwd_fm.m_index,
                   m_parent_lb, m_parent_rb, m_rev_lb, m_rev_rb, m_index->fwd_fm.m_index.comp2char[c],
                   _m_fwd_lb, _m_fwd_rb, _m_rev_lb, _m_rev_rb))
        {
            ++c;
        }

        if (c != m_index->fwd_fm.m_index.sigma)
        {
            m_fwd_lb = _m_fwd_lb;
            m_fwd_rb = _m_fwd_rb;
            m_rev_lb = _m_rev_lb;
            m_rev_rb = _m_rev_rb;

            m_last_char = c;

            return true;
        }
        return false;
    }

    /*!\brief Tries to replace the leftmost character of the query by the next lexicographically larger character such
     *        that the query is found in the text.
     *        \cond DEV
     *            Moves the iterator to the right sibling of the current suffix tree node. It would be equivalent to
     *            going up an edge and going down that edge with the smallest character that is larger than the
     *            previous searched character. Calling cycle_*() on an iterator pointing to the root node is undefined
     *            behaviour!
     *        \endcond
     * \returns `true` if there exists a query in the text where the leftmost character of the query is
     *          lexicographically larger than the current leftmost character of the query.
     *
     * Example:
     *
     * \snippet test/snippet/index/bi_fm_index_iterator.cpp cycle
     *
     * ### Complexity
     *
     * O(SIGMA) * O(T_BACKWARD_SEARCH). It scans linearly over the alphabet starting from the leftmost character
     * until it finds the query with a larger leftmost character.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    bool cycle_front() noexcept
    {
        // cycle_front() can only be used into the same direction that has previously been used for down(...)
        #ifndef NDEBUG
            assert(!fwd_iter_last_used);
        #endif

        assert(m_index != nullptr && query_length() > 0);

        sdsl_char_type c = m_last_char + 1;
        size_type _m_fwd_lb, _m_fwd_rb, _m_rev_lb, _m_rev_rb;

        while (c < m_index->rev_fm.m_index.sigma && !bidirectional_search_cycle(m_index->rev_fm.m_index,
                   m_parent_lb, m_parent_rb, m_fwd_lb, m_fwd_rb, m_index->rev_fm.m_index.comp2char[c],
                   _m_rev_lb, _m_rev_rb, _m_fwd_lb, _m_fwd_rb))
        {
            ++c;
        }

        if (c != m_index->rev_fm.m_index.sigma)
        {
            m_fwd_lb = _m_fwd_lb;
            m_fwd_rb = _m_fwd_rb;
            m_rev_lb = _m_rev_lb;
            m_rev_rb = _m_rev_rb;

            m_last_char = c;

            return true;
        }
        return false;
    }


    /*!\brief Outputs the rightmost respectively leftmost character depending on whether extend_right() or extend_left()
     *        has been called last.
     * \returns Rightmost or leftmost character.
     *
     * Example:
     *
     * \snippet test/snippet/index/bi_fm_index_iterator.cpp cycle
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    typename index_t::char_type last_char() noexcept
    {
        assert(m_index != nullptr && query_length() > 0);

        typename index_t::char_type c;
        c.assign_rank(m_index->fwd_fm.m_index.comp2char[m_last_char] - 1); // text is not allowed to contain ranks of 0
        return c;
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
    size_type query_length() const noexcept
    {
        assert(m_index != nullptr);
        // depth == 0 -> root node
        assert(m_depth != 0 ||
            (m_fwd_lb == m_rev_lb && m_fwd_rb == m_rev_rb && m_fwd_lb == 0 && m_fwd_rb == m_index->size() - 1));

        return m_depth;
    }

    /*!\brief Returns a unidirectional seqan3::fm_index_iterator on the original text. query() on the returned
     *        unidirectional index iterator will be equal to query() on the bidirectional index iterator.
     *        cycle_back() and last_char() will be undefined behavior if the last extension on the bidirectional
     *        FM index has been to the left. The behavior will be well-defined after the first extension to the right
     *        on the unidirectional index.
     * \returns Returns a unidirectional seqan3::fm_index_iterator on the index of the original text.
     *
     * Example:
     *
     * \snippet test/snippet/index/bi_fm_index_iterator.cpp get_fwd_iterator
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    fwd_iterator get_fwd_iterator() const noexcept
    {
        assert(m_index != nullptr);

        fwd_iterator it(m_index->fwd_fm);
        it.m_parent_lb = m_parent_lb;
        it.m_parent_rb = m_parent_rb;
        it.m_node = {m_fwd_lb, m_fwd_rb, m_depth, m_last_char};

        #ifndef NDEBUG
            if (!fwd_iter_last_used)
            {
                // invalidate parent range
                it.m_parent_lb = 1;
                it.m_parent_rb = 0;
            }
        #endif

        return it;
    }

    /*!\brief Returns a unidirectional seqan3::fm_index_iterator on the reversed text. query() on the returned
     *        unidirectional index iterator will be equal to reversing query() on the bidirectional index iterator.
     *        Note that because of the text being reversed, extend_right() resp. cycle_back()
     *        correspond to extend_left() resp. cycle_front() on the bidirectional index iterator.
     *        Furthermore cycle_back() and last_char() will be undefined behavior if the last extension on the
     *        bidirectional FM index has been to the left. The behavior will be well-defined after the first
     *        extension to the right on the unidirectional index.
     * \returns Returns a unidirectional seqan3::fm_index_iterator on the index of the reversed text.
     *
     * Example:
     *
     * \snippet test/snippet/index/bi_fm_index_iterator.cpp get_rev_iterator
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    rev_iterator get_rev_iterator() const noexcept
    {
        assert(m_index != nullptr);

        rev_iterator it(m_index->rev_fm);
        it.m_parent_lb = m_parent_lb;
        it.m_parent_rb = m_parent_rb;
        it.m_node = {m_rev_lb, m_rev_rb, m_depth, m_last_char};

        #ifndef NDEBUG
            if (fwd_iter_last_used)
            {
                // invalidate parent range
                it.m_parent_lb = 1;
                it.m_parent_rb = 0;
            }
        #endif

        return it;
    }

    /*!\brief Returns the searched query.
     *        \cond DEV
     *            Returns the concatenation of all edges from the root node to the iterators current node.
     *        \endcond
     *
     * ### Complexity
     *
     * O(SAMPLING_RATE * T_BACKWARD_SEARCH) + query_length().
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    auto query() const noexcept
    {
        assert(m_index != nullptr && m_index->text != nullptr);

        size_type const query_begin = offset() - m_index->fwd_fm.m_index[m_fwd_lb];
        return *m_index->text | ranges::view::slice(query_begin, query_begin + query_length());
    }

    //!\copydoc query()
    auto operator*() const noexcept
    {
       assert(m_index != nullptr && m_index->text != nullptr);

       return query();
    }

    /*!\brief Counts the number of occurrences of the searched query in the text.
     * \returns Number of occurrences of the searched query in the text.
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
        assert(m_index != nullptr && 1 + m_fwd_rb - m_fwd_lb == 1 + m_rev_rb - m_rev_lb);

        return 1 + m_fwd_rb - m_fwd_lb;
    }

    /*!\brief Locates the occurrences of the searched query in the text.
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
        assert(m_index != nullptr);

        std::vector<size_type> occ(count());
        for (size_type i = 0; i < occ.size(); ++i)
        {
            occ[i] = offset() - m_index->fwd_fm.m_index[m_fwd_lb + i];
        }
        return occ;
    }

    /*!\brief Locates the occurrences of the searched query in the text on demand, i.e. a ranges::view is returned
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
        assert(m_index != nullptr);

        size_type const _offset = offset();
        return ranges::view::iota(m_fwd_lb, m_fwd_lb + count())
               | ranges::view::transform([*this, _offset] (auto sa_pos) {
                   return _offset - m_index->fwd_fm.m_index[sa_pos];
               });
    }

};

//!\}

} // namespace seqan3
