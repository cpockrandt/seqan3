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
 * \brief Provides the seqan3::fm_index_iterator for searching in the unidirectional seqan3::fm_index.
 */

#pragma once

#include <array>
#include <type_traits>

#include <sdsl/suffix_trees.hpp>

#include <range/v3/view/iota.hpp>
#include <range/v3/view/slice.hpp>

#include <seqan3/alphabet/all.hpp>
#include <seqan3/core/metafunction/range.hpp>
#include <seqan3/search/fm_index/detail/csa_alphabet_strategy.hpp>
#include <seqan3/search/fm_index/detail/fm_index_iterator.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/std/view/transform.hpp>

namespace seqan3
{

/*!\addtogroup submodule_fm_index
 * \{
 */

/*!\brief The SeqAn FM Index Iterator.
 * \implements seqan3::fm_index_iterator_concept
 * \tparam index_t The type of the underlying index; must model seqan3::fm_index_concept.
 * \details
 *
 * The iterator's interface provides searching a string from left to right in the indexed text.
 * All methods modifying the iterator (e.g. extending by a character with extend_right()) return a `bool` value whether
 * the operation was successful or not. In case of an unsuccessful operation the iterator remains unmodified, i.e. an
 * iterator can never be in an invalid state except default constructed iterators that are always invalid.
 *
 * \cond DEV
 *     The behaviour is equivalent to a suffix tree with the space and time efficiency of the underlying pure FM index.
 *     The iterator traverses the implicit suffix tree beginning at the root node.
 *     The implicit suffix tree is not compacted, i.e. going down an edge using extend_right(char) will increase the
 *     query by only one character.
 * \endcond
 *
 * The asymptotic running times for using the iterator depend on the SDSL index configuration. To determine the exact
 * running times, you have to additionally look up the running times of the used traits (configuration).
 */
template <typename index_t>
class fm_index_iterator
{

public:

    /*!\name Member types
     * \{
     */
    //!\brief Type of the index.
    using index_type = index_t;
    //!\brief Type for representing positions in the indexed text.
    using size_type = typename index_type::size_type;
    //!\}

protected:
    //!\privatesection

    /*!\name Member types
     * \{
     */
    //!\brief Type of the representation of a suffix tree node.
    using node_type = detail::fm_index_iterator_node<index_t>;
    //!\brief Type of the representation of characters in the underlying SDSL index.
    using sdsl_char_type = typename index_type::sdsl_char_type;
    //!\}

    //!\brief Type of the underlying FM index.
    index_type const * index;
    //!\brief Left suffix array interval of the parent node. Needed for right().
    size_type parent_lb;
    //!\brief Right suffix array interval of the parent node. Needed for right().
    size_type parent_rb;
    //!\brief Underlying index from the SDSL.
    node_type node;

    template <typename _index_t>
    friend class bi_fm_index_iterator;

    //!\brief Helper function to recompute text positions since the indexed text is reversed.
    size_type offset() const noexcept
    {
        return index->m_index.size() - query_length() - 1; // since the string is reversed during construction
    }

    //!\brief Optimized backward search without alphabet mapping
    template <detail::sdsl_index_concept csa_t>
    bool backward_search(csa_t const & csa, size_type const l, size_type const r, sdsl_char_type const c,
                         size_type & l_res, size_type & r_res) const noexcept
    {
        assert(l <= r && r < csa.size());

        if constexpr(std::Same<typename csa_t::alphabet_type, sdsl::plain_byte_alphabet>)
        {
            size_type const c_begin = csa.C[c];
            if (r + 1 - l == csa.size()) // [[unlikely]]
            {
                l_res = c_begin;
                r_res = csa.C[c + 1] - 1;
            }
            else
            {
                l_res = c_begin + csa.bwt.rank(l, c);		  // count c in bwt[0..l-1]
                r_res = c_begin + csa.bwt.rank(r + 1, c) - 1; // count c in bwt[0..r]
            }
            assert(r_res + 1 - l_res >= 0);
            return r_res >= l_res;
        }
        else
        {
            size_type const cc = csa.char2comp[c];
            if (cc == 0 && c > 0) // [[unlikely]]
                return false;

            size_type const c_begin = csa.C[cc];
            if (l == 0 && r + 1 == csa.size()) // [[unlikely]]
            {
                l_res = c_begin;
                r_res = csa.C[cc + 1] - 1;
                return true;
            }
            else
            {
                l_res = c_begin + csa.bwt.rank(l, c);		  // count c in bwt[0..l-1]
                r_res = c_begin + csa.bwt.rank(r + 1, c) - 1; // count c in bwt[0..r]
                assert(r_res + 1 - l_res >= 0);
                return r_res >= l_res;
            }
        }
    }

public:

    /*!\name Constructors and destructor
     * \{
     */
    //!\brief Default constructor. Accessing member functions on a default constructed object is undefined behavior.
    //        Default construction is necessary to make this class semi-regular and e.g., to allow construction of
    //        std::array of iterators.
    fm_index_iterator() noexcept = default;
    fm_index_iterator(fm_index_iterator const &) noexcept = default;
    fm_index_iterator & operator=(fm_index_iterator const &) noexcept = default;
    fm_index_iterator(fm_index_iterator &&) noexcept = default;
    fm_index_iterator & operator=(fm_index_iterator &&) noexcept = default;

    fm_index_iterator(index_t const & _index) noexcept : index(&_index), node({0, _index.m_index.size() - 1, 0, 0})
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
    bool operator==(fm_index_iterator const & rhs) const noexcept
    {
        assert(index != nullptr);
        assert(node != rhs.node ||
            (query_length() == 0 || (parent_lb == rhs.parent_lb && parent_rb == rhs.parent_rb)));

        // position in the implicit suffix tree is defined by the SA interval and depth.
        // No need to compare parent intervals
        return node == rhs.node;
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
    bool operator!=(fm_index_iterator const & rhs) const noexcept
    {
        assert(index != nullptr);

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
     * \f$O(\Sigma) * O(T_{BACKWARD\_SEARCH})\f$
     *
     * It scans linearly over the alphabet until it finds the smallest character that is represented by an edge.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    bool extend_right() noexcept
    {
        // TODO: specialize extend_right() and cycle_back() for EPR-dictionaries
        // store all iterators at once in a private std::array of iterators
        assert(index != nullptr);

        sdsl_char_type c = 1; // NOTE: start with 0 or 1 depending on implicit_sentintel
        size_type _lb, _rb;
        while (c < index->m_index.sigma &&
               !backward_search(index->m_index, node.lb, node.rb, index->m_index.comp2char[c], _lb, _rb))
        {
            ++c;
        }

        if (c != index->m_index.sigma)
        {
            parent_lb = node.lb;
            parent_rb = node.rb;
            node = {_lb, _rb, node.depth + 1, c};
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
     * \f$O(T_{BACKWARD\_SEARCH})\f$
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
        assert(index != nullptr);

        size_type _lb, _rb;

        sdsl_char_type c_char = to_rank(c) + 1;

        if (backward_search(index->m_index, node.lb, node.rb, c_char, _lb, _rb))
        {
            parent_lb = node.lb;
            parent_rb = node.rb;
            node = {_lb, _rb, node.depth + 1, c_char};
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
     * \f$|seq| * O(T_{BACKWARD\_SEARCH})\f$
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
        auto first = seq.begin();
        auto last = seq.end();

        assert(index != nullptr && first != last); // range must not be empty!

        size_type _lb = node.lb, _rb = node.rb;
        size_type _parent_lb = node.lb, _parent_rb = node.rb;

        sdsl_char_type c;

        for (auto it = first; it != last; ++it)
        {
            c = to_rank(*it) + 1;

            _parent_lb = _lb;
            _parent_rb = _rb;
            if (!backward_search(index->m_index, _parent_lb, _parent_rb, c, _lb, _rb))
                return false;
        }

        parent_lb = _parent_lb;
        parent_rb = _parent_rb;
        node = {_lb, _rb, last - first + node.depth, c};
        return true;
    }

    /*!\brief Tries to replace the rightmost character of the query by the next lexicographically larger character such
     *        that the query is found in the text.
     *        \cond DEV
     *            Moves the iterator to the right sibling of the current suffix tree node. It would be equivalent to
     *            going up an edge and going down that edge with the smallest character that is larger than the
     *            previous searched character. Calling cycle_back() on an iterator pointing to the root node is
     *            undefined behaviour!
     *        \endcond
     * \returns `true` if there exists a query in the text where the rightmost character of the query is
     *          lexicographically larger than the current rightmost character of the query.
     *
     * Example:
     *
     * \snippet test/snippet/index/fm_index_iterator.cpp cycle
     *
     * ### Complexity
     *
     * \f$O(\Sigma) * O(T_{BACKWARD\_SEARCH})\f$
     *
     * It scans linearly over the alphabet starting from the rightmost character until it finds the query with a larger
     * rightmost character.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    bool cycle_back() noexcept
    {
        // parent_lb > parent_rb --> invalid interval
        assert(index != nullptr && query_length() > 0 && parent_lb <= parent_rb);

        sdsl_char_type c = node.last_char + 1;
        size_type _lb, _rb;

        while (c < index->m_index.sigma &&
               !backward_search(index->m_index, parent_lb, parent_rb, index->m_index.comp2char[c], _lb, _rb))
        {
            ++c;
        }

        if (c != index->m_index.sigma)
        {
            node = {_lb, _rb, node.depth, c};
            return true;
        }
        return false;
    }

    /*!\brief Outputs the rightmost character.
     * \returns Rightmost character.
     *
     * Example:
     *
     * \snippet test/snippet/index/fm_index_iterator.cpp cycle
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
        // parent_lb > parent_rb --> invalid interval
        assert(index != nullptr && query_length() > 0 && parent_lb <= parent_rb);

        typename index_t::char_type c;
        c.assign_rank(index->m_index.comp2char[node.last_char] - 1); // text is not allowed to contain ranks of 0
        return c;
    }

    /*!\brief Returns the length of the searched query.
     *        \cond DEV
     *            Returns the depth of the iterator node in the implicit suffix tree.
     *        \endcond
     * \returns Length of query.
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
        assert(index != nullptr);
        assert(node.depth != 0 || (node.lb == 0 && node.rb == index->size() - 1)); // depth == 0 -> root node

        return node.depth;
    }

    /*!\brief Returns the searched query.
     *        \cond DEV
     *            Returns the concatenation of all edges from the root node to the iterators current node.
     *        \endcond
     * \returns Searched query.
     *
     * ### Complexity
     *
     * \f$O(SAMPLING\_RATE * T_{BACKWARD\_SEARCH}) + query\_length()\f$
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    auto query() const noexcept
    {
        assert(index != nullptr && index->text != nullptr);

        size_type const query_begin = offset() - index->m_index[node.lb];
        return *index->text | ranges::view::slice(query_begin, query_begin + query_length());
    }

    //!\copydoc query()
    auto operator*() const noexcept
    {
       assert(index != nullptr && index->text != nullptr);

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
        assert(index != nullptr);

        return 1 + node.rb - node.lb;
    }

    /*!\brief Locates the occurrences of the searched query in the text.
     * \returns Positions in the text.
     *
     * ### Complexity
     *
     * \f$count() * O(T_{BACKWARD\_SEARCH} * SAMPLING\_RATE)\f$
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
            occ[i] = offset() - index->m_index[node.lb + i];
        }
        return occ;
    }

    /*!\brief Locates the occurrences of the searched query in the text on demand, i.e. a ranges::view is returned and
     *        every position is located once it is accessed.
     * \returns Positions in the text.
     *
     * ### Complexity
     *
     * \f$count() * O(T_{BACKWARD\_SEARCH} * SAMPLING\_RATE)\f$
     *
     * ### Exceptions
     *
     * Strong exception guarantee (no data is modified in case an exception is thrown).
     */
    auto lazy_locate() const
    {
        assert(index != nullptr);

        size_type const _offset = offset();
        return ranges::view::iota(node.lb, node.lb + count())
               | view::transform([*this, _offset] (auto sa_pos) { return _offset - index->m_index[sa_pos]; });
    }

};

//!\}

} // namespace seqan3
