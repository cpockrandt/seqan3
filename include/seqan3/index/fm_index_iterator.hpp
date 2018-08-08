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

// TODO: unify terminology (pattern, sequence)

/*!\brief The undirectional FM index iterator.
 * \ingroup fm_index
 * \implements seqan3::fm_index_iterator_concept
 * \tparam index_t The type of the underlying index; must satisfy seqan3::fm_index_concept.
 * \details
 *
 * The iterator's interface and behaviour is equivalent to a suffix tree with the space and time efficiency of the
 * underlying pure FM index. The implicit suffix tree is not compacted, i.e. going down an edge will increase the
 * searched sequence by only one character.
 *
 * The iterator traverses the implicit suffix tree beginning at the root node. All methods
 * modifying the iterator (e.g. going down an edge) return a `bool` value whether the operation was successful or not.
 * In case of an unsuccessful operation the iterator remains unmodified, i.e. an iterator can never be in an
 * invalid state except default constructed iterators that are always invalid. E.g., going down an edge while the
 * iterator is in a leaf node, will return `false` and will remain in that leaf node.
 *
 * The asymptotic running times for using the iterator depend on the SDSL index configuration. To determine the exact
 * running times, you have to additionally look up the running times of the used traits (configuration).
 */
template <typename index_t>
class fm_index_iterator
{

public:
    //!\brief Type of the index.
    using index_type = index_t;
    //!\brief Type for representing positions in the indexed text.
    using size_type = typename index_type::size_type;

protected:
    //!\privatesection

    //!\brief Type of the representation of a suffix tree node.
    using node_type = detail::fm_index_iterator_node<index_t>;

    //!\brief Type of the representation of characters in the underlying SDSL index.
    using sdsl_char_type = typename index_type::sdsl_char_type;

    //!\brief Type of the underlying FM index.
    index_type const * index;
    //!\brief Left suffix array range of the parent node. Needed for right().
    size_type parent_lb;
    //!\brief Right suffix array range of the parent node. Needed for right().
    size_type parent_rb;
    //!\brief Underlying index from the SDSL.
    node_type node;

    //!\brief Helper function to recompute text positions since the indexed text is reversed.
    size_type offset() const noexcept
    {
        return index->m_index.size() - depth() - 1; // since the string is reversed during construction
    }

    //!\publicsection

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
     * Constant. (\todo: will change when comparing indices, not yet supported by SDSL)
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    bool operator==(fm_index_iterator const & rhs) const noexcept
    {
        assert(index != nullptr);
        assert(node != rhs.node || (depth() == 0 || (parent_lb == rhs.parent_lb && parent_rb == rhs.parent_rb)));

        // position in the implicit suffix tree is defined by the SA range and depth. No need to compare parent ranges
        return node == rhs.node;
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
    bool operator!=(fm_index_iterator const & rhs) const noexcept
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
        assert(index != nullptr);

        sdsl_char_type c = 1; // NOTE: start with 0 or 1 depending on implicit_sentintel
        size_type _lb, _rb;
        while (c < index->m_index.sigma &&
               !sdsl::backward_search(index->m_index, node.lb, node.rb, index->m_index.comp2char[c], _lb, _rb))
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
        assert(index != nullptr);

        size_type _lb, _rb;

        sdsl_char_type c_char = to_rank(c) + 1;

        if (sdsl::backward_search(index->m_index, node.lb, node.rb, c_char, _lb, _rb))
        {
            parent_lb = node.lb;
            parent_rb = node.rb;
            node = {_lb, _rb, node.depth + 1, c_char};
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
        auto first = pattern.cbegin();
        auto last = pattern.cend();

        assert(index != nullptr && first != last); // range must not be empty!

        size_type _lb = node.lb, _rb = node.rb;
        size_type _parent_lb = node.lb, _parent_rb = node.rb;

        sdsl_char_type c;

        for (auto it = first; it != last; ++it)
        {
            c = to_rank(*it) + 1;

            _parent_lb = _lb;
            _parent_rb = _rb;
            if (!sdsl::backward_search(index->m_index, _parent_lb, _parent_rb, c, _lb, _rb))
                return false;
        }

        parent_lb = _parent_lb;
        parent_rb = _parent_rb;
        node = {_lb, _rb, node.depth + last - first, c};
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
        assert(index != nullptr && depth() > 0 && parent_lb <= parent_rb); // parent_lb > parent_rb --> invalid range

        sdsl_char_type c = node.last_char + 1;
        size_type _lb, _rb;

        while (c < index->m_index.sigma &&
               !sdsl::backward_search(index->m_index, parent_lb, parent_rb, index->m_index.comp2char[c], _lb, _rb))
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
    std::array<fm_index_iterator, alphabet_size_v<typename index_type::char_type>> children() const noexcept
    {
        assert(index != nullptr);

        std::array<fm_index_iterator, alphabet_size_v<typename index_type::char_type>> result;

        sdsl_char_type c = 1; // NOTE: start with 0 or 1 depending on implicit_sentintel
        size_type _lb, _rb;

        uint8_t i = 0;
        while (c < index->m_index.sigma)
        {
            // TODO: this will be implemented much more efficiently in the future (performed rank queries are almost the
            // same). Rank information for different characters are located in the same cache line.
            if (sdsl::backward_search(index->m_index, node.lb, node.rb, index->m_index.comp2char[c], _lb, _rb))
            {
                result[i].parent_lb = node.lb;
                result[i].parent_rb = node.rb;
                result[i].index = index;
                result[i].node = {_lb, _rb, node.depth + 1, c};
                ++i;
            }
            ++c;
        }

        // fill the rest with iterators pointing to root
        while (i < index_type::char_type::value_size)
        {
            result[i++] = fm_index_iterator(*this->index);
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
        assert(node.depth != 0 || (node.lb == 0 && node.rb == index->size() - 1)); // depth == 0 -> root node

        return node.depth;
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

        size_type const pattern_begin = offset() - index->m_index[node.lb];
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
        assert(index != nullptr);

        return 1 + node.rb - node.lb;
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
            occ[i] = offset() - index->m_index[node.lb + i];
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
        return ranges::view::iota(node.lb, node.lb + count())
               | ranges::view::transform([*this, _offset] (auto sa_pos) { return _offset - index->m_index[sa_pos]; });
    }

};

//!\}

} // namespace seqan3
