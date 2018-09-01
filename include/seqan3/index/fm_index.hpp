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
 * \brief Provides the unidirectional seqan3::fm_index.
 */

#pragma once

#include <sdsl/suffix_trees.hpp>

#include <range/v3/view/reverse.hpp>
#include <range/v3/view/transform.hpp>

#include <seqan3/index/concept.hpp>
#include <seqan3/index/fm_index_iterator.hpp>
#include <seqan3/index/detail/fm_index_iterator.hpp>
#include <seqan3/index/detail/csa_alphabet_strategy.hpp>
#include <seqan3/core/metafunction/range.hpp>

namespace seqan3
{

//!\cond
// forward declarations
template <typename index_t>
class fm_index_iterator;

template <typename index_t>
class bi_fm_index_iterator;
//!\endcond

/*!\addtogroup index
 * \{
 */

/*!\brief The default FM Index Configuration.
 * \ingroup fm_index
 *
 * \details
 *
 * ### Running time / Space consumption
 *
 * SAMPLING_RATE = 16
 * SIGMA: alphabet_size<char_type> where char_type is the seqan3 alphabet type (e.g. dna4 has an alphabet size of 4)
 * T_BACKWARD_SEARCH: O(log SIGMA)
 *
 * \todo Asymptotic space consumption:
 *
 */
struct fm_index_default_traits
{
    //!\brief Type of the underlying SDSL index.
    using sdsl_index_type = sdsl::csa_wt<
        sdsl::wt_blcd<
            sdsl::bit_vector,
            sdsl::rank_support_v<>,
            sdsl::select_support_scan<>,
            sdsl::select_support_scan<0>
        >,
        16,
        10000000,
        sdsl::sa_order_sa_sampling<>,
        sdsl::isa_sampling<>,
        sdsl::plain_byte_alphabet
    >;
};

/*!\brief The SeqAn FM Index.
 * \ingroup fm_index
 * \implements seqan3::fm_index_concept
 * \tparam text_t The type of the text to be indexed; must satisfy std::ranges::ForwardRange.
 * \tparam fm_index_traits The traits determining the implementation of the underlying SDSL index;
                           must satisfy seqan3::fm_index_traits_concept.
 * \details
 *
 * The seqan3::fm_index is a fast and space-efficient string index to search strings and collections of strings.
 *
 * ### General information
 *
 * Here is a short example on how to build an index and search a pattern using an iterator. Please note that there is a
 * very powerful search module with a high-level interface \todo seqan3::search that encapsulates the use of iterators.
 *
 * \snippet test/snippet/index/fm_index.cpp all
 *
 * Here is an example using a collection of strings (e.g. a genome with multiple chromosomes or a protein database):
 *
 * Coming soon. Stay tuned!
 *
 * ### Choosing an index implementation
 *
 * \todo The underlying implementation of the FM Index (Rank data structure, sampling rates, etc.) can be specified ...
 */
template <std::ranges::ForwardRange text_t, fm_index_traits_concept fm_index_traits = fm_index_default_traits>
//!\cond
    requires alphabet_concept<innermost_value_type_t<text_t>>
        && std::is_same_v<typename underlying_rank<innermost_value_type_t<text_t>>::type, uint8_t>
//!\endcond
class fm_index
{
protected:
    //!\privatesection

    //!\brief The type of the underlying SDSL index.
    using sdsl_index_type = typename fm_index_traits::sdsl_index_type;
    /*!\brief The type of the reduced alphabet type. (The reduced alphabet might be smaller than the original alphabet
     *        in case not all possible characters occur in the indexed text.)
     */
    using sdsl_char_type = typename sdsl_index_type::alphabet_type::char_type;

    //!\brief Underlying index from the SDSL.
    sdsl_index_type m_index;
    //!\brief Pointer to the indexed text.
    text_t const * text = nullptr;

    //!\publicsection

public:
    //!\brief The type of the indexed text.
    using text_type = text_t;
    //!\brief The type of the underlying character of text_type.
    using char_type = innermost_value_type_t<text_t>;
    //!\brief Type for representing positions in the indexed text.
    using size_type = typename sdsl_index_type::size_type;

    //!\brief The type of the (unidirectional) iterator.
    using iterator_type = fm_index_iterator<fm_index<text_t, fm_index_traits>>;

    template <typename bi_fm_index_t>
    friend class bi_fm_index_iterator;

    template <typename fm_index_t>
    friend class fm_index_iterator;

    template <typename fm_index_t>
    friend class detail::fm_index_iterator_node;

    /*!\name Constructors and destructor
     * \{
     */
    fm_index() = default;
    fm_index(fm_index const &) = default;
    fm_index & operator=(fm_index const &) = default;
    fm_index(fm_index &&) = default;
    fm_index & operator=(fm_index &&) = default;
    ~fm_index() = default;
    //!\}

    /*!\brief Constructor that immediately constructs the index given a range.
     *        The range cannot be an rvalue (i.e. a temporary object).
     * \tparam text_t The type of range to construct from; must satisfy std::ranges::ForwardRange.
     * \param[in] text The text to construct from.
     *
     * ### Complexity
     *
     * \todo At least linear.
     *
     * ### Exceptions
     *
     * No guarantees.
     */
    fm_index(text_t const & text)
    {
        construct(text);
    }

    //!\overload
    fm_index(text_t &&) = delete;

    //!\overload
    fm_index(text_t const &&) = delete;

    /*!\brief Constructs the index given a range. The range cannot be an rvalue (i.e. a temporary object).
     *        \todo Poorely implemented with regard to the memory peak due to not matching interfaces with the SDSL.
     * \tparam text_t The type of range to construct from; must satisfy std::ranges::ForwardRange.
     * \param[in] text The text to construct from.
     *
     * ### Complexity
     *
     * \todo At least linear.
     *
     * ### Exceptions
     *
     * No guarantees.
     */
    void construct(text_t const & text)
    {
        assert(text.begin() != text.end()); // text must not be empty
        this->text = &text;
        // TODO:
        // * check what happens in sdsl when constructed twice!
        // * choose between in-memory/external and construction algorithms
        // * sdsl construction currently only works for int_vector, std::string and char *, not ranges in general
        sdsl::int_vector<8> tmp_text(text.size());

        // uint8_t largest_char = 0;
        for (auto it = text.begin(); it != text.end(); it++)
        {
            // largest_char = std::max(largest_char, static_cast<uint8_t>(to_rank(*it) + 1));
            tmp_text[text.end() - it - 1] = to_rank(*it) + 1; // reverse and increase rank by one
        }
        sdsl::construct_im(m_index, tmp_text, 0);

        // TODO: would be nice but doesn't work since it's private and the public member references are const
        // m_index.m_C.resize(largest_char);
        // m_index.m_C.shrink_to_fit();
        // m_index.m_sigma = largest_char;
    }

    //!\overload
    void construct(text_t &&) = delete;

    //!\overload
    void construct(text_t const &&) = delete;

    /*!\brief Returns the length of the indexed text including sentinel characters.
     * \returns Returns the length of the indexed text including sentinel characters.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    size_type size() const noexcept
    {
        return m_index.size();
    }

    /*!\brief Checks whether the index is empty.
     * \returns `true` if the index is empty, `false` otherwise.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    bool empty() const noexcept
    {
        return size() == 0;
    }

    // operator== not implemented by sdsl indices yet
    // bool operator==(fm_index const & rhs) const noexcept
    // {
    //     return m_index == rhs.m_index;
    // }

    // operator== not implemented by sdsl indices yet
    // bool operator!=(fm_index const & rhs) const noexcept
    // {
    //     return !(*this == rhs);
    // }

    /*!\brief Returns a seqan3::fm_index_iterator on the index that can be used for searching.
     *        \cond DEV
     *            Iterator is pointing to the root node of the implicit suffix tree.
     *        \endcond
     * \returns Returns a (unidirectional) seqan3::fm_index_iterator on the index.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    iterator_type begin() const noexcept
    {
        return iterator_type(*this);
    }

    /*!\brief Loads the index from disk. Temporary function until cereal is supported.
     *        \todo cereal
     * \returns `true` if the index was successfully loaded from disk.
     *
     * ### Complexity
     *
     * Linear.
     *
     * ### Exceptions
     *
     * No guarantees.
     */
    bool load(std::string const & path)
    {
        return sdsl::load_from_file(m_index, path);
    }

    /*!\brief Stores the index to disk. Temporary function until cereal is supported.
     *        \todo cereal
     * \returns `true` if the index was successfully stored to disk.
     *
     * ### Complexity
     *
     * Linear.
     *
     * ### Exceptions
     *
     * No guarantees.
     */
    bool store(std::string const & path) const
    {
        return sdsl::store_to_file(m_index, path);
    }

};

//!\}

} // namespace seqan3
