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
#include <seqan3/index/fm_index.hpp>
#include <seqan3/index/bi_fm_index_iterator.hpp>
#include <seqan3/core/metafunction/range.hpp>

namespace seqan3
{

/*!\addtogroup index
 * \{
 */

/*!\brief The default Bidirectional FM Index Configuration.
 * \ingroup bi_fm_index
 *
 * \details
 *
 * \todo write me
 *
 */
struct bi_fm_index_default_traits
{
    // TODO: seems a bit over-engineered to me and complicated for users. it would be nicer to have:
    // using sdsl_index_type = ...; and
    // using rev_sdsl_index_type = ...;
    // but we need unidirectional trait objects
    //!\brief Type of the underlying forward SDSL index.
    using fm_index_traits = fm_index_default_traits;

    //!\brief Type of the underlying reverse SDSL index.
    using rev_fm_index_traits = fm_index_default_traits; // trait object without sampling.
    // but not good for creating unidirectional iterators on the reverse index ...
};

/*!\brief The SeqAn Bidirectional FM Index.
 * \ingroup bi_fm_index
 * \implements seqan3::bi_fm_index_concept
 * \tparam text_t The type of the text to be indexed; must satisfy std::ranges::ForwardRange.
 * \tparam bi_fm_index_traits The traits determining the implementation of the underlying SDSL indices;
                              must satisfy seqan3::bi_fm_index_traits_concept.
 * \details
 *
 * \todo
 */
template <std::ranges::ForwardRange text_t, bi_fm_index_traits_concept bi_fm_index_traits = bi_fm_index_default_traits>
//!\cond
    requires alphabet_concept<innermost_value_type_t<text_t>>
        && std::is_same_v<typename underlying_rank<innermost_value_type_t<text_t>>::type, uint8_t>
//!\endcond
class bi_fm_index
{
protected:
    //!\privatesection
    //!\brief Pointer to the indexed text.
    text_t const * text = nullptr;

    //!\publicsection

public:
    //!\brief The type of the forward indexed text.
    using text_type = text_t;
    //!\brief The type of the forward indexed text.
    using rev_text_type = decltype(ranges::view::reverse(*text));


    using index_traits = bi_fm_index_traits;

protected:
    //!\privatesection

    //!\brief The type of the underlying forward SDSL index.
    using sdsl_index_type = typename bi_fm_index_traits::fm_index_traits::sdsl_index_type;

    //!\brief The type of the underlying reverse SDSL index.
    using rev_sdsl_index_type = typename bi_fm_index_traits::rev_fm_index_traits::sdsl_index_type;

    /*!\brief The type of the reduced alphabet type. (The reduced alphabet might be smaller than the original alphabet
     *        in case not all possible characters occur in the indexed text.)
     */
    using sdsl_char_type = typename sdsl_index_type::alphabet_type::char_type;

    //!\brief \todo
    rev_text_type rev_text;

    using fm_index_type = fm_index<text_t, typename bi_fm_index_traits::fm_index_traits>;
    using rev_fm_index_type = fm_index<rev_text_type, typename bi_fm_index_traits::rev_fm_index_traits>;

    fm_index_type fwd_fm;
    rev_fm_index_type rev_fm;

    //!\publicsection

public:
    //!\brief The type of the underlying character of text_type.
    using char_type = innermost_value_type_t<text_t>;
    //!\brief Type for representing positions in the indexed text.
    using size_type = typename sdsl_index_type::size_type;

    //!\brief The type of the (unidirectional) iterator.
    using iterator_type = bi_fm_index_iterator<bi_fm_index<text_t, bi_fm_index_traits>>;
    using fwd_iterator_type = fm_index_iterator<fm_index_type>;
    using rev_iterator_type = fm_index_iterator<rev_fm_index_type>;

    // friend class bi_fm_index_iterator<fm_index<text_t, bi_fm_index_traits>>;
    template <typename index_t>
    friend class bi_fm_index_iterator/*<bi_fm_index<text_t, bi_fm_index_traits>>*/;
    template <typename index_t>
    friend class fm_index_iterator/*<fm_index_type>*/;
    // friend class fm_index_iterator<rev_fm_index_type>;
    //friend class detail::fm_index_iterator_node<fm_index<text_t, fm_index_traits>>;

    // static constexpr bool is_bidirectional = false;
    // static constexpr uint8_t dimensions = DIMENSIONS;

    /*!\name Constructors and destructor
     * \{
     */
    bi_fm_index() = default;
    bi_fm_index(bi_fm_index const &) = default;
    bi_fm_index & operator=(bi_fm_index const &) = default;
    bi_fm_index(bi_fm_index &&) = default;
    bi_fm_index & operator=(bi_fm_index &&) = default;
    ~bi_fm_index() = default;
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
    bi_fm_index(text_t const & text)
    {
        construct(text);
    }

    //!\overload
    bi_fm_index(text_t &&) = delete;

    //!\overload
    bi_fm_index(text_t const &&) = delete;

    /*!\brief Constructs the index given a range. The range cannot be an rvalue (i.e. a temporary object).
     *        \todo Poorely implemented with regard to the memory peak due to not matching interfaces with the SDSL
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
        rev_text = ranges::view::reverse(text);
        fwd_fm.construct(text);
        rev_fm.construct(rev_text);
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
        return fwd_fm.size();
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

    /*!\brief Returns an iterator on the index.
     * \returns Returns a (unidirectional) iterator on the index pointing to the root node of the implicit suffix tree.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    iterator_type root() const noexcept
    {
        return iterator_type(*this);
    }

    fwd_iterator_type fwd_root() const noexcept
    {
       return fwd_iterator_type(fwd_fm);
    }

    rev_iterator_type rev_root() const noexcept
    {
       return rev_iterator_type(rev_fm);
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
        return fwd_fm.load(path + ".fwd") && rev_fm.load(path + ".rev");
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
        return fwd_fm.store(path + ".fwd") && rev_fm.store(path + ".rev");
    }

};

//!\}

} // namespace seqan3
