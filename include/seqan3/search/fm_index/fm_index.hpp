// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 * \brief Provides the unidirectional seqan3::fm_index.
 */

#pragma once

#include <sdsl/suffix_trees.hpp>

#include <seqan3/core/type_traits/range.hpp>
#include <seqan3/std/filesystem>
#include <seqan3/range/shortcuts.hpp>
#include <seqan3/range/view/to_rank.hpp>
#include <seqan3/search/fm_index/concept.hpp>
#include <seqan3/search/fm_index/detail/csa_alphabet_strategy.hpp>
#include <seqan3/search/fm_index/detail/fm_index_cursor.hpp>
#include <seqan3/search/fm_index/fm_index_cursor.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/ranges>

namespace seqan3
{

//!\cond
// forward declarations
template <typename index_t>
class fm_index_cursor;

template <typename index_t>
class bi_fm_index_cursor;
//!\endcond

/*!\addtogroup submodule_fm_index
 * \{
 */

/*!\brief The FM Index Configuration using a Wavelet Tree.
 *
 * \details
 *
 * ### Running time / Space consumption
 *
 * \f$SAMPLING\_RATE = 16\f$
 * \f$\Sigma\f$: alphabet_size<char_type> where char_type is the seqan3 alphabet type (e.g. dna4 has an alphabet size
 *               of 4).
 *
 * For an index over a text collection a delimiter is added inbetween the texts. This causes sigma to increase by 1.
 * \attention For any alphabet, the symbol with rank 255 is not allowed to occur in the text. Addtionally,
 *            rank 254 cannot occur when indexing text collections.
 *
 * \if DEV
 * Rank 255 is generally not allowed. When constructing the index, we increase the rank of every letter by 1, since the
 * SDSL uses 0 as a sentinel. Incrementing 255 by 1 will cause the rank (uint8_t) to overflow to 0.
 * For text collections, 255 is reserved as delimiter between the individual texts. Therefore the letter with rank 254
 * cannot occur in the text.
 *
 * This index will only work for byte alphabets, i.e. alphabets with a size <= 256.
 * When switching to bigger alphabets, the static_asserts, the delimiter choice and the transform view for the
 * construction need to be adjusted. Additionally, all occurrences of uint8_t should be double checked to make sure
 * they also apply to bigger alphabets.
 * \endif
 *
 * \f$T_{BACKWARD\_SEARCH}: O(\log \Sigma)\f$
 *
 * \if DEV \todo Asymptotic space consumption: \endif
 *
 */
using sdsl_wt_index_type =
    sdsl::csa_wt<sdsl::wt_blcd<sdsl::bit_vector,
                               sdsl::rank_support_v<>,
                               sdsl::select_support_scan<>,
                               sdsl::select_support_scan<0>>,
                 16,
                 10000000,
                 sdsl::sa_order_sa_sampling<>,
                 sdsl::isa_sampling<>,
                 sdsl::plain_byte_alphabet>;

using sdsl_epr_index_type =
 sdsl::csa_wt<sdsl::wt_pc_epr<11>,
              16,
              10000000,
              sdsl::sa_order_sa_sampling<>,
              sdsl::isa_sampling<>,
              sdsl::plain_byte_alphabet>;

/*!\brief The default FM Index Configuration.
 * \attention The default might be changed in a future release. If you rely on a stable API and on-disk-format,
 *            please hard-code your sdsl_index_type to a concrete type.
 */
using default_sdsl_index_type = sdsl_wt_index_type;

//!\brief The possible text layouts (single, collection) the seqan3::fm_index and seqan3::bi_fm_index can support.
enum text_layout : bool
{
    //!\brief The text is a single range.
    single,
    //!\brief The text is a range of ranges.
    collection
};

//!\cond
SEQAN3_DEPRECATED_310
void fm_index_deprecation(bool);

template <typename t>
void fm_index_deprecation(t);
//!\endcond

/*!\brief The SeqAn FM Index.
 * \implements seqan3::FmIndex
 * \tparam is_collection    Indicates whether this index works on a text collection or a single text.
 *                          See seqan3::text_layout.
 * \tparam sdsl_index_type_ The type of the underlying SDSL index, must model seqan3::SdslIndex.
 * \details
 *
 * The seqan3::fm_index is a fast and space-efficient string index to search strings and collections of strings.
 *
 * ### General information
 *
 * Here is a short example on how to build an index and search a pattern using an cursor. Please note that there is a
 * very powerful search module with a high-level interface seqan3::search that encapsulates the use of cursors.
 *
 * \include test/snippet/search/fm_index.cpp
 *
 * \attention When building an index for a **single text** over any alphabet, the symbol with rank 255 is reserved
 *            and may not occur in the text.
 *
 * Here is an example using a collection of strings (e.g. a genome with multiple chromosomes or a protein database):
 *
 * \include test/snippet/search/fm_index_collection.cpp
 *
 * \attention When building an index for a **text collection** over any alphabet, the symbols with rank 254 and 255
 *            are reserved and may not be used in the text.
 *
 * \if DEV
 * ### Choosing an index implementation
 *
 * The underlying implementation of the FM Index (rank data structure, sampling rates, etc.) can be specified by
 * passing a new SDSL index type as second template parameter:
 *
 * \todo Link to SDSL documentation or write our own once SDSL3 documentation is available somewhere....
 *
 * \endif
 *
 * \deprecated Use seqan3::text_layout to indicate single texts and text collections. The use of bool is deprecated.
 */
template <auto is_collection = text_layout::single, detail::SdslIndex sdsl_index_type_ = default_sdsl_index_type>
class fm_index
{
protected:
    //!\brief The alphabet size of the text.
    size_t sigma{0};
    //!\brief Indicates whether index is built over a collection.
    static constexpr bool is_collection_{is_collection};

    /*!\name Member types
     * \{
     */
    //!\brief The type of the underlying SDSL index.
    using sdsl_index_type = sdsl_index_type_;
    /*!\brief The type of the reduced alphabet type. (The reduced alphabet might be smaller than the original alphabet
     *        in case not all possible characters occur in the indexed text.)
     */
    using sdsl_char_type = typename sdsl_index_type::alphabet_type::char_type;
    //!\brief The type of the alphabet size of the underlying SDSL index.
    using sdsl_sigma_type = typename sdsl_index_type::alphabet_type::sigma_type;
    //!\}

    //!\brief Underlying index from the SDSL.
    sdsl_index_type index;

    //!\brief Bitvector storing begin positions for collections.
    sdsl::sd_vector<> text_begin;
    //!\brief Select support for text_begin.
    sdsl::select_support_sd<1> text_begin_ss;
    //!\brief Rank support for text_begin.
    sdsl::rank_support_sd<1> text_begin_rs;

    //!\cond
    using unused_t [[maybe_unused]] = decltype(fm_index_deprecation(is_collection));
    //!\endcond

public:
    /*!\name Member types
     * \{
     */
    //!\brief Type for representing positions in the indexed text.
    using size_type = typename sdsl_index_type::size_type;
    //!\brief The type of the (unidirectional) cursor.
    using cursor_type = fm_index_cursor<fm_index<is_collection, sdsl_index_type>>;
    //!\}

    template <typename bi_fm_index_t>
    friend class bi_fm_index_cursor;

    template <typename fm_index_t>
    friend class fm_index_cursor;

    template <typename fm_index_t>
    friend class detail::fm_index_cursor_node;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    fm_index() = default;              //!< Defaulted.

    fm_index(fm_index const & rhs) :   //!< When copy constructing, also update internal data structures.
        sigma{rhs.sigma}, index{rhs.index}, text_begin{rhs.text_begin}, text_begin_ss{rhs.text_begin_ss},
        text_begin_rs{rhs.text_begin_rs}
    {
        text_begin_ss.set_vector(&text_begin);
        text_begin_rs.set_vector(&text_begin);
    }

    fm_index(fm_index && rhs) :        //!< When move constructing, also update internal data structures.
        sigma{std::move(rhs.sigma)}, index{std::move(rhs.index)}, text_begin{std::move(rhs.text_begin)},
        text_begin_ss{std::move(rhs.text_begin_ss)}, text_begin_rs{std::move(rhs.text_begin_rs)}
    {
        text_begin_ss.set_vector(&text_begin);
        text_begin_rs.set_vector(&text_begin);
    }

    fm_index & operator=(fm_index rhs) //!< When copy/move assigning, also update internal data structures.
    {
        index = std::move(rhs.index);
        sigma = std::move(rhs.sigma);
        text_begin = std::move(rhs.text_begin);
        text_begin_ss = std::move(rhs.text_begin_ss);
        text_begin_rs = std::move(rhs.text_begin_rs);

        text_begin_ss.set_vector(&text_begin);
        text_begin_rs.set_vector(&text_begin);

        return *this;
    }

    ~fm_index() = default;             //!< Defaulted.

    /*!\brief Constructor that immediately constructs the index given a range. The range cannot be empty.
     * \tparam text_t The type of range to construct from; must model std::ranges::BidirectionalRange.
     * \param[in] text The text to construct from.
     *
     * ### Complexity
     *
     * \if DEV \todo \endif At least linear.
     */
    template <std::ranges::Range text_t>
    fm_index(text_t && text)
    {
        construct(std::forward<text_t>(text));
    }
    //!\}

    /*!\brief Constructs the index given a range.
              The range cannot be an rvalue (i.e. a temporary object) and has to be non-empty.
     * \tparam text_t The type of range to construct from; must model std::ranges::BidirectionalRange.
     * \param[in] text The text to construct from.
     *
     * \details
     * \if DEV
     * \todo This has to be better implemented with regard to the memory peak due to not matching interfaces
     *       with the SDSL.
     * \endif
     *
     * ### Complexity
     *
     * \if DEV \todo \endif At least linear.
     *
     * ### Exceptions
     *
     * No guarantee. \if DEV \todo Ensure strong exception guarantee. \endif
     */
    template <std::ranges::Range text_t>
    void construct(text_t && text)
        //!\cond
        requires !is_collection_
        //!\endcond
    {
        static_assert(std::ranges::BidirectionalRange<text_t>, "The text must model BidirectionalRange.");
        static_assert(alphabet_size<innermost_value_type_t<text_t>> <= 256, "The alphabet is too big.");
        static_assert(dimension_v<text_t> == 1, "The input cannot be a text collection.");

        // text must not be empty
        if (std::ranges::begin(text) == std::ranges::end(text))
            throw std::invalid_argument("The text that is indexed cannot be empty.");

        constexpr auto cexpr_sigma = alphabet_size<innermost_value_type_t<text_t>>;
        sigma = cexpr_sigma;
        // TODO:
        // * check what happens in sdsl when constructed twice!
        // * choose between in-memory/external and construction algorithms
        // * sdsl construction currently only works for int_vector, std::string and char *, not ranges in general
        // uint8_t largest_char = 0;
        sdsl::int_vector<8> tmp_text(text.size());

        std::ranges::copy(text
                          | view::to_rank
                          | std::view::transform([] (uint8_t const r)
                          {
                              if constexpr (cexpr_sigma == 256)
                              {
                                  if (r == 255)
                                      throw std::out_of_range("The input text cannot be indexed, because for full"
                                                              "character alphabets the last one/two values are reserved"
                                                              "(single sequence/collection).");
                              }
                              return r + 1;
                          })
                          | std::view::reverse,
                          seqan3::begin(tmp_text)); // reverse and increase rank by one

        sdsl::construct_im(index, tmp_text, 0);

        // TODO: would be nice but doesn't work since it's private and the public member references are const
        // index.m_C.resize(largest_char);
        // index.m_C.shrink_to_fit();
        // index.m_sigma = largest_char;
    }

    //!\overload
    template <std::ranges::Range text_t>
    void construct(text_t && text)
        //!\cond
        requires is_collection_
        //!\endcond
    {
        static_assert(std::ranges::BidirectionalRange<text_t>, "The text collection must model BidirectionalRange.");
        static_assert(std::ranges::BidirectionalRange<reference_t<text_t>>,
                      "The elements of the text collection must model BidirectionalRange.");
        static_assert(alphabet_size<innermost_value_type_t<text_t>> <= 256, "The alphabet is too big.");
        static_assert(dimension_v<text_t> == 2, "The input must be a text collection.");

        // text collection must not be empty
        if (std::ranges::begin(text) == std::ranges::end(text))
            throw std::invalid_argument("The text that is indexed cannot be empty.");

        size_t text_size{0}; // text size including delimiters

        // there must be at least one non-empty text in the collection
        bool all_empty = true;

        for (auto && t : text)
        {
            if (std::ranges::begin(t) != std::ranges::end(t))
            {
                all_empty = false;
            }
            text_size += 1 + t.size(); // text size and delimiter (sum will be 1 for empty texts)
        }

        if (all_empty)
            throw std::invalid_argument("A text collection that only contains empty texts cannot be indexed.");

        constexpr auto cexpr_sigma = alphabet_size<innermost_value_type_t<text_t>>;
        sigma = cexpr_sigma;

        // bitvector where 1 marks the begin position of a single text from the collection in the concatenated text
        sdsl::bit_vector pos(text_size, 0);
        size_t prefix_sum{0};

        for (auto && t : text)
        {
            pos[prefix_sum] = 1;
            prefix_sum += t.size() + 1;
        }

        text_begin    = sdsl::sd_vector(pos);
        text_begin_ss = sdsl::select_support_sd<1>(&text_begin);
        text_begin_rs = sdsl::rank_support_sd<1>(&text_begin);

        sdsl::int_vector<8> tmp_text(text_size - 1); // last text in collection needs no delimiter

        constexpr uint8_t delimiter = cexpr_sigma >= 255 ? 255 : cexpr_sigma + 1;

        std::vector<uint8_t> tmp = text
                                   | view::deep{view::to_rank}
                                   | view::deep
                                   {
                                       std::view::transform([] (uint8_t const r)
                                       {
                                           if constexpr (cexpr_sigma >= 255)
                                           {
                                               if (r >= 254)
                                                   throw std::out_of_range("The input text cannot be indexed, because"
                                                                           " for full character alphabets the last one/"
                                                                           "two values are reserved (single sequence/"
                                                                           "collection).");
                                           }
                                           return r + 1;
                                       })
                                   }
                                   | std::view::join(delimiter);

        std::ranges::copy((tmp | std::view::reverse), seqan3::begin(tmp_text));

        //!\if DEV \todo Replace with this once this does not cause debug builds to exceed max memory on travis \endif
        // std::ranges::copy(text
        //                   | view::deep{view::to_rank}
        //                   | view::deep{std::view::transform([] (uint8_t const r) { return r + 1; })} // increase rank
        //                   | view::deep{std::view::reverse}
        //                   | std::view::reverse
        //                   | std::view::join(delimiter), // join with delimiter
        //                   seqan3::begin(tmp_text));

        sdsl::construct_im(index, tmp_text, 0);
    }

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
        return index.size();
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

    /*!\brief Compares two indices.
     * \returns `true` if the indices are equal, false otherwise.
     *
     * ### Complexity
     *
     * Linear.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    bool operator==(fm_index const & rhs) const noexcept
    {
        // (void) rhs;
        return (index == rhs.index) && (text_begin == rhs.text_begin);
    }

    /*!\brief Compares two indices.
     * \returns `true` if the indices are unequal, false otherwise.
     *
     * ### Complexity
     *
     * Linear.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    bool operator!=(fm_index const & rhs) const noexcept
    {
        return !(*this == rhs);
    }

    /*!\brief Returns a seqan3::fm_index_cursor on the index that can be used for searching.
     *        \if DEV
     *            Cursor is pointing to the root node of the implicit suffix tree.
     *        \endif
     * \returns Returns a (unidirectional) seqan3::fm_index_cursor on the index.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    cursor_type begin() const noexcept
    {
        return {*this};
    }

    /*!\cond DEV
     * \brief Serialisation support function.
     * \tparam archive_t Type of `archive`; must satisfy seqan3::CerealArchive.
     * \param archive The archive being serialised from/to.
     *
     * \attention These functions are never called directly, see \ref serialisation for more details.
     */
    template <CerealArchive archive_t>
    void CEREAL_SERIALIZE_FUNCTION_NAME(archive_t & archive)
    {
        archive(index);
        archive(text_begin);
        archive(text_begin_ss);
        text_begin_ss.set_vector(&text_begin);
        archive(text_begin_rs);
        text_begin_rs.set_vector(&text_begin);
        archive(sigma);
        bool tmp = is_collection_;
        archive(tmp);
        assert(tmp == is_collection_);
    }
    //!\endcond

};

/*!\name Template argument type deduction guides
 * \{
 */
//! \brief Deduces the dimensions of the text.
template <std::ranges::Range text_t>
fm_index(text_t &&) -> fm_index<text_layout{dimension_v<text_t> != 1}>;
//!\}

//!\}
} // namespace seqan3
