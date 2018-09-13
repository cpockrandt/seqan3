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
 * \brief Provides the concepts for seqan3::fm_index and seqan3::bi_fm_index and its traits and iterators.
 */

#pragma once

#include <type_traits>

#include <sdsl/suffix_arrays.hpp>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/metafunction/range.hpp>
#include <seqan3/range/concept.hpp>

namespace seqan3::detail
{

/*!\addtogroup submodule_fm_index
 * \{
 */

 // ============================================================================
 //  sdsl_index_concept
 // ============================================================================

/*!\interface seqan3::detail::sdsl_index_concept <>
 * \brief Concept for SDSL FM indices (which are called compressed suffix arrays in the SDSL).
 */
//!\cond
template <typename t>
concept sdsl_index_concept = requires (t sdsl_index)
{
    typename t::size_type;

    { sdsl_index.size() } -> typename t::size_type;
    { sdsl_index[0] }; // suffix array access
    { sdsl_index.comp2char[0] } -> uint8_t;
    { sdsl_index.char2comp[0] } -> uint8_t;
    { sdsl_index.sigma };
    { sdsl_index.C[0] };

    requires requires (t sdsl_index, typename t::char_type const c, typename t::size_type const lb,
                                     typename t::size_type const rb, sdsl::int_vector<8> const text)
    {
        { sdsl_index.bwt.rank(lb, c) };
        { sdsl_index.wavelet_tree.lex_count(lb, rb, c) };
        { sdsl::construct_im(sdsl_index, text, 0) };
    };
};
//!\endcond
/*!\name Requirements for seqan3::detail::sdsl_index_concept
 * \relates seqan3::detail::sdsl_index_concept
 * \brief The SDSL index must support the following interface to work with SeqAn3 FM indices.
 * \{
 *
 * \typedef typename t::size_type size_type
 * \memberof seqan3::detail::sdsl_index_concept
 * \brief Type for representing the size of the indexed text.
 *
 * \todo Write me.
 *
 * \}
 */

//!\}

} // namespace seqan3::detail

namespace seqan3
{

/*!\addtogroup submodule_fm_index
 * \{
 */

// ============================================================================
//  fm_index_traits_concept
// ============================================================================

/*!\interface seqan3::fm_index_traits_concept <>
 * \brief Concept for unidirectional FM Index traits.
 *
 * The traits object must contain an index type of the SDSL namespace.
 */
//!\cond
template <typename t>
concept fm_index_traits_concept = requires (t v)
{
    typename t::sdsl_index_type;

    requires detail::sdsl_index_concept<typename t::sdsl_index_type>;
};
//!\endcond
/*!\name Requirements for seqan3::fm_index_traits_concept
 * \relates seqan3::fm_index_traits_concept
 * \brief The SDSL index must support the following interface to work with SeqAn3.
 * \{
 *
 * \typedef typename t::sdsl_index_type sdsl_index_type
 * \memberof seqan3::fm_index_traits_concept
 * \brief Declares the type of the underlying SDSL index. Must satisfy the seqan3::detail::sdsl_index_concept.
 *
 * \}
 */

// ============================================================================
//  fm_index_concept
// ============================================================================

/*!\interface seqan3::fm_index_concept <>
 * \brief Concept for unidirectional FM indices.
 *
 * This concept defines the interface for unidirectional FM indices.
 */
//!\cond
template <typename t>
concept fm_index_concept = std::Semiregular<t> && requires (t index)
{
    typename t::text_type;
    typename t::char_type;
    typename t::size_type;
    typename t::iterator_type;

    // NOTE: circular dependency
    // requires fm_index_iterator_concept<typename t::iterator_type>;

    requires requires (t index, std::vector<dna4> const text)
    {
        { t(text) };
        { index.construct(text) } -> void;
    };

    { index.begin() } -> typename t::iterator_type;

    { index.size()  } -> typename t::size_type;
    { index.empty() } -> bool;

    { index.load(std::string{})  } -> bool;
    { index.store(std::string{}) } -> bool;
};
//!\endcond
/*!\name Requirements for seqan3::fm_index_concept
 * \relates seqan3::fm_index_concept
 * \brief You can expect these member types and member functions on all types that satisfy seqan3::fm_index_concept.
 * \{
 *
 * \typedef typename t::text_type text_type
 * \memberof seqan3::fm_index_concept
 * \brief Type of the indexed text.
 *
 * \typedef typename t::char_type char_type
 * \memberof seqan3::fm_index_concept
 * \brief Type of the underlying character of text_type.
 *
 * \typedef typename t::size_type size_type
 * \memberof seqan3::fm_index_concept
 * \brief Type for representing the size of the indexed text.
 *
 * \typedef typename t::iterator_type iterator_type
 * \memberof seqan3::fm_index_concept
 * \brief Type of the unidirectional FM index iterator.
 *
 * \todo Write me!
 *
 * \}
 */

// ============================================================================
//  fm_index_iterator_concept
// ============================================================================

/*!\interface seqan3::fm_index_iterator_concept <>
 * \brief Concept for unidirectional FM index iterators.
 *
 * This concept defines the interface for iterators for unidirectional FM indices.
 */
//!\cond
template <typename t>
concept fm_index_iterator_concept = std::Semiregular<t> && requires (t it)
{
    typename t::index_type;
    typename t::size_type;

    requires fm_index_concept<typename t::index_type>;

    requires requires (typename t::index_type const index) { { t(index) } };

    requires requires (t it, typename t::index_type::char_type const c,
                             std::vector<typename t::index_type::char_type> const seq)
    {
        { it.extend_right()    } -> bool;
        { it.extend_right(c)   } -> bool;
        { it.extend_right(seq) } -> bool;
        { it.cycle_back()      } -> bool;
    };

    { it.last_char()    } -> typename t::index_type::char_type;
    { it.query_length() } -> typename t::size_type;
    { it.query()        } -> auto;
    { *it               } -> auto;
    { it.count()        } -> typename t::size_type;
    { it.locate()       } -> std::vector<typename t::size_type>;
    { it.lazy_locate()  } -> auto;
};
//!\endcond
/*!\name Requirements for seqan3::fm_index_iterator_concept
 * \relates seqan3::fm_index_iterator_concept
 * \brief You can expect these member types and member functions on all types that satisfy seqan3::fm_index_iterator_concept.
 * \{
 *
 * \typedef typename t::index_type index_type
 * \memberof seqan3::fm_index_iterator_concept
 * \brief Type of the underlying SeqAn FM index (not the underlying SDSL index).
 *
 * \typedef typename t::size_type size_type
 * \memberof seqan3::fm_index_iterator_concept
 * \brief Type for representing the size of the indexed text.
 *
 * \todo Write me!
 *
 * \}
 */

// ============================================================================
//  bi_fm_index_traits_concept
// ============================================================================

/*!\interface seqan3::bi_fm_index_traits_concept <>
 * \brief Concept for bidirectional FM Index traits.
 *
 * The traits object must contain two unidirectional FM Index traits.
 */
//!\cond
template <typename t>
concept bi_fm_index_traits_concept = requires (t v)
{
    requires fm_index_traits_concept<typename t::fm_index_traits>;
    requires fm_index_traits_concept<typename t::rev_fm_index_traits>;

    requires std::Same<typename t::fm_index_traits::sdsl_index_type::size_type,
                       typename t::rev_fm_index_traits::sdsl_index_type::size_type>;
};
//!\endcond
/*!\name Requirements for seqan3::bi_fm_index_traits_concept
 * \relates seqan3::bi_fm_index_traits_concept
 * \brief The bidirectional FM index traits must provide the following types:
 * \{
 *
 * \typedef typename t::fm_index_traits fm_index_traits
 * \memberof seqan3::bi_fm_index_traits_concept
 * \brief Declares the type of the underlying unidirectional FM index on the original text.
 *        Must satisfy seqan3::fm_index_traits_concept.
 *
 * \typedef typename t::rev_fm_index_traits rev_fm_index_traits
 * \memberof seqan3::bi_fm_index_traits_concept
 * \brief Declares the type of the underlying unidirectional FM index on the reversed text.
 *        Must satisfy seqan3::fm_index_traits_concept.
 *
 * \}
 */

// ============================================================================
//  bi_fm_index_concept
// ============================================================================

/*!\interface seqan3::bi_fm_index_concept <>
 * \brief Concept for bidirectional FM indices.
 *
 * This concept defines the interface for bidirectional FM indices.
 */
//!\cond
template <typename t>
concept bi_fm_index_concept = fm_index_concept<t> && requires (t index)
{
    typename t::iterator_type; // already required by fm_index_concept but has a different documentation
    typename t::fwd_iterator_type;
    typename t::rev_iterator_type;

    // NOTE: circular dependency
    // requires bi_fm_index_iterator_concept<typename t::iterator_type>;

    { index.fwd_begin() } -> typename t::fwd_iterator_type;
    { index.rev_begin() } -> typename t::rev_iterator_type;
};
//!\endcond
/*!\name Requirements for seqan3::bi_fm_index_concept
 * \relates seqan3::bi_fm_index_concept
 * \brief You can expect these member types and member functions on all types that satisfy seqan3::bi_fm_index_concept.
 * \{
 *
 * \typedef typename t::iterator_type iterator_type
 * \memberof seqan3::bi_fm_index_concept
 * \brief Type of the bidirectional FM index iterator.
 *
 * \typedef typename t::fwd_iterator_type fwd_iterator_type
 * \memberof seqan3::bi_fm_index_concept
 * \brief Type of the unidirectional FM index iterator based on the unidirectional FM index on the original text.
 *
 * \typedef typename t::rev_iterator_type rev_iterator_type
 * \memberof seqan3::bi_fm_index_concept
 * \brief Type of the unidirectional FM index iterator based on the unidirectional FM index on the reversed text.
 *
 * \todo Write me!
 *
 * \}
 */

// ============================================================================
//  bi_fm_index_iterator_concept
// ============================================================================

/*!\interface seqan3::bi_fm_index_iterator_concept <>
 * \brief Concept for bidirectional FM index iterators.
 *
 * This concept defines the interface for iterators for bidirectional FM indices.
 */
//!\cond
template <typename t>
concept bi_fm_index_iterator_concept = fm_index_iterator_concept<t> && requires (t it)
{
    requires bi_fm_index_concept<typename t::index_type>;

    requires requires (typename t::index_type const index) { { t(index) } };

    requires requires (t it, typename t::index_type::char_type const c,
                             std::vector<typename t::index_type::char_type> const seq)
    {
        { it.extend_left()    } -> bool;
        { it.extend_left(c)   } -> bool;
        { it.extend_left(seq) } -> bool;
        { it.cycle_front()    } -> bool;
    };

};
//!\endcond
/*!\name Requirements for seqan3::bi_fm_index_iterator_concept
 * \relates seqan3::bi_fm_index_iterator_concept
 * \brief You can expect these member types and member functions on all types that satisfy
 *        seqan3::fm_index_iterator_concept.
 * \{
 *
 * \typedef typename t::index_type index_type
 * \memberof seqan3::bi_fm_index_iterator_concept
 * \brief Type of the underlying SeqAn FM index (not the underlying SDSL index).
 *
 * \typedef typename t::size_type size_type
 * \memberof seqan3::bi_fm_index_iterator_concept
 * \brief Type for representing the size of the indexed text.
 *
 * \todo Write me!
 *
 * \}
 */

//!\}

} // namespace seqan3
