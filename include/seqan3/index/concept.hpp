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
 * \brief Provides the concepts for the seqan3::fm_index and its traits and iterators.
 */

#pragma once

#include <type_traits>

#include <sdsl/suffix_arrays.hpp>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/metafunction/range.hpp>
#include <seqan3/range/concept.hpp>

namespace seqan3
{

/*!\addtogroup index
 * \{
 */

/*!\interface seqan3::fm_index_traits_concept <>
 * \brief Concept for FM Index traits.
 *
 * The traits object must contain an index type of the SDSL namespace.
 * \todo:
 * * this concept will be documented later once it is finalised.
 * * sdsl index requirements (maybe move them into a separate concept)
 */
//!\cond
template <typename t>
concept fm_index_traits_concept = requires (t v, typename t::sdsl_index_type::size_type lb,
                                                      typename t::sdsl_index_type::size_type rb,
                                                      typename t::sdsl_index_type::char_type c,
                                                      typename t::sdsl_index_type m_index)
{
    typename t::sdsl_index_type;

    // sdsl index requirements:
    typename t::sdsl_index_type::size_type;
    (typename t::sdsl_index_type{}).size();
    (typename t::sdsl_index_type{})[0]; // suffix array access
    (typename t::sdsl_index_type{}).comp2char[0];
    (typename t::sdsl_index_type{}).char2comp[0];
    (typename t::sdsl_index_type{}).sigma;
    (typename t::sdsl_index_type{}).C[0];
    (typename t::sdsl_index_type{}).bwt.rank(lb, c);
    (typename t::sdsl_index_type{}).wavelet_tree.lex_count(lb, rb, c);
    { sdsl::construct_im(m_index, sdsl::int_vector<8> {}, 0) } -> void;
};
//!\endcond

/*!\name Requirements for seqan3::fm_index_traits_concept
 * \relates seqan3::fm_index_traits_concept
 * \brief The SDSL index must support the following interface to work with SeqAn3.
 * \{
 */

/*!\typedef typename t::sdsl_index_type sdsl_index_type
 * \memberof seqan3::fm_index_traits_concept
 * \brief Declares the type of the underlying SDSL index.
 */

//!\}



/*!\interface seqan3::fm_index_concept <>
 * \brief Concept for FM indices.
 *
 * This concept defines the interface for (unidirectional) FM indices.
 */
//!\cond
template <typename t>
concept fm_index_concept = std::Semiregular<t> && requires (t v)
{
    typename t::text_type;
    typename t::char_type;
    typename t::size_type;
    typename t::iterator_type;

    // NOTE: circular dependency
    // requires fm_index_iterator_concept<typename t::iterator_type>;

    requires requires (t index, std::vector<dna4> const & text) { { t(text) } };
    requires requires (t index, std::vector<dna4> const & text) { { index.construct(text) } -> void; };

    { v.root() } -> typename t::iterator_type;

    { v.size()  } -> typename t::size_type;
    { v.empty() } -> bool;

    { v.load(std::string{})  } -> bool;
    { v.store(std::string{}) } -> bool;
};
//!\endcond

/*!\name Requirements for seqan3::fm_index_concept
 * \relates seqan3::fm_index_concept
 * \brief You can expect these member types and member functions on all types that satisfy seqan3::fm_index_concept.
 * \{
 */

/*!\typedef typename t::char_type char_type
 * \memberof seqan3::fm_index_concept
 * \brief Type of the underlying character of text_type.
 */

/*!\typedef typename t::size_type size_type
 * \memberof seqan3::fm_index_concept
 * \brief Type for representing positions in the indexed text.
 */

/*!\typedef typename t::iterator_type iterator_type
 * \memberof seqan3::fm_index_concept
 * \brief Type of the iterator.
 */

//!\todo Write me!

//!\}



/*!\interface seqan3::fm_index_iterator_concept <>
 * \brief Concept for FM index iterators.
 *
 * This concept defines the interface for iterators for (unidirectional) FM indices.
 */
//!\cond
template <typename t>
concept fm_index_iterator_concept = std::Semiregular<t> && requires (t it)
{
    typename t::index_type;
    typename t::size_type;

    requires fm_index_concept<typename t::index_type>;

    requires requires (typename t::index_type const & index) { { t(index) } };

    { it.down()                                                 } -> bool;
    { it.down(typename t::index_type::char_type{})              } -> bool;
    { it.down(std::vector<typename t::index_type::char_type>{}) } -> bool;
    { it.right()                                                } -> bool;

    { it.children() } -> std::array<t, alphabet_size_v<typename t::index_type::char_type>>;

    { it.depth()       } -> typename t::size_type;
    { it.path_label()  } -> auto;
    { *it              } -> auto;
    { it.count()       } -> typename t::size_type;
    { it.locate()      } -> std::vector<typename t::size_type>;
    { it.lazy_locate() } -> auto;
};
//!\endcond

/*!\name Requirements for seqan3::fm_index_iterator_concept
 * \relates seqan3::fm_index_iterator_concept
 * \brief You can expect these member types and member functions on all types that satisfy seqan3::fm_index_iterator_concept.
 * \{
 */

/*!\typedef typename t::index_type index_type
 * \memberof seqan3::fm_index_iterator_concept
 * \brief Type of the underlying SeqAn FM index wrapper (not the underlying SDSL index).
 */

/*!\typedef typename t::size_type size_type
 * \memberof seqan3::fm_index_iterator_concept
 * \brief Type for representing positions in the indexed text.
 */

//!\todo Write me!

//!\}




template <typename t>
concept bi_fm_index_traits_concept = requires (t v)
{
    requires fm_index_traits_concept<typename t::fm_index_traits>;
    requires fm_index_traits_concept<typename t::rev_fm_index_traits>;

    requires std::is_same_v<typename t::fm_index_traits::sdsl_index_type::size_type,
                            typename t::rev_fm_index_traits::sdsl_index_type::size_type>;
};



template <typename t>
concept bi_fm_index_concept = std::Semiregular<t> && requires (t v)
{
    typename t::text_type;
    typename t::char_type;
    typename t::size_type;
    typename t::iterator_type;
    typename t::fwd_iterator_type;
    typename t::rev_iterator_type;

    // NOTE: circular dependency
    // requires bi_fm_index_iterator_concept<typename t::iterator_type>;

    // (bool)t::is_bidirectional;

    requires requires (t index, std::vector<dna4> const & text) { { t(text) } };
    requires requires (t index, std::vector<dna4> const & text) { { index.construct(text) } -> void; };

    { v.root() } -> typename t::iterator_type;
    { v.fwd_root() } -> typename t::fwd_iterator_type;
    { v.rev_root() } -> typename t::rev_iterator_type;

    { v.size()  } -> typename t::size_type;
    { v.empty() } -> bool;

    { v.load(std::string{})  } -> bool;
    { v.store(std::string{}) } -> bool;
};

template <typename t>
concept bi_fm_index_iterator_concept = fm_index_iterator_concept<t> && requires (t it)
{
    requires bi_fm_index_concept<typename t::index_type>;

    requires requires (typename t::index_type const & index) { { t(index) } };

    { it.down_rev()                                                 } -> bool;
    { it.down_rev(typename t::index_type::char_type{})              } -> bool;
    { it.down_rev(std::vector<typename t::index_type::char_type>{}) } -> bool;
    { it.right_rev()                                                } -> bool;

    { it.children_rev() } -> std::array<t, alphabet_size_v<typename t::index_type::char_type>>;

    // TODO: should we offer _rev() methods?
    // { it.path_label()  } -> auto;
    // { it.locate()      } -> std::vector<typename t::size_type>;
    // { it.lazy_locate() } -> auto;
};






//!\}

} // namespace seqan3
