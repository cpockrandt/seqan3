// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2017, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2017, Knut Reinert & MPI Molekulare Genetik
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
// Author: Christopher Pockrandt <github-noreply AT cpockrandt.de>
// ============================================================================

#pragma once

#include <vector>

namespace seqan3
{

template <typename t>
concept bool index_tree_iter_node_concept = requires (t n)
{
    typename t::index_type; // TODO

    typename t::text_pos_type;
    typename t::label_iterator_type; // view of range library?

    // requires index_tree_iter_concept<typename t::iter_type>;

    { n.edge_label()    } -> typename t::label_iterator_type;
    // { n.path_label()     } -> typename t::label_iterator_type;
    { n.depth()          } -> typename t::index_type::text_size_type;
    // { n.left_boundary()  } -> typename t::index_type::text_size_type;
    // { n.right_boundary() } -> typename t::index_type::text_size_type;

    /*
    left bounary, right boundary, depth,
    global: count_occurences
    */

    // { it.count()  } -> typename t::text_size_type; // TODO: rename. don't use verbs since there is no effect
    // { it.locate() } -> std::vector<typename t::text_pos_type>;
    // ...
};

template <typename t>
concept bool index_tree_iter_concept = requires (t it)
{
    typename t::index_type;
    // typename t::iter_node_type;

    //requires index_tree_concept<typename t::index_type>;

    // TODO: constructors
    // requires requires (t it, typename t::index_type const & index) { it(index) };

    { it.go_down()                                                 } -> bool; //t &;
    { it.go_down(typename t::index_type::char_type{})              } -> bool; //t &;
    { it.go_down(std::vector<typename t::index_type::char_type>{}) } -> bool; //t &;
    { it.go_right()                                                } -> bool; //t &;
    { it.is_leaf()                                                 } -> bool;
    { it.is_root()                                                 } -> bool;
    // { ++it } -> t&;
    // { it++ } -> t;
    { *it  } -> typename t::iter_node_type const &;
};

template <typename t>
concept bool bi_index_tree_iter_concept = requires (t it)
{
    requires index_tree_iter_concept<t>;
    requires t::is_bidirectional;

    { it.go_down_inv()    } -> bool;
    { it.go_right_inv()   } -> bool;

    // { it.edge_label_inv() } -> typename t::edge_label_type;
    // { it.is_leaf_inv()    } -> bool;
};

template <typename t>
concept bool index_tree_concept = requires(t tree)
{
    typename t::tree_iter_type;
    typename t::trie_iter_type;

    typename t::iter_node_type;

    typename t::char_type;
    typename t::text_size_type;

    // requires alphabet_concept<typename t::char_type>;
    requires index_tree_iter_concept<typename t::tree_iter_type>;
    requires index_tree_iter_concept<typename t::trie_iter_type>;

    requires index_tree_iter_node_concept<typename t::iter_node_type>;

    (uint8_t)t::dimensions; // number of nested containers
    (bool)t::is_bidirectional;

    // TODO: constructors?

    { tree.tree_root() } -> typename t::tree_iter_type;
    { tree.trie_root() } -> typename t::trie_iter_type;

    // template <typename construction_policy_t>
    // { tree.construct(typename t::text_type const &, construction_policy_t const &) } -> void;  // throws exception
    // TODO
    // requires requires (t tree, typename t::text_type const & text) { { tree.construct(text) } -> void; }; // throws exception
};

} // namespace seqan3

#ifndef NDEBUG

#include <fm_index.hpp>

static_assert(index_tree_concept<fm_index<char, 1> >);
// static_assert(index_tree_iter_concept<index_tree_iter<fm_index<char, 1> > >);

#endif
