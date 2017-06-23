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

#include <sdsl/suffix_trees.hpp>

namespace seqan3
{

template <typename index_t>
class index_trie_iter;

template <typename index_t>
class index_trie_iter_node
{
public:
    using index_type = index_t;
    using text_pos_type = uint32_t; // depends on the dimensions!
    using label_iterator_type = void; // TODO: view of range library
    friend class index_trie_iter<index_type>;

    index_trie_iter_node() = default;

    index_trie_iter_node(index_type const & _index): index(_index)
    {
        sdsl_node = index.sdsl_index.root();
    }

    label_iterator_type edge_label()
    {
        // ...
    }

    // label_iterator_type path_label()
    // {
    //     // ...
    // }

    typename index_type::text_size_type depth()
    {
        return index.sdsl_index.depth(sdsl_node);
    }

    typename index_type::text_size_type lb()
    {
        return index.sdsl_index.lb(sdsl_node);
    }

    typename index_type::text_size_type rb()
    {
        return index.sdsl_index.rb(sdsl_node);
    }

    // std::pair<typename index_t::text_size_type, typename index_t::text_size_type> range;
    // typename tree_iter_type::edge_label_type last_char; // TODO: rename
    // typename tree_iter_type::search_depth_size_type path_label_length;

    // bool operator!=(const tree_iter_node<tree_iter_type> & other)
    // {
    //     return
    //         this->range == other.range &&
    //         this->last_char == other.last_char &&
    //         this->path_label_length == other.path_label_length;
    // }

protected:
    index_type index;
    typename index_type::sdsl_index_type::node_type sdsl_node;
};

template <typename index_t>
class index_trie_iter
{
public:
    using index_type = index_t;
    using iter_node_type = typename index_type::iter_node_type;

    index_trie_iter(index_type const & _index) : index(_index), node(_index)
    {
        node.sdsl_node = index.sdsl_index.root();
    }

    // rule of six constructors
    index_trie_iter() = delete;
    index_trie_iter(index_trie_iter const &) = default;
    index_trie_iter & operator=(index_trie_iter const &) = default;
    index_trie_iter(index_trie_iter &&) = default;
    index_trie_iter & operator=(index_trie_iter &&) = default;

    bool go_down()
    {
        if (!is_leaf())
        {
            node.sdsl_node = index.sdsl_index.leftmost_leaf(node.sdsl_node);
            std::cout << "down: 1" << std::endl;
            return true;
        }
        std::cout << "down: 0" << std::endl;
        return false;
    }

    bool go_up()
    {
        if (!is_root())
        {
            std::cout << "up: 1" << std::endl;
            node.sdsl_node = index.sdsl_index.parent(node.sdsl_node);
            return true;
        }
        std::cout << "up: 0" << std::endl;
        return false;
    }

    bool go_down(typename index_type::char_type const & c)
    {
        node.sdsl_node = index.sdsl_index.child(node.sdsl_node, c);
        return true;
    }

    bool go_down(std::vector<typename index_type::char_type> const & pattern)
    {
        // ...
        return true;
    }

    bool go_right()
    {
        auto sibling = index.sdsl_index.sibling(node.sdsl_node);
        if (sibling != index.sdsl_index.root())
        {
            node.sdsl_node = sibling;
            std::cout << "right: 1" << std::endl;
            return true;
        }
        std::cout << "right: 0" << std::endl;
        return false;
    }

    bool is_leaf()
    {
        // ...
        return index.sdsl_index.is_leaf(node.sdsl_node);
    }

    bool is_root()
    {
        return node.sdsl_node == index.sdsl_index.root();
    }

    iter_node_type const & operator*()
    {
        return node;
    }

protected:
    index_type index;
    std::pair<typename index_type::text_size_type, typename index_type::text_size_type> parent_range;
    iter_node_type node;
};

} // namespace seqan3
