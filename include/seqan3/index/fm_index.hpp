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

#include <sdsl/suffix_trees.hpp>

namespace seqan3
{

template <typename char_t, uint8_t DIMENSIONS>
class fm_index
{
public:
    using char_type = char_t;
    using text_size_type = uint32_t;

    using tree_iter_type = index_tree_iter<fm_index<char_type, DIMENSIONS> >;
    using trie_iter_type = index_tree_iter<fm_index<char_type, DIMENSIONS> >;
    using iter_node_type = index_tree_iter_node<fm_index<char_type, DIMENSIONS> >;
    friend class index_tree_iter<fm_index<char_type, DIMENSIONS> >;
    friend class index_tree_iter_node<fm_index<char_type, DIMENSIONS> >;

    static const bool is_bidirectional = false;
    static const uint8_t dimensions = DIMENSIONS;

    // rule of six constructors
    fm_index() = default;
    fm_index(fm_index const &) = default;
    fm_index & operator=(fm_index const &) = default;
    fm_index(fm_index &&) = default;
    fm_index & operator=(fm_index &&) = default;

    tree_iter_type tree_root()
    {
        tree_iter_type iter(*this);
        return iter;
    }

    trie_iter_type trie_root()
    {
        trie_iter_type iter(*this);
        return iter;
    }

    template <typename text_type> // TODO: text container_type and lowest value_type = char_t? Hannes ...
    void construct(text_type & text)
    {
        sdsl::construct_im(sdsl_index, text, 1);
    }

protected:
    // typedef sdsl::wt_blcd<sdsl::bit_vector, sdsl::rank_support_v<>, sdsl::select_support_scan<>, sdsl::select_support_scan<0> > _wt_type;
    // typedef sdsl::csa_wt<_wt_type, 32, 1000000, sdsl::sa_order_sa_sampling<>, sdsl::isa_sampling<>/*, wt_alphabet_trait*/ > _csa_type;
    typedef sdsl::cst_sct3<> sdsl_index_type;
    sdsl_index_type sdsl_index;
};

} // namespace seqan3
