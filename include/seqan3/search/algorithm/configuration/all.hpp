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
 * \brief Meta-header for the \link search_configuration search configuration module \endlink.
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 */

#pragma once

/*!\defgroup search_configuration Configuration
 * \brief Data structures and utility functions for configuring search algorithm.
 * \ingroup search
 *
 * \todo Write detailed landing page.
 */

#include <seqan3/search/algorithm/configuration/max_error.hpp>
#include <seqan3/search/algorithm/configuration/max_error_rate.hpp>
#include <seqan3/search/algorithm/configuration/mode.hpp>
#include <seqan3/search/algorithm/configuration/output_index_iterator.hpp>
#include <seqan3/search/algorithm/configuration/output_text_position.hpp>
#include <seqan3/search/algorithm/configuration/utility.hpp>

/*!\namespace seqan3::search_cfg
 * \brief A special sub namespace for the search configurations.
 */
