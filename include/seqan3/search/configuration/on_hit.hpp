// // ============================================================================
// //                 SeqAn - The Library for Sequence Analysis
// // ============================================================================
// //
// // Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
// // Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
// // All rights reserved.
// //
// // Redistribution and use in source and binary forms, with or without
// // modification, are permitted provided that the following conditions are met:
// //
// //     * Redistributions of source code must retain the above copyright
// //       notice, this list of conditions and the following disclaimer.
// //     * Redistributions in binary form must reproduce the above copyright
// //       notice, this list of conditions and the following disclaimer in the
// //       documentation and/or other materials provided with the distribution.
// //     * Neither the name of Knut Reinert or the FU Berlin nor the names of
// //       its contributors may be used to endorse or promote products derived
// //       from this software without specific prior written permission.
// //
// // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// // AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// // IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// // ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// // FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// // DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// // SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// // CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// // LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// // OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// // DAMAGE.
// //
// // ============================================================================
//
// /*!\file
//  * \brief Provides gap configurations.
//  * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
//  */
//
// #pragma once
//
// #include <seqan3/search/configuration/utility.hpp>
// #include <seqan3/core/algorithm/all.hpp>
// #include <seqan3/core/metafunction/basic.hpp>
// #include <seqan3/core/metafunction/template_inspection.hpp>
//
// namespace seqan3::detail
// {
// /*!\brief A configuration element for the maximum number of errors.
//  * \ingroup configuration
//  */
// template <typename delegate_t>
// struct search_config_on_hit
// {
//     //!\brief The actual value.
//     delegate_t value;
// };
//
// /*!\brief The on_hit adaptor enabling pipe notation.
//  * \ingroup configuration
//  */
// template <template <typename ...> typename delegate_t>
// struct search_config_on_hit_adaptor : public configuration_fn_base<search_config_on_hit_adaptor<delegate_t>>
// {
//
//     /*!\brief Adds to the configuration a on_hit configuration element.
//      * \param[in] cfg  The configuration to be extended.
//      * \param[in] nbr The number of maximum errors used for the algorithm. (TODO: maximum vs maximal?)
//      * \returns A new configuration containing the on_hit configuration element.
//      */
//     template <typename configuration_t, typename delegate2_t>
//     //!\cond
//         requires is_algorithm_configuration_v<remove_cvref_t<configuration_t>>
//     //!\endcond
//     constexpr auto invoke(configuration_t && cfg, delegate2_t const nbr) const
//     {
//         static_assert(is_valid_search_configuration_v<search_cfg::id::on_hit, remove_cvref_t<configuration_t>>,
//                       SEQAN3_INVALID_CONFIG(search_cfg::id::on_hit));
//
//         search_config_on_hit<delegate2_t> tmp{nbr};
//         return std::forward<configuration_t>(cfg).push_front(std::move(tmp));
//     }
// };
//
// //!\brief Helper template meta-function associated with detail::search_config_on_hit.
// //!\ingroup configuration
// template <>
// struct on_search_config<search_cfg::id::on_hit>
// {
//     //!\brief Type alias used by meta::find_if
//     template <config_element_concept t>
//     using invoke = typename is_type_specialisation_of<t, search_config_on_hit>::type;
// };
//
// //!\brief Mapping from the detail::search_config_on_hit type to it's corresponding seqan3::search_cfg::id.
// //!\ingroup configuration
// template <>
// struct search_config_type_to_id<search_config_on_hit>
// {
//     //!\brief The associated seqan3::search_cfg::id.
//     static constexpr search_cfg::id value = search_cfg::id::on_hit;
// };
// } // namespace seqan3::detail
//
// namespace seqan3::search_cfg
// {
// /*!\brief A configuration adaptor for linear on_hit.
//  * \ingroup configuration
//  */
// inline constexpr detail::search_config_on_hit_adaptor on_hit; // TODO: all singular
//
// // inline constexpr detail::search_config_gap_adaptor<seqan3::gap_affine> gap_affine;
// } // namespace seqan3::search_cfg
