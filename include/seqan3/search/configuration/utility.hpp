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
 * \brief Provides functionality to access get function by enum values.
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 */

#pragma once

#include <meta/meta.hpp>

#include <seqan3/core/algorithm/configuration.hpp>
#include <seqan3/core/detail/reflection.hpp>

namespace seqan3::search_cfg
{

/*!\brief Specifies an id for every configuration element.
 * \ingroup configuration
 *
 * \details
 *
 * The seqan3::search_cfg::id is used to identify a specific search configuration element independent of
 * it's concrete type and position within the seqan3::search_cfg::search_configuration object.
 * Thus one can access the value of the corresponding configuration element via the special get interface.
 *
 * ### Example
 *
 * ```cpp
 * search_cfg::search_configuration cfg = search_cfg::max_total_errors(3);
 * auto max_total_errors = get<search_cfg::id::max_total_errors>(cfg);  // max_total_errors = 3;
 * ```
 */
enum struct id : uint8_t
{
    //!\brief Identifier for max_errors configuration.
    max_deletion_error,
    max_deletion_error_rate,
    max_insertion_error,
    max_insertion_error_rate,
    max_substitution_error,
    max_substitution_error_rate,
    max_total_error,
    max_total_error_rate,
    on_hit,
    return_index_iterator,
    return_text_position,
    strategy_all,
    strategy_all_best,
    strategy_best,
    strategy_strata,
    //!\cond
    // ATTENTION: Must always be the last item; will be used to determine the number of ids.
    SIZE
    //!\endcond
};

} // namespace seqan3::search_cfg

namespace seqan3::detail
{
// ----------------------------------------------------------------------------
// SEQAN3_INVALID_CONFIG
// ----------------------------------------------------------------------------

//!\cond
// A common error message for testing valid configuration combinations.
#define SEQAN3_INVALID_CONFIG(e) "Configuration error: The configuration <" #e "> is not combinable with "  \
                                 "one of the previous config elements. Please see the documentation to get "\
                                 "more information about which configurations can be combined."
//!\endcond

// ----------------------------------------------------------------------------
// on_search_config
// ----------------------------------------------------------------------------

/*!\brief Checks if a specific type corresponds to the given seqan3::search_cfg::id.
 * \ingroup configuration
 * \tparam e The id to test for.
 *
 * \details
 *
 * This class is an invocable trait meta-function which is used in combination with the meta::find_if type function.
 * It is used to find the configuration element type within the detail::configuration, which associated with the
 * seqan3::search_cfg::id given as non-type template parameter.
 * To use this, implementors need to overload this class with the specific seqan3::search_cfg::id and provide
 * a type alias called `invoke`, which resolves to a std::bool_constant (e.g. std::true_type).
 *
 * ### Example
 *
 * ```cpp
 * template <>
 * struct on_search_config<search_cfg::id::max_total_errors>
 * {
 *     template <typename t>
 *     using invoke = typename std::is_same<t, search_config_max_total_errors>::type;
 * };
 * ```
 *
 * In this example `search_config_max_total_errors` is a class representing the configuration element for maximum numbers
 * of total errors (i.e. across all error types).
 */
template <search_cfg::id e>
struct on_search_config
{
    /*!\brief A template alias for the invocation. Defaults to std::false_type.
     * \tparam t The type of the configuration element to test.
     */
    template <config_element_concept t>
    using invoke = std::false_type;
};

// ----------------------------------------------------------------------------
// search_config_type_to_id
// ----------------------------------------------------------------------------

/*!\brief Maps the given configuration element type to it's associated seqan3::search_cfg::id.
 * \ingroup configuration
 * \tparam config_element_t The type to get the mapped seqan3::search_cfg::id for.
 * \see seqan3::detail::search_config_type_to_id_v
 */
template <config_element_concept config_element_t>
struct search_config_type_to_id
{
    //!\brief The mapped seqan3::search_cfg::id. Defaults to seqan3::search_cfg::id::SIZE.
    static constexpr search_cfg::id value = search_cfg::id::SIZE;
};

//!\brief Helper variable template for seqan3::detail::search_config_type_to_id.
//!\ingroup configuration
template <config_element_concept config_element_t>
inline constexpr search_cfg::id search_config_type_to_id_v = search_config_type_to_id<config_element_t>::value;

// ----------------------------------------------------------------------------
// search_config_validation_matrix
// ----------------------------------------------------------------------------

/*!\brief Validation matrix to check how search configuration elements can be combined.
 * \ingroup configuration
 *
 * \details
 *
 * This matrix is used to check if the specified search configurations can be combined with each other.
 * A cell value `true`, indicates that the corresponding seqan3::search_cfg::id in the current column can be combined
 * with the associated seqan3::search_cfg::id in the current row. The size of the matrix is determined by the enum value
 * SIZE of seqan3::search_cfg::id.
 */
inline constexpr std::array<std::array<bool, static_cast<uint8_t>(search_cfg::id::SIZE)>,
                            static_cast<uint8_t>(search_cfg::id::SIZE)> search_config_validation_matrix =
{
    {
        // max_deletion_error, max_deletion_error_rate, max_insertion_error, max_insertion_error_rate,
        // max_substitution_error, max_substitution_error_rate, max_total_error, max_total_error_rate, on_hit,
        // return_index_iterator, return_text_position, strategy_all, strategy_all_best, strategy_best, strategy_strata
        { 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1 },
        { 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0 }
    }
};

/*!\brief Determines the first type in reverse order of the given detail::configuration that is not combinable with
 *        the given configuration identified by seqan3::search_cfg::id.
 * \ingroup configuration
 *
 * \tparam query           The seqan3::search_cfg::id to test for valid combination.
 * \tparam configuration_t The detail::configuration type.
 *
 * \returns The first type that is not combinable with the specified seqan3::search_cfg::id.
 * \see
 */
template <search_cfg::id query, typename configuration_t>
//!\cond
    requires is_algorithm_configuration_v<remove_cvref_t<configuration_t>>
//!\endcond
struct invalid_search_configuration
{
protected:

    /*!\brief A type trait meta-function that checks if the given type is compatible with `query`.
     * \tparam target_t The type to check.
     */
    template <typename target_t>
    struct invoke
    {
        //!\brief Column index of `query`.
        static constexpr uint8_t pos1 = static_cast<uint8_t>(query);
        //!\brief Row index of the seqan3::search_cfg::id associated with `target_t`.
        static constexpr uint8_t pos2 = static_cast<uint8_t>(search_config_type_to_id_v<remove_cvref_t<target_t>>);

        static_assert(pos1 != static_cast<uint8_t>(search_cfg::id::SIZE),
                      "Unknown search_cfg::id! "
                      "Did you use an unknown config or did you forget to specialize search_config_type_to_id?");
        static_assert(pos2 != static_cast<uint8_t>(search_cfg::id::SIZE),
                      "Unknown search_cfg::id! "
                      "Did you use an unknown config or did you forget to specialize search_config_type_to_id?");

        //!\brief Type alias resolving to `std::false_type` if `target_t` can be combined with query, `std::true_type` otherwise.
        using type = std::conditional_t<search_config_validation_matrix[pos1][pos2], std::false_type, std::true_type>;
    };

    //!\brief seqan3::type_list alias for the types contained in `configuration_t`.
    using target_list_t = tuple_type_list_t<typename remove_cvref_t<configuration_t>::base_type>;
    //!\brief seqan3::type_list alias for tail list returned by meta::find_if.
    using tail_list_t   = meta::find_if<meta::reverse<target_list_t>, meta::quote_trait<invoke>>;
    //!\brief Defines a list of void element if the tail_list_t is empty.
    using final_list_t  = std::conditional_t<meta::empty<tail_list_t>::value, type_list<void>, tail_list_t>;

public:

    //!\brief The type alias for the type that is not compatible with query or void, if no conflicts are detected.
    using type = meta::front<final_list_t>;
};

//!\brief Helper template definition for seqan3::detail::invalid_search_configuration.
//!\ingroup configuration
template <search_cfg::id query, typename configuration_t>
//!\cond
    requires is_algorithm_configuration_v<remove_cvref_t<configuration_t>>
//!\endcond
using invalid_search_configuration_t = typename invalid_search_configuration<query, configuration_t>::type;

/*!\brief Returns whether the given `query` is compatible with `configuration_t`.
 * \ingroup configuration
 *
 * \tparam query           The seqan3::search_cfg::id to check.
 * \tparam configuration_t The seqan3::detail::configuration type to test against.
 *
 * \returns std::true_type, if `query` is compatible, `std::false_type` otherwise.
 * \see seqan3::detail::is_valid_search_configuration_v
 */
template <search_cfg::id query, typename configuration_t>
//!\cond
    requires is_algorithm_configuration_v<remove_cvref_t<configuration_t>>
//!\endcond
struct is_valid_search_configuration :
    public std::is_same<invalid_search_configuration_t<query, configuration_t>, void>
{};

//!\brief Helper variable template for seqan3::detail::is_valid_search_configuration.
//!\ingroup configuration
template <search_cfg::id query, typename configuration_t>
//!\cond
    requires is_algorithm_configuration_v<remove_cvref_t<configuration_t>>
//!\endcond
inline constexpr bool is_valid_search_configuration_v =
    is_valid_search_configuration<query, configuration_t>::value;

} // namespace seqan3::detail

namespace seqan3
{

// TODO: my own implementation to check whether a cfg element is in the configuration. This will be replaced by rrahns implementation

//!\cond
template <search_cfg::id e, typename ... cfg_elements_t>
constexpr bool contains(detail::configuration<cfg_elements_t...> & cfg) noexcept
{
    using type_list_t = detail::tuple_type_list_t<typename detail::configuration<cfg_elements_t...>::base_type>;

    constexpr size_t pos = meta::size<type_list_t>::value -
                           meta::size<meta::find_if<type_list_t, detail::on_search_config<e>>>::value;

    return pos < meta::size<type_list_t>::value;
}

template <search_cfg::id e, typename ... cfg_elements_t>
constexpr bool contains(detail::configuration<cfg_elements_t...> const & cfg) noexcept
{
    using type_list_t = detail::tuple_type_list_t<typename detail::configuration<cfg_elements_t...>::base_type>;

    constexpr size_t pos = meta::size<type_list_t>::value -
                           meta::size<meta::find_if<type_list_t, detail::on_search_config<e>>>::value;

    return pos < meta::size<type_list_t>::value;
}

template <search_cfg::id e, typename ... cfg_elements_t>
constexpr bool contains(detail::configuration<cfg_elements_t...> && cfg) noexcept
{
    using type_list_t = detail::tuple_type_list_t<typename detail::configuration<cfg_elements_t...>::base_type>;

    constexpr size_t pos = meta::size<type_list_t>::value -
                           meta::size<meta::find_if<type_list_t, detail::on_search_config<e>>>::value;

    return pos < meta::size<type_list_t>::value;
}

template <search_cfg::id e, typename ... cfg_elements_t>
constexpr bool contains(detail::configuration<cfg_elements_t...> const && cfg) noexcept
{
    using type_list_t = detail::tuple_type_list_t<typename detail::configuration<cfg_elements_t...>::base_type>;

    constexpr size_t pos = meta::size<type_list_t>::value -
                           meta::size<meta::find_if<type_list_t, detail::on_search_config<e>>>::value;

    return pos < meta::size<type_list_t>::value;
}
//!\endcond

/*!\name Enum-based get-interface
 * \ingroup configuration
 * \brief Provides a special overload for the get interface using seqan3::search_cfg::id as identifier.
 * \{
 */

/*!\brief Access the value of the search configuration element identified by seqan3::search_cfg::id.
 * \relates seqan3::detail::configuration
 *
 * \tparam e                 The seqan3::search_cfg::id to get the stored configuration value for.
 * \tparam cfg_elements_t... The configuration elements stored in detail::search_configuration.
 *
 * \param[in] cfg            The configuration for which to get the configuration element.
 * \returns The stored configuration value associated with seqan3::search_cfg::id `e`.
 */
template <search_cfg::id e, typename ... cfg_elements_t>
constexpr auto & get(detail::configuration<cfg_elements_t...> & cfg) noexcept
{
    using type_list_t = detail::tuple_type_list_t<typename detail::configuration<cfg_elements_t...>::base_type>;

    constexpr size_t pos = meta::size<type_list_t>::value -
                           meta::size<meta::find_if<type_list_t, detail::on_search_config<e>>>::value;

    static_assert(static_cast<int8_t>(e) < static_cast<int8_t>(search_cfg::id::SIZE),
                  "Did you forget to update search_cfg::id::SIZE value?");
    static_assert(pos < meta::size<type_list_t>::value,
                  "The specified config element is not contained in the configuration.");

    return get<pos>(cfg).value;
}

/*!\brief Access the value of the search configuration element identified by seqan3::search_cfg::id.
 * \relates seqan3::detail::configuration
 *
 * \tparam e                 The seqan3::search_cfg::id to get the stored configuration value for.
 * \tparam cfg_elements_t... The configuration elements stored in detail::search_configuration.
 *
 * \param[in] cfg            The configuration for which to get the configuration element.
 * \returns The stored configuration value associated with seqan3::search_cfg::id `e`.
 */
template <search_cfg::id e, typename ... cfg_elements_t>
constexpr auto const & get(detail::configuration<cfg_elements_t...> const & cfg) noexcept
{
    using type_list_t = detail::tuple_type_list_t<typename detail::configuration<cfg_elements_t...>::base_type>;
    constexpr size_t pos = meta::size<type_list_t>::value -
                           meta::size<meta::find_if<type_list_t, detail::on_search_config<e>>>::value;

    static_assert(static_cast<int8_t>(e) < static_cast<int8_t>(search_cfg::id::SIZE),
                  "Did you forget to update search_cfg::id::SIZE value?");
    static_assert(pos < meta::size<type_list_t>::value,
                  "The specified config element is not contained in the configuration.");

    return get<pos>(cfg).value;
}

/*!\brief Access the value of the search configuration element identified by seqan3::search_cfg::id.
 * \relates seqan3::detail::configuration
 *
 * \tparam e                 The seqan3::search_cfg::id to get the stored configuration value for.
 * \tparam cfg_elements_t... The configuration elements stored in detail::search_configuration.
 *
 * \param[in] cfg            The configuration for which to get the configuration element.
 * \returns The stored configuration value associated with seqan3::search_cfg::id `e`.
 */
template <search_cfg::id e, typename ... cfg_elements_t>
constexpr auto && get(detail::configuration<cfg_elements_t...> && cfg) noexcept
{
    using type_list_t = detail::tuple_type_list_t<typename detail::configuration<cfg_elements_t...>::base_type>;
    constexpr size_t pos = meta::size<type_list_t>::value -
                           meta::size<meta::find_if<type_list_t, detail::on_search_config<e>>>::value;

    static_assert(static_cast<int8_t>(e) < static_cast<int8_t>(search_cfg::id::SIZE),
                  "Did you forget to update search_cfg::id::SIZE value?");
    static_assert(pos < meta::size<type_list_t>::value,
                  "The specified config element is not contained in the configuration.");

    return get<pos>(std::move(cfg)).value;
}

/*!\brief Access the value of the search configuration element identified by seqan3::search_cfg::id.
 * \relates seqan3::detail::configuration
 *
 * \tparam e                 The seqan3::search_cfg::id to get the stored configuration value for.
 * \tparam cfg_elements_t... The configuration elements stored in detail::search_configuration.
 *
 * \param[in] cfg            The configuration for which to get the configuration element.
 * \returns The stored configuration value associated with seqan3::search_cfg::id `e`.
 */
template <search_cfg::id e, typename ... cfg_elements_t>
constexpr auto const && get(detail::configuration<cfg_elements_t...> const && cfg) noexcept
{
    using type_list_t = detail::tuple_type_list_t<typename detail::configuration<cfg_elements_t...>::base_type>;
    constexpr size_t pos = meta::size<type_list_t>::value -
                           meta::size<meta::find_if<type_list_t, detail::on_search_config<e>>>::value;

    static_assert(static_cast<int8_t>(e) < static_cast<int8_t>(search_cfg::id::SIZE),
                  "Did you forget to update search_cfg::id::SIZE value?");
    static_assert(pos < meta::size<type_list_t>::value,
                  "The specified config element is not contained in the configuration.");

    // TODO: Remove outer move once std::get is fixed for gcc 7.
    return std::move(get<pos>(std::move(cfg)).value);
}
//!\}
} // namespace seqan3
