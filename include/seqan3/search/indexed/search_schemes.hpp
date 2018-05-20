#pragma once

#include <seqan3/std/concept/iterator.hpp>
#include <sdsl/suffix_trees.hpp>

namespace seqan3
{

/*template <typename cst_t, typename search_parameters_t>
requires search_parameters_concept<search_parameters_t>
inline auto search(cst_t const & cst, typename cst_t::string_type const & pattern, search_params const & params)
{
    std::vector<...> v;
    auto callback = []() {
        v.add(...);
    }
    search_and_then(cst, pattern, params, callback);
    return v;
}*/
// TODO: filter duplicates in params??
/*
    cst_fwd
    first -- requires random_access_iterator_concept
    last
    search_parameters_t (error, distance_metric)
    callback
    TODO: wie sicherstellen, dass gleiches alphabet genutzt wird?
*/ // TODO: search_parameters struct ist aber nicht gut um den user zu constexpr metric/error zu bringen!
template <typename cst_t, typename iter_t, typename search_parameters_t, typename callback_f>
requires random_access_iterator_concept<iter_t> && search_parameters_concept<search_parameters_t> // TODO const iterator?
inline void search_and_then(cst_t const & cst, iter_t first, iter_t last, search_parameters_t const & params, callback_f && callback)
{
    // assert(params.min_errors <= 255);
    assert(params.max_errors <= 255);
    assert(first != last);
    detail::search_backtracking(cst, first, last, params, callback);
}

/*
    cst_fwd
    cst_bwd
    first -- requires random_access_iterator_concept
    last
    search_parameters_t (error, distance_metric)
*/


// template <typename cst_t, typename search_parameters_t, typename callback_f>
// requires search_parameters_concept<search_parameters_t>
// inline void search_and_then(cst_t const & cst, typename cst_t::string_type const & pattern, search_parameters_t const & params, callback_f && callback)
// {
//     auto & scheme = seqan3::detail::schemes[params.min_errors][params.max_errors];
//     seqan3::detail::computeBlocklength(scheme, pattern.size());
//
//     bi_fm_index_iter<cst_t> it(cst);
//
//     for (seqan3::detail::Search const & s : scheme)
//         seqan3::detail::search(callback, it, pattern, s, params);
// }

} //
