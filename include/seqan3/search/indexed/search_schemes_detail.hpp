#pragma once

#include <sdsl/suffix_trees.hpp>

namespace seqan3::detail
{

    template <typename cst_t, typename iter_t, typename search_parameters_t, typename callback_f>
    // TODO: concepts überall wiederholen?
    // requires random_access_iterator_concept<iter_t> && search_parameters_concept<search_parameters_t> // TODO const iterator?
    inline void _search_backtracking(cst_t const & cst, iter_t first, iter_t const last, search_parameters_t const & params, callback_f && callback, uint64_t l, uint64_t r, uint8_t errors) // TODO: iterator als referenz übergeben?
    {
        // Exact case.
        if (errors == params.max_errors)
        {
            if (sdsl::backward_search(cst, l, r, first, last, l, r))
                callback(l, r, errors);
        }
        // Approximate case.
        else // if (errors < threshold)
        {
            // Base case.
            if (first == last)
            {
                callback(l, r, errors);
            }
            // Recursive case.
            else
            {
                // Insertion.
                if (params.metric == search_parameters_metric::levenshtein) // TODO: if constexpr?
                    _search_backtracking(cst, first, last - 1, params, callback, l, r, errors + 1);

                uint64_t l_parent = l, r_parent = r;
                for (typename cst_t::comp_char_type c = 1; c < cst.sigma; ++c)
                {
                    auto cc = cst.comp2char[c]; // TODO: use of comp2char sucks since it will be converted back to c in backward_search
                    // std::cout << "loop c = " << (unsigned)c << ", cc = " << (unsigned)cc << "\n";
                    if (!sdsl::backward_search(cst, l_parent, r_parent, cc, l, r))
                        continue;

                    // Match / Mismatch.
                    const bool delta = *(last - 1) != cc; // TODO: check
                    _search_backtracking(cst, first, last - 1, params, callback, l, r, errors + delta);

                    // Deletion.
                    if (params.metric == search_parameters_metric::levenshtein)
                        _search_backtracking(cst, first, last, params, callback, l, r, errors + 1);
                }
            }
        }
    }

    template <typename cst_t, typename iter_t, typename search_parameters_t, typename callback_f>
    // requires random_access_iterator_concept<iter_t> && search_parameters_concept<search_parameters_t>
    inline void search_backtracking(cst_t const & cst, iter_t first, iter_t last, search_parameters_t const & params, callback_f && callback)
    {
        uint8_t errors = 0;
        uint64_t l = 0, r = cst.size() - 1;
        detail::_search_backtracking(cst, first, last, params, callback, l, r, errors);
    }

}
