#pragma once

namespace seqan3
{

enum struct search_parameters_metric : uint8_t
{
    hamming,
    levenshtein
};

struct search_parameters
{
    // static constexpr uint8_t min_errors = 0u;
    static constexpr uint8_t max_errors = 0u;
    static constexpr search_parameters_metric metric {search_parameters_metric::hamming};
    static constexpr bool output_alignments = false;
};

template <typename t>
concept bool search_parameters_concept = requires(t c)
{
    // c.min_errors; // TODO: some functions might not allow min_errors. how to model that?
    c.max_errors;
    c.metric;

    // requires std::is_same_v<decltype(c.metric), search_parameters_metric> == true;
};

}
