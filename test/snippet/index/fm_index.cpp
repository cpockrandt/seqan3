//! [all]
#include <vector>
#include <iostream>
#include <seqan3/search/fm_index/all.hpp>

using namespace seqan3;
using namespace seqan3::literal;

int main()
{
std::vector<dna4> genome{"ATCGATCGAAGGCTAGCTAGCTAAGGGA"_dna4};
fm_index index{genome}; // build the index

auto it = index.begin(); // create an iterator
it.extend_right("AAGG"_dna4); // search the pattern "AAGG"
std::cout << "Number of hits: " << it.count() << '\n'; // outputs: 2
std::cout << "Positions in the genome: ";
for (auto const & pos : it.locate()) // outputs: 8, 22
    std::cout << pos << ' ';
std::cout << '\n';
return 0;
}
//! [all]
