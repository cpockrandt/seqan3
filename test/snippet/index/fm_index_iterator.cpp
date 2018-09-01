#include <vector>
#include <iostream>
#include <seqan3/index/all.hpp>

using namespace seqan3;
using namespace seqan3::literal;

int main()
{
    //! [cycle]
    std::vector<dna4> genome{"AATAATAAC"_dna4};
    fm_index<std::vector<dna4>> index{genome}; // build the index

    // TODO(comments): outputs [A, A, G] instead of AAG

    auto it = index.begin(); // create an iterator
    // it.cycle_back(); // cycle_back on begin() is undefined behaviour!
    it.extend_right("AAC"_dna4); // search the sequence "AAC"
    std::cout << it.query() << '\n'; // outputs "AAC"
    std::cout << it.last_char() << '\n'; // outputs 'C'

    it.cycle_back(); // search the sequence "AAT"
    std::cout << it.query() << '\n'; // outputs "AAT"
    std::cout << it.last_char() << '\n'; // outputs 'T'

    it.cycle_back(); // iterator does not change because the rightmost character is already the largest dna4 character.
    std::cout << it.query() << '\n'; // outputs "AAT"
    std::cout << it.last_char() << '\n'; // outputs 'T'
    //! [cycle]

    return 0;
}
