#include <vector>
#include <iostream>
#include <seqan3/search/fm_index/all.hpp>

using namespace seqan3;
using namespace seqan3::literal;

int main()
{

// TODO(comments): outputs [A, A, G] instead of AAG

//! [cycle]
std::vector<dna4> genome{"AATAATAAC"_dna4};
fm_index index{genome};              // build the index

auto it = index.begin();             // create an iterator
// it.cycle_back();                  // cycle_back on begin() is undefined behaviour!
it.extend_right("AAC"_dna4);         // search the sequence "AAC"
std::cout << it.query() << '\n';     // outputs "AAC"
std::cout << it.last_char() << '\n'; // outputs 'C'

it.cycle_back();                     // search the sequence "AAT"
std::cout << it.query() << '\n';     // outputs "AAT"
std::cout << it.last_char() << '\n'; // outputs 'T'

it.cycle_back();                     // "it" doesn't change because the rightmost char is already the largest dna4 char.
std::cout << it.query() << '\n';     // outputs "AAT"
std::cout << it.last_char() << '\n'; // outputs 'T'
//! [cycle]

return 0;
}
