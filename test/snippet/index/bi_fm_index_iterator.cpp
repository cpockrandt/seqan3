#include <vector>
#include <iostream>
#include <seqan3/search/fm_index/all.hpp>

using namespace seqan3;
using namespace seqan3::literal;

int main()
{

{
std::cout << "Example cycle_back() and cycle_front()\n";
//! [cycle]
std::vector<dna4> genome{"GAATTAATGAAC"_dna4};
bi_fm_index index{genome}; // build the index

// TODO(comments): outputs [A, A, G] instead of AAG

auto it = index.begin(); // create an iterator
// it.cycle_back(); // cycle_back / cycle_front on begin() is undefined behaviour!
it.extend_right("AAC"_dna4); // search the sequence "AAC"
std::cout << it.query() << '\n'; // outputs "AAC"
std::cout << it.last_char() << '\n'; // outputs 'C'

// it.cycle_front(); // undefined behaviour! only cycle_back() is allowed after extend_right()
it.cycle_back(); // search the sequence "AAT"
std::cout << it.query() << '\n'; // outputs "AAT"
std::cout << it.last_char() << '\n'; // outputs 'T'

it.extend_left(dna4::G); // search the sequence "GAAT"
std::cout << it.query() << '\n'; // outputs "GAAC"
std::cout << it.last_char() << '\n'; // outputs 'G'

// it.cycle_back(); // undefined behaviour! only cycle_front() is allowed after extend_left()
it.cycle_front(); // search the sequence "TAAT"
std::cout << it.query() << '\n'; // outputs "TAAT"
std::cout << it.last_char() << '\n'; // outputs 'T'

it.cycle_front(); // search the sequence "TAAT"
std::cout << it.query() << '\n'; // outputs "TAAT"
std::cout << it.last_char() << '\n'; // outputs 'T'
//! [cycle]
}

{
std::cout << "Example to_fwd_iterator()\n";
//! [to_fwd_iterator]
std::vector<dna4> genome{"GAATTAACGAAC"_dna4};
bi_fm_index index{genome}; // build the index

auto it = index.begin(); // create an iterator
it.extend_left("AAC"_dna4); // search the sequence "AAC"
std::cout << it.query() << '\n'; // outputs "AAC"
auto uni_it = it.to_fwd_iterator(); // unidirectional iterator on the text "GAATTAACGAAC"
std::cout << uni_it.query() << '\n'; // outputs "CAA"
// Undefined behaviour! Cannot be called on the forward iterator if the last extension on the bidirectional
// iterator was to the left:
// it.cycle_back();
// std::cout << it.last_char() << '\n';

uni_it.extend_right(dna4::G); // search the sequence "AACG"
std::cout << uni_it.query() << '\n'; // outputs "AACG"
std::cout << uni_it.last_char() << '\n'; // outputs 'G'
uni_it.cycle_back(); // returns false since there is no sequence "AACT" in the text.
//! [to_fwd_iterator]
}

{
std::cout << "Example to_rev_iterator()\n";
//! [to_rev_iterator]
std::vector<dna4> genome{"GAATTAACGAAC"_dna4};
bi_fm_index index{genome}; // build the index

auto it = index.begin(); // create an iterator
it.extend_right("AAC"_dna4); // search the sequence "AAC"
std::cout << it.query() << '\n'; // outputs "AAC"
auto uni_it = it.to_rev_iterator(); // unidirectional iterator on the text "CAAGCAATTAAG"
std::cout << uni_it.query() << '\n'; // outputs "CAA"
// Undefined behaviour! Cannot be called on the reversed iterator if the last extension on the bidirectional
// iterator was to the right:
// it.cycle_back();
// std::cout << it.last_char() << '\n';

uni_it.extend_right(dna4::G); // search the sequence "CAAG"
std::cout << uni_it.query() << '\n'; // outputs "CAAG"
std::cout << uni_it.last_char() << '\n'; // outputs 'G'
uni_it.cycle_back(); // search the sequence "CAAT"
//! [to_rev_iterator]
}

    return 0;
}
