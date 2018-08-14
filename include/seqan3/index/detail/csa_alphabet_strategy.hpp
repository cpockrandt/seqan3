#pragma once

#include <sdsl/int_vector.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/rank_support.hpp>
#include <sdsl/select_support.hpp>
#include <sdsl/sdsl_concepts.hpp>
#include <sdsl/config.hpp>
#include <string>

// This file extends sdsl/csa_alphabet_stragy.hpp

namespace sdsl {

    //!\brief Byte alphabet that does no mapping of char_type to comp_char_type and vice versa.
    //        This is recommended when the underlying text uses the entire alphabet and not just a small subset.
    class plain_byte_alphabet {
    public:
        class mapping_wrapper;

        typedef int_vector<>::size_type size_type;
        typedef mapping_wrapper char2comp_type;
        typedef mapping_wrapper comp2char_type;
        typedef int_vector<64> C_type;
        typedef uint16_t sigma_type;
        typedef uint8_t char_type;
        typedef uint8_t comp_char_type;
        typedef std::string string_type;
        typedef byte_alphabet_tag alphabet_category;
        enum { int_width = 8 };

        //! Helper class for the char2comp and comp2char mapping
        class mapping_wrapper {
        public:
            mapping_wrapper() {}

            constexpr char_type operator[](char_type const c) const noexcept { return c; }
        };

        const char2comp_type char2comp;
        const comp2char_type comp2char;
        const C_type &C;
        const sigma_type &sigma;

    private:
        C_type m_C;            // Cumulative counts for the compact alphabet [0..sigma].
        sigma_type m_sigma;        // Effective size of the alphabet.

    public:
        //! Default constructor
        plain_byte_alphabet()
                : C(m_C), sigma(m_sigma), m_sigma(0) {
        }

        //! Construct from a byte-stream
        /*!
         *  \param text_buf Byte stream.
         *  \param len      Length of the byte stream.
         */
        plain_byte_alphabet(int_vector_buffer<8> &text_buf, int_vector_size_type len)
                : C(m_C), sigma(m_sigma) {
            m_sigma = 0;
            if (0 == len or 0 == text_buf.size()) return;
            assert(len <= text_buf.size());
            // initialize vectors
            m_C = int_vector<64>(257, 0);
            // count occurrences of each symbol
            for (size_type i = 0; i < len; ++i) {
                ++m_C[text_buf[i]];
            }
            assert(1 == m_C[0]); // null-byte should occur exactly once
            m_sigma = 255;
            for (int i = 0; i < 256; ++i)
                if (m_C[i]) {
                    m_sigma = i + 1;
                    // m_C[m_sigma]	   = m_C[i];
            //        ++m_sigma;
                }
            //m_C.resize(m_sigma + 1);
            for (int i = (int) 256; i > 0; --i)
                m_C[i] = m_C[i - 1];
            m_C[0] = 0;
            for (int i = 1; i <= (int) 256; ++i)
                m_C[i] += m_C[i - 1];
            assert(C[sigma] == len);
        }

        plain_byte_alphabet(const plain_byte_alphabet &strat)
                : C(m_C), sigma(m_sigma), m_C(strat.m_C), m_sigma(strat.m_sigma) {
        }

        plain_byte_alphabet(plain_byte_alphabet &&strat)
                : C(m_C), sigma(m_sigma), m_C(std::move(strat.m_C)), m_sigma(strat.m_sigma) {}

        plain_byte_alphabet &operator=(const plain_byte_alphabet &strat) {
            if (this != &strat) {
                plain_byte_alphabet tmp(strat);
                *this = std::move(tmp);
            }
            return *this;
        }

        plain_byte_alphabet &operator=(plain_byte_alphabet &&strat) {
            if (this != &strat) {
                m_C = std::move(strat.m_C);
                m_sigma = std::move(strat.m_sigma);
            }
            return *this;
        }

        size_type
        serialize(std::ostream &out, structure_tree_node *v, std::string name = "") const {
            structure_tree_node *child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += m_C.serialize(out, child, "m_C");
            written_bytes += write_member(m_sigma, out, child, "m_sigma");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        void load(std::istream &in) {
            m_C.load(in);
            read_member(m_sigma, in);
        }
    };

}
