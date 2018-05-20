#pragma once

// #include <sdsl/suffix_trees.hpp>
//
// namespace seqan3::detail
// {
//     inline void _getErrorDistributions(vector<uint8_t> l, vector<uint8_t> u,
//                                        vector<vector<uint8_t> > & res,
//                                        uint8_t const e = 0)
//     {
//         if (l.size() == 0)
//         {
//             vector<uint8_t> _v;
//             res.push_back(_v);
//             return;
//         }
//
//         uint8_t l1 = l[0];
//         uint8_t u1 = u[0];
//         l.erase(l.begin());
//         u.erase(u.begin());
//
//         for (uint8_t i = max(e, l1); i <= u1; ++i)
//         {
//             vector<vector<uint8_t> > _res;
//             _getErrorDistributions(l, u, _res, i);
//             for (vector<uint8_t> &_resElem : _res)
//             {
//                 _resElem.insert(_resElem.begin(), i - e);
//                 res.push_back(_resElem);
//             }
//         }
//     }
//
//     // Compute all possible error distributions given a search. The result is in the
//     // same order as the search (s.pi)
//     // test [] [] _ = [[]]
//     // test (l:ls) (u:us) e = [[a-e] ++ x | a <- [(min l e)..u], x <- test ls us (e+a)]
//     inline void getErrorDistributions(Search const & s,
//                                       vector<vector<uint8_t> > & res)
//     {
//         vector<uint8_t> l = s.l;
//         vector<uint8_t> u = s.u;
//         _getErrorDistributions(l, u, res, 0u);
//     }
//
//     template <typename T>
//     inline void orderVector(Search const & s, vector<T> & v)
//     {
//         vector<T> v_tmp = v;
//         // Reorder blocks s.t. they are in a sequential order (from left to right)
//         for (uint8_t i = 0; i < s.pi.size(); ++i)
//         {
//             uint8_t index = s.pi[i] - 1;
//             v[index] = v_tmp[i];
//         }
//     }
//
//     inline void getOrderedSearch(Search const & s, Search & os)
//     {
//         // Reorder blocks s.t. they are in a sequential order (from left to right)
//         // Blocklength is stored as absolute values instead of cumulative values.
//         for (uint8_t i = 0; i < s.pi.size(); ++i)
//         {
//             uint8_t index = s.pi[i] - 1;
//             os.pi[index] = s.pi[i];
//             os.l[index] = s.l[i];
//             os.u[index] = s.u[i];
//             os.blocklength[index] = s.blocklength[i] - ((i > 0) ? s.blocklength[i - 1] : 0);
//         }
//     }
//
//     template <typename TDelegate, typename TText, typename TIndex, typename TIndexSpec, typename TText2, typename TDirection>
//     inline void _searchDeletion(TDelegate & delegate,
//                          Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
//                          TText2 const & needle,
//                          unsigned needleLeftIt,
//                          unsigned needleRightIt,
//                          uint8_t const errors,
//                          Search const & s,
//                          TDirection const /**/,
//                          uint8_t const blockIndex
//                          #ifdef ENABLE_DEBUG_MACRO
//                             , std::string str
//                          #endif
//     )
//     {
//         uint8_t maxErrorsLeftInBlock = s.u[blockIndex] - errors;
//         uint8_t minErrorsLeftInBlock = (s.l[blockIndex] > errors) ? (s.l[blockIndex] - errors) : 0;
//
//         if (minErrorsLeftInBlock == 0)
//         {
//             uint8_t const blockIndex2 = min((size_t) blockIndex + 1, s.u.size() - 1);
//             bool const goToRight2 = s.pi[blockIndex2] > s.pi[blockIndex2 - 1];
//
//             if (goToRight2)
//                 _search(delegate, iter, needle, needleLeftIt, needleRightIt, errors, s, Rev(), blockIndex2, true
//                     #ifdef ENABLE_DEBUG_MACRO
//                         , str
//                     #endif
//                 );
//             else
//                 _search(delegate, iter, needle, needleLeftIt, needleRightIt, errors, s, Fwd(), blockIndex2, true
//                     #ifdef ENABLE_DEBUG_MACRO
//                         , str
//                     #endif
//                 );
//         }
//
//         if (maxErrorsLeftInBlock > 0)
//         {
//             if (goDown(iter, TDirection()))
//             {
//                 do {
//                     _searchDeletion(delegate, iter, needle, needleLeftIt, needleRightIt, errors + 1, s, TDirection(), blockIndex
//                         #ifdef ENABLE_DEBUG_MACRO
//                             , is_same<TDirection, Rev>::value ? str + "-" : "-" + str
//                         #endif
//                     );
//                 } while (goRight(iter, TDirection()));
//             }
//         }
//     }
//
//     template <typename TDelegate, typename TText, typename TIndex, typename TIndexSpec, typename TText2, typename TDir>
//     inline void _searchChildren(TDelegate & delegate,
//                         Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
//                         TText2 const & needle,
//                         unsigned const needleLeftIt,
//                         unsigned const needleRightIt,
//                         uint8_t const errors,
//                         Search const & s,
//                         TDir const /**/,
//                         uint8_t const blockIndex,
//                         bool const indels,
//                         auto minErrorsLeftInBlock
//                         #ifdef ENABLE_DEBUG_MACRO
//                             , std::string str
//                         #endif
//     )
//     {
//         bool goToRight = is_same<TDir, Rev>::value;
//         //cout << "vor goDown: " << representative(iter.fwdIter) << endl;
//         if (goDown(iter, TDir()))
//         {
//             //cout << "nac goDown: " << representative(iter.fwdIter) << endl;
//             unsigned charsLeft = s.blocklength[blockIndex] - (needleRightIt - needleLeftIt - 1/*(needleLeftIt != needleRightIt)*/);
//             do {
//                 //if (!(needleLeftIt == 0 && needleRightIt == length(needle) + 1))
//                 //{
//                 bool delta = !ordEqual(parentEdgeLabel(iter, TDir()), needle[goToRight ? needleRightIt - 1 : needleLeftIt - 1]);
//
//                 // TODO: this is not optimal yet! we have more edges than in the theoretical model,
//                 // since we go down an edge before we check whether it can even work out!
//                 if (!indels && minErrorsLeftInBlock > 0 && charsLeft - 1 < minErrorsLeftInBlock - delta)
//                     continue;
//
//                 auto needleLeftIt2 = needleLeftIt - !goToRight;
//                 auto needleRightIt2 = needleRightIt + goToRight;
//
//                 #ifdef ENABLE_DEBUG_MACRO
//                     CharString c = (CharString) parentEdgeLabel(iter, TDir());
//                 #endif
//
//                 if (needleRightIt - needleLeftIt == s.blocklength[blockIndex])
//                 {
//                     // leave the possibility for one or multiple deletions! therefore, don't change direction, etc!
//                     if (indels)
//                         _searchDeletion(delegate, iter, needle, needleLeftIt2, needleRightIt2, errors + delta, s, TDir(), blockIndex
//                             #ifdef ENABLE_DEBUG_MACRO
//                                 , goToRight ? str + toCString(c) : toCString(c) + str
//                             #endif
//                         );
//                     else
//                     {
//                         uint8_t blockIndex2 = min((size_t) blockIndex + 1, s.u.size() - 1);
//                         bool goToRight2 = s.pi[blockIndex2] > s.pi[blockIndex2 - 1];
//                         if (goToRight2)
//                             _search(delegate, iter, needle, needleLeftIt2, needleRightIt2, errors + delta, s, Rev(), blockIndex2, indels
//                                 #ifdef ENABLE_DEBUG_MACRO
//                                     , goToRight ? str + toCString(c) : toCString(c) + str
//                                 #endif
//                             );
//                         else
//                             _search(delegate, iter, needle, needleLeftIt2, needleRightIt2, errors + delta, s, Fwd(), blockIndex2, indels
//                                 #ifdef ENABLE_DEBUG_MACRO
//                                     , goToRight ? str + toCString(c) : toCString(c) + str
//                                 #endif
//                             );
//                     }
//                 }
//                 else
//                     _search(delegate, iter, needle, needleLeftIt2, needleRightIt2, errors + delta, s, TDir(), blockIndex, indels
//                         #ifdef ENABLE_DEBUG_MACRO
//                             , goToRight ? str + toCString(c) : toCString(c) + str
//                         #endif
//                     );
//                 //}
//
//                 // Deletion
//                 if (indels)
//                 {
//                     _search(delegate, iter, needle, needleLeftIt, needleRightIt, errors + 1, s, TDir(), blockIndex, indels
//                         #ifdef ENABLE_DEBUG_MACRO
//                             , goToRight ? str + "-" : "-" + str
//                         #endif
//                     );
//                 }
//             } while (goRight(iter, TDir()));
//         }
//     }
//
//     template <typename TDelegate, typename TText, typename TIndex, typename TIndexSpec, typename TText2, typename TDir>
//     inline void _searchExact(TDelegate & delegate,
//                         Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
//                         TText2 const & needle,
//                         unsigned const needleLeftIt,
//                         unsigned const needleRightIt,
//                         uint8_t const errors,
//                         Search const & s,
//                         TDir const /**/,
//                         uint8_t const blockIndex,
//                         bool const indels
//                         #ifdef ENABLE_DEBUG_MACRO
//                             , std::string str
//                         #endif
//     )
//     {
//         bool goToRight2 = s.pi[blockIndex + 1] > s.pi[blockIndex]; // TODO: segfault? value doesn't matter for last block, but we should code it better anyway
//         if (is_same<TDir, Rev>::value)
//         {
//             unsigned infixPosLeft = needleRightIt /*+ (needleLeftIt == needleRightIt)*/ - 1;
//             unsigned infixPosRight = needleLeftIt + s.blocklength[blockIndex] - 1;
//
//             //if (infixPosLeft > infixPosRight)
//             //    exit(23);
//
//             if (!goDown(iter, infix(needle, infixPosLeft, infixPosRight + 1), TDir()))
//                 return;
//
//             if (goToRight2)
//                 _search(delegate, iter, needle, needleLeftIt, infixPosRight + 2, errors, s, Rev(), min((size_t) blockIndex + 1, s.u.size() - 1), indels
//                     #ifdef ENABLE_DEBUG_MACRO
//                         , str + toCString((CharString) infix(needle, infixPosLeft, infixPosRight + 1))
//                     #endif
//                 );
//             else
//                 _search(delegate, iter, needle, needleLeftIt, infixPosRight + 2, errors, s, Fwd(), min((size_t) blockIndex + 1, s.u.size() - 1), indels
//                     #ifdef ENABLE_DEBUG_MACRO
//                         , str + toCString((CharString) infix(needle, infixPosLeft, infixPosRight + 1))
//                     #endif
//                 );
//         }
//         else
//         {
//             // has to be signed, otherwise we run into troubles when checking for -1 >= 0u
//             signed infixPosLeft = needleRightIt - s.blocklength[blockIndex] - 1;
//             signed infixPosRight = needleLeftIt /*- (needleLeftIt == needleRightIt)*/ - 1;
//
//             //if (infixPosLeft > infixPosRight)
//             //    exit(23);
//
//             #ifdef ENABLE_DEBUG_MACRO
//                 auto str2 = str;
//             #endif
//
//             while (infixPosRight >= infixPosLeft)
//             {
//                 if (!goDown(iter, needle[infixPosRight], TDir()))
//                     return;
//                 #ifdef ENABLE_DEBUG_MACRO
//                     str2 = toCString((CharString) needle[infixPosRight]) + str2;
//                 #endif
//                 --infixPosRight;
//                 // --needleLeftIt;
//             }
//             if (goToRight2)
//                 _search(delegate, iter, needle, /*needleLeftIt*/infixPosLeft, needleRightIt, errors, s, Rev(), min((size_t)blockIndex + 1, s.u.size() - 1), indels
//                     #ifdef ENABLE_DEBUG_MACRO
//                         , str2
//                     #endif
//                 );
//             else
//                 _search(delegate, iter, needle, /*needleLeftIt*/infixPosLeft, needleRightIt, errors, s, Fwd(), min((size_t)blockIndex + 1, s.u.size() - 1), indels
//                     #ifdef ENABLE_DEBUG_MACRO
//                         , str2
//                     #endif
//                 );
//         }
//     }
//
//     // TODO(knut): prefer paths with less errors over paths with more errors
//     template <typename TDelegate, typename TText, typename TIndex, typename TIndexSpec, typename TText2, typename TDir>
//     inline void _search(TDelegate & delegate,
//                         Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
//                         TText2 const & needle,
//                         unsigned const needleLeftIt,
//                         unsigned const needleRightIt,
//                         uint8_t const errors,
//                         Search const & s,
//                         TDir const /**/,
//                         uint8_t const blockIndex,
//                         bool const indels
//                         #ifdef ENABLE_DEBUG_MACRO
//                             , std::string str
//                         #endif
//     )
//     {
//         uint8_t maxErrorsLeftInBlock = s.u[blockIndex] - errors;
//         uint8_t minErrorsLeftInBlock = (s.l[blockIndex] > errors) ? (s.l[blockIndex] - errors) : 0;
//
//         // Done.
//         if (minErrorsLeftInBlock == 0 && needleLeftIt == 0 && needleRightIt == length(needle) + 1) // NOTE: switch to iterator syntax
//         {
//             delegate(iter
//                 #ifdef ENABLE_DEBUG_MACRO
//                     , str
//                 #endif
//             );
//             //if (maxErrorsLeftInBlock == 0)
//                 //return;
//         }
//
//         //if (needleLeftIt == 0 && needleRightIt == length(needle) + 1)
//         //    return;
//
//         // Exact search in current block.
//         else if (maxErrorsLeftInBlock == 0 && /*!(blockIndex > 0 && */needleRightIt - needleLeftIt - 1 /*==*/ != s.blocklength[blockIndex]/*)*/) // TODO: why blockIndex > 0? vermutlich wegen needleRightIt == nedleLeftIt
//         {
//             _searchExact(delegate, iter, needle, needleLeftIt, needleRightIt, errors, s, TDir(), blockIndex, indels
//                 #ifdef ENABLE_DEBUG_MACRO
//                     , str
//                 #endif
//             );
//         }
//         // Approximate search in current block.
//         else // if (s.blocklength[blockIndex] - (needleRightIt - needleLeftIt - (needleLeftIt != needleRightIt)) >= minErrorsLeftInBlock)
//         {
//             // Insertion
//             if (indels) // && !(needleLeftIt == 0 && needleRightIt == length(needle) + 1) && charsLeft > 0 // no insertions at the end of the block when we didnt increment the blockid (for deletions)
//             {
//                 bool const goToRight = is_same<TDir, Rev>::value;
//                 auto const needleLeftIt2 = needleLeftIt - !goToRight;
//                 auto const needleRightIt2 = needleRightIt + goToRight;
//
//                 if (needleRightIt - needleLeftIt == s.blocklength[blockIndex]) //
//                     // TODO: check!
//                     // leave the possibility for one or multiple deletions! therefore, don't change direction, etc!
//                     _searchDeletion(delegate, iter, needle, needleLeftIt2, needleRightIt2, errors + 1, s, TDir(), blockIndex
//                         #ifdef ENABLE_DEBUG_MACRO
//                             , goToRight ? str + "+" : "+" + str
//                         #endif
//                     );
//                 else
//                     _search(delegate, iter, needle, needleLeftIt2, needleRightIt2, errors + 1, s, TDir(), blockIndex, indels
//                         #ifdef ENABLE_DEBUG_MACRO
//                             , goToRight ? str + "+" : "+" + str
//                         #endif
//                     );
//             }
//
//             _searchChildren(delegate, iter, needle, needleLeftIt, needleRightIt, errors, s, TDir(), blockIndex, indels, minErrorsLeftInBlock
//                 #ifdef ENABLE_DEBUG_MACRO
//                     , str
//                 #endif
//             );
//         }
//     }
//
//     template <typename TDelegate, typename TText, typename TIndex, typename TIndexSpec, typename TText2>
//     inline void search(TDelegate & delegate,
//                        Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > it,
//                        TText2 const & needle, Search const & s, bool const indels,
//                        bool const debug)
//     {
//         if (s.initialDirection)
//         {
//             _search(delegate, it, needle, s.startPos, s.startPos + 1, 0, s, Rev(), 0, indels
//                 #ifdef ENABLE_DEBUG_MACRO
//                     , ""
//                 #endif
//             );
//         }
//         else
//         {
//             _search(delegate, it, needle, s.startPos - 1, s.startPos, 0, s, Fwd(), 0, indels
//                 #ifdef ENABLE_DEBUG_MACRO
//                     , ""
//                 #endif
//             );
//         }
//     }
//
//     template <typename TChar>
//     inline void _countSearch(unsigned const needleLength,
//                             unsigned const needleLeftIt,
//                             unsigned const needleRightIt,
//                             uint8_t const errors, Search const & s,
//                             bool const goToRight, uint8_t const blockIndex,
//                             bool const debug, unsigned & count, bool const indels)
//     {
//         // Done.
//         if (needleLeftIt == 0 && needleRightIt == needleLength + 1)
//         {
//             return;
//         }
//
//         // Exact search in current block.
//         uint8_t maxErrorsLeftInBlock = s.u[blockIndex] - errors;
//         uint8_t minErrorsLeftInBlock = (s.l[blockIndex] > errors) ? (s.l[blockIndex] - errors) : 0;
//         if (maxErrorsLeftInBlock == 0)
//         {
//             signed infixPosLeft, infixPosRight;
//
//             if (goToRight)
//             {
//                 infixPosLeft = needleRightIt + /*(needleLeftIt == needleRightIt)*/ - 1;
//                 infixPosRight = needleLeftIt + s.blocklength[blockIndex] - 1;
//             }
//             else
//             {
//                 infixPosLeft = needleRightIt - s.blocklength[blockIndex] - 1;
//                 infixPosRight = needleLeftIt - 1;
//             }
//
//             bool goToRight2 = s.pi[blockIndex + 1] > s.pi[blockIndex];
//             count += infixPosRight - infixPosLeft + 1;
//             if (goToRight)
//                 return _countSearch<TChar>(needleLength, needleLeftIt, infixPosRight + 2, errors, s, goToRight2, blockIndex + 1, debug, count, indels);
//             else
//                 return _countSearch<TChar>(needleLength, infixPosLeft, needleRightIt, errors, s, goToRight2, blockIndex + 1, debug, count, indels);
//         }
//         // Approximate search in current block.
//         else
//         {
//             unsigned charsLeft = s.blocklength[blockIndex] - (needleRightIt - needleLeftIt - 1/*(needleLeftIt != needleRightIt)*/);
//             // if (charsLeft < minErrorsLeftInBlock)
//             //     return;
//
//             // Insertion
//             if (indels)
//             {
//                 uint8_t blockIndex2 = blockIndex;
//                 bool goToRight2 = goToRight;
//                 if (needleRightIt - needleLeftIt == s.blocklength[blockIndex])
//                 {
//                     ++blockIndex2;
//                     goToRight2 = s.pi[blockIndex2] > s.pi[blockIndex];
//                 }
//
//                 _countSearch<TChar>(needleLength, needleLeftIt - !goToRight, needleRightIt + goToRight, errors + 1, s, goToRight2, blockIndex2, debug, count, indels);
//             }
//
//             for (unsigned i = 0; i < ValueSize<TChar>::VALUE; ++i)
//             {
//                 bool delta = (i != 0);
//
//                 // Deletion
//                 if (indels)
//                 {
//                     // TODO: the following check is probably obsolete?
//                     /*if (minErrorsLeftInBlock > 0 && charsLeft < minErrorsLeftInBlock - 1)
//                     {
//                         std::cout << "WAAAAAaaaaaa" << std::endl;
//                         continue;
//                     }*/
//
//                     _countSearch<TChar>(needleLength, needleLeftIt, needleRightIt, errors + 1, s, goToRight, blockIndex, debug, count, indels);
//                 }
//
//                 if (minErrorsLeftInBlock > 0 && charsLeft - 1 < minErrorsLeftInBlock - delta)
//                     continue;
//
//                 uint8_t blockIndex2 = blockIndex;
//                 bool goToRight2 = goToRight;
//                 if (needleRightIt - needleLeftIt == s.blocklength[blockIndex])
//                 {
//                     ++blockIndex2;
//                     goToRight2 = s.pi[blockIndex2] > s.pi[blockIndex];
//                 }
//
//                 ++count;
//                 if (goToRight)
//                     _countSearch<TChar>(needleLength, needleLeftIt, needleRightIt + 1, errors + delta, s, goToRight2, blockIndex2, debug, count, indels);
//                 else
//                     _countSearch<TChar>(needleLength, needleLeftIt - 1, needleRightIt, errors + delta, s, goToRight2, blockIndex2, debug, count, indels);
//             }
//         }
//     }
//
//     template <typename TChar>
//     inline unsigned countSearch(unsigned int needleLength, Search const & s,
//                                 bool const indels, bool const debug = false)
//     {
//         unsigned count = 0;
//         _countSearch<TChar>(needleLength, s.startPos - !s.initialDirection, s.startPos + s.initialDirection, 0, s, s.initialDirection /* true <=> Rev() */, 0, debug, count, indels);
//         return count;
//     }
//
//     template <typename TChar>
//     inline unsigned countSearchScheme(unsigned int needleLength, SearchScheme const & ss,
//                                 bool const indels, bool const debug = false)
//     {
//         unsigned count = 0;
//         for (Search const & s : ss)
//             count += countSearch<TChar>(needleLength, s, indels, debug);
//         return count;
//     }
//
//     //////////////////////////////////////////////////////////////////////////////
//
//     inline void checkSearchScheme(SearchScheme const & ss)
//     {
//         // max number of errors
//         uint8_t errors = 0u;
//
//         unsigned blocks = ss[0].pi.size();
//
//         // Get all ordered error distributions
//         vector<vector<uint8_t> > distrs;
//         for (Search const & s : ss)
//         {
//             errors = max(errors, s.u.back());
//             vector<vector<uint8_t> > _distrs;
//             getErrorDistributions(s, _distrs);
//             for (vector<uint8_t > & distr : _distrs)
//             {
//                 orderVector(s, distr);
//             }
//             distrs.insert(distrs.end(), _distrs.begin(), _distrs.end());
//
//             if (s.pi.size() != blocks)
//             {
//                 cerr << "ERROR: Different number of blocks across searches!" << endl;
//                 exit(15);
//             }
//         }
//
//         Search trivialSearch;
//         trivialSearch.pi.resize(blocks);
//         trivialSearch.l.resize(blocks);
//         trivialSearch.u.resize(blocks);
//         for (unsigned i = 0; i < blocks; ++i)
//         {
//             trivialSearch.pi[i] = i + 1;
//             trivialSearch.l[i] = 0;
//             trivialSearch.u[i] = errors;
//         }
//
//         vector<vector<uint8_t> > trivial_distrs;
//         getErrorDistributions(trivialSearch, trivial_distrs);
//
//         sort(distrs.begin(), distrs.end());
//         sort(trivial_distrs.begin(), trivial_distrs.end());
//
//         /*if (distrs != trivial_distrs)
//         {
//             cerr << "ERROR: SearchScheme seems to be broken (e.g. does not cover all possible error configurations)" << endl;
//             exit(15);
//         }*/
//     }
//
//     inline void _readSearchScheme_parsLine(string const & str, vector<uint8_t> & v)
//     {
//         istringstream ss(str);
//         while (ss)
//         {
//             string s;
//             if (!getline(ss, s, ','))
//                 break;
//             v.push_back(stoul(s));
//         }
//     }
//
//     inline void readSearchScheme(CharString const & path, SearchScheme & ss)
//     {
//         ifstream infile(toCString(path));
//
//         while (infile)
//         {
//             string pi, l, u;
//             if (!getline(infile, pi) || !getline(infile, l) || !getline(infile, u))
//                 break;
//
//             Search s;
//             _readSearchScheme_parsLine(pi, s.pi);
//             _readSearchScheme_parsLine(l, s.l);
//             _readSearchScheme_parsLine(u, s.u);
//             s.blocklength.assign(s.pi.size(), 0);
//             ss.push_back(s);
//         }
//
//         checkSearchScheme(ss);
//     }
//
//     template <typename TChar>
//     inline void computeBlocklength(SearchScheme & ss, uint16_t const needleLength, bool optimal = false)
//     {
//         // count edges for all blocklengths and choose theoretically optimal one!
//         if (optimal)
//         {
//             unsigned maxErrors = 0, countEdges, countEdgesOptimal = (unsigned int) -1; // maximum value of unsigned int
//             vector<uint8_t> blocklength(ss[0].pi.size()), blocklengthOptimal(ss[0].pi.size());
//             for (Search const & s : ss)
//             {
//                 maxErrors = max(maxErrors, (unsigned) s.u.back());
//             }
//
//             for (unsigned i = 1; i < ss[0].pi.size(); ++i)
//             {
//                 blocklength[i-1] = i * maxErrors;
//             }
//             back(blocklength) = needleLength;
//
//             while (true)
//             {
//                 vector<uint8_t> _blocklengthAbsolute = blocklength;
//                 for (unsigned i = _blocklengthAbsolute.size(); i > 0; --i)
//                 {
//                     _blocklengthAbsolute[i] -= _blocklengthAbsolute[i - 1];
//                 }
//
//                 setBlockLength(ss, _blocklengthAbsolute);
//                 initSearchScheme(ss);
//                 countEdges = countSearchScheme<TChar>(needleLength, ss, false); // TODO: document that indels are always ignored!
//                 if (countEdges < countEdgesOptimal)
//                 {
//                     countEdgesOptimal = countEdges;
//                     blocklengthOptimal = _blocklengthAbsolute;
//                 }
//
//                 // compute next blockpos
//                 signed i;
//                 for (i = blocklength.size() - 2; i >= 0; --i)
//                 {
//                     // find rightmost element in blockpos that we can increment
//                     if (blocklength[i] < blocklength[i + 1] - maxErrors)
//                     {
//                         ++blocklength[i];
//                         // reset all elements after pos i
//                         for (unsigned j = i + 1; j < blocklength.size() - 1; ++j)
//                         {
//                             blocklength[j] = blocklength[j - 1] + maxErrors;
//                         }
//                         break;
//                     }
//                 }
//                 if (i < 0)
//                     break;
//             }
//
//             setBlockLength(ss, blocklengthOptimal); // TODO: print optimal block lengths!
//         }
//         else
//         {
//             uint8_t blocks = ss[0].pi.size();
//             uint8_t blocklength = needleLength / blocks;
//             uint16_t rest = needleLength - blocks * blocklength;
//             vector<uint8_t> blocklengths;
//             for (unsigned i = 0; i < blocks; ++i)
//             {
//                 if (i < rest)
//                     blocklengths.push_back(blocklength + 1);
//                 else
//                     blocklengths.push_back(blocklength);
//             };
//
//             // TODO: remove the following check
//             if (accumulate(blocklengths.begin(), blocklengths.end(), 0) != needleLength)
//             {
//                 cerr << "ERROR: blocklength was not initialized correctly!" << endl;
//                 exit(33);
//             }
//
//             setBlockLength(ss, blocklengths);
//         }
//         initSearchScheme(ss);
//     }
//
//     inline void printInputFormatSearchScheme()
//     {
//         cout << "Format of input files for search schemes (for 2 errors):" << endl;
//
//         cout << "1,2,3" << endl;
//         cout << "0,0,2" << endl;
//         cout << "0,1,2" << endl;
//
//         cout << "3,2,1" << endl;
//         cout << "0,0,0" << endl;
//         cout << "0,2,2" << endl;
//
//         cout << "2,3,1" << endl;
//         cout << "0,1,1" << endl;
//         cout << "0,1,2" << endl;
//     }
// }
