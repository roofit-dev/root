// Authors: Roel Aaij, Patrick Bos, Netherlands eScience Center / NIKHEF 2015-2021

#ifndef ZEROMQ_FUNCTIONS_H
#define ZEROMQ_FUNCTIONS_H 1

#include <cstddef> // std::size_t

namespace ZMQ {

template <class T>
std::size_t defaultSizeOf(const T &)
{
   return sizeof(T);
}

std::size_t stringLength(const char &cs);

} // namespace ZMQ

#endif // ZEROMQ_FUNCTIONS_H
