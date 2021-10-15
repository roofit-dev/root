// Authors: Roel Aaij, Patrick Bos, Netherlands eScience Center / NIKHEF 2015-2021

#include "RooFit_ZMQ/functions.h"

#include <cstring>

namespace ZMQ {
std::size_t stringLength(const char &cs)
{
   return strlen(&cs);
}
} // namespace ZMQ
