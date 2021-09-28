#include "RooFit_ZMQ/functions.h"

#define ZMQ_BUILD_DRAFT_API
#include <zmq.hpp>
#include "RooFit_ZMQ/ZeroMQSvc.h"

namespace ZMQ {
size_t stringLength(const char &cs)
{
   return strlen(&cs);
}
} // namespace ZMQ
