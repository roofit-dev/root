#include "RooFit_ZMQ/ppoll.h"

namespace ZMQ {

// This function can throw, so wrap in try-catch!
int ppoll(zmq_pollitem_t *items_, size_t nitems_, long timeout_, const sigset_t *sigmask_)
{
   int rc = zmq_ppoll(items_, static_cast<int>(nitems_), timeout_, sigmask_);
   if (rc < 0)
      throw ppoll_error_t();
   return rc;
}

// This function can throw, so wrap in try-catch!
int ppoll(std::vector<zmq_pollitem_t> &items, long timeout_, const sigset_t *sigmask_)
{
   return ppoll(items.data(), items.size(), timeout_, sigmask_);
}

} // namespace ZMQ
