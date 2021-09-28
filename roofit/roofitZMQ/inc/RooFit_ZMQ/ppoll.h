#ifndef ROOT_ROOFIT_ZMQ_ppoll
#define ROOT_ROOFIT_ZMQ_ppoll

#define ZMQ_BUILD_DRAFT_API
#include <zmq.hpp>
#include <vector>

namespace ZMQ {

int ppoll(zmq_pollitem_t *items_, size_t nitems_, long timeout_, const sigset_t *sigmask_);
int ppoll(std::vector<zmq_pollitem_t> &items, long timeout_, const sigset_t *sigmask_);
class ppoll_error_t : public zmq::error_t {
};

} // namespace ZMQ

#endif // ROOT_ROOFIT_ZMQ_ppoll
