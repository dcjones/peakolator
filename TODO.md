

The priority queue should be over interval bounds, not intervals.

This means we need a regular old queue (or stack) for intervals.


* Overflow is a real concern. Sums should be computed in `uint64_t` if it
  doesn't dramatically affect speed.

* Start writing the algorithm proper.

* We need a means of 

