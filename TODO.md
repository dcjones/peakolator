
* We should support a x_max parameter that controls how large the value of any
  given x_i may be.

* Overflow is a real concern. Sums should be computed in `uint64_t` if it
  doesn't dramatically affect speed.

* Start writing the algorithm proper.


