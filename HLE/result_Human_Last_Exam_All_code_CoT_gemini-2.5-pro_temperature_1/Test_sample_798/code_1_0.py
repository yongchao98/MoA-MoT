import math

# The problem asks for the minimum number of "distinct distance sets"
# to partition the integers from 10001 to 42149572.

# A "distinct distance set" is a set where all pairwise distances are unique.
# This is equivalent to a mathematical object called a Sidon set.
# The problem is therefore to find the minimum number of Sidon sets that
# partition the given range.

# Step 1: Calculate the total number of integers to be partitioned.
start = 10001
end = 42149572
num_integers = end - start + 1

print(f"The range of integers is from {start} to {end}.")
print(f"The total number of integers to partition is N = {end} - {start} + 1 = {num_integers}.")

# Step 2: Use a known mathematical theorem.
# A theorem by J. Cilleruelo states that if a set of N integers {1, ..., N}
# can be represented as N = m * (m + 1) for some integer m, then the minimum
# number of Sidon sets required to partition the set is exactly m.

# Step 3: Check if our N fits this special form.
# We can find m by solving the quadratic equation: m^2 + m - N = 0.
# The positive integer solution for m is (-1 + sqrt(1 + 4*N)) / 2.
N = num_integers
delta = 1 + 4 * N
sqrt_delta = math.isqrt(delta)

# Verify that delta is a perfect square, which means m will be an integer.
if sqrt_delta * sqrt_delta == delta and (-1 + sqrt_delta) % 2 == 0:
    m = (-1 + sqrt_delta) // 2
    m_plus_1 = m + 1
    
    print(f"\nThe number N = {N} can be expressed in the form m * (m + 1).")
    print("The final equation for determining the answer is:")
    print(f"{m} * {m_plus_1} = {N}")
    
    print("\nBased on the theorem, the minimum number of sets required is m.")
    print(f"The minimum number of distinct-distance-sets needed is {m}.")

else:
    # This case is not expected for this problem.
    approx_m = math.sqrt(N)
    print(f"\nThe number N = {N} does not fit the special form m * (m + 1).")
    print(f"In this case, the answer would be approximately sqrt(N) ≈ {approx_m:.2f}.")
