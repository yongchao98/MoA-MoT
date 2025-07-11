import math

# The problem asks us to represent N = 10^100 as a sum of terms of the form:
# T = tet(2, a) * pow(2, b)
# where tet(2, a) = 2**(2**(a-1)) and the constraint is tet(2, a) > pow(2, b).

# The number to be represented. Python handles large integers automatically.
N = 10**100

# Part 1: Find the count of sequences.
# A greedy decomposition of N is equivalent to expressing N in binary. Each term
# in the sum corresponds to a '1' in the binary representation of N.
# Therefore, the number of terms is the population count (number of set bits) of N.
# Python's .bit_count() method (3.10+) or bin(N).count('1') can find this.
try:
    count = N.bit_count()
except AttributeError:
    # Fallback for Python versions older than 3.10
    count = bin(N).count('1')

# Part 2: Find a1 and b1 for the largest term in the sequence.
# The largest term, T1, corresponds to the most significant bit of N.
# The term is T1 = 2**p, where p is the position of this most significant bit.
# p can be calculated as floor(log2(N)).
p = N.bit_length() - 1

# We need to express T1 = 2**p in the form tet(2, a1) * pow(2, b1).
# 2**p = (2**(2**(a1-1))) * (2**b1)
# This gives the equation for the exponent: p = 2**(a1-1) + b1
# with the constraint b1 < 2**(a1-1).

# This equation has a unique solution for a1 and b1.
# 2**(a1-1) must be the largest power of 2 less than or equal to p.
# This means a1-1 = floor(log2(p)).
# In Python, for a positive integer x, floor(log2(x)) + 1 is equivalent to x.bit_length().
a1 = p.bit_length()

# Now we can calculate b1 from the equation.
# b1 = p - 2**(a1-1)
power_of_2_for_a1 = 1 << (a1 - 1)
b1 = p - power_of_2_for_a1

# The problem asks for the count, a1, and b1.
# The phrase "output each number in the final equation" can be interpreted as showing
# the components of the result. We will print the required values as specified.
# For example, for the largest term, the equation for its exponent 'p' is p = 2**(a1-1) + b1.
# Let's see the values:
# p = 332
# a1 = 9, so 2**(a1-1) = 2**8 = 256
# b1 = 332 - 256 = 76
# Final equation for p: 332 = 256 + 76
# The requested output format is "count a1 b1".

print(f"{count} {a1} {b1}")

# The full decomposition of the number N is:
# N = sum( tet(2, a_i) * pow(2, b_i) ) for i=1 to count
# Our method corresponds to a greedy choice at each step, which is equivalent to
# N = sum(2**p_i) for all p_i where the bit at position p_i in N is 1.
# Each 2**p_i term is then uniquely represented as tet(2, a_i) * pow(2, b_i).

<<<168 9 76>>>