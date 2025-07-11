import math

# Question 1: How many lines of code have complier errors?
# Line 3: `unsigned long long` is not a defined type.
# Line 4: `scanf` uses "%d" for `unsigned long` pointers (should be "%n").
# Line 9: `printf` uses "%d" for an `unsigned long` value (should be "%n").
# Total lines with errors:
q1_errors = 3

# Question 2: What is the total of memory (in D) used for variables?
# The most efficient solution requires 3 variables to store n, m, and a.
# n and m can be up to 4*10^9, which requires the `unsigned long` type (12D).
# `a` should also be `unsigned long` for calculations.
# Total memory = 3 variables * 12D/variable
q2_memory = 3 * 12

# Question 3: Can this problem be solved without branching statements?
# Yes, by using ceiling division formula: (numerator + denominator - 1) / denominator
q3_branching = 'Y'

# Question 4: What is the expected output m = n = 4*10^9 and a = 1.
n = 4 * 10**9
m = 4 * 10**9
a = 1
# Using the ceiling division formula for integer arithmetic
# math.ceil can also be used here as Python handles large integers.
na = math.ceil(n / a)
ma = math.ceil(m / a)
q4_output = na * ma

# Print the final answer in the specified format x:y:z:t
print(f"{q1_errors}:{q2_memory}:{q3_branching}:{q4_output}")
