import math

# This script will calculate the answers to the four questions and print them
# in the specified format.

# --- Question 1: How many lines of code have compiler errors? ---
# The following lines in the original C code have errors based on XVM rules:
# 1. `unsigned long long n = 0, m = 0, a = 0, na = 0, ma = 0;`
#    - Reason: The data type `unsigned long long` does not exist in XVM.
# 2. `scanf("%d %d %d", &n, &m, &a);`
#    - Reason: The format specifier `%d` is for a 1D `digit` type, which is incorrect for the variables.
# 3. `printf("%d",na * ma);`
#    - Reason: The format specifier `%d` is for a 1D `digit` type, which is incorrect for the result.
# Total lines with errors:
answer1 = 3

# --- Question 2: Memory usage in an optimized version ---
# An optimized version would use the fewest variables and no branches.
# We need variables for n, m, and a. We can reuse n and m to store the results.
# n: unsigned long (to hold up to 4*10^9) -> 12D
# m: unsigned long (to hold up to 4*10^9) -> 12D
# a: unsigned int (to hold up to 10^6) -> 6D
# Total memory = 12 + 12 + 6
answer2 = 30

# --- Question 3: Can this problem be solved without branching statements? ---
# Yes. The logic `ceil(n/a)` can be implemented using integer arithmetic
# as `(n + a - 1) / a`, which avoids using an `if` statement.
answer3 = "Y"

# --- Question 4: What is the expected output for n=m=4*10^9 and a=1? ---
# The number of flagstones is ceil(n/a) * ceil(m/a).
n = 4 * 10**9
m = 4 * 10**9
a = 1
# Python's integer division `//` and math.ceil can be used.
# Or we can use the integer arithmetic trick: (numerator + denominator - 1) // denominator
flagstones_n = (n + a - 1) // a
flagstones_m = (m + a - 1) // a
answer4 = flagstones_n * flagstones_m

# --- Final Output ---
# Print the answers in the format x:y:z:t
print(f"{answer1}:{answer2}:{answer3}:{answer4}")
