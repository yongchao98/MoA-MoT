# Parameters from the problem
p = 3
k = 3

# The inequality is p**k - (p-1)*n > 0
# We need to find the largest integer n that satisfies this.
# This is equivalent to n < p**k / (p-1)

# Calculate p^k and p-1
p_k = p**k
p_minus_1 = p - 1

# Calculate the upper bound for n
n_upper_bound = p_k / p_minus_1

# The largest integer n is the floor of this upper bound.
# Note: For strict inequality, if the bound is an integer, we must subtract 1.
# e.g., if n < 14.0, largest n is 13. if n < 13.5, largest n is 13.
# math.floor works for positive numbers, but for robustness with the inequality,
# it's better to use int() or careful logic.
# int(13.5) = 13
# Let's say we have N > a*n. Then n < N/a. The largest integer n is floor(N/a) if N/a is not integer,
# and N/a - 1 if N/a is an integer. Or more simply ceil(N/a) - 1.
import math
n_largest = math.ceil(n_upper_bound) - 1

print(f"The ring is Z/{p**k}.")
print(f"The stability condition is {p}^{k} - ({p}-1)*n > 0.")
print(f"{p_k} - {p_minus_1}*n > 0")
print(f"{p_k} > {p_minus_1}*n")
print(f"{n_upper_bound} > n")
print(f"The largest natural number n satisfying this is {n_largest}.")
<<<13>>>