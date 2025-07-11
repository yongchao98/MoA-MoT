import math

# We will calculate the value of f(n) / (n * log2(n)) for a large n
# to numerically verify the limit we found analytically.
n = 100000

# According to a known theorem, f(n) can be calculated as the sum of
# the number of binary '1's for all integers from 1 to n.
# f(n) = sum_{k=1 to n} s_2(k)
# The function s_2(k) can be implemented in Python using bin(k).count('1').
f_n = sum(bin(k).count('1') for k in range(1, n + 1))

# Now we calculate the terms in the denominator of the limit expression.
log2_n = math.log2(n)
denominator = n * log2_n

# Finally, we calculate the ratio.
ratio = f_n / denominator

# As requested, we print each number in the final equation.
# The equation is: f(n) / (n * log_2(n))
print("To find the limit, we compute the expression for a large value of n.")
print(f"Let's use n = {n}.")
print(f"The corresponding value of f(n) is sum_{k=1..n} s_2(k) = {f_n}.")
print(f"The value of log_2(n) is approximately {log2_n:.6f}.")
print("\nPlugging these into the expression f(n) / (n * log_2(n)), we get:")
print(f"Final Equation: {f_n} / ({n} * {log2_n:.6f})")
print(f"= {f_n} / {denominator:.6f}")
print(f"= {ratio:.8f}")
print("\nThis numerical result is very close to the analytical limit, which is 0.5.")