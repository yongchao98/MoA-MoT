import math

# The limit of the sequence is given by the greatest common divisor (GCD)
# of the polynomial values P(k) for all integers k.
# Let's denote this GCD by d.
# From the analysis, d has the form 2^a * 3^b * 5^c.
# The exponents are determined by the minimum p-adic valuations of P(k).

# v_2(d) = min(v_2(P(k))) for all integers k
a = 10

# v_3(d) = min(v_3(P(k))) for all integers k
b = 2

# v_5(d) = min(v_5(P(k))) for all integers k
c = 1

# Calculate the final result
val_a = 2**a
val_b = 3**b
val_c = 5**c

result = val_a * val_b * val_c

# The final equation is d = 2^10 * 3^2 * 5^1
print(f"{val_a} * {val_b} * {val_c} = {result}")
