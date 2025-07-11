import math

# Parameters from the problem description
k1 = 10**3000
k2 = 10**500
lambda1 = 1.0
lambda2 = 0.5
h_norm = 1000.0

# Asymptotic norm calculations based on the hypothesis
# For n -> +inf, we use the parameters k2, lambda2
lim_sup_norm = k2 * h_norm * (lambda2 / (1 - lambda2))

# For n -> -inf, we use the parameters k1, lambda1.
# The standard summation diverges for lambda1=1. We assume the effective multiplier is 1.
lim_inf_norm = k1 * h_norm

# Values to be used in the final formula
val1_num = lim_sup_norm
val2_num = lim_inf_norm

# The constants in the expression
c1 = 100
c2 = 10
denominator = 3.0

# Calculate the two terms in the expression
term1 = c1 * math.log10(val1_num / denominator)
term2 = c2 * math.log10(val2_num / denominator)

# The final result
result = term1 + term2

# We can also compute this by hand using properties of logarithms to verify
# term1 = 100 * (log10(k2*h_norm) - log10(3)) = 100 * (503 - log10(3))
# term2 = 10 * (log10(k1*h_norm) - log10(3)) = 10 * (3003 - log10(3))
# result = 100*503 - 100*log10(3) + 10*3003 - 10*log10(3) = 80330 - 110*log10(3)

print(f"k1 = 10^{math.log10(k1)}")
print(f"k2 = 10^{math.log10(k2)}")
print(f"lambda1 = {lambda1}")
print(f"lambda2 = {lambda2}")
print(f"|||h||| = {h_norm}")
print(f"lim_sup ||x_n|| as n->+inf = {val1_num:.1e}")
print(f"lim_inf ||x_n|| as n->-inf = {val2_num:.1e}")
print("\nFinal Expression:")
print(f"{c1}*log10((1/ {denominator}) * {val1_num:.1e}) + {c2}*log10((1/{denominator}) * {val2_num:.1e})")
print(f"\nResult: {result}")