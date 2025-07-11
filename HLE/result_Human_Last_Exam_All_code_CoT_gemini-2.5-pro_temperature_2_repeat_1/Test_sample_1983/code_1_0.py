import math

# Step 1: Define the given parameters.
k1 = 10**3000
k2 = 10**500
h_norm = 1000
lambda2 = 0.5

# Step 2: Address the contradiction and assume a corrected value for lambda1.
# The given condition lambda2 = 0.5 * lambda1 = 0.5 implies lambda1 = 1.
# This contradicts the definition of discrete dichotomy in the reference paper, which requires lambda1 > 1.
# We assume a typo and that the intended relation was lambda1 = 1 / lambda2.
lambda1 = 1 / lambda2  # This gives lambda1 = 2, which satisfies lambda1 > 1.

print(f"Given parameters:")
print(f"k1 = 10^{math.log10(k1)}")
print(f"k2 = 10^{math.log10(k2)}")
print(f"|||h||| = {h_norm}")
print(f"lambda2 = {lambda2}")
print(f"Assumed lambda1 = {lambda1}\n")

# Step 3: Calculate the bounds on the solution's norm based on Theorem 3.2 from the reference.
# Let's denote X_plus = lim sup ||x_n|| as n -> +inf and X_minus = lim inf ||x_n|| as n -> -inf
X_plus = (k1 * h_norm) / (lambda1 - 1)
X_minus = (k2 * h_norm) / (1 - lambda2)

print(f"Calculating the bounds on the norm:")
print(f"X_plus = lim sup ||x_n|| = (k1 * |||h|||) / (lambda1 - 1)")
print(f"       = (10^3000 * 1000) / (2 - 1) = 10^{int(math.log10(X_plus))}")
print(f"X_minus = lim inf ||x_n|| = (k2 * |||h|||) / (1 - lambda2)")
# Use logarithms for safe calculation and printing
X_minus_log10 = math.log10(k2) + math.log10(h_norm) - math.log10(1-lambda2)
print(f"        = (10^500 * 1000) / (1 - 0.5) = 2 * 10^{int(X_minus_log10)}\n")


# Step 4: Calculate the two terms in the expression.
# The full expression is: 100 * log10(X_plus/3) + 10 * log10(X_minus/3)
lim_sup_term = math.log10(X_plus / 3)
lim_inf_term = math.log10(X_minus / 3)

# To fulfill the requirement "output each number in the final equation"
# We can represent the final equation with the computed terms.
final_eq_str = f"100 * ({lim_sup_term:.4f}) + 10 * ({lim_inf_term:.4f})"
print("Final Equation to calculate: 100 * lim_sup_term + 10 * lim_inf_term")
print(f"lim_sup_term = log10(X_plus/3) = log10(10^3003 / 3) = {lim_sup_term:.4f}")
print(f"lim_inf_term = log10(X_minus/3) = log10(2*10^503 / 3) = {lim_inf_term:.4f}\n")


# Step 5: Compute the final result.
result = 100 * lim_sup_term + 10 * lim_inf_term

print(f"The calculation is: {final_eq_str}")
print(f"Final Result: {result}")
