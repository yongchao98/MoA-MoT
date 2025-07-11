import math

# Given parameters
# k1 and k2 are very large, so we will work with their logarithms.
# k1 = 10^3000 -> log10(k1) = 3000
log10_k1 = 3000
# k2 = 10^500 -> log10(k2) = 500
log10_k2 = 500
lambda1 = 0.5
lambda2 = 0.5
# |||h||| = 1000 -> log10(|||h|||) = 3
log10_h_norm = 3

# Step 1: Calculate the value of L_plus = lim_sup ||x_n|| as n -> +inf
# L_plus = (k2 / (1 - lambda2)) * |||h|||
# log10(L_plus) = log10(k2) - log10(1 - lambda2) + log10(|||h|||)
log10_L_plus = log10_k2 - math.log10(1 - lambda2) + log10_h_norm
# We can print the full expression for L_plus for clarity
# 1 - lambda2 = 0.5
# L_plus = 10^500 / 0.5 * 1000 = 2 * 10^500 * 10^3 = 2e503
# The coefficients are C_L_plus = 1 / (1-lambda2) * 1 = 2
log10_C_L_plus = -math.log10(1-lambda2)

# Step 2: Calculate the value of L_minus_minus = lim_inf ||x_n|| as n -> -inf
# We assume lim_inf = lim_sup as n -> -inf
# L_minus_minus = (k1 / lambda1) * |||h|||
# log10(L_minus_minus) = log10(k1) - log10(lambda1) + log10(|||h|||)
log10_L_minus_minus = log10_k1 - math.log10(lambda1) + log10_h_norm
# We can print the full expression for L_minus_minus for clarity
# lambda1 = 0.5
# L_minus_minus = 10^3000 / 0.5 * 1000 = 2 * 10^3000 * 10^3 = 2e3003
# The coefficients are C_L_minus = 1 / lambda1 * 1 = 2
log10_C_L_minus = -math.log10(lambda1)


# Step 3: Calculate the two terms of the final expression.
# The expression is: 100 * log10(L_plus / 3) + 10 * log10(L_minus_minus / 3)

# First term: 100 * log10(L_plus / 3) = 100 * (log10(L_plus) - log10(3))
term1 = 100 * (log10_L_plus - math.log10(3))

# Second term: 10 * log10(L_minus_minus / 3) = 10 * (log10(L_minus_minus) - log10(3))
term2 = 10 * (log10_L_minus_minus - math.log10(3))

# Step 4: Calculate the final result
result = term1 + term2

# Print the results step-by-step
print("This script calculates the value of the expression based on the properties of a difference equation admitting a discrete dichotomy.")
print("\nParameters:")
print(f"k1 = 10^{log10_k1}, k2 = 10^{log10_k2}, lambda1 = {lambda1}, lambda2 = {lambda2}, |||h||| = 10^{log10_h_norm}")

print("\nStep 1: Calculate L_plus = sup lim ||x_n|| as n -> +infinity")
# Using the derived form log10(2) + 503
print(f"L_plus = (k2 / (1 - lambda2)) * |||h||| = (10^{log10_k2} / {1 - lambda2}) * 10^{log10_h_norm} = 2 * 10^{log10_k2 + log10_h_norm}")
print(f"L_plus = 2e{log10_k2 + log10_h_norm}")
log10_L_plus_val = math.log10(2) + log10_k2 + log10_h_norm
print(f"log10(L_plus) = {log10_L_plus_val:.4f}")

print("\nStep 2: Calculate L_minus_minus = inf lim ||x_n|| as n -> -infinity")
# Using the derived form log10(2) + 3003
print(f"L_minus_minus = (k1 / lambda1) * |||h||| = (10^{log10_k1} / {lambda1}) * 10^{log10_h_norm} = 2 * 10^{log10_k1 + log10_h_norm}")
print(f"L_minus_minus = 2e{log10_k1 + log10_h_norm}")
log10_L_minus_minus_val = math.log10(2) + log10_k1 + log10_h_norm
print(f"log10(L_minus_minus) = {log10_L_minus_minus_val:.4f}")

print("\nStep 3: Calculate the final expression: 100 * log10(L_plus / 3) + 10 * log10(L_minus_minus / 3)")
print(f"First term = 100 * (log10(L_plus) - log10(3)) = 100 * ({log10_L_plus_val:.4f} - {math.log10(3):.4f}) = {term1:.4f}")
print(f"Second term = 10 * (log10(L_minus_minus) - log10(3)) = 10 * ({log10_L_minus_minus_val:.4f} - {math.log10(3):.4f}) = {term2:.4f}")

print("\nFinal Calculation:")
print(f"Result = {term1:.4f} + {term2:.4f}")
print(f"Result = {result:.4f}")