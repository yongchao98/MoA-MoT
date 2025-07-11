import math

# Step 1: Define the given constants from the problem description.
# The value k1 is given as a power of 10, so we store its exponent.
log10_k1 = 3000
# The value k2 is given as a power of 10, so we store its exponent.
log10_k2 = 500
# The norm of h.
h_norm = 1000

# Step 2: Formulate the limits based on the working hypothesis.
# We assume the asymptotic norm as n -> +inf is k1 * |||h|||.
# We take the log10 of this value.
log10_L_sup = log10_k1 + math.log10(h_norm)

# We assume the asymptotic norm as n -> -inf is k2 * |||h|||.
# We take the log10 of this value.
log10_L_inf = log10_k2 + math.log10(h_norm)

# Step 3: Calculate the two terms in the expression.
# The first term is 100 * lim_sup[log10(1/3 * ||x_n||)]
# This is equivalent to 100 * log10(1/3 * L_sup)
term1_val = 100 * (math.log10(1/3) + log10_L_sup)

# The second term is 10 * lim_inf[log10(1/3 * ||x_n||)]
# This is equivalent to 10 * log10(1/3 * L_inf)
term2_val = 10 * (math.log10(1/3) + log10_L_inf)

# Step 4: Calculate the final result.
final_result = term1_val + term2_val

# Step 5: Print the breakdown of the calculation for clarity.
print("Based on the hypothesis that the asymptotic norms are proportional to the dichotomy constants and the forcing term norm:")
print(f"Let overline_lim ||x_n|| = k1 * |||h||| = 10^{log10_k1} * {h_norm} = 10^{log10_L_sup}")
print(f"Let underline_lim ||x_n|| = k2 * |||h||| = 10^{log10_k2} * {h_norm} = 10^{log10_L_inf}")
print("\nThe expression to calculate is:")
print(f"100 * log10(1/3 * (10^{log10_L_sup:.0f})) + 10 * log10(1/3 * (10^{log10_L_inf:.0f}))")
print("\nBreaking it down:")
print(f"Term 1 = 100 * (log10(1/3) + {log10_L_sup:.0f}) = 100 * ({math.log10(1/3):.4f} + {log10_L_sup:.0f}) = {term1_val:.4f}")
print(f"Term 2 = 10 * (log10(1/3) + {log10_L_inf:.0f}) = 10 * ({math.log10(1/3):.4f} + {log10_L_inf:.0f}) = {term2_val:.4f}")
print(f"\nTotal = {term1_val:.4f} + {term2_val:.4f} = {final_result:.4f}")
print(f"\nThe final computed value is: {final_result}")
