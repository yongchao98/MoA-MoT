import math

# --- Parameters of the function field ---
# q: the size of the finite field of constants (must be a prime power, e.g., 2, 3, 4, 5, 7, 8, 9, ...)
# g: the genus of the function field (a non-negative integer)
# h: the ideal class number of the affine ring R (a positive integer)
# Let's use some example values.
q = 2
g = 1
h = 1

# --- Constants ---
# gamma: The Euler-Mascheroni constant
gamma = 0.57721566490153286060651209008240243104215933593992

# --- Calculation Step-by-Step ---
# The formula for the limit is: e^(-gamma) * (q^g * log(q)) / h

# 1. Calculate e^(-gamma)
e_neg_gamma = math.exp(-gamma)

# 2. Calculate q^g
q_pow_g = q**g

# 3. Calculate log(q) (natural logarithm)
log_q = math.log(q)

# 4. The value of h is given
h_val = h

# 5. Compute the final result
result = (e_neg_gamma * q_pow_g * log_q) / h_val

# --- Output the results ---
print(f"For the given parameters (q={q}, g={g}, h={h}):")
print("-" * 30)
print(f"The term e^(-gamma) is: {e_neg_gamma}")
print(f"The term q^g is: {q_pow_g}")
print(f"The term log(q) is: {log_q}")
print(f"The term h is: {h_val}")
print("-" * 30)
print(f"The final equation is: ({e_neg_gamma} * {q_pow_g} * {log_q}) / {h_val}")
print(f"The value of the limit is: {result}")