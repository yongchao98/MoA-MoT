import math

# Please replace these values with the specific parameters of your function field and ring.
# q: size of the constant field (a prime power)
q = 2
# g: genus of the function field
g = 1
# h: ideal class number of the affine ring R
h = 1

# The Euler-Mascheroni constant
gamma = 0.57721566490153286060651209008240243104215933593992

# Calculate the components of the formula
e_neg_gamma = math.exp(-gamma)
q_pow_g = q**g
log_q = math.log(q)

# Calculate the final result
result = (e_neg_gamma * q_pow_g * log_q) / h

# Output the equation with the calculated values
print(f"Using the values q={q}, g={g}, h={h}:")
print(f"The formula is (e^-gamma * q^g * log(q)) / h")
print("We calculate each part:")
print(f"  e^-gamma       = {e_neg_gamma}")
print(f"  q^g            = {q_pow_g}")
print(f"  log(q)         = {log_q}")
print(f"  h              = {h}")
print("\nThe final equation is:")
print(f"  ({e_neg_gamma} * {q_pow_g} * {log_q}) / {h} = {result}")

print("\nFinal Result:")
print(result)