import math

# The derivation leads to the variance being a rational number p/q.
# E[Y] = 1/2
# E[Y^2] = 12/35
# Var(Y) = E[Y^2] - (E[Y])^2 = 12/35 - (1/2)^2 = 12/35 - 1/4
p_var = 48 - 35
q_var = 140

# We can simplify the fraction if needed, but 13 is prime.
common_divisor = math.gcd(p_var, q_var)
p_final = p_var // common_divisor
q_final = q_var // common_divisor

print("The variance of Y can be expressed as a fraction.")
print(f"Let the variance be Var(Y) = E[Y^2] - (E[Y])^2.")
print(f"E[Y^2] = 12/35")
print(f"E[Y] = 1/2")
print(f"Var(Y) = 12/35 - (1/2)^2 = 12/35 - 1/4")
print(f"Var(Y) = (48 - 35) / 140 = {p_var}/{q_var}")
print(f"The final simplified fraction is {p_final}/{q_final}.")