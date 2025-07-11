import math

# The goal is to find the minimum q such that the probability of Alice winning, p_A, is at least 0.95.
# We solve the equation for q when p_A is exactly 0.95.
p_A = 0.95

# The equation relating q and p_A is: q_0 = [ (1 - (1 - p_A)^(1/3)) / p_A^3 ]^(1/4)
# Let's calculate the components of this equation.

# The value of p_A
val_p_A = p_A
print(f"The value for p_A is: {val_p_A}")

# The value of the term (1 - p_A)
val_one_minus_p_A = 1 - val_p_A
print(f"The value for (1 - p_A) is: {val_one_minus_p_A:.2f}")

# The value of the term (1 - p_A)^(1/3)
val_cbrt = val_one_minus_p_A**(1/3)
print(f"The value for (1 - p_A)^(1/3) is: {val_cbrt}")

# The value of the numerator 1 - (1 - p_A)^(1/3)
val_numerator = 1 - val_cbrt
print(f"The value for the numerator 1 - (1 - p_A)^(1/3) is: {val_numerator}")

# The value of the denominator p_A^3
val_denominator = val_p_A**3
print(f"The value for the denominator p_A^3 is: {val_denominator}")

# The value of q^4
val_q_fourth = val_numerator / val_denominator
print(f"The value for q^4 is: {val_q_fourth}")

# The value of q_0
q0 = val_q_fourth**(1/4)
print(f"The minimum value of q, q_0, is: {q0}")

# The problem asks for floor(100 * q_0)
result = math.floor(100 * q0)
print(f"\nThe final result floor(100 * q_0) is: {result}")
