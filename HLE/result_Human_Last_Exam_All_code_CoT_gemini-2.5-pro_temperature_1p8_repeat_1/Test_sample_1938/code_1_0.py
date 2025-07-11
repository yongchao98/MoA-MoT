import math

# Let p_A be Alice's winning probability. We are given that p_A >= 0.95.
# To find the minimum q, we set p_A = 0.95.
p_A = 0.95

# The relationship between p_A and q is given by the equation:
# p_A = 1 - (1 - q^4 * p_A^3)^3
# We solve this equation for q.
# q = ((1 - (1 - p_A)^(1/3)) / p_A^3)^(1/4)

# We will now compute the values step-by-step as per the formula.

print(f"Let Alice's winning probability p_A = {p_A}")

# Step 1: Calculate p_A^3
p_A_cubed = p_A**3
print(f"p_A^3 = {p_A}^3 = {p_A_cubed}")

# Step 2: Calculate 1 - p_A
one_minus_p_A = 1 - p_A
print(f"1 - p_A = 1 - {p_A} = {one_minus_p_A}")

# Step 3: Calculate (1 - p_A)^(1/3)
one_minus_p_A_cubert = one_minus_p_A**(1/3)
print(f"(1 - p_A)^(1/3) = {one_minus_p_A}^(1/3) = {one_minus_p_A_cubert}")

# Step 4: Calculate the numerator 1 - (1 - p_A)^(1/3)
numerator = 1 - one_minus_p_A_cubert
print(f"Numerator = 1 - {one_minus_p_A_cubert} = {numerator}")

# Step 5: Calculate q^4
q_pow_4 = numerator / p_A_cubed
print(f"q^4 = Numerator / p_A^3 = {numerator} / {p_A_cubed} = {q_pow_4}")

# Step 6: Calculate q
q0 = q_pow_4**(1/4)
print(f"q = (q^4)^(1/4) = {q_pow_4}^(1/4) = {q0}")

# Step 7: Calculate 100 * q
val = 100 * q0
print(f"100 * q = 100 * {q0} = {val}")

# Step 8: Calculate floor(100 * q)
result = math.floor(val)
print(f"The floor of 100 * q is {result}")
print("\nFinal Answer:")
print(result)
