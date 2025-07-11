import math

# Alice's winning probability, p_A
p_A = 0.95

# The equation we need to solve for q is:
# q = ((1 - (1 - p_A)**(1/3)) / p_A**3)**(1/4)
# We will calculate the value step-by-step, showing the numbers in the equation.

print(f"We are solving for q with Alice's winning probability p_A = {p_A}")
print("-" * 20)

# Numerator calculation
# Step 1: 1 - p_A
val1 = 1 - p_A
print(f"1 - p_A = {val1}")

# Step 2: (1 - p_A)^(1/3)
val2 = val1**(1/3)
print(f"(1 - p_A)^(1/3) = {val2}")

# Step 3: 1 - (1 - p_A)^(1/3)
numerator = 1 - val2
print(f"Numerator = 1 - (1 - p_A)^(1/3) = {numerator}")
print("-" * 20)

# Denominator calculation
# Step 4: p_A^3
denominator = p_A**3
print(f"Denominator = p_A^3 = {denominator}")
print("-" * 20)

# q^4 calculation
# Step 5: q^4 = Numerator / Denominator
q_fourth = numerator / denominator
print(f"q^4 = Numerator / Denominator = {q_fourth}")
print("-" * 20)

# q calculation
# Step 6: q = (q^4)^(1/4)
q0 = q_fourth**(1/4)
print(f"q = (q^4)^(1/4) = {q0}")
print("-" * 20)

# Final answer calculation
# Step 7: floor(100 * q)
final_answer = math.floor(100 * q0)
print(f"The final answer is floor(100 * q) = {final_answer}")

print(f"\nFinal result: {final_answer}")