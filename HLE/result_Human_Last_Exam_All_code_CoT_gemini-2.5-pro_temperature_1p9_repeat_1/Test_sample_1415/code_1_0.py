import math

# Step 1: Calculate the values of F₃(x) at x=0 and x=1
# The function is F₃(x) = ln(1 + e^(sin(x)))

# At x=0: F₃(0) = ln(1 + e^sin(0)) = ln(1 + e^0) = ln(2)
f3_at_0 = math.log(2)

# At x=1: F₃(1) = ln(1 + e^sin(1))
f3_at_1 = math.log(1 + math.exp(math.sin(1)))

# Step 2: Evaluate the definite integral V
# ∫ (F₃'(x)/F₃(x)) dx from 0 to 1 = [ln(F₃(x))] from 0 to 1
# V = ln(F₃(1)) - ln(F₃(0)) = ln(F₃(1) / F₃(0))
V = math.log(f3_at_1 / f3_at_0)

# Step 3: Calculate the final required value
# Find the closest integer to 10000 * V
final_result = round(10000 * V)

# Print the final equation with all its numerical components
print("The definite integral ∫ (F₃'(x)/F₃(x)) dx from 0 to 1 evaluates to ln(F₃(1) / F₃(0)).")
print(f"F₃(0) = {f3_at_0}")
print(f"F₃(1) = {f3_at_1}")
print("The final equation for the integral's value V is:")
print(f"V = ln({f3_at_1} / {f3_at_0})")
print(f"V ≈ {V}")
print("\nWe need to find the closest integer to 10000 * V:")
print(f"10000 * V ≈ {10000 * V}")
print(f"The closest integer is: {final_result}")
