import math

# Step 1: Calculate the values of F_3(x) at the integral limits x=1 and x=0.

# Calculation for x = 1
f3_1 = math.log(1 + math.exp(math.sin(1)))

# Calculation for x = 0
# F_3(0) = ln(1 + e^(sin(0))) = ln(1 + e^0) = ln(1 + 1) = ln(2)
f3_0 = math.log(2)

# Step 2: The definite integral V is ln(F_3(1)) - ln(F_3(0)).
# Using logarithm properties, V = ln(F_3(1) / F_3(0)).
V = math.log(f3_1 / f3_0)

# Step 3: The final requested value is the closest integer to 10000 * V.
final_value = 10000 * V
rounded_value = round(final_value)

# Step 4: Print the numbers involved in the final equation as requested.
print("--- Integral Evaluation ---")
print(f"The integral is evaluated as ln(F₃(1)) - ln(F₃(0)).")
print("The components of the final equation are:")
print(f"F₃(1) = ln(1 + e^(sin(1))) = {f3_1}")
print(f"F₃(0) = ln(2) = {f3_0}")
print(f"Value of the integral, V = ln(F₃(1) / F₃(0)) = {V}")
print("\n--- Final Answer Calculation ---")
print(f"The target value is 10000 * V = {final_value}")
print(f"The closest integer is: {rounded_value}")
