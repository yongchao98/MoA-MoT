import math

# This script calculates the value of the definite integral and the final requested value.

# 1. Define the function F3(x) based on the problem description.
# F₀(x) = x
# F₁(x) = sin(x)
# F₂(x) = e^(sin(x))
# F₃(x) = ln(1 + e^(sin(x)))
def F3(x):
    """Calculates the value of F₃(x) = ln(1 + e^(sin(x)))"""
    return math.log(1 + math.exp(math.sin(x)))

# 2. The integral of F₃'(x)/F₃(x) from 0 to 1 is V = ln(F₃(1)) - ln(F₃(0)).
# Calculate F₃(x) at the integration limits x=1 and x=0.
F3_at_1 = F3(1)
F3_at_0 = F3(0)

# 3. Calculate the value of the integral V.
V = math.log(F3_at_1) - math.log(F3_at_0)

# 4. Calculate the final value as requested: 10000 * V.
final_value = 10000 * V

# 5. Find the closest integer to the final value.
result = round(final_value)

# 6. Print the steps of the calculation, including the final equation with numbers.
print("The integral is V = ln(F₃(1)) - ln(F₃(0))")
print(f"F₃(1) = ln(1 + e^(sin(1))) = {F3_at_1}")
print(f"F₃(0) = ln(1 + e^(sin(0))) = ln(2) = {F3_at_0}")
print("The final equation with the numerical values is:")
print(f"V = ln({F3_at_1}) - ln({F3_at_0})")
print(f"V = {math.log(F3_at_1)} - ({math.log(F3_at_0)})")
print(f"V = {V}")
print("\nThe final requested value is the closest integer to 10000 * V.")
print(f"10000 * V = {final_value}")
print(f"The closest integer is: {result}")