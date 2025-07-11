import math

# The integral of F₃'(x)/F₃(x) from 0 to 1 is V = ln(F₃(1)) - ln(F₃(0)).
# We need to calculate F₃(0) and F₃(1) first.

# Calculate F₃(0)
# F₃(0) = ln(1 + e^(sin(0))) = ln(1 + e^0) = ln(2)
F3_0 = math.log(2)

# Calculate F₃(1)
# F₃(1) = ln(1 + e^(sin(1)))
# Note: math.sin() takes the angle in radians.
F3_1 = math.log(1 + math.exp(math.sin(1)))

# Calculate the value of the integral V
V = math.log(F3_1) - math.log(F3_0)

# Calculate the final result: the closest integer to 10000 * V
final_result = round(10000 * V)

# As requested, output the numbers in the final equation for V
print("The integral V is calculated as ln(F₃(1)) - ln(F₃(0)).")
print(f"F₃(0) = ln(2) ≈ {F3_0}")
print(f"F₃(1) = ln(1 + e^(sin(1))) ≈ {F3_1}")
print("\nThe final equation for V with these numbers is:")
print(f"V = ln({F3_1}) - ln({F3_0})")
print(f"V ≈ {V}")
print("\nThe value of 10000 * V is approximately:")
print(10000 * V)
print("\nThe closest integer to 10000 * V is:")
print(final_result)
