import math

# Step 1: Define the function F_3(x)
# F_3(x) = ln(1 + e^(sin(x)))

# Step 2: The definite integral V = ∫ (F_3'(x)/F_3(x)) dx from 0 to 1
# is given by ln(F_3(1)) - ln(F_3(0)).

# Step 3: Calculate the values of F_3(0) and F_3(1).
# These are the numbers in the final equation for the integral's value.
# F_3(0) = ln(1 + e^sin(0)) = ln(1 + e^0) = ln(2)
f3_at_0 = math.log(2)

# F_3(1) = ln(1 + e^sin(1))
# Note: math.sin() takes the angle in radians.
f3_at_1 = math.log(1 + math.exp(math.sin(1)))

print(f"The equation for the integral V is: ln(F_3(1)) - ln(F_3(0))")
print(f"Value of F_3(1): {f3_at_1}")
print(f"Value of F_3(0): {f3_at_0}")

# Step 4: Calculate the value of the integral V
V = math.log(f3_at_1) - math.log(f3_at_0)
print(f"Value of V: {V}")

# Step 5: Calculate 10000 * V and find the closest integer
final_value = 10000 * V
closest_integer = round(final_value)

print(f"The result of 10000 * V is: {final_value}")
print(f"The closest integer is: {closest_integer}")

print("<<<{}>>>".format(closest_integer))