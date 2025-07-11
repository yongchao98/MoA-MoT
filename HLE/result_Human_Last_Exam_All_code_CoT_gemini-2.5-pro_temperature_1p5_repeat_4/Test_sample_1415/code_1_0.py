import math

# The problem is to evaluate the definite integral of F₃'(x)/F₃(x) from 0 to 1.
# The integral of f'(x)/f(x) dx is ln|f(x)|.
# So, ∫ (F₃'(x)/F₃(x)) dx from 0 to 1 = ln|F₃(1)| - ln|F₃(0)|.
# Let's call the value of this integral V.

# Step 1: Define and calculate F₃(0)
# F₀(0) = 0
# F₁(0) = sin(0) = 0
# F₂(0) = e⁰ = 1
# F₃(0) = ln(1 + 1) = ln(2)
f3_0 = math.log(2)

# Step 2: Define and calculate F₃(1)
# F₀(1) = 1
# F₁(1) = sin(1)
# F₂(1) = e^(sin(1))
# F₃(1) = ln(1 + e^(sin(1)))
f3_1 = math.log(1 + math.exp(math.sin(1)))

# Step 3: Evaluate V = ln(F₃(1)) - ln(F₃(0))
# The arguments to ln, F₃(0) and F₃(1), are both positive, so we can ignore the absolute value.
V = math.log(f3_1) - math.log(f3_0)

# Step 4: Calculate 10000 * V and find the closest integer
final_value = 10000 * V
closest_integer = round(final_value)

# Print the evaluation steps with the numbers in the final equation
print("The value of the integral V is calculated as: ln(F₃(1)) - ln(F₃(0))")
print("-----------------------------------------------------------------")
print(f"Calculating F₃(0):")
print(f"F₃(0) = ln(2) ≈ {f3_0}")
print("\nCalculating F₃(1):")
print(f"F₃(1) = ln(1 + e^(sin(1))) ≈ {f3_1}")
print("\nCalculating V:")
print(f"V = ln({f3_1}) - ln({f3_0})")
print(f"V ≈ {V}")
print("-----------------------------------------------------------------")
print("\nCalculating the final result:")
print(f"10000 * V ≈ 10000 * {V} = {final_value}")
print(f"\nThe closest integer to 10000 * V is: {closest_integer}")

print(f"\n<<<{closest_integer}>>>")