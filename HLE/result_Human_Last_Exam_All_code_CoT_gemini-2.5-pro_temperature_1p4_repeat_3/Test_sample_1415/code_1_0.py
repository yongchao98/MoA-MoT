import math

# Step 1: Define F₃(x) = ln(1 + e^(sin(x)))
def f3(x):
    return math.log(1 + math.exp(math.sin(x)))

# The integral ∫ (F₃'(x)/F₃(x)) dx from 0 to 1 is V = ln(F₃(1)) - ln(F₃(0)).

# Step 2: Calculate F₃(1) and F₃(0)
val_f3_1 = f3(1)
val_f3_0 = f3(0) # This simplifies to math.log(2)

# Step 3: Calculate the value of the integral V
# V = ln(F₃(1)) - ln(F₃(0))
V = math.log(val_f3_1) - math.log(val_f3_0)

# Step 4: Calculate the final result required by the user
result = 10000 * V
closest_integer = round(result)

# Step 5: Print the outputs showing the equation and numbers
print("The value of the definite integral is V = ln(F₃(1)) - ln(F₃(0)).")
print("\nFirst, we evaluate the terms in the equation:")
print(f"F₃(1) = ln(1 + e^(sin(1))) = {val_f3_1}")
print(f"F₃(0) = ln(2) = {val_f3_0}")

print("\nPlugging these into the equation for V:")
# We use the equation V = ln(F₃(1)) - ln(F₃(0))
print(f"V = ln({val_f3_1}) - ln({val_f3_0})")
print(f"The numerical value is V ≈ {V}")

print("\nFinally, we calculate the requested value:")
print(f"10000 * V ≈ 10000 * {V} = {result}")
print(f"The closest integer to 10000 * V is: {closest_integer}")