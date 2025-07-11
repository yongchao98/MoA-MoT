import math

# Step 1: Define the functions to calculate F3(x) at x=1 and x=0
def calculate_f3(x):
    """Calculates F₃(x) = ln(1 + e^(sin(x)))"""
    f0 = x
    f1 = math.sin(f0)
    f2 = math.exp(f1)
    f3 = math.log(1 + f2)
    return f3

# Step 2: Calculate F3(1) and F3(0)
f3_at_1 = calculate_f3(1)
f3_at_0 = calculate_f3(0) # This will be ln(2)

# Step 3: Calculate the value of the integral V
# V = ln(F3(1)) - ln(F3(0))
V = math.log(f3_at_1) - math.log(f3_at_0)

# Step 4: Calculate 10000 * V and round to the nearest integer
final_result = 10000 * V
rounded_result = round(final_result)

# Step 5: Print the results as requested
print("The integral is of the form ∫(f'(x)/f(x)) dx, which evaluates to ln(f(x)).")
print("Using the Fundamental Theorem of Calculus, the definite integral is ln(F₃(1)) - ln(F₃(0)).")
print("\nFirst, we calculate the numbers for the final equation:")
print(f"F₃(1) = {f3_at_1}")
print(f"F₃(0) = {f3_at_0}")
print("\nNow, we substitute these numbers into the equation for the integral V:")
print(f"V = ln({f3_at_1}) - ln({f3_at_0})")
print(f"V ≈ {V}")
print(f"\nFinally, we calculate the required value:")
print(f"10000 * V ≈ {final_result}")
print(f"The closest integer to 10000*V is {rounded_result}")
