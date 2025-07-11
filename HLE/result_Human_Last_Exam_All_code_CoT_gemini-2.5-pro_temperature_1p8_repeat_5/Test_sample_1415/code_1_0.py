import math

# Define the function F3(x) = ln(1 + e^(sin(x)))
def F3(x):
    """Calculates the value of F3(x)."""
    return math.log(1 + math.exp(math.sin(x)))

# Evaluate F3 at the bounds of the integral
F3_at_1 = F3(1)
F3_at_0 = F3(0)

# The integral V = ∫ (F₃'(x)/F₃(x)) dx from 0 to 1 equals ln(F₃(1)) - ln(F₃(0))
V = math.log(F3_at_1) - math.log(F3_at_0)

# Calculate the final value to report
value_to_report = 10000 * V
closest_integer = round(value_to_report)

# Print the final equation with the computed numbers
print(f"The integral evaluates to the expression: ln(F₃(1)) - ln(F₃(0))")
print(f"The value of F₃(1) is: {F3_at_1}")
print(f"The value of F₃(0) is: {F3_at_0}")
print(f"So, the value of the integral V is: ln({F3_at_1}) - ln({F3_at_0}) = {V}")
print(f"We need to find the closest integer to 10000 * V.")
print(f"10000 * V = {value_to_report}")
print(f"The closest integer is: {closest_integer}")

<<<5490>>>