import math

# The problem is to evaluate the definite integral of F₃'(x)/F₃(x) from 0 to 1.
# The antiderivative of a function in the form f'(x)/f(x) is ln(f(x)).
# Therefore, the value of the integral, V, is ln(F₃(1)) - ln(F₃(0)).

# Define the functions to find F₃(x)
def f3(x):
    """Calculates the value of F₃(x) = ln(1 + e^(sin(x)))"""
    f0_x = x
    f1_x = math.sin(f0_x)
    f2_x = math.exp(f1_x)
    f3_x = math.log(1 + f2_x)
    return f3_x

# Calculate F₃(1) and F₃(0)
f3_at_1 = f3(1)
f3_at_0 = f3(0)

# The value of the integral V is ln(F₃(1)) - ln(F₃(0))
V = math.log(f3_at_1) - math.log(f3_at_0)

# The final requested value is 10000 * V, rounded to the nearest integer
final_answer = round(10000 * V)

# Print the step-by-step evaluation of the final equation
print(f"The value of the integral is given by the equation: V = ln(F₃(1)) - ln(F₃(0))")
print(f"First, we evaluate the terms:")
print(f"F₃(1) = ln(1 + e^(sin(1))) = {f3_at_1}")
print(f"F₃(0) = ln(1 + e^(sin(0))) = {f3_at_0}")
print(f"Now, we substitute these values back into the equation for V:")
print(f"V = ln({f3_at_1}) - ln({f3_at_0})")
print(f"V = {V}")
print(f"Finally, we calculate 10000 * V and round to the closest integer:")
print(f"10000 * {V} = {10000 * V}")
print(f"The closest integer is: {final_answer}")
<<<5488>>>