import math

# Define the sequence of functions
def f0(x):
    """F₀(x) = x"""
    return x

def f1(x):
    """F₁(x) = sin(F₀(x))"""
    return math.sin(f0(x))

def f2(x):
    """F₂(x) = e^(F₁(x))"""
    return math.exp(f1(x))

def f3(x):
    """F₃(x) = ln(1 + F₂(x))"""
    return math.log(1 + f2(x))

# The integral of F₃'(x)/F₃(x) from 0 to 1 is V = ln(F₃(1)) - ln(F₃(0)).
# First, let's calculate the values of F₃(x) at the boundaries x=1 and x=0.
f3_at_1 = f3(1)
f3_at_0 = f3(0)

print(f"The integral is V = ln(F₃(1)) - ln(F₃(0))")
print("------------------------------------------")

# Print the components of the equation
print(f"F₃(1) = ln(1 + e^(sin(1))) = {f3_at_1}")
print(f"F₃(0) = ln(1 + e^(sin(0))) = {f3_at_0}")
print("------------------------------------------")

# Calculate V
# The arguments to math.log(), f3_at_1 and f3_at_0, are both positive.
ln_f3_at_1 = math.log(f3_at_1)
ln_f3_at_0 = math.log(f3_at_0)
V = ln_f3_at_1 - ln_f3_at_0

# Print the evaluation of V
print(f"V = ln({f3_at_1}) - ln({f3_at_0})")
print(f"V = {ln_f3_at_1} - ({ln_f3_at_0})")
print(f"V = {V}")
print("------------------------------------------")


# Final step: calculate 10000 * V and find the closest integer
result = 10000 * V
closest_integer = round(result)

print(f"The value to report is the closest integer to 10000 * V.")
print(f"10000 * {V} = {result}")
print(f"The closest integer is: {closest_integer}")

print(f"<<<{closest_integer}>>>")