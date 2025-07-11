import math

# Calculate the values for the equation
pi_val = math.pi
log2_val = math.log(2)

# Calculate the denominator 2 * log(2)
denominator = 2 * log2_val

# Calculate the final result: pi / (2 * log(2))
result = pi_val / denominator

# The final answer is the expected maximum earthquake magnitude.
# The following print statements show the formula and the calculation steps as requested.
print("The expected maximum earthquake magnitude E[M_N] is given by the equation:")
print(f"E[M_N] = pi / (2 * log(2))")
print("\nSubstituting the values:")
print(f"E[M_N] = {pi_val} / (2 * {log2_val})")
print(f"E[M_N] = {pi_val} / {denominator}")
print("\nFinal Result:")
print(f"E[M_N] = {result}")