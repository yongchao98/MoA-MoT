# The problem specifies an integer k >= 2.
# Since a specific value for k is not provided, we will use k=2 as an example
# to demonstrate the calculation. The code can be easily adapted for any other k.
k = 2

# Based on the mathematical derivation, the limit is given by the formula 1 - 1/(2*k).
result = 1 - 1 / (2 * k)

# The prompt asks to output each number in the final equation,
# so we format the output to show the calculation.
print(f"For the example case where k = {k}:")
print(f"The calculation is: 1 - 1 / (2 * {k})")
print(f"Result = {result}")