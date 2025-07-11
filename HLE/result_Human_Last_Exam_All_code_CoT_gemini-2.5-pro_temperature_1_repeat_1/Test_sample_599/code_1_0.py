# This script calculates the 50th segmented number.
# As derived in the plan, the sequence of segmented numbers is the sequence of powers of 2.
# The n-th segmented number is 2^(n-1).
# Therefore, the 50th segmented number is 2^(50-1) = 2^49.

# Define the base and the exponent for the calculation.
base = 2
exponent = 49

# Calculate the result. Python handles large integers automatically.
result = base ** exponent

# Print the final equation with all its components, as requested.
print(f"{base} ** {exponent} = {result}")