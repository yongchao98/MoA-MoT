# The problem asks for the 50th segmented number.
n = 50

# As determined by the logic, the nth segmented number is 2 to the power of (n-1).
# The sequence of segmented numbers are the powers of 2: 2^0, 2^1, 2^2, ...
# So, for the 50th element, the exponent is 50 - 1.
exponent = n - 1

# Calculate the final result. Python handles large integers automatically.
result = 2**exponent

# Print the final result in an equation format as requested.
# This shows the values used in the calculation.
print(f"The {n}th segmented number is calculated as 2^{exponent}, which equals:")
print(result)