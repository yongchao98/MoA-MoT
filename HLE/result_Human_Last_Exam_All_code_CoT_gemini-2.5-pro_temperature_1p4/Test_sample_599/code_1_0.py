# The problem asks for the 50th segmented number.
n = 50

# Based on the analysis, the nth segmented number is 2^(n-1).
# We calculate the exponent first.
exponent = n - 1

# Calculate the final result.
# Python's ** operator is used for exponentiation and handles large integers automatically.
result = 2**exponent

# Print the final result in the format of an equation,
# showing each number used in the calculation as requested.
print(f"The 50th segmented number is given by the formula 2^(n-1).")
print(f"For n={n}, the equation is 2^({n}-1) = {result}")
