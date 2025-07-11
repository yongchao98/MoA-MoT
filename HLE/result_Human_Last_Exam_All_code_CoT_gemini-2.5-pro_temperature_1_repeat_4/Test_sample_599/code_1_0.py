# The term to find in the sequence of segmented numbers.
term_number = 50

# As determined by the analysis, the segmented numbers are powers of 2.
# The n-th segmented number is given by the formula 2^(n-1).
# For the 50th term, we calculate 2^(50-1).
base = 2
exponent = term_number - 1

# Calculate the result. Python handles large integers automatically.
result = base ** exponent

# Print the final equation with all numbers involved.
print(f"The {term_number}th segmented number = {base}^{exponent} = {result}")