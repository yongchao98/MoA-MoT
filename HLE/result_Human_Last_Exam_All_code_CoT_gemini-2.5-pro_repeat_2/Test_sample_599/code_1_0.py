# The Nth segmented number is given by the formula 2^(N-1).
# We need to find the 50th element of this sequence.

# Define the inputs for our equation
n = 50
base = 2
exponent = n - 1

# Calculate the result
result = base ** exponent

# The final equation is 2^(50-1) = result.
# We print all the numbers from this equation as requested.
print(f"The {n}th segmented number is {base}^({n}-1) which is equal to {result}")