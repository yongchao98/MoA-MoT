# Let's calculate the maximal number of complex zeros for N=5.
N = 5

# The formula derived is N * 2^(N-1).
# First, calculate the exponent.
exponent = N - 1

# Then, calculate the value of 2 raised to that exponent.
power_of_two = 2**exponent

# Finally, multiply by N to get the result.
max_zeros = N * power_of_two

# Print the calculation step-by-step, showing each number.
print(f"For N = {N}, the maximal number of complex zeros is:")
print(f"{N} * 2^({N} - 1) = {N} * 2^{exponent} = {N} * {power_of_two} = {max_zeros}")