# Define the value of z
z = 0.2

# The PDF of Z is f(z) = 6 * z * (1 - z)
# We need to calculate this for z = 0.2
val1 = 6
val2 = z
val3 = 1 - z

# Calculate the final result
result = val1 * val2 * val3

# Print the equation with the numbers and the final result
print(f"{val1} * {val2} * {val3:.1f} = {result}")