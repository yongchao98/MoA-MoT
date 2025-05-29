import math

# Calculate the square root of 0.6
sqrt_value = math.sqrt(0.6)

# Calculate the two possible values of k
k1 = sqrt_value
k2 = -sqrt_value

# Round the results to four decimal places
k1_rounded = round(k1, 4)
k2_rounded = round(k2, 4)

# Print the results
print(f"{k2_rounded}, {k1_rounded}")