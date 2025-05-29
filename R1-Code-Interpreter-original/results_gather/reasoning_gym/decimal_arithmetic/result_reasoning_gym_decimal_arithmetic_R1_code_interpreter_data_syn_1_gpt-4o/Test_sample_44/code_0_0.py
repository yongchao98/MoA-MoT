# Define the values
numerator = 0.865
denominator = 4.615
multiplier = 8.024
initial_value = 6.384

# Perform the calculations
division_result = numerator / denominator
multiplication_result = division_result * multiplier
final_result = initial_value - multiplication_result

# Print the result rounded to 12 significant digits
print(f"{final_result:.12g}")