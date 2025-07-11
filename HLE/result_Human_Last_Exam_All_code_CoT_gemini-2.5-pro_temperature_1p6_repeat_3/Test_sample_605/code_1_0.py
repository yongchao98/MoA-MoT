# 1. Define the weights for the ambient space.
weights = [22, 29, 49, 50, 75]

# 2. Calculate the sum of the weights.
sum_of_weights = sum(weights)

# 3. Determine the degree 'd' based on the Calabi-Yau condition.
# For a Calabi-Yau variety, the degree 'd' must equal the sum of the weights.
d = sum_of_weights

# 4. Define and calculate the Crawley-Nordstrom invariant.
# The formula is: c_N = (sum of weights) - 2*d
cn_invariant = sum_of_weights - 2 * d

# 5. Print the calculation steps and the final result.
print(f"The weights are: {weights}")
print("For the link to be Calabi-Yau, the degree 'd' must equal the sum of the weights.")
print(f"Sum of weights = {' + '.join(map(str, weights))} = {sum_of_weights}")
print(f"Therefore, the degree d = {d}.")
print("\nThe Crawley-Nordström invariant (c_N) is calculated as: c_N = (sum of weights) - 2 * d")
print("\nSubstituting the values into the formula:")
# Final output showing the full equation as requested
print(f"c_N = {sum_of_weights} - 2 * {d}")
print(f"c_N = {sum_of_weights} - {2 * d}")
print(f"c_N = {cn_invariant}")