# Define the numerical constants present in the derived inequality.
# The relationship is k_k,inf <= 1 * V^0.5
coefficient = 1.0
exponent = 0.5

# Print the final equation representing the upper bound.
# This fulfills the requirement to output each number in the equation.
print("The upper bound for the maximum norm (k_k,inf) in relation to the covolume (V) is given by the inequality:")
print(f"k_k,inf <= {coefficient} * V^{exponent}")
