import math

# Define the weights and the number of variables
weights = [22, 29, 49, 50, 75]
n = len(weights)

# The degree 'd' of a Calabi-Yau hypersurface is the sum of its weights.
# We confirmed this matches the degree of most terms in the polynomial.
d = sum(weights)

# Calculate the numerator of the Crawley-Nordstrom invariant: d^(n-1)
numerator = d**(n - 1)

# Calculate the denominator: product of all weights
denominator = 1
for w in weights:
    denominator *= w

# Calculate the final invariant
invariant = numerator / denominator

# To present the full equation, we need the weights as strings
weights_str = " * ".join(map(str, weights))

# Print the final result in the requested format
print("The Crawley-Nordstrom Invariant is calculated as d^(n-1) / (product of weights).")
print(f"Given weights: ({', '.join(map(str, weights))})")
print(f"Number of variables (n): {n}")
print(f"Degree (d = sum of weights): {d}")
print("\nFinal Calculation:")
print(f"CN Invariant = {d}^{n-1} / ({weights_str})")
print(f"CN Invariant = {numerator} / {denominator}")
print(f"CN Invariant = {invariant}")

# Simplify the fraction for a more precise representation
# Use math.gcd to find the greatest common divisor and simplify the fraction.
common_divisor = math.gcd(numerator, denominator)
simplified_numerator = numerator // common_divisor
simplified_denominator = denominator // common_divisor

print(f"\nAs a simplified fraction:")
print(f"CN Invariant = {simplified_numerator} / {simplified_denominator}")

# Return the final numerical answer in the specified format
# The problem asks for the numerical invariant itself.
final_answer = invariant
print(f"\n<<< {final_answer} >>>")