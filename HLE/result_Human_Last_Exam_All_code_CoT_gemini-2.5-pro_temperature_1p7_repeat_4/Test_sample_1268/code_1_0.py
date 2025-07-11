# The problem seeks the relationship between the "maximum norm" (k_k,∞) and
# the covolume (V) for a quadratic field K = Q(sqrt(N)).
#
# We interpret k_k,∞ as the Minkowski bound, M_K, which is a fundamental
# quantity in algebraic number theory. The covolume V is the square root
# of the field's discriminant, V = sqrt(|d_K|).
#
# The formula for the Minkowski bound is:
# M_K = (4/π)^s * (n! / n^n) * sqrt(|d_K|)
#
# For a real quadratic field, the degree n=2 and the number of complex
# embeddings s=0. Plugging these values in:
# M_K = (4/π)^0 * (2! / 2^2) * sqrt(|d_K|)
# M_K = 1 * (2 / 4) * sqrt(|d_K|)
# M_K = (1/2) * sqrt(|d_K|)
#
# Since we defined k_k,∞ = M_K and V = sqrt(|d_K|), the relationship is:
# k_k,∞ = (1/2) * V
#
# The following code prints this final relationship.

# Define the components of the equation as strings and numbers
var_k = "k_{k,∞}"
equality_symbol = "="
fraction_numerator = 1
fraction_denominator = 2
var_V = "V"

# Print the final equation as requested, showing each number
print(f"The upper bound relationship is given by the equation:")
print(f"{var_k} {equality_symbol} ({fraction_numerator}/{fraction_denominator}) * {var_V}")