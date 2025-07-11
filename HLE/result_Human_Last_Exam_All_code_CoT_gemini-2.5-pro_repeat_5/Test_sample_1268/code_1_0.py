import math

# The problem is to find the upper bound for k_k,inf in relation to the covolume V.
# This corresponds to finding the constant C in the inequality M_K <= C * V, where M_K
# is the Minkowski bound for a quadratic number field K = Q(sqrt(N)).

# The general formulas are:
# Minkowski bound M_K = (4/pi)^r2 * (n! / n^n) * sqrt(|d_K|)
# Covolume V = 2^(-r2) * sqrt(|d_K|)
# The constant C is their ratio: C = M_K / V = ((4/pi)^r2 * (n! / n^n)) / (2^(-r2))

# For any quadratic field, the degree n is 2.
n = 2
print(f"For a quadratic field, the degree n = {n}.\n")

# --- Case 1: Real Quadratic Field (e.g., N > 0) ---
# For a real quadratic field, there are r1 = 2 real embeddings and r2 = 0 complex embeddings.
r1_real = 2
r2_real = 0
print("--- Case 1: Real Quadratic Field (N > 0) ---")
print(f"Here, r1 = {r1_real} and r2 = {r2_real}.")

# Calculate the constant C for the real case
c_real_numerator = (math.factorial(n) / (n**n)) # (4/pi)^0 is 1
c_real_denominator = (2**(-r2_real)) # 2^0 is 1
c_real = c_real_numerator / c_real_denominator
print(f"The constant C_real = ({n}! / {n}^{n}) / (2^{(-{r2_real})}) = {c_real}")
print(f"For real quadratic fields, the relationship is: k_k,inf <= {c_real} * V\n")


# --- Case 2: Imaginary Quadratic Field (e.g., N < 0) ---
# For an imaginary quadratic field, there are r1 = 0 real embeddings and r2 = 1 pair of complex embeddings.
r1_imag = 0
r2_imag = 1
print("--- Case 2: Imaginary Quadratic Field (N < 0) ---")
print(f"Here, r1 = {r1_imag} and r2 = {r2_imag}.")

# Calculate the constant C for the imaginary case
c_imag_numerator = ((4 / math.pi)**r2_imag * (math.factorial(n) / (n**n)))
c_imag_denominator = (2**(-r2_imag))
c_imag = c_imag_numerator / c_imag_denominator
print(f"The constant C_imaginary = ((4/pi)^{r2_imag} * ({n}! / {n}^{n})) / (2^{(-{r2_imag})}) = (4/pi) * (1/2) / (1/2) = 4/pi")
print(f"Numerically, C_imaginary is approximately {c_imag:.5f}")
print(f"For imaginary quadratic fields, the relationship is: k_k,inf <= (4/pi) * V\n")


# --- General Upper Bound for all Quadratic Fields ---
# The general upper bound must hold for both cases, so we take the larger of the two constants.
print("--- General Upper Bound ---")
print(f"To find a single bound for all cases, we compare the constants: {c_real} (real) and {c_imag:.5f} (imaginary).")
print("The larger constant provides the universal upper bound.")

# The final equation, showing each number as requested.
final_numerator = 4
final_denominator_symbol = "pi"

print("\nThe final upper bound relationship is:")
print(f"k_k,inf <= ({final_numerator} / {final_denominator_symbol}) * V")
