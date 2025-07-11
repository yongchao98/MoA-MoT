import math

# Parameters for IEEE 754 64-bit floating-point numbers
n = 64
m = 52
e = 11
B = 2**(e - 1) - 1 # Bias B is 1023

print(f"Using IEEE 754 64-bit float as an example: m = {m}, e = {e}, B = {B}")
print("-" * 20)

# Step 1: Calculate the total number of representable normalized numbers in (0, 1].
# This includes numbers where E is from 1 to B-1 (with any mantissa),
# plus the number 1.0 (where E=B and M=0).
num_values = (B - 1) * (2**m) + 1

# Step 2: Calculate the theoretical minimum number of bits required (entropy).
# This is log base 2 of the number of possible values.
bits_required = math.log2(num_values)

print(f"Number of representable values in (0, 1]: N = ({B} - 1) * 2^{m} + 1")
print(f"Minimum bits required = log2(N) = {bits_required:.4f} bits")
print("-" * 20)

# Step 3: Evaluate the formula from choice A.
choice_a_val = m + math.log2(B)
print("Evaluating Choice A: m + log_2(B)")
print(f"The calculation is: {m} + log2({B}) = {choice_a_val:.4f}")
print("-" * 20)

# Step 4: Compare the results.
print("Comparing the theoretical value with Choice A's formula:")
print(f"log2(N) = {bits_required:.4f}")
print(f"m + log2(B) = {choice_a_val:.4f}")
print("\nThe values are almost identical. The number of bits required to pick uniformly from the set of representable numbers in [0,1] is best described by the formula in Choice A.")
