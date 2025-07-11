import math

# Step 1 & 2: Establish base numbers and naming patterns from the problem description.
total_records = 720

# 'L' is defined as one-sixth of the total records.
# This value is referred to as both "Octavius-lost records" and "lost records" in the constraints.
lost_records = total_records // 6
L = lost_records

# 'N' is the count of distinct naming patterns mentioned.
# 1. "Octavius" (lost)
# 2. "Imperator Caesar Divi Filius"
# 3. "Augustus" and "Caesar" (dual-named)
# 4. "Octavianus" or "Augustus" (single-variant)
# 5. "Imperator Caesar Augustus" (full imperial title)
# 6. "Caesar" only
num_patterns = 6
N = num_patterns

# Other explicitly mentioned record counts from the 'if' clause and setup
caesar_only_records = 80
icdf_records = (total_records - L) // 3


# Step 3: Use the mathematical constraints from the "if" clause to determine S, D, and F.
# These constraints override the narrative proportions (1/4, 2/5, etc.).

# Calculate 'S' based on its constraint: S is a perfect square whose root is L / 4.
s_root = L / 4
single_variant_docs = s_root ** 2
S = int(single_variant_docs)

# For 'D', we have two conditions:
# 1. From the constraint F = sqrt(D) + 3*L, F must be an integer, so D must be a perfect square.
# 2. The sum of all found numbers (L, D, S, F, and others mentioned) must be divisible by N.
# (L + D + S + F + caesar_only_records + icdf_records) % N == 0
# Substituting F = sqrt(D) + 3*L:
# (L + D + S + (sqrt(D) + 3*L) + caesar_only_records + icdf_records) % N == 0
# (D + sqrt(D) + 4*L + S + caesar_only_records + icdf_records) % N == 0
#
# Let's find the value of D. We will test perfect squares.
# The narrative mentioned D was 1/4 of the remaining 600 records (i.e., 150).
# This suggests D should be a perfect square near 150. The closest is 169 (13*13).
# Let's verify if D = 169 fits our divisibility rule.
test_D = 169
test_F = math.sqrt(test_D) + 3 * L
total_sum = L + test_D + S + test_F + caesar_only_records + icdf_records
if total_sum % N == 0:
    dual_named_docs = test_D
    D = int(dual_named_docs)
else:
    # Fallback in case the logic is wrong, though 169 is the intended value.
    # This block should not be reached with the correct logic.
    D = 0 # Placeholder for error

# Calculate 'F' based on its constraint, now that D is known.
full_title_docs = math.sqrt(D) + 3 * L
F = int(full_title_docs)

# Step 4 & 5: Perform the final calculation.
# Product = a * b * c * d
product = L * D * S * F

# Division by the number of patterns
divided_by_patterns = product / N

# Final answer is floor(result / 1000)
final_answer = math.floor(divided_by_patterns / 1000)

# Output the numbers used in the final equation as requested.
print(f"a) Lost records (L): {L}")
print(f"b) Dual-named documents (D): {D}")
print(f"c) Single-variant documents (S): {S}")
print(f"d) Full imperial title documents (F): {F}")
print(f"Divisor (Distinct naming patterns, N): {N}")
print(f"---")
print(f"Product (L * D * S * F) = {L} * {D} * {S} * {F} = {product}")
print(f"Product / N = {product} / {N} = {divided_by_patterns}")
print(f"Answer = floor({divided_by_patterns} / 1000)")
print(f"\nFinal Answer: {final_answer}")
