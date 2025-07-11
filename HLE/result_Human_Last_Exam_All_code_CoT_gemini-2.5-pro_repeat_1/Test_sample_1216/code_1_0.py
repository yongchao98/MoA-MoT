import math

# Step 1-5: Based on the logical derivation in the plan, we establish the values.
# L must be a multiple of 60. The most likely candidate is 1/6 * 720 = 120.
lost_records = 120  # (a)
total_records = 720
readable_records = total_records - lost_records

# Following the structure where 'icdf' records are calculated and removed first.
icdf_records = readable_records / 3
remaining_pool = readable_records - icdf_records

# Calculate the values based on the remaining pool.
dual_named_docs = remaining_pool / 4  # (b)
not_dual_named = remaining_pool - dual_named_docs
single_variant_docs = (2/5) * not_dual_named  # (c)
remaining_after_single = not_dual_named - single_variant_docs
full_title_docs = remaining_after_single / 2  # (d)

# Ensure all calculated values are integers for the final calculation.
a = int(lost_records)
b = int(dual_named_docs)
c = int(single_variant_docs)
d = int(full_title_docs)

print(f"Lost records (a): {a}")
print(f"Documents with dual naming (b): {b}")
print(f"Single-variant documents (c): {c}")
print(f"Full imperial title documents (d): {d}")

# Step 6: Perform the final calculation as requested.
product = a * b * c * d
print(f"\nProduct ({a} * {b} * {c} * {d}) = {product}")

# There are 6 distinct naming patterns mentioned:
# 1. "Octavius" (lost)
# 2. "Imperator Caesar Divi Filius"
# 3. "Augustus" and "Caesar" (dual)
# 4. "Octavianus" or "Augustus" (single)
# 5. "Imperator Caesar Augustus" (full)
# 6. "Only Caesar"
distinct_patterns = 6
print(f"Number of distinct naming patterns: {distinct_patterns}")

division_result = product / distinct_patterns
print(f"Product / Distinct Patterns = {division_result}")

final_answer = math.floor(division_result / 1000)
print(f"\nFinal Answer (floor(result / 1000)): {final_answer}")

print(f"\nEquation steps:")
print(f"Product = {a} * {b} * {c} * {d} = {product}")
print(f"Result = {product} / {distinct_patterns} = {int(division_result)}")
print(f"Final = floor({int(division_result)} / 1000) = {final_answer}")