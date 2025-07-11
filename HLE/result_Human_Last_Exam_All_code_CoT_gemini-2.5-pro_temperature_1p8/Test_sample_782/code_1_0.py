# Cost of primitive operations in terms of field multiplications
cost_convert_affine_to_extended = 1
cost_madd = 10  # Mixed addition
cost_dbl = 8    # Doubling

# Step 1: Compute A + C where C = -B. A and C are affine.
# Convert A to extended (1M) and then add C (affine) using madd (10M).
cost_step1 = cost_convert_affine_to_extended + cost_madd

# Step 2: Double the result of Step 1.
cost_step2 = cost_dbl

# Step 3: Add C (affine) to the result of Step 2 using madd.
cost_step3 = cost_madd

# Total cost is the sum of all steps
total_cost = cost_step1 + cost_step2 + cost_step3

print(f"The total cost is the sum of computing (A+C), doubling the result, and adding C back.")
print(f"Cost of (A+C): {cost_step1}M")
print(f"Cost of doubling: {cost_step2}M")
print(f"Cost of final addition: {cost_step3}M")
print(f"Total Cost Equation: {cost_step1} + {cost_step2} + {cost_step3} = {total_cost}")