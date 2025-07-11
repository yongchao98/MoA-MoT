# This script calculates the total cost in multiplications.

# Cost of Step 1: Compute C = A - B
# 1M (to convert A to extended) + 9M (for mixed addition)
cost_step1 = 10

# Cost of Step 2: Compute D = 2C
# 7M (for doubling in extended coordinates)
cost_step2 = 7

# Cost of Step 3: Compute Result = D - B
# 9M (for mixed addition)
cost_step3 = 9

# Calculate the total cost
total_cost = cost_step1 + cost_step2 + cost_step3

# Print the final equation showing each component of the cost
print(f"The cost of the most efficient computation is broken down as follows:")
print(f"Cost(A-B) + Cost(2*C) + Cost(D-B) = Total Cost")
print(f"{cost_step1}M + {cost_step2}M + {cost_step3}M = {total_cost}M")