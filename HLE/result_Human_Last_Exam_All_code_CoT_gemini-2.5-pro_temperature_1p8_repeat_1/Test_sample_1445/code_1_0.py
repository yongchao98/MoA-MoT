# Step 1: Calculate the cost to transform S_C ("0101...") to S_B ("111...").
# S_C has 50 blocks of '0's. These must be removed.
cost_remove_zeros_from_sc = 50
# After removing the 50 '0's (total 50 digits), 50 '1's remain.
# We need to make the sequence a block of 100 '1's.
# This requires inserting 50 '1's into the existing merged block of '1's.
cost_insert_ones_for_sb = 1
# Total cost for S_C -> S_B
cost_sc_to_sb = cost_remove_zeros_from_sc + cost_insert_ones_for_sb

# Step 2: Calculate the cost to transform S_B ("111...") to S_D ("1010...").
# S_B has 100 '1's. S_D has 50 '1's and 50 '0's.
# We need to remove 50 '1's. This can be done as a single block.
cost_remove_ones_from_sb = 1
# Now we have a sequence of 50 '1's. We need to insert 50 '0's to form the '1010...' pattern.
# Each '0' splits a block of '1's, so they must be separate operations.
cost_insert_zeros_for_sd = 50
# Total cost for S_B -> S_D
cost_sb_to_sd = cost_remove_ones_from_sb + cost_insert_zeros_for_sd

# Step 3: Calculate the total worst-case cost, n.
# This is the sum of the costs for the path S_C -> S_B -> S_D.
n = cost_sc_to_sb + cost_sb_to_sd

print(f"{cost_sc_to_sb} + {cost_sb_to_sd} = {n}")