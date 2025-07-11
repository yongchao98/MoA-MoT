import math

def factorial(n):
    """Calculates the factorial of a non-negative integer."""
    if n < 0:
        return 0
    if n == 0:
        return 1
    return math.factorial(n)

def combinations(n, k):
    """Calculates combinations (n choose k)."""
    if k < 0 or k > n:
        return 0
    return math.comb(n,k)

# Step 1: Calculate the number of ways to choose the rowers
# Choose 2 male scientists from 6: C(6, 2)
# Choose 1 male mathematician from 2: C(2, 1)
num_male_scientists = 6
num_male_mathematicians = 2
ways_to_choose_s_rowers = combinations(num_male_scientists, 2)
ways_to_choose_m_rowers = combinations(num_male_mathematicians, 1)
num_rower_choices = ways_to_choose_s_rowers * ways_to_choose_m_rowers

# Step 2: Calculate permutations for the block of regular classicists
# There are 4 regular classicists, so they can be arranged in 4! ways.
perms_creg = factorial(4)

# Step 3: Calculate permutations for the two main scenarios

# --- Scenario 1: Block_SM is not next to Cassie (neighbors are Ethicists) ---
# There are 4 ways to arrange the blocks and interface people in this configuration.
# Block_SM can be S-M or M-S (2 orientations).
# Internal perms for an unconstrained Block_SM: (10! * 2! * 3!)
# where 10! is for other scientists, 2! for scientist rowers, 3! for other mathematicians.
num_block_arrangements_1 = 4
perms_sm_unconstrained = 2 * (factorial(10) * factorial(2) * factorial(3))
total_ways_scenario1 = num_block_arrangements_1 * perms_sm_unconstrained * perms_creg

# --- Scenario 2: Block_SM is next to Cassie ---
# There are 8 ways to arrange the blocks and interface people in this configuration.
# One end of Block_SM MUST be female. We need the sum of permutations for an S...M block
# with a female at the M-end, and an M...S block with a female at the S-end.
# This corresponds to the two possible orientations of the block when seated next to Cassie.
num_block_arrangements_2 = 8

# Perms for S...M block with female M_end: (10! * 2!) for S-part * (2 * 2!) for M-part
perms_sm_female_m_end = (factorial(10) * factorial(2)) * (2 * factorial(2))

# Perms for M...S block with female S_end: (3!) for M-part * (6 * 9! * 2!) for S-part
perms_sm_female_s_end = factorial(3) * (6 * factorial(9) * factorial(2))

# Total permutations for a constrained Block_SM
perms_sm_constrained = perms_sm_female_m_end + perms_sm_female_s_end
total_ways_scenario2 = num_block_arrangements_2 * perms_sm_constrained * perms_creg

# Final Calculation: Combine all parts
total_ways = num_rower_choices * (total_ways_scenario1 + total_ways_scenario2)

# Print the final equation with all its components
print("The total number of arrangements is calculated as follows:")
print("Ways = (Ways to Choose Rowers) * [ (Ways for Scenario 1) + (Ways for Scenario 2) ]\n")

print("Ways to Choose Rowers = C(6, 2) * C(2, 1)")
print(f"  = {ways_to_choose_s_rowers} * {ways_to_choose_m_rowers} = {num_rower_choices}\n")

print("Scenario 1: Science/Math block (SM) is between two Ethicists.")
print("Ways = (Block Layouts) * (Internal Perms of SM) * (Internal Perms of Classicists)")
print(f"  = {num_block_arrangements_1} * {perms_sm_unconstrained} * {perms_creg}")
print(f"  = {total_ways_scenario1}\n")

print("Scenario 2: Science/Math block (SM) is next to Cassie.")
print("Ways = (Block Layouts) * (Internal Perms of SM, constrained) * (Internal Perms of Classicists)")
print(f"  = {num_block_arrangements_2} * ({perms_sm_female_m_end} + {perms_sm_female_s_end}) * {perms_creg}")
print(f"  = {num_block_arrangements_2} * {perms_sm_constrained} * {perms_creg}")
print(f"  = {total_ways_scenario2}\n")

print("Total Ways = (Ways to Choose Rowers) * [ (Ways for Scenario 1) + (Ways for Scenario 2) ]")
print(f"Total Ways = {num_rower_choices} * [ {total_ways_scenario1} + {total_ways_scenario2} ]")
print(f"Total Ways = {num_rower_choices} * {total_ways_scenario1 + total_ways_scenario2}")
print(f"Total Ways = {total_ways}")

<<<821348311040000>>>