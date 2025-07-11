import math

def nCr(n, r):
    """Calculates the number of combinations (n choose r)."""
    return math.comb(n, r)

def factorial(n):
    """Calculates the factorial of a number."""
    return math.factorial(n)

# Step 1: Define group compositions
num_s_male, num_s_female = 6, 6
num_m_male, num_m_female = 2, 2
num_rower_s, num_rower_m = 2, 1
num_ethicists = 2
num_classicists = 5

# The rowing team members are chosen from the male academics
other_s_male = num_s_male - num_rower_s
other_s_female = num_s_female
other_m_male = num_m_male - num_rower_m
other_m_female = num_m_female

# Number of ways to choose the specific individuals for the rowing team
ways_to_choose_rowers = nCr(num_s_male, num_rower_s) * nCr(num_m_male, num_rower_m)
print(f"Number of ways to choose the rowing team members: C({num_s_male},{num_rower_s}) * C({num_m_male},{num_rower_m}) = {ways_to_choose_rowers}")

# Step 2: Calculate ways to arrange the 7 "outsiders" based on A-Block ends
# The outsiders are 5 classicists (4 normal, 1 Cassie) and 2 ethicists.
# We are arranging them in a line of 7 seats between the two ends of the A-Block.
# Let the ends of the line be Left_End and Right_End.
# Normal classicists cannot sit at Left_End or Right_End.
# Cassie can only sit there if the adjacent academic is female.
# Ethicists can sit anywhere.

# Case O1: Both ends are occupied by Ethicists.
# Place 2 ethicists at the 2 ends (2! ways). Arrange 5 classicists in the middle (5! ways).
ways_O_E_E = factorial(2) * factorial(5)
print(f"Ways to seat outsiders with Ethicists at both ends: 2! * 5! = {ways_O_E_E}")

# Case O2: Cassie at one end, Ethicist at the other.
# Place Cassie at an end (1 way), place one of 2 Ethicists at the other end (2 ways).
# Arrange remaining 5 people (4 classicists, 1 ethicist) in middle (5! ways).
ways_O_C_E = 1 * 2 * factorial(5)
print(f"Ways to seat outsiders with Cassie at one end and an Ethicist at the other: 2 * 5! = {ways_O_C_E}")

# Total ways for outsiders depend on the A-Block end genders:
# Ends are Male-Male: only E-E arrangement is possible. Total_O_MM = ways_O_E_E
# Ends are Female-Male: Can have E-E or C-E (Cassie next to Female). Total_O_FM = ways_O_E_E + ways_O_C_E
# Ends are Male-Female: Can have E-E or E-C. Total_O_MF = ways_O_E_E + ways_O_C_E
# Ends are Female-Female: Can have E-E, C-E, or E-C. Total_O_FF = ways_O_E_E + ways_O_C_E + ways_O_C_E

# Step 3: Calculate ways for A-Block arrangements for one orientation [S...M]
# The A-Block is [other_scientists (10)]-[rower_scientists (2)]-[rower_mathematician (1)]-[other_mathematicians (3)]
# The ends are one of the 10 other scientists and one of the 3 other mathematicians.

# Ways to arrange the 2 scientist rowers within their slot
ways_permute_s_rowers = factorial(2)

# Calculate internal permutations for each end-gender case.
# A(LS, RM): Left Scientist is Male, Right Mathematician is Male
ways_A_LM_RM = ways_permute_s_rowers * (other_s_male * factorial(9)) * (other_m_male * factorial(2))
# A(LS, RF): Left Scientist is Male, Right Mathematician is Female
ways_A_LM_RF = ways_permute_s_rowers * (other_s_male * factorial(9)) * (other_m_female * factorial(2))
# A(LF, RM): Left Scientist is Female, Right Mathematician is Male
ways_A_LF_RM = ways_permute_s_rowers * (other_s_female * factorial(9)) * (other_m_male * factorial(2))
# A(LF, RF): Left Scientist is Female, Right Mathematician is Female
ways_A_LF_RF = ways_permute_s_rowers * (other_s_female * factorial(9)) * (other_m_female * factorial(2))

# Combine with outsider arrangements
# For SM orientation (L=Scientist, R=Mathematician)
total_SM = 0
total_SM += ways_A_LM_RM * ways_O_E_E
total_SM += ways_A_LM_RF * (ways_O_E_E + ways_O_C_E) # Right end is female, Cassie can sit there
total_SM += ways_A_LF_RM * (ways_O_E_E + ways_O_C_E) # Left end is female, Cassie can sit there
total_SM += ways_A_LF_RF * (ways_O_E_E + ways_O_C_E + ways_O_C_E)

# Step 4: Repeat for the other orientation [M...S]
# A(LM, RS): Left Mathematician is Male, Right Scientist is Male
ways_A_LM_RS = ways_permute_s_rowers * (other_m_male * factorial(2)) * (other_s_male * factorial(9))
# A(LM, RF): Left Mathematician is Male, Right Scientist is Female
ways_A_LM_RF_2 = ways_permute_s_rowers * (other_m_male * factorial(2)) * (other_s_female * factorial(9))
# A(LF, RS): Left Mathematician is Female, Right Scientist is Male
ways_A_LF_RS = ways_permute_s_rowers * (other_m_female * factorial(2)) * (other_s_male * factorial(9))
# A(LF, RF): Left Mathematician is Female, Right Scientist is Female
ways_A_LF_RF_2 = ways_permute_s_rowers * (other_m_female * factorial(2)) * (other_s_female * factorial(9))

# For MS orientation (L=Mathematician, R=Scientist)
total_MS = 0
total_MS += ways_A_LM_RS * ways_O_E_E
total_MS += ways_A_LM_RF_2 * (ways_O_E_E + ways_O_C_E)
total_MS += ways_A_LF_RS * (ways_O_E_E + ways_O_C_E)
total_MS += ways_A_LF_RF_2 * (ways_O_E_E + ways_O_C_E + ways_O_C_E)

# The total number of arrangements is the sum for both orientations, multiplied by
# the number of ways to choose the rowers in the first place.
final_answer = ways_to_choose_rowers * (total_SM + total_MS)

print("\n--- Detailed Calculation ---")
print(f"Total internal arrangements for the S...M block, weighted by outsider perms = {ways_to_choose_rowers} * ({total_SM})")
print(f"Total internal arrangements for the M...S block, weighted by outsider perms = {ways_to_choose_rowers} * ({total_MS})")

# Let's break down one term for clarity
# Term 1: S...M block, ends are Female Scientist and Female Mathematician
# Ways to arrange A-block: ways_A_LF_RF = 2! * (6 * 9!) * (2 * 2!) = 17,418,240
# Ways to arrange outsiders: ways_O_E_E + ways_O_C_E + ways_O_C_E = 240 + 240 + 240 = 720
# Contribution = ways_to_choose_rowers * ways_A_LF_RF * 720
# = 30 * 17418240 * 720
# Our total_SM already aggregates this.

print("\n--- Final Answer ---")
print("Final Equation: (ways_to_choose_rowers) * (total_ways_for_SM_orientation + total_ways_for_MS_orientation)")
print(f"Final Calculation: {ways_to_choose_rowers} * ({total_SM} + {total_MS})")
print(f"The total number of ways to arrange the table is: {final_answer}")