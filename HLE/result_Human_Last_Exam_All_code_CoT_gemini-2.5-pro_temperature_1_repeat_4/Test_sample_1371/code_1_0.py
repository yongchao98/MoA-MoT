import math

# This script calculates the number of seating arrangements based on a set of complex constraints.

# Step 1: Define key numbers and pre-calculate factorials.
num_s_total = 12
num_m_total = 4
num_e_total = 2
num_c_total = 5

num_s_male = 6
num_s_female = 6
num_m_male = 2
num_m_female = 2
# For the 3 non-rower mathematicians, there are 1 male and 2 females.
# For the 10 non-rower scientists, there are 4 males and 6 females.

f10 = math.factorial(10)
f9 = math.factorial(9)
f5 = math.factorial(5)
f3 = math.factorial(3)
f2 = math.factorial(2)

# Step 2: Calculate the number of internal arrangements for the combined SM_BLOCK.
# The block is formed because all Scientists must sit together, all Mathematicians must sit together,
# and a rowing team (2 male S, 1 male M) must sit together, forcing the groups to be adjacent.

# First, find the number of ways to select and arrange the rowers at the boundary.
# Choose 2 male scientist rowers from 6.
ways_to_choose_s_rowers = math.comb(num_s_male, 2)
# Choose 1 male mathematician rower from 2.
ways_to_choose_m_rowers = math.comb(num_m_male, 1)
# The 2 scientist rowers can be arranged in 2! ways.
rower_boundary_arrangements = ways_to_choose_s_rowers * ways_to_choose_m_rowers * f2

# Now, calculate the different types of internal arrangements based on who is at the ends.
# W_total: Total internal arrangements.
# This involves arranging the 10 non-rower scientists (10!) and 3 non-rower mathematicians (3!).
W_total = rower_boundary_arrangements * f10 * f3

# W_S_female: Arrangements where the Scientist at the free end is female.
# There are 6 female scientists available for the end position.
# Arrange the remaining 9 scientists (9!) and 3 mathematicians (3!).
arrangements_s_female_end = num_s_female * f9
W_S_female = rower_boundary_arrangements * arrangements_s_female_end * f3

# W_M_female: Arrangements where the Mathematician at the free end is female.
# There are 2 female mathematicians available for the end position.
# Arrange the 10 scientists (10!) and the remaining 2 mathematicians (2!).
arrangements_m_female_end = num_m_female * f2
W_M_female = rower_boundary_arrangements * f10 * arrangements_m_female_end

# Step 3: Calculate external arrangements around the fixed SM_BLOCK.
# The people not neighboring the block (always 5 of them) can be arranged in 5! ways.
external_person_arrangements = f5

# Case 1: The two Ethicists are the neighbors of the SM_BLOCK.
# The 2 Ethicists can be arranged in 2! ways.
ways_case1 = f2 * external_person_arrangements * W_total

# Case 2: One Ethicist and Cassie are the neighbors.
# Cassie's placement depends on the gender of the end person of the block.
# Choose which of the 2 Ethicists is the neighbor.
# Subcase 2a: Cassie is at the Scientist end (requires female Scientist).
ways_case2a = 2 * external_person_arrangements * W_S_female
# Subcase 2b: Cassie is at the Mathematician end (requires female Mathematician).
ways_case2b = 2 * external_person_arrangements * W_M_female

# Step 4: Sum all valid cases to get the total number of arrangements.
total_arrangements = ways_case1 + ways_case2a + ways_case2b

# To display the final equation clearly, we can factor out common terms.
# The common factor for placing neighbors and remaining people is 2 * 5! = 240 (or 2! * 5! = 240).
# The final equation simplifies to: 240 * (W_total + W_S_female + W_M_female)
# We can further simplify the W terms in relation to 10!
# W_total = 360 * 10!
# W_S_female = 216 * 10!
# W_M_female = 240 * 10!
# So, the sum of weighted internal arrangements is (360 + 216 + 240) = 816.

final_coeff_1 = 240
final_coeff_2 = int((W_total + W_S_female + W_M_female) / f10)

print("The total number of arrangements is found by summing the possibilities for who can sit next to the main block of scientists and mathematicians.")
print("\nFinal Equation:")
print(f"Total Arrangements = (Neighbor Arrangements) * (Sum of Weighted Internal Block Arrangements)")
print(f"Total Arrangements = {final_coeff_1} * ({int(W_total/f10)}*10! + {int(W_S_female/f10)}*10! + {int(W_M_female/f10)}*10!)")
print(f"Total Arrangements = {final_coeff_1} * {final_coeff_2} * 10!")
print(f"Total Arrangements = {final_coeff_1} * {final_coeff_2} * {f10}")
print(f"Total Arrangements = {total_arrangements}")
<<<710893568000>>>