import math

def combinations(n, k):
    """Helper function to calculate combinations C(n, k)"""
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

# Define group sizes and gender distribution
s_men, s_women = 6, 6
m_men, m_women = 2, 2

print("This program calculates the number of seating arrangements based on a set of complex constraints.")
print("The total number of ways is the sum of arrangements from two main scenarios.\n")

# --- Scenario A: Classicists are NOT next to Scientists/Mathematicians ---
# This forces the block arrangement: [SM]-E1-[C]-E2, where E1 and E2 are the two distinct ethicists.
# There are 2 such circular block arrangements ([SM]-E1-[C]-E2 and [SM]-E2-[C]-E1).
blocks_A = 2

# Calculate internal arrangements for the [SM] super-block, considering the rowing team constraint.
# The rowing team (2 male scientists, 1 male mathematician) forces S and M blocks to be adjacent,
# with the rowers at the boundary.
ways_choose_s_rowers = combinations(s_men, 2)
ways_choose_m_rower = combinations(m_men, 1)
order_sm_blocks = 2  # S-M or M-S
perms_s_block_A = math.factorial(10) * math.factorial(2) # 10 other scientists, 2 rowers arrange themselves
perms_m_block_A = math.factorial(3) # 3 other mathematicians

internal_SM = ways_choose_s_rowers * ways_choose_m_rower * order_sm_blocks * perms_s_block_A * perms_m_block_A

# Calculate internal arrangements for the [C] block (5 people).
internal_C = math.factorial(5)

# Total ways for Scenario A
ways_A = internal_SM * internal_C * blocks_A
print("--- Scenario A: Classicists fully separated from Scientists/Mathematicians ---")
print(f"Number of block arrangements: {blocks_A}")
print(f"Internal [SM] arrangements = C({s_men},2) * C({m_men},1) * {order_sm_blocks} * ({10}! * {2}!) * {3}! = {internal_SM}")
print(f"Internal [C] arrangements = {5}! = {internal_C}")
print(f"Total Ways for Scenario A = {internal_SM} * {internal_C} * {blocks_A} = {ways_A}\n")


# --- Scenario B: Cassie sits next to a female Scientist or Mathematician ---
# This allows the block arrangement: [SM]-[C]-E-E. There are 4 such circular arrangements.
blocks_B = 4

# Case B1: Cassie is next to a female Mathematician (...S-M-[C]...).
choices_B1 = combinations(s_men, 2) * combinations(m_men, 1) * combinations(m_women, 1)
perms_C_B = math.factorial(4) # Cassie is fixed at the boundary
perms_S_B1 = math.factorial(10) * math.factorial(2) # Rowers at S-M boundary
perms_M_B1 = math.factorial(2) # Rower at one end, female at other; 2 in middle
ways_B1 = choices_B1 * perms_C_B * perms_S_B1 * perms_M_B1

# Case B2: Cassie is next to a female Scientist (...M-S-[C]...).
choices_B2 = combinations(s_men, 2) * combinations(s_women, 1) * combinations(m_men, 1)
perms_S_B2 = math.factorial(9) * math.factorial(2) # Female at one end, rowers at other; 9 in middle
perms_M_B2 = math.factorial(3) # Rower at S-M boundary
ways_B2 = choices_B2 * perms_C_B * perms_S_B2 * perms_M_B2

# Total ways for Scenario B
ways_B = (ways_B1 + ways_B2) * blocks_B
print("--- Scenario B: Cassie sits next to a Scientist/Mathematician ---")
print(f"Number of block arrangements: {blocks_B}")
print(f"Ways for sub-case (Cassie-Math) = C({s_men},2)*C({m_men},1)*C({m_women},1) * {4}! * ({10}!*{2}!) * {2}! = {ways_B1}")
print(f"Ways for sub-case (Cassie-Sci) = C({s_men},2)*C({s_women},1)*C({m_men},1) * {4}! * ({9}!*{2}!) * {3}! = {ways_B2}")
print(f"Total Ways for Scenario B = ({ways_B1} + {ways_B2}) * {blocks_B} = {ways_B}\n")

# --- Total Number of Arrangements ---
total_ways = ways_A + ways_B
print("--- Final Calculation ---")
print(f"Total Ways = Ways Scenario A + Ways Scenario B")
print(f"Total Ways = {ways_A} + {ways_B} = {total_ways}")
print(f"\nThe total number of ways to arrange the table is {total_ways}.")

<<<216576 * 3628800>>>