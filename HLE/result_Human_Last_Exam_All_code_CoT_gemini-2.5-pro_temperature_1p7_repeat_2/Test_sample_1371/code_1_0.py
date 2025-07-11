import math

def combinations(n, k):
    """Calculates n C k"""
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

def permutations(n, k):
    """Calculates n P k"""
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // math.factorial(n - k)

# Step 1 & 2: Define groups and structure
# Scientists: 12 total -> 6 men, 6 women
# Mathematicians: 4 total -> 2 men, 2 women
# Classicists: 5 total -> 2 men, 3 women (Cassie is a woman)
# Ethicists: 2 total -> 1 man, 1 woman

# The structure is a large SM-block (16 people) and a CE-group (7 people)
# We fix the SM-block, making the arrangement linear for the CE-group.

# Step 3: Calculate internal arrangements of Block SM
# Choose rowers: 1 male M, 2 male S
m_men_total = 2
s_men_total = 6
ways_to_choose_rowers = combinations(m_men_total, 1) * combinations(s_men_total, 2)

# Non-rower members
m_non_rowers_m = m_men_total - 1
m_non_rowers_f = 2
s_non_rowers_m = s_men_total - 2
s_non_rowers_f = 6
s_non_rowers_total = s_non_rowers_m + s_non_rowers_f

# Internal arrangement of the 2-man scientist rower team
s_rower_perms = math.factorial(2)

# Calculate BSM arrangement coefficient for each case
# The arrangements of the non-rowers within the blocks are needed.
# For M block (3 non-rowers): arranging them gives 3! ways.
# For S block (10 non-rowers): arranging them gives 10! ways.
# But we need to specify the gender of the end person.

# BSM Case A: M-end is female, S-end is female
ways_m_f_end = m_non_rowers_f * math.factorial(m_non_rowers_m + m_non_rowers_f - 1)
ways_s_f_end = s_non_rowers_f * s_rower_perms * math.factorial(s_non_rowers_total - 1)
bsm_coef_A = ways_to_choose_rowers * ways_m_f_end * (s_non_rowers_f * s_rower_perms) # Coef for 9!

# BSM Case B: M-end is male, S-end is female
ways_m_m_end = m_non_rowers_m * math.factorial(m_non_rowers_m + m_non_rowers_f - 1)
bsm_coef_B = ways_to_choose_rowers * ways_m_m_end * (s_non_rowers_f * s_rower_perms)

# BSM Case C: M-end is female, S-end is male
ways_s_m_end = s_non_rowers_m * s_rower_perms * math.factorial(s_non_rowers_total - 1)
bsm_coef_C = ways_to_choose_rowers * ways_m_f_end * (s_non_rowers_m * s_rower_perms)

# BSM Case D: M-end is male, S-end is male
bsm_coef_D = ways_to_choose_rowers * ways_m_m_end * (s_non_rowers_m * s_rower_perms)


# Step 4: Calculate arrangements of the CE group
# People who can be at the ends: Cassie (Ca), 2 Ethicists (E). Total 3.
# People who must be in the middle: 4 "normal" Classicists (CN).
# Total CE people to arrange in the middle: 5
ways_ce_middle = math.factorial(5)

# CE Case A: S_end(f), M_end(f). Ends can be anyone from {Ca, E1, E2}.
ce_coef_A = permutations(3, 2) * ways_ce_middle

# CE Case B: S_end(f), M_end(m). S1 from {Ca, E1, E2}, S7 from {E1, E2}.
# Num ways to choose ends = 2 (for S7) * 2 (for S1) = 4
ce_coef_B = 4 * ways_ce_middle

# CE Case C: S_end(m), M_end(f). Symmetric to Case B.
ce_coef_C = 4 * ways_ce_middle

# CE Case D: S_end(m), M_end(m). Ends must be {E1, E2}.
ce_coef_D = permutations(2, 2) * ways_ce_middle

# Step 5: Combine and sum results
# BSM coefficient calculation (part that multiplies 9!)
s_block_factor = math.factorial(s_non_rowers_total - 1) # This is 9!

# Let's write the coefficients that will be multiplied by 9!
c_bsm_A = ways_to_choose_rowers * ways_m_f_end * (s_non_rowers_f * s_rower_perms)
c_bsm_B = ways_to_choose_rowers * ways_m_m_end * (s_non_rowers_f * s_rower_perms)
c_bsm_C = ways_to_choose_rowers * ways_m_f_end * (s_non_rowers_m * s_rower_perms)
c_bsm_D = ways_to_choose_rowers * ways_m_m_end * (s_non_rowers_m * s_rower_perms)

# Calculate final product terms for the sum
term_A = (c_bsm_A * ce_coef_A)
term_B = (c_bsm_B * ce_coef_B)
term_C = (c_bsm_C * ce_coef_C)
term_D = (c_bsm_D * ce_coef_D)

total_coef_sum = term_A + term_B + term_C + term_D

total_arrangements = total_coef_sum * s_block_factor

print("The calculation finds the total arrangements by summing four cases based on the gender of the people at the ends of the Scientist-Mathematician block.")
print("Total Arrangements = (Arrangements for Case FF + FM + MF + MM)")
print("Each case = (Ways to form SM block for case) * (Ways to arrange CE group for case)")
print(f"Ways to arrange the 9 non-end scientists in their block = 9! = {s_block_factor}")
print("\nThe equation is structured as: 9! * [ (BSM_coef_A * CE_coef_A) + ... ]")
print("\nCase A (ends are F, F):")
print(f"  BSM coefficient part = C(2,1)*C(6,2) * (ways for M-block with F end) * (ways for S-block with F end) = {c_bsm_A}")
print(f"  CE arrangements = (ways to pick 2 of 3 people for ends) * 5! = {ce_coef_A}")
print("\nCase B (ends are F, M):")
print(f"  BSM coefficient part = {c_bsm_B}")
print(f"  CE arrangements = (4 ways to pick ends) * 5! = {ce_coef_B}")
print("\nCase C (ends are M, F):")
print(f"  BSM coefficient part = {c_bsm_C}")
print(f"  CE arrangements = (4 ways to pick ends) * 5! = {ce_coef_C}")
print("\nCase D (ends are M, M):")
print(f"  BSM coefficient part = {c_bsm_D}")
print(f"  CE arrangements = (ways to pick 2 of 2 people for ends) * 5! = {ce_coef_D}")

print("\nCombining these parts for the final calculation:")
final_eq_str = (
    f"Total = {s_block_factor} * [({c_bsm_A} * {ce_coef_A}) + ({c_bsm_B} * {ce_coef_B}) + "
    f"({c_bsm_C} * {ce_coef_C}) + ({c_bsm_D} * {ce_coef_D})]"
)
print(final_eq_str)

final_eq_str_2 = (
    f"Total = {s_block_factor} * [{term_A} + {term_B} + "
    f"{term_C} + {term_D}]"
)
print(final_eq_str_2)
print(f"Total = {s_block_factor} * {total_coef_sum}")
print(f"Total = {total_arrangements}")
print(f"\nFinal Answer: {total_arrangements}")
>>>710612736000