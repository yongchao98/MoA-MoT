import math

def solve_seating_arrangement():
    """
    Calculates the number of ways to arrange the group at a circular table
    based on the given constraints.
    """

    # Step 1: Define the number of people in each category.
    # Scientists: 6 Male (2 rowers), 6 Female
    # Mathematicians: 2 Male (1 rower), 2 Female
    # Ethicists: 1 Male, 1 Female
    # Classicists: 2 Male, 3 Female (including Cassie)
    
    num_f_sci_non_rower = 6
    num_m_sci_non_rower = 4
    num_f_math_non_rower = 2
    num_m_math_non_rower = 1

    # Step 2: Calculate the number of internal arrangements for the SM block for each
    # of the four end-gender combinations. The block is (10 non-rower Sci) - (3 rowers) - (3 non-rower Math).
    # The internal permutations of the non-end people are 9! for scientists, 2! for mathematicians,
    # and 3! for the rowing team.

    # Base permutations for the groups within the block
    perm_9_sci = math.factorial(9)
    perm_2_math = math.factorial(2)
    perm_3_rowers = math.factorial(3)

    # Case 1: Female Scientist end, Female Mathematician end
    ways_sm_ff = (num_f_sci_non_rower * perm_9_sci) * (num_f_math_non_rower * perm_2_math) * perm_3_rowers
    
    # Case 2: Female Scientist end, Male Mathematician end
    ways_sm_fm = (num_f_sci_non_rower * perm_9_sci) * (num_m_math_non_rower * perm_2_math) * perm_3_rowers

    # Case 3: Male Scientist end, Female Mathematician end
    ways_sm_mf = (num_m_sci_non_rower * perm_9_sci) * (num_f_math_non_rower * perm_2_math) * perm_3_rowers

    # Case 4: Male Scientist end, Male Mathematician end
    ways_sm_mm = (num_m_sci_non_rower * perm_9_sci) * (num_m_math_non_rower * perm_2_math) * perm_3_rowers

    # Step 3: Calculate the number of ways to arrange the people outside the SM block for each case.
    # The people available to sit next to the SM block are Cassie (C) and 2 Ethicists (E1, E2).
    # The remaining 5 people (4 regular classicists + whoever is left) are arranged in 5! ways in the gap.
    perm_5_gap = math.factorial(5)

    # Case 1 (F-F ends): C, E1, E2 can sit on either side. P(3,2) ways.
    outer_1 = math.perm(3, 2) * perm_5_gap

    # Case 2 (F-M ends): C,E1,E2 can sit by F end. E1,E2 can sit by M end. (1*2 + 2*1) = 4 ways.
    outer_2 = 4 * perm_5_gap

    # Case 3 (M-F ends): Symmetric to Case 2. 4 ways.
    outer_3 = 4 * perm_5_gap

    # Case 4 (M-M ends): E1, E2 can sit on either side. P(2,2) ways.
    outer_4 = math.perm(2, 2) * perm_5_gap

    # Step 4: Calculate total for each case and sum them up.
    total_1 = ways_sm_ff * outer_1
    total_2 = ways_sm_fm * outer_2
    total_3 = ways_sm_mf * outer_3
    total_4 = ways_sm_mm * outer_4

    grand_total = total_1 + total_2 + total_3 + total_4

    # Step 5: Print the breakdown of the calculation.
    # The calculation can be simplified to: (5! * 9! * 3!) * (sum of coefficients)
    # The coefficients are derived from (num_at_end * perm_of_others * outer_perm_factor)
    # Let's calculate the central coefficient:
    # Coeff = (6*2*6) + (6*1*4) + (4*2*4) + (4*1*2) = 72 + 24 + 32 + 8 = 136. Wait, my math is off.
    # Let's re-calculate the coefficient from the formula:
    # Total = 5! * [ (ways_sm_ff*P(3,2)) + (ways_sm_fm*4) + (ways_sm_mf*4) + (ways_sm_mm*P(2,2)) ] / 5!
    # Total = 5! * 9! * 2! * 3! * [ (6*2*6) + (6*1*4) + (4*2*4) + (4*1*2) ]
    # Total = 5! * 9! * 2 * 6 * [ 72 + 24 + 32 + 8 ] = 5! * 9! * 12 * 136. Still not matching my manual calc.
    
    # Let's use the manual calculation which was more robust.
    # Total = (9! * 3! * 5!) * [ (24 * 6) + (12 * 4) + (16 * 4) + (8 * 2) ]
    # The coefficients 24, 12, 16, 8 come from (ways_sm_case / (9! * 3!))
    # Let's verify: ways_sm_ff / (9!*3!) = (6*9!*2*2!*3!)/(9!*3!) = 6*2*2 = 24. Correct.
    
    c1 = (num_f_sci_non_rower * num_f_math_non_rower * perm_2_math) * math.perm(3,2) # 6 * 2 * 2 * 6 = 144
    c2 = (num_f_sci_non_rower * num_m_math_non_rower * perm_2_math) * 4 # 6 * 1 * 2 * 4 = 48
    c3 = (num_m_sci_non_rower * num_f_math_non_rower * perm_2_math) * 4 # 4 * 2 * 2 * 4 = 64
    c4 = (num_m_sci_non_rower * num_m_math_non_rower * perm_2_math) * math.perm(2,2) # 4 * 1 * 2 * 2 = 16
    
    total_coefficient = c1 + c2 + c3 + c4 # 144 + 48 + 64 + 16 = 272

    f5 = math.factorial(5)
    f9 = math.factorial(9)
    f3 = math.factorial(3)
    
    final_result = f5 * f9 * f3 * total_coefficient

    print("The total number of arrangements is calculated by summing four cases based on the gender of the people at the ends of the Scientist-Mathematician (SM) block.")
    print("The formula simplifies to: (Arrangements of people in the 'safe' gap) * (Internal arrangements of non-end members of SM block) * (A coefficient representing the sum of all case permutations)")
    print("\nThis is: 5! * 9! * 3! * Coefficient\n")
    print(f"The coefficient is the sum of possibilities for each case: {int(c1)} + {int(c2)} + {int(c3)} + {int(c4)} = {int(total_coefficient)}\n")
    print("Calculating the final number:")
    print(f"5! = {f5}")
    print(f"9! = {f9}")
    print(f"3! = {f3}")
    print(f"Coefficient = {int(total_coefficient)}\n")
    print(f"Total Ways = {f5} * {f9} * {f3} * {int(total_coefficient)}")
    print(f"Total Ways = {final_result}")
    
    return final_result

result = solve_seating_arrangement()
print(f"\n<<<{result}>>>")