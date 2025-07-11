import math

def solve_seating_arrangement():
    """
    Calculates the number of ways to arrange the members at a circular table
    based on a set of complex constraints.
    """
    
    # Plan Explanation:
    # 1. The problem is broken down into two main blocks: a STEM block (Scientists + Mathematicians)
    #    and a Humanities block (Classicists + Ethicists).
    # 2. We calculate arrangements within the STEM block, forced into a specific structure by the
    #    rowing team constraint.
    # 3. We then analyze the two people at the ends of this STEM block. The seating arrangement
    #    of the Humanities block depends on the gender of these two people.
    # 4. We use a case-based approach based on the gender of the STEM block's end-members.
    # 5. The total arrangements are the sum of possibilities from all cases.

    # Group Composition:
    # Scientists (S): 12 total (6M, 6W). 2M are rowers. -> Other S = 10 people (4M, 6W).
    # Mathematicians (M): 4 total (2M, 2W). 1M is a rower. -> Other M = 3 people (1M, 2W).
    # Humanities (H): 7 people -> 2 Ethicists, 4 non-Cassie Classicists, 1 Cassie (Classicist).

    # --- Step 1: Calculate permutations for components of the STEM block ---
    
    # The STEM block has a structure: (10 'other' scientists)-(2 scientist rowers)-(1 math rower)-(3 'other' mathematicians).
    # The ends of this entire block are one of the 'other' scientists and one of the 'other' mathematicians.

    # Number of ways to arrange the scientist portion with a WOMAN at the outer end.
    # (Choose 1 of 6 women for the end) * (arrange other 9 scientists) * (arrange the 2 rowers)
    n_s_end_W = 6 * math.factorial(9) * math.factorial(2)
    
    # Number of ways to arrange the scientist portion with a MAN at the outer end.
    # (Choose 1 of 4 men for the end) * (arrange other 9 scientists) * (arrange the 2 rowers)
    n_s_end_M = 4 * math.factorial(9) * math.factorial(2)

    # Number of ways to arrange the mathematician portion with a WOMAN at the outer end.
    # (Choose 1 of 2 women for the end) * (arrange other 2 mathematicians)
    n_m_end_W = 2 * math.factorial(2)

    # Number of ways to arrange the mathematician portion with a MAN at the outer end.
    # (Choose 1 of 1 man for the end) * (arrange other 2 mathematicians)
    n_m_end_M = 1 * math.factorial(2)
    
    # --- Step 2: Calculate permutations for the Humanities block based on cases ---

    # There are 7 people in the Humanities block. The 5 in the middle can be arranged in 5! ways.
    # The arrangement possibilities of the two end-members change based on the case.
    perms_middle_H = math.factorial(5)

    # Case A: STEM block ends are (Woman, Woman).
    # Humanities ends can be (Cassie, Ethicist), (Ethicist, Cassie), or (Ethicist, Ethicist). This gives 6 permutations for the end-pair.
    ways_H_A = 6 * perms_middle_H
    
    # Case B: STEM block ends are (Woman, Man).
    # One end (next to Woman) can be Cassie or Ethicist. The other end (next to Man) can only be an Ethicist. This gives 4 valid permutations for the end-pair.
    ways_H_B = 4 * perms_middle_H

    # Case C: STEM block ends are (Man, Woman). Symmetric to Case B.
    ways_H_C = 4 * perms_middle_H
    
    # Case D: STEM block ends are (Man, Man).
    # Both ends must be Ethicists. This gives 2 permutations for the end-pair.
    ways_H_D = 2 * perms_middle_H

    # --- Step 3: Combine cases for the final answer ---

    # Total for Case A = (arrangements for S-end=W) * (arrangements for M-end=W) * (arrangements for H-block in this case)
    total_A = n_s_end_W * n_m_end_W * ways_H_A

    # Total for Case B = (arrangements for S-end=W) * (arrangements for M-end=M) * (arrangements for H-block in this case)
    total_B = n_s_end_W * n_m_end_M * ways_H_B

    # Total for Case C = (arrangements for S-end=M) * (arrangements for M-end=W) * (arrangements for H-block in this case)
    total_C = n_s_end_M * n_m_end_W * ways_H_C
    
    # Total for Case D = (arrangements for S-end=M) * (arrangements for M-end=M) * (arrangements for H-block in this case)
    total_D = n_s_end_M * n_m_end_M * ways_H_D

    total_arrangements = total_A + total_B + total_C + total_D

    # For clarity, let's show the simplified equation Total = C * 9! * 5!
    # The coefficient C is the sum of the products of the non-factorial parts.
    c_s_W = 6 * 2 # 12
    c_s_M = 4 * 2 # 8
    c_m_W = 2 * 2 # 4
    c_m_M = 1 * 2 # 2
    c_h_A = 6
    c_h_B = 4
    c_h_C = 4
    c_h_D = 2
    
    term_A = c_s_W * c_m_W * c_h_A # 12 * 4 * 6 = 288
    term_B = c_s_W * c_m_M * c_h_B # 12 * 2 * 4 = 96
    term_C = c_s_M * c_m_W * c_h_C # 8 * 4 * 4 = 128
    term_D = c_s_M * c_m_M * c_h_D # 8 * 2 * 2 = 32
    
    coefficient = term_A + term_B + term_C + term_D

    print("The final calculation can be simplified to the form: C * 9! * 5!")
    print(f"The coefficient C is the sum of the possibilities from four distinct cases:")
    print(f"C = ({c_s_W}* {c_m_W}* {c_h_A}) + ({c_s_W}* {c_m_M}* {c_h_B}) + ({c_s_M}* {c_m_W}* {c_h_C}) + ({c_s_M}* {c_m_M}* {c_h_D})")
    print(f"C = {term_A} + {term_B} + {term_C} + {term_D} = {coefficient}")
    print("\nThus, the total number of arrangements is:")
    print(f"Total = {coefficient} * 9! * 5!")
    f9 = math.factorial(9)
    f5 = math.factorial(5)
    print(f"Total = {coefficient} * {f9} * {f5}")
    print(f"\nFinal Answer: {total_arrangements}")


solve_seating_arrangement()
<<<22186368000>>>