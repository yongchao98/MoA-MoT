import math

def solve_seating_arrangement():
    """
    Calculates the number of ways to arrange the attendees at a circular table
    based on a complex set of constraints.
    """

    # --- Step 1: Define group compositions ---
    # From the problem statement, we deduce the following:
    # Scientists: 10 non-rowers (4m, 6f), 2 male rowers
    # Mathematicians: 3 non-rowers (1m, 2f), 1 male rower
    # Ethicists: 1m, 1f
    # Classicists: 4 non-Cassie (2m, 2f), Cassie (f)

    # --- Step 2: Calculate internal arrangements for the Scientist-Mathematician (SM) block ---
    # We calculate the number of ways to arrange the SM block for the "S-M" order,
    # broken down by the gender of the people at the ends of the block.
    # A symmetric calculation exists for the "M-S" order, so we will multiply by 2 later.

    # Arrangements of the Scientist part: [end-person, 9 others, 2 rowers]
    # - 6 choices for a female end-person, 4 for a male.
    # - 9! ways to arrange the 9 other non-rower scientists.
    # - 2! ways to arrange the 2 male rowers.
    ways_s_part_f_end = 6 * math.factorial(9) * math.factorial(2)
    ways_s_part_m_end = 4 * math.factorial(9) * math.factorial(2)

    # Arrangements of the Mathematician part: [1 rower, 2 others, end-person]
    # - 2 choices for a female end-person, 1 for a male.
    # - 2! ways to arrange the 2 other non-rower mathematicians.
    ways_m_part_f_end = 2 * math.factorial(2)
    ways_m_part_m_end = 1 * math.factorial(2)

    # Combine to get arrangements for each end-gender case
    W_FF = ways_s_part_f_end * ways_m_part_f_end
    W_FM = ways_s_part_f_end * ways_m_part_m_end # S-end is F, M-end is M
    W_MF = ways_s_part_m_end * ways_m_part_f_end # S-end is M, M-end is F
    W_MM = ways_s_part_m_end * ways_m_part_m_end

    # --- Step 3: Calculate external arrangements for the 7 other people ---
    # The 7 people are 2 Ethicists and 5 Classicists (including Cassie).
    # The 4 non-Cassie classicists cannot sit next to the SM block.
    fact_5 = math.factorial(5)

    # Case 1: SM block ends are Female-Female.
    # The 2 adjacent seats can be filled by 2 Ethicists or Cassie (3 people).
    # Ways to choose and place 2 out of 3 people: P(3,2) = 6.
    # The remaining 5 people can be arranged in 5! ways.
    Ext_FF = 3 * 2 * fact_5

    # Case 2: SM block ends are Female-Male.
    # Seat next to F-end: Ethicists or Cassie (3 choices).
    # Seat next to M-end: Ethicists only (2 choices).
    # Total ways to fill adjacent seats: (choose Ethicist for M-end)* (choose from rest for F-end) = 2 * 2 = 4
    Ext_FM = 4 * fact_5
    
    # Case 3: SM block ends are Male-Female. Symmetric to the F-M case.
    Ext_MF = 4 * fact_5

    # Case 4: SM block ends are Male-Male.
    # The 2 adjacent seats must be filled by the 2 Ethicists. 2! ways.
    # The remaining 5 people (all classicists) are arranged in 5! ways.
    Ext_MM = math.factorial(2) * fact_5

    # --- Step 4: Combine internal and external arrangements ---
    total_sm_order_ways = (W_FF * Ext_FF) + (W_FM * Ext_FM) + (W_MF * Ext_MF) + (W_MM * Ext_MM)

    # --- Step 5: Final calculation ---
    # Multiply by 2 for the two possible SM block orders (S-M and M-S)
    final_answer = 2 * total_sm_order_ways

    # --- Step 6: Print the detailed breakdown of the equation ---
    print("The final calculation is structured by considering the gender of the people at the ends of the combined Scientist-Mathematician (SM) block.")
    print("Total ways = 2 * ( Sum of [Ways for SM internal] * [Ways for external] for each case )\n")
    
    print("Let's break down the final equation:")
    print("Total = 2 * ( W(F,F)*Ext(F,F) + W(F,M)*Ext(F,M) + W(M,F)*Ext(M,F) + W(M,M)*Ext(M,M) )")
    
    fact_9_val = math.factorial(9)
    # The numbers here (48, 24, etc.) are derived from the component calculations above.
    # e.g., W(F,F) = (6 * 9! * 2!) * (2 * 2!) = 48 * 9!
    print(f"Total = 2 * ( ({int(W_FF/fact_9_val)} * {fact_9_val})*{Ext_FF} + ({int(W_FM/fact_9_val)} * {fact_9_val})*{Ext_FM} + ({int(W_MF/fact_9_val)} * {fact_9_val})*{Ext_MF} + ({int(W_MM/fact_9_val)} * {fact_9_val})*{Ext_MM} )")
    
    sum_of_prods = (W_FF/fact_9_val * Ext_FF) + (W_FM/fact_9_val * Ext_FM) + (W_MF/fact_9_val * Ext_MF) + (W_MM/fact_9_val * Ext_MM)
    print(f"Total = 2 * {fact_9_val} * ( {int(W_FF/fact_9_val)}*{Ext_FF} + {int(W_FM/fact_9_val)}*{Ext_FM} + {int(W_MF/fact_9_val)}*{Ext_MF} + {int(W_MM/fact_9_val)}*{Ext_MM} )")
    
    prod1 = int(W_FF/fact_9_val * Ext_FF)
    prod2 = int(W_FM/fact_9_val * Ext_FM)
    prod3 = int(W_MF/fact_9_val * Ext_MF)
    prod4 = int(W_MM/fact_9_val * Ext_MM)
    print(f"Total = 2 * {fact_9_val} * ( {prod1} + {prod2} + {prod3} + {prod4} )")

    inner_sum = int(sum_of_prods)
    print(f"Total = 2 * {fact_9_val} * {inner_sum}")
    
    final_calc_str = f"Total = {2 * fact_9_val * inner_sum}"
    print(final_calc_str)
    
    return final_answer

final_result = solve_seating_arrangement()
print(f"\nFinal Answer: {final_result}")
<<<47318092800>>>