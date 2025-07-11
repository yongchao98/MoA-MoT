import math

def solve_seating_arrangement():
    """
    Calculates the number of ways to arrange guests at a circular table
    based on a complex set of constraints.
    """
    # --- Step 0: Define constants and helper function ---
    f = math.factorial

    # Total members per group
    num_total_sci = 12
    num_total_math = 4
    num_ethicists = 2
    num_classicists = 5
    
    # Rowing team details (all male)
    num_rower_sci = 2
    num_rower_math = 1

    # Non-rower counts and gender breakdown
    # Scientists: 12 total (6M, 6F). 2M are rowers. -> 10 non-rowers (4M, 6F)
    # Mathematicians: 4 total (2M, 2F). 1M is a rower. -> 3 non-rowers (1M, 2F)
    num_non_rower_sci = num_total_sci - num_rower_sci
    num_non_rower_math = num_total_math - num_rower_math
    num_female_non_rower_sci = 6
    num_female_non_rower_math = 2

    print("This script calculates the total number of seating arrangements based on the problem's constraints.")
    print("The final result is the sum of three distinct cases for who sits next to the Scientist/Mathematician block.")
    print("-" * 50)

    # --- Step 1: Case 1: Ethicist - (SM_Block) - Ethicist ---
    print("Case 1: The Scientist-Mathematician block (SM_Block) is shielded by two Ethicists.")
    
    # Internal arrangement of the SM_Block: [10 non-rower Sci]-[2 Sci rowers]-[1 Math rower]-[3 non-rower Math]
    ways_other_sci = f(num_non_rower_sci)
    ways_sci_rowers = f(num_rower_sci)
    ways_other_math = f(num_non_rower_math)
    ways_sm_internal = ways_other_sci * ways_sci_rowers * ways_other_math
    
    # Ways to place the two ethicists on the flanks
    ways_place_ethicists = f(num_ethicists)
    
    # Ways to arrange the new [E-SM-E] block with the 5 classicists in a circle (6 items total)
    ways_circular_arrangement = f((1 + num_classicists) - 1)
    
    N1 = ways_sm_internal * ways_place_ethicists * ways_circular_arrangement
    print(f"Internal SM_Block ways = f({num_non_rower_sci}) * f({num_rower_sci}) * f({num_non_rower_math}) = {ways_sm_internal:,}")
    print(f"Equation for Case 1: ({ways_other_sci} * {ways_sci_rowers} * {ways_other_math}) * {ways_place_ethicists} * {ways_circular_arrangement}")
    print(f"Total for Case 1 = {N1:,}\n")

    # --- Step 2: Case 2: Cassie - (SM_Block with female Scientist at end) - Ethicist ---
    print("Case 2: Cassie shields the Scientist-end, an Ethicist shields the Mathematician-end.")
    
    # Internal ways for an SM_Block with a female scientist at the designated end
    ways_sci_part_fS_end = num_female_non_rower_sci * f(num_non_rower_sci - 1)
    ways_sm_fS_end = ways_sci_part_fS_end * ways_sci_rowers * ways_other_math
    
    # Ways to choose which of the 2 ethicists takes the available spot
    ways_choose_ethicist = num_ethicists

    # Ways to arrange the [Ca-SM-E] block, the other ethicist, and the 4 other classicists (6 items total)
    ways_circular_arrangement_2 = f((1 + 1 + (num_classicists - 1)) - 1)

    N2 = ways_sm_fS_end * ways_choose_ethicist * ways_circular_arrangement_2
    print(f"Internal SM_Block (Female Sci end) ways = ({num_female_non_rower_sci} * f({num_non_rower_sci-1})) * f({num_rower_sci}) * f({num_non_rower_math}) = {ways_sm_fS_end:,}")
    print(f"Equation for Case 2: {ways_sm_fS_end} * {ways_choose_ethicist} * {ways_circular_arrangement_2}")
    print(f"Total for Case 2 = {N2:,}\n")
    
    # --- Step 3: Case 3: Ethicist - (SM_Block with female Mathematician at end) - Cassie ---
    print("Case 3: An Ethicist shields the Scientist-end, Cassie shields the Mathematician-end.")
    
    # Internal ways for an SM_Block with a female mathematician at the designated end
    ways_math_part_fM_end = num_female_non_rower_math * f(num_non_rower_math - 1)
    ways_sm_fM_end = ways_other_sci * ways_sci_rowers * ways_math_part_fM_end

    # Ways to choose which of the 2 ethicists takes the available spot
    ways_choose_ethicist_3 = num_ethicists

    # Ways to arrange the [E-SM-Ca] block, the other ethicist, and the 4 other classicists (6 items total)
    ways_circular_arrangement_3 = f((1 + 1 + (num_classicists - 1)) - 1)
    
    N3 = ways_sm_fM_end * ways_choose_ethicist_3 * ways_circular_arrangement_3
    print(f"Internal SM_Block (Female Math end) ways = f({num_non_rower_sci}) * f({num_rower_sci}) * ({num_female_non_rower_math} * f({num_non_rower_math - 1})) = {ways_sm_fM_end:,}")
    print(f"Equation for Case 3: {ways_sm_fM_end} * {ways_choose_ethicist_3} * {ways_circular_arrangement_3}")
    print(f"Total for Case 3 = {N3:,}\n")

    # --- Step 4: Final Calculation ---
    print("-" * 50)
    total_arrangements = N1 + N2 + N3
    print("The total number of arrangements is the sum of all cases.")
    print(f"Final Equation: Case 1 + Case 2 + Case 3")
    print(f"Total = {N1:,} + {N2:,} + {N3:,}")
    print(f"Total Arrangements = {total_arrangements:,}")
    print("-" * 50)
    print("Final equation expressed with factorials:")
    print(f"Total = [ (f(10)*f(2)*f(3)) * f(2) * f(5) ] + [ (6*f(9)*f(2)*f(3)) * 2 * f(5) ] + [ (f(10)*f(2)*(2*f(2))) * 2 * f(5) ]")
    print(f"Total = (2880 * f(10)) + (1728 * f(10)) + (1920 * f(10)) = 6528 * f(10)")
    final_f10 = f(10)
    final_result = 6528 * final_f10
    print(f"Total = 6528 * {final_f10} = {final_result:,}")
    

if __name__ == '__main__':
    solve_seating_arrangement()