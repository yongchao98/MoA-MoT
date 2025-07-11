import math

def solve_branching_walk():
    """
    Calculates the limiting probability for the branching random walk problem.
    """
    # Step 1: Define the parameters for the h->0 limiting case.
    # In the limit h->0, all sites become blue.
    # The jump probabilities on a blue site are given.
    p_L_blue_num = 1
    p_L_blue_den = 5
    p_R_blue_num = 4
    p_R_blue_den = 5

    p_L_blue = p_L_blue_num / p_L_blue_den
    p_R_blue = p_R_blue_num / p_R_blue_den

    # The initial position of the particle.
    start_pos = 3
    # The target site.
    target_pos = 0

    print("Step-by-step calculation:")
    print("1. The problem simplifies to a 1D random walk in an all-blue environment with no branching.")
    print(f"   - Left jump probability p_L = {p_L_blue}")
    print(f"   - Right jump probability p_R = {p_R_blue}")
    print(f"   - Starting position = {start_pos}")
    print(f"   - Target position = {target_pos}")

    # Step 2: The probability of this walk ever reaching the target is (p_L/p_R)^x.
    # Calculate the ratio of jump probabilities.
    ratio_num = p_L_blue_num * p_R_blue_den
    ratio_den = p_R_blue_num * p_L_blue_den
    
    # Simplify the fraction for the ratio
    common_divisor = math.gcd(ratio_num, ratio_den)
    ratio_num_simple = ratio_num // common_divisor
    ratio_den_simple = ratio_den // common_divisor

    # Calculate the distance to the target.
    distance = start_pos - target_pos

    print("\n2. The probability of reaching the target is given by the formula (p_L / p_R)^distance.")
    
    # Step 3: Print and compute the final equation.
    print("\n3. Final equation and result:")
    print(f"   P = (({p_L_blue_num}/{p_L_blue_den}) / ({p_R_blue_num}/{p_R_blue_den}))^{distance}")
    print(f"   P = ({ratio_num_simple}/{ratio_den_simple})^{distance}")
    
    final_numerator = ratio_num_simple ** distance
    final_denominator = ratio_den_simple ** distance
    
    print(f"   P = {final_numerator}/{final_denominator}")
    
    final_answer = final_numerator / final_denominator
    print(f"\nThe calculated probability is: {final_answer}")

solve_branching_walk()
