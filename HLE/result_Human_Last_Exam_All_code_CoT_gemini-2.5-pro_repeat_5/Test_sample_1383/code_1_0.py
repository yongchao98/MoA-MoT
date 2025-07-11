def solve_betting_problem():
    """
    This script provides a step-by-step symbolic solution to the bike race betting problem.
    It calculates the achieved doubling rate (W) based on incorrect probabilities,
    and the decrease from the optimal rate (ΔW = W* - W).
    """

    print("This script calculates the doubling rate W you will achieve and the decrease from the optimal rate, ΔW.")
    print("The final answers are provided in terms of natural logs and fractions.")
    print("-" * 70)

    # --- Part 1: Calculate the achieved doubling rate W ---
    print("Part 1: Calculating the achieved doubling rate W")
    print("You bet based on incorrect probabilities q = [1/4, 1/2, 1/8, 1/8].")
    print("The actual growth rate W depends on the true probabilities p = [1/2, 1/4, 1/8, 1/8].")
    print("The odds are o = [4, 3, 7, 7].")
    print("The formula is: W = Σ p_i * ln(q_i * o_i)\n")
    
    print("W = p_1*ln(q_1*o_1) + p_2*ln(q_2*o_2) + p_3*ln(q_3*o_3) + p_4*ln(q_4*o_4)")
    print("W = (1/2) * ln((1/4) * 4) + (1/4) * ln((1/2) * 3) + (1/8) * ln((1/8) * 7) + (1/8) * ln((1/8) * 7)")
    print("W = (1/2) * ln(1) + (1/4) * ln(3/2) + (1/8) * ln(7/8) + (1/8) * ln(7/8)")
    print("Since ln(1) = 0:")
    print("W = 0 + (1/4) * ln(3/2) + (2/8) * ln(7/8)")
    print("\nThe final expression for W is:")
    print("W = (1/4) * ln(3/2) + (1/4) * ln(7/8)")
    print("-" * 70)

    # --- Part 2: Calculate the optimal doubling rate W* ---
    print("Part 2: Calculating the optimal doubling rate W* for reference")
    print("The optimal rate is achieved by betting based on the true probabilities p.")
    print("The formula is: W* = Σ p_i * ln(p_i * o_i)\n")

    print("W* = p_1*ln(p_1*o_1) + p_2*ln(p_2*o_2) + p_3*ln(p_3*o_3) + p_4*ln(p_4*o_4)")
    print("W* = (1/2) * ln((1/2) * 4) + (1/4) * ln((1/4) * 3) + (1/8) * ln((1/8) * 7) + (1/8) * ln((1/8) * 7)")
    print("W* = (1/2) * ln(2) + (1/4) * ln(3/4) + (2/8) * ln(7/8)")
    print("\nThe final expression for W* is:")
    print("W* = (1/2) * ln(2) + (1/4) * ln(3/4) + (1/4) * ln(7/8)")
    print("-" * 70)

    # --- Part 3: Calculate the decrease in doubling rate ΔW ---
    print("Part 3: Calculating the decrease in doubling rate ΔW = W* - W\n")
    print("ΔW = W* - W")
    print("ΔW = [ (1/2)ln(2) + (1/4)ln(3/4) + (1/4)ln(7/8) ] - [ (1/4)ln(3/2) + (1/4)ln(7/8) ]")
    print("The (1/4)ln(7/8) terms cancel out:")
    print("ΔW = (1/2)ln(2) + (1/4)ln(3/4) - (1/4)ln(3/2)")
    print("Using log properties ln(a/b) = ln(a) - ln(b) and ln(x^y) = y*ln(x):")
    print("ΔW = (1/2)ln(2) + (1/4)(ln(3) - ln(4)) - (1/4)(ln(3) - ln(2))")
    print("ΔW = (1/2)ln(2) + (1/4)ln(3) - (1/4)ln(2^2) - (1/4)ln(3) + (1/4)ln(2)")
    print("ΔW = (1/2)ln(2) + (1/4)ln(3) - (2/4)ln(2) - (1/4)ln(3) + (1/4)ln(2)")
    print("ΔW = (1/2)ln(2) - (1/2)ln(2) + (1/4)ln(2)")
    print("\nThe final expression for ΔW is:")
    print("ΔW = (1/4) * ln(2)")
    print("-" * 70)

if __name__ == '__main__':
    solve_betting_problem()