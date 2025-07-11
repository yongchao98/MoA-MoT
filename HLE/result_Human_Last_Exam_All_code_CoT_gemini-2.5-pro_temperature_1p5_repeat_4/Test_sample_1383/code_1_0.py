def solve_betting_problem():
    """
    This function calculates and prints the derivation for the achieved growth rate (W)
    and the decrease in growth rate (ΔW) for a betting scenario.
    """

    # p: True probabilities [p1, p2, p3, p4]
    # q: Incorrectly believed probabilities [q1, q2, q3, q4]
    # o: Bookmaker odds (e.g., 4-for-1)
    # d: Decimal odds (multiplier) [d1, d2, d3, d4]
    
    p = ["1/2", "1/4", "1/8", "1/8"]
    q = ["1/4", "1/2", "1/8", "1/8"]
    d = [4, 3, 7, 7]

    print("--- Calculating the Optimal Growth Rate (W*) ---")
    print("The optimal growth rate W* is achieved by betting fractions equal to the true probabilities p_i.")
    print("W* = p1*log(p1*d1) + p2*log(p2*d2) + p3*log(p3*d3) + p4*log(p4*d4)")
    print(f"W* = ({p[0]})*log(({p[0]})*{d[0]}) + ({p[1]})*log(({p[1]})*{d[1]}) + ({p[2]})*log(({p[2]})*{d[2]}) + ({p[3]})*log(({p[3]})*{d[3]})")
    print("W* = (1/2)*log(2) + (1/4)*log(3/4) + (1/8)*log(7/8) + (1/8)*log(7/8)")
    print("W* = (1/2)*log(2) + (1/4)*log(3/4) + (1/4)*log(7/8)\n")

    print("--- Calculating the Achieved Growth Rate (W) with Incorrect Probabilities ---")
    print("The achieved rate W results from using betting fractions equal to the incorrect probabilities q_i.")
    print("The expectation is still taken over the true probabilities p_i.")
    print("W = p1*log(q1*d1) + p2*log(q2*d2) + p3*log(q3*d3) + p4*log(q4*d4)")
    print(f"W = ({p[0]})*log(({q[0]})*{d[0]}) + ({p[1]})*log(({q[1]})*{d[1]}) + ({p[2]})*log(({q[2]})*{d[2]}) + ({p[3]})*log(({q[3]})*{d[3]})")
    print("W = (1/2)*log(1) + (1/4)*log(3/2) + (1/8)*log(7/8) + (1/8)*log(7/8)")
    print("W = (1/4)*log(3/2) + (1/4)*log(7/8)\n")

    print("--- Calculating the Decrease in Growth Rate (ΔW) ---")
    print("ΔW is the difference between the optimal rate W* and the achieved rate W.")
    print("ΔW = W* - W")
    print("This is also the Kullback-Leibler divergence D_KL(p || q) = sum(p_i * log(p_i / q_i)).")
    print("ΔW = p1*log(p1/q1) + p2*log(p2/q2) + p3*log(p3/q3) + p4*log(p4/q4)")
    print(f"ΔW = ({p[0]})*log(({p[0]})/({q[0]})) + ({p[1]})*log(({p[1]})/({q[1]})) + ({p[2]})*log(({p[2]})/({q[2]})) + ({p[3]})*log(({p[3]})/({q[3]}))")
    print("ΔW = (1/2)*log( (1/2)/(1/4) ) + (1/4)*log( (1/4)/(1/2) ) + (1/8)*log( (1/8)/(1/8) ) + (1/8)*log( (1/8)/(1/8) )")
    print("ΔW = (1/2)*log(2) + (1/4)*log(1/2) + (1/8)*log(1) + (1/8)*log(1)")
    print("Since log(1/2) = -log(2) and log(1) = 0:")
    print("ΔW = (1/2)*log(2) - (1/4)*log(2)")
    print("ΔW = (1/4)*log(2)\n")

    print("The achieved doubling rate is W = (1/4)*log(3/2) + (1/4)*log(7/8).")
    print("The decrease in the doubling rate is ΔW = (1/4)*log(2).")


solve_betting_problem()
<<<W = (1/4)*log(3/2) + (1/4)*log(7/8), ΔW = (1/4)*log(2)>>>