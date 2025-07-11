def solve_betting_problem():
    """
    This function calculates and prints the symbolic expressions for the
    achieved growth rate (W) and the decrease in growth rate (ΔW).
    """

    # Given data
    # p: true probabilities
    # q: incorrect probabilities
    # o: odds
    p = ["1/2", "1/4", "1/8", "1/8"]
    q = ["1/4", "1/2", "1/8", "1/8"]
    o = [4, 3, 7, 7]

    print("Step 1: Calculate the achieved growth rate W.")
    print("The formula is W = sum(p_i * ln(q_i * o_i)).")
    print("W = p1*ln(q1*o1) + p2*ln(q2*o2) + p3*ln(q3*o3) + p4*ln(q4*o4)")
    print(f"W = ({p[0]})*ln(({q[0]})*{o[0]}) + ({p[1]})*ln(({q[1]})*{o[1]}) + ({p[2]})*ln(({q[2]})*{o[2]}) + ({p[3]})*ln(({q[3]})*{o[3]})")
    print("W = (1/2)*ln(1) + (1/4)*ln(3/2) + (1/8)*ln(7/8) + (1/8)*ln(7/8)")
    print("Since ln(1) = 0, the first term vanishes.")
    print("W = (1/4)*ln(3/2) + (1/4)*ln(7/8)")
    print("-" * 30)

    print("Step 2: Calculate the decrease in growth rate ΔW = W* - W.")
    print("First, the formula for the optimal rate W* is W* = sum(p_i * ln(p_i * o_i)).")
    print(f"W* = ({p[0]})*ln(({p[0]})*{o[0]}) + ({p[1]})*ln(({p[1]})*{o[1]}) + ({p[2]})*ln(({p[2]})*{o[2]}) + ({p[3]})*ln(({p[3]})*{o[3]})")
    print("W* = (1/2)*ln(2) + (1/4)*ln(3/4) + (1/4)*ln(7/8)")
    print("\nNow, we compute ΔW = W* - W:")
    print("ΔW = [(1/2)*ln(2) + (1/4)*ln(3/4) + (1/4)*ln(7/8)] - [(1/4)*ln(3/2) + (1/4)*ln(7/8)]")
    print("The term (1/4)*ln(7/8) cancels out.")
    print("ΔW = (1/2)*ln(2) + (1/4)*ln(3/4) - (1/4)*ln(3/2)")
    print("Using log properties ln(a/b) = ln(a) - ln(b):")
    print("ΔW = (1/2)*ln(2) + (1/4)*(ln(3) - ln(4)) - (1/4)*(ln(3) - ln(2))")
    print("ΔW = (1/2)*ln(2) + (1/4)*ln(3) - (1/4)*ln(4) - (1/4)*ln(3) + (1/4)*ln(2)")
    print("The term (1/4)*ln(3) cancels out. Since ln(4) = ln(2^2) = 2*ln(2):")
    print("ΔW = (1/2)*ln(2) - (1/4)*(2*ln(2)) + (1/4)*ln(2)")
    print("ΔW = (1/2)*ln(2) - (1/2)*ln(2) + (1/4)*ln(2)")
    print("ΔW = (1/4)*ln(2)")
    print("-" * 30)

solve_betting_problem()