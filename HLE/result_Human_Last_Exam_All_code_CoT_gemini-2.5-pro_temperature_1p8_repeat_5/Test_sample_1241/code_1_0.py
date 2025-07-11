import sys
# It seems there is a typo in the original equation P1', since sum of P'i should be 0. We'll use the corrected one
# P1'(t) = lambda_01*P0(t) - (lambda_10+lambda_12)*P1(t) + lambda_21*P2(t) + lambda_31*P3(t)

def solve_markov_chain():
    """
    Solves for the stationary distribution of the given Kolmogorov-Chepmen system
    and calculates P0(inf) + P1(inf).
    """
    # Given transition rates
    lambda_01 = 0.019
    lambda_10 = 0.65
    lambda_12 = 0.4
    lambda_21 = 0.392
    lambda_23 = 0.008
    lambda_31 = 0.008

    print("Step 1: Define the system parameters.")
    print(f"lambda_01 = {lambda_01}")
    print(f"lambda_10 = {lambda_10}")
    print(f"lambda_12 = {lambda_12}")
    print(f"lambda_21 = {lambda_21}")
    print(f"lambda_23 = {lambda_23}")
    print(f"lambda_31 = {lambda_31}\n")

    print("Step 2: Set up the steady-state equations (P_i' = 0).")
    print("Based on the standard form of the Kolmogorov-Chapman equations, we have:")
    print("0 = -lambda_01*P0 + lambda_10*P1")
    print("0 = - (lambda_21 + lambda_23)*P2 + lambda_12*P1")
    print("0 = -lambda_31*P3 + lambda_23*P2")
    print("And the normalization condition: P0 + P1 + P2 + P3 = 1\n")

    print("Step 3: Solve for P0, P2, and P3 in terms of P1.")
    # From the first equation: P0 = (lambda_10 / lambda_01) * P1
    p0_coeff = lambda_10 / lambda_01
    print(f"From the first equation, P0 = ({lambda_10} / {lambda_01}) * P1 = {p0_coeff:.4f} * P1")
    
    # From the third equation: P2 = (lambda_12 / (lambda_21 + lambda_23)) * P1
    p2_coeff = lambda_12 / (lambda_21 + lambda_23)
    print(f"From the third equation, P2 = ({lambda_12} / ({lambda_21} + {lambda_23})) * P1 = {p2_coeff:.4f} * P1")
    
    # From the fourth equation: P3 = (lambda_23 / lambda_31) * P2
    p3_coeff_vs_p2 = lambda_23 / lambda_31
    p3_coeff = p3_coeff_vs_p2 * p2_coeff
    print(f"From the fourth equation, P3 = ({lambda_23} / {lambda_31}) * P2 = {p3_coeff_vs_p2:.4f} * P2 = {p3_coeff:.4f} * P1\n")
    
    print("Step 4: Substitute into the normalization equation to find P1.")
    # P0 + P1 + P2 + P3 = 1  => p0_coeff*P1 + P1 + p2_coeff*P1 + p3_coeff*P1 = 1
    total_coeff = p0_coeff + 1 + p2_coeff + p3_coeff
    print(f"({p0_coeff:.4f} + 1 + {p2_coeff:.4f} + {p3_coeff:.4f}) * P1 = 1")
    print(f"{total_coeff:.4f} * P1 = 1")
    
    p1 = 1 / total_coeff
    print(f"P1 = {p1}\n")

    print("Step 5: Calculate P0 using the value of P1.")
    p0 = p0_coeff * p1
    print(f"P0 = {p0_coeff:.4f} * P1 = {p0}\n")
    
    print("Step 6: Calculate the final sum P0 + P1.")
    result = p0 + p1
    # This print statement fulfills the "output each number in the final equation" requirement.
    print(f"The value of P0(+infinity) + P1(+infinity) is:")
    print(f"{p0} + {p1} = {result}")

solve_markov_chain()