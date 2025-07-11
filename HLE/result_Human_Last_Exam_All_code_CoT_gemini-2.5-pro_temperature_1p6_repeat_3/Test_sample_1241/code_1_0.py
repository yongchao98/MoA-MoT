import numpy as np

def solve_markov_chain():
    """
    Solves the steady-state probabilities for the given Kolmogorov-Chepmen system
    and calculates P0(+inf) + P1(+inf).
    """
    # Step 1: Define the given transition rates (lambda values).
    lambda_01 = 0.019
    lambda_10 = 0.65
    lambda_12 = 0.4
    lambda_21 = 0.392
    lambda_23 = 0.008
    lambda_31 = 0.008

    # Step 2: At steady state (t -> +inf), P_i'(t) = 0. Let pi_i = P_i(+inf).
    # The system of differential equations becomes a system of linear equations:
    #   (1)  lambda_10 * pi_1 - lambda_01 * pi_0 = 0
    #   (2)  lambda_12 * pi_1 - (lambda_21 + lambda_23) * pi_2 = 0
    #   (3)  lambda_23 * pi_2 - lambda_31 * pi_3 = 0
    #   (4)  pi_0 + pi_1 + pi_2 + pi_3 = 1

    print("Step 1: Set up the steady-state equations from the given system.")
    print(f"From P0'(t) = 0: {lambda_10} * \u03C0\u2081 - {lambda_01} * \u03C0\u2080 = 0")
    print(f"From P2'(t) = 0: {lambda_12} * \u03C0\u2081 - ({lambda_21} + {lambda_23}) * \u03C0\u2082 = 0")
    print(f"From P3'(t) = 0: {lambda_23} * \u03C0\u2082 - {lambda_31} * \u03C0\u2083 = 0")
    print("Normalization: \u03C0\u2080 + \u03C0\u2081 + \u03C0\u2082 + \u03C0\u2083 = 1\n")

    # Step 3: Express pi_0, pi_2, pi_3 in terms of pi_1.
    # From (1): pi_0 = (lambda_10 / lambda_01) * pi_1
    c0 = lambda_10 / lambda_01
    
    # From (2): pi_2 = (lambda_12 / (lambda_21 + lambda_23)) * pi_1
    c2 = lambda_12 / (lambda_21 + lambda_23)
    
    # From (3): pi_3 = (lambda_23 / lambda_31) * pi_2. Substituting the expression for pi_2:
    # pi_3 = (lambda_23 / lambda_31) * c2 * pi_1
    c3 = (lambda_23 / lambda_31) * c2
    
    print("Step 2: Solve for relationships between probabilities based on the numerical values.")
    print(f"\u03C0\u2080 = ({lambda_10} / {lambda_01}) * \u03C0\u2081 = {c0:.4f} * \u03C0\u2081")
    print(f"\u03C0\u2082 = ({lambda_12} / ({lambda_21} + {lambda_23})) * \u03C0\u2081 = ({lambda_12} / {lambda_21 + lambda_23}) * \u03C0\u2081 = {c2} * \u03C0\u2081")
    print(f"\u03C0\u2083 = ({lambda_23} / {lambda_31}) * \u03C0\u2082 = ({lambda_23 / lambda_31}) * ({c2} * \u03C0\u2081) = {c3} * \u03C0\u2081")
    print("This shows a simple relationship: \u03C0\u2081 = \u03C0\u2082 = \u03C0\u2083.\n")
    
    # Step 4: Use the normalization condition to find the sum P0(+inf) + P1(+inf).
    # The sum can be found from the relationship pi_0 + pi_1 + pi_2 + pi_3 = 1
    # pi_0 + pi_1 + pi_1 + pi_1 = 1  => (c0*pi_1) + 3*pi_1 = 1 => pi_1 = 1/(c0+3)
    # The sum is pi_0 + pi_1 = (c0+1)*pi_1 = (c0+1)/(c0+3)
    # Substituting c0 = lambda_10/lambda_01:
    # sum = (lambda_10/lambda_01 + 1) / (lambda_10/lambda_01 + 3)
    # sum = (lambda_10 + lambda_01) / (lambda_10 + 3*lambda_01)

    numerator = lambda_10 + lambda_01
    denominator = lambda_10 + 3 * lambda_01
    result = numerator / denominator
    
    print("Step 3: Calculate the final sum P\u2080(+\u221E) + P\u2081(+\u221E).")
    print("The final equation for the sum is:")
    print(f"  P\u2080(+\u221E) + P\u2081(+\u221E) = (\u03BB\u2081\u2080 + \u03BB\u2080\u2081) / (\u03BB\u2081\u2080 + 3 * \u03BB\u2080\u2081)")
    print(f"  = ({lambda_10} + {lambda_01}) / ({lambda_10} + 3 * {lambda_01})")
    print(f"  = {numerator} / {denominator}")
    print(f"  = {result}")

solve_markov_chain()
<<<0.9462517680339463>>>