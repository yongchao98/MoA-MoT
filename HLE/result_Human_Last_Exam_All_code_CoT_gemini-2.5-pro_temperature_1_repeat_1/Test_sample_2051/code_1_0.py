def solve_moduli_volume_question():
    """
    Calculates the degree of the polynomial Z_g,n and formats the final answer.
    """
    # (a) Continuity Analysis
    # The function Z is piecewise polynomial. As a geometric volume, it must be
    # continuous across the boundaries of the polynomial pieces.
    answer_a = "Yes"

    # (b) Degree Calculation
    # The degree of the polynomial Z_{g,n} is given by the formula 2 * (3g - 3 + n).
    
    # Parameters for the specific case
    g = 0
    n_plus = 3
    n_minus = 1

    # Total number of boundaries
    n = n_plus + n_minus
    
    # Calculation of the degree
    # We will print the step-by-step evaluation of the formula as requested.
    print("Calculating the degree for (g, n) = (0, 4):")
    print(f"Degree = 2 * (3 * g - 3 + n)")
    print(f"Degree = 2 * (3 * {g} - 3 + {n})")
    
    term_in_parentheses = 3 * g - 3 + n
    print(f"Degree = 2 * ({term_in_parentheses})")
    
    degree = 2 * term_in_parentheses
    print(f"Degree = {degree}")

    # Format the final answer for both parts
    answer_b = degree
    final_answer = f"(a) {answer_a}; (b) {answer_b}"
    
    print("\nFinal formatted answer:")
    print(final_answer)

solve_moduli_volume_question()
