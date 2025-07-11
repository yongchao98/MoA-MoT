def solve_moduli_volume_problem():
    """
    Solves the two-part question about the volume of the moduli space
    of oriented metric ribbon graphs.
    """
    
    # Part (a) reasoning:
    # The piecewise polynomiality of the volume Z implies continuity.
    # This is because Z represents a geometric volume, and the polynomial
    # pieces are known to match up at the boundaries between cells in the
    # parameter space.
    answer_a = "Yes"

    # Part (b) calculation:
    # The degree of the polynomial Z is given by the dimension of the
    # moduli space M_{g,n}, which is 3g - 3 + n.
    g = 0
    n_plus = 3
    n_minus = 1

    # Calculate the total number of boundaries n
    n = n_plus + n_minus
    
    # Calculate the degree of the polynomial
    degree = 3 * g - 3 + n
    
    print("Calculating the degree for part (b):")
    print(f"Genus g = {g}")
    print(f"Positive boundaries n_+ = {n_plus}")
    print(f"Negative boundaries n_- = {n_minus}")
    print(f"Total boundaries n = n_+ + n_- = {n_plus} + {n_minus} = {n}")
    print("The degree is calculated using the formula: 3*g - 3 + n")
    print("Final calculation:")
    print(f"3 * {g} - 3 + {n} = {degree}")

    answer_b = degree
    
    # Format the final answer as requested.
    final_answer = f"<<<(a) {answer_a}; (b) {answer_b}>>>"
    print("\nFinal Answer:")
    print(final_answer)

solve_moduli_volume_problem()