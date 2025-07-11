import sys

def solve_moduli_volume_properties():
    """
    Solves the two-part question about the properties of Z_g,n.
    
    Part (a) addresses the continuity of the function.
    Part (b) calculates the degree of the polynomial for a specific case.
    """
    
    # Part (a): Reasoning for continuity
    # The function Z_g,n is known to be continuous from the theory of moduli spaces.
    answer_a = "Yes"
    
    # Part (b): Calculation of the polynomial degree
    # Given parameters
    g = 0
    n_plus = 3
    n_minus = 1
    
    # Total number of boundaries
    n = n_plus + n_minus
    
    # Formula for the degree of the polynomial Z_g,n
    # Degree = 6g - 6 + 2n
    degree = 6 * g - 6 + 2 * n
    
    # Print the final combined answer as requested.
    print(f"The solution to the question is:")
    print(f"(a) {answer_a}")
    print(f"(b) The degree is determined by the formula: Degree = 6g - 6 + 2n.")
    print(f"    For g = {g} and n = n_+ + n_- = {n_plus} + {n_minus} = {n}, the calculation is:")
    # The final equation with numbers is printed below
    print(f"    Degree = 6 * {g} - 6 + 2 * {n} = {degree}")

solve_moduli_volume_properties()