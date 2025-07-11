import sys

def solve_moduli_volume_problem():
    """
    Solves the two-part question about the volume of the moduli space
    of oriented metric ribbon graphs.
    """
    
    # Part (a): Continuity
    # The function Z is a volume of a geometric object parametrized by lengths L.
    # Such volumes are continuous functions. While 'piecewise polynomial' alone
    # doesn't imply continuity, the geometric context does. The polynomial
    # pieces must match at the cell boundaries.
    answer_a = "Yes"
    
    # Part (b): Degree of the polynomial
    # The parameters for the specific case are:
    g = 0
    n_plus = 3
    n_minus = 1
    
    # The degree of the polynomial Z_{g, n_+, n_-} is given by the formula:
    # Degree = 4*g - 4 + n_+ + n_-
    # This formula is valid for (g, n_+, n_-) != (0,1,1) and (0,2,0).
    # Our case (0, 3, 1) is not one of these exceptions.
    
    degree = 4 * g - 4 + n_plus + n_minus
    
    # Print the answers and the calculation for clarity.
    print("Part (a): Does the property of piecewise polynomiality of Z imply that it is continuous?")
    print(f"Answer: {answer_a}. In the context of geometric volumes, the function is continuous. The polynomial pieces must join continuously at the boundaries between cells.\n")
    
    print("Part (b): For g = 0, n_+ = 3, and n_- = 1, determine the degree of the polynomial Z.")
    print("Using the formula: Degree = 4*g - 4 + n_+ + n_-")
    print(f"Calculation: Degree = 4 * {g} - 4 + {n_plus} + {n_minus} = {degree}")
    print(f"Answer: The degree of the polynomial is {degree}.\n")
    
    # Final answer in the required format
    final_answer = f"<<<[{answer_a}]; [{degree}]>>>"
    
    # This check is to prevent the final answer string from being printed
    # when the script is imported, as per good practice.
    if __name__ == "__main__":
        print(final_answer, file=sys.stdout)

if __name__ == "__main__":
    solve_moduli_volume_problem()