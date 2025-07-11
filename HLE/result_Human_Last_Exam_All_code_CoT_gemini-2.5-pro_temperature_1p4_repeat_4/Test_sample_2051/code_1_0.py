import sys

def solve_moduli_volume_properties():
    """
    Solves the two-part question about the properties of the volume of the moduli space
    of oriented metric ribbon graphs, Z_{g, n_+, n_-}.
    """
    
    # --- Part (a): Continuity ---
    # The function Z_{g, n_+, n_-} represents a geometric volume. While a generic
    # piecewise polynomial function is not guaranteed to be continuous, in the context
    # of volumes of moduli spaces, continuity is a known property. The polynomial
    # pieces defined on different cells of the parameter space are known to agree on
    # the boundaries between these cells, ensuring the overall function is continuous.
    answer_a = "Yes"

    # --- Part (b): Degree of the Polynomial ---
    # The degree 'd' of the volume polynomial Z_{g, n_+, n_-} is given by the
    # formula: d = 3g - 3 + n, where n = n_+ + n_- is the total number of boundaries.
    
    # Given parameters for the specific case:
    g = 0       # Genus
    n_plus = 3  # Number of positively oriented boundaries
    n_minus = 1 # Number of negatively oriented boundaries

    # Calculate the total number of boundaries, n
    n = n_plus + n_minus
    
    # Calculate the degree of the polynomial using the formula
    degree = 3 * g - 3 + n

    # --- Output the Explanation and Results ---
    print("--- Analysis and Solution ---")
    print("\n(a) Does the property of piecewise polynomiality imply continuity?")
    print("Answer: Yes.")
    print("Explanation: The function Z represents a geometric volume. For such functions derived from moduli spaces, the polynomial pieces defined on different regions are known to join continuously at their boundaries. Therefore, the function is continuous.")

    print("\n(b) What is the degree of the polynomial Z_{0,3,1}?")
    print("Explanation: The degree 'd' is calculated using the formula d = 3g - 3 + n, where n is the total number of boundaries.")
    
    # Print the step-by-step calculation as required
    print(f"The parameters are: g = {g}, n_+ = {n_plus}, n_- = {n_minus}.")
    print(f"The total number of boundaries is n = {n_plus} + {n_minus} = {n}.")
    print("The calculation for the degree 'd' is:")
    print(f"d = 3*g - 3 + n = 3*{g} - 3 + {n} = {degree}")
    print(f"Answer: The degree is {degree}.")
    
    # Format the final answer string as specified in the problem.
    answer_b = degree
    final_answer_string = f"<<<(a) [{answer_a}]; (b) [{answer_b}]>>>"
    
    # Writing to stdout which is standard for such tasks.
    sys.stdout.write(f"\n{final_answer_string}\n")

if __name__ == '__main__':
    solve_moduli_volume_properties()
