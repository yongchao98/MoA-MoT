import sys

def solve_feynman_diagram_properties():
    """
    This function solves two problems related to 3-loop, 4-point scalar phi^3 theory diagrams.
    1. It determines the number of distinct planar graphs without vertex corrections.
    2. It calculates the power of the leading divergent term in the epsilon expansion of the integral.
    """

    # --- Part 1: Number of distinct planar graphs ---
    
    print("--- Part 1: How many distinct planar graphs? ---")
    print("The problem asks for the number of distinct 3-loop, 4-point planar graphs in scalar phi^3 theory.")
    print("The condition 'excluding the diagrams with vertex corrections' means we should count 'primitive' graphs.")
    print("Primitive graphs are those that cannot be constructed by inserting 1-loop self-energy or vertex corrections onto a lower-loop graph.")
    print("The number of such graph topologies is a known combinatorial result in physics.")
    print("For a 4-point function at 3-loop order, there are 3 distinct primitive planar graph topologies.")
    
    num_graphs = 3
    print(f"\nResult for Part 1: The number of distinct planar graphs is {num_graphs}.")
    print("-" * 20)

    # --- Part 2: Power of the leading divergent term ---

    print("\n--- Part 2: What is the power of the leading divergence? ---")
    print("The diagrams are for massless particles with on-shell external legs (p_i^2 = 0) in d = 4 - 2*epsilon dimensions.")
    print("By power counting, the superficial degree of divergence shows these diagrams are Ultraviolet (UV) finite.")
    print("Thus, all divergences are of Infrared (IR) origin, arising from soft or collinear loop momenta.")
    print("For a planar L-loop amplitude with massless on-shell particles, the leading IR divergence has the form 1/epsilon^(2L).")

    # The number of loops L is given.
    L = 3
    
    # The formula for the power of the leading divergent term is k = 2 * L.
    power_numerator = 2
    power_loop_factor = L
    leading_power = power_numerator * power_loop_factor

    print(f"\nThe number of loops is L = {L}.")
    print("The power 'k' of the leading divergent term (1/epsilon^k) is calculated as:")
    
    # As requested, printing each number in the final equation.
    print(f"k = {power_numerator} * {power_loop_factor} = {leading_power}")
    
    print(f"\nResult for Part 2: The number of the power of the leading divergent term is {leading_power}.")
    print("-" * 20)
    
    # Final answer in the specified format
    # For question 1, the answer is an integer. For question 2, it is also an integer.
    # Since the request asks for a single return value, I will return the result of the second part,
    # as it's the one involving a calculation.
    return leading_power

if __name__ == '__main__':
    solve_feynman_diagram_properties()
    # The final combined answer format is ambiguous for two questions.
    # The prompt asks for <<<answer content>>> at the end. I will provide the answers to both questions.
    # For question 1, the answer is 3. For question 2, the answer is 6.
    # Let's provide them as a comma-separated string.
    sys.stdout.write("\n<<<3, 6>>>\n")
