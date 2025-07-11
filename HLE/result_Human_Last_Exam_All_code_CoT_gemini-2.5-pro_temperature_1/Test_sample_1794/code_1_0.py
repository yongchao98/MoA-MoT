def solve_feynman_questions():
    """
    This function calculates and prints the answers to the two questions
    based on the principles of quantum field theory.
    """

    # --- Question 1: Number of distinct planar graphs ---
    print("--- Question 1: Number of Diagrams ---")
    
    # Total number of planar 3-loop 4-point graphs in phi^3 theory
    total_graphs = 6
    
    # Number of graphs with 1-loop vertex correction subgraphs
    vertex_correction_graphs = 1
    
    # The problem asks to exclude diagrams with vertex corrections.
    # We subtract the number of graphs with vertex corrections from the total.
    distinct_graphs = total_graphs - vertex_correction_graphs
    
    print(f"The total number of planar 3-loop 4-point graphs is {total_graphs}.")
    print(f"The number of these graphs that have vertex correction subgraphs is {vertex_correction_graphs}.")
    print("The final count is obtained by the equation:")
    print(f"{total_graphs} - {vertex_correction_graphs} = {distinct_graphs}")
    print(f"So, there are {distinct_graphs} distinct planar graphs when excluding those with vertex corrections.")
    print("-" * 35)
    print()

    # --- Question 2: Power of the leading divergent term ---
    print("--- Question 2: Power of Divergence ---")
    
    # The number of loops in the diagrams
    L = 3
    
    # For massless on-shell scattering amplitudes, the maximal power of the 1/epsilon
    # pole at L-loop order is given by the formula 2*L due to infrared divergences.
    leading_power = 2 * L
    
    print(f"The diagrams are at L = {L} loop order.")
    print("The power of the leading IR divergence is given by the equation:")
    print(f"2 * L = 2 * {L} = {leading_power}")
    print(f"So, the power of the leading divergent term is {leading_power}.")
    print("-" * 35)
    
    # Store the final numerical answers
    answer1 = distinct_graphs
    answer2 = leading_power
    
    # The final answer format is not requested for this problem, but we can print it for clarity
    # print(f"\nFinal numerical answers: Q1 = {answer1}, Q2 = {answer2}")


solve_feynman_questions()