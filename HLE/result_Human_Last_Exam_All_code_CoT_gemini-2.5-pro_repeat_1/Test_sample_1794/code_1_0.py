import sys

def solve_feynman_queries():
    """
    This function calculates and explains the solutions to the two questions
    about 3-loop diagrams in massless phi^3 theory.
    """

    # --- Part 1: Counting distinct planar graphs ---

    print("--- Part 1: Number of Distinct Planar Graphs ---")
    print("The task is to count the number of 3-loop, 4-point, planar Feynman diagrams in phi^3 theory, excluding those with vertex corrections.")
    print("\nInterpretation:")
    print("In d=4 dimensions, primitive 1-loop self-energy graphs are UV divergent (logarithmically), while primitive 1-loop vertex graphs are UV finite.")
    print("This suggests a physical reason for the specific wording: we should exclude graphs with 1PI 3-point subgraphs (vertex corrections) but include those with 1PI 2-point subgraphs (self-energy corrections).")
    print("\nCounting Procedure:")
    print("The total number of graphs is the sum of graphs from three categories:")
    
    # Category 1: Primitive 3-loop graphs
    num_primitive_3_loop = 6
    print(f"1. Primitive 3-loop graphs (no self-energy or vertex corrections): {num_primitive_3_loop}")

    # Category 2: 1-loop self-energy insertions on 2-loop graphs
    # There are 3 distinct 2-loop 1PI planar graphs: double-box, penta-box, and box-with-bubble.
    # We count the number of non-equivalent propagators to insert a 1-loop self-energy bubble.
    num_insert_on_double_box = 3
    num_insert_on_penta_box = 4 # Based on known symmetries
    num_insert_on_box_bubble = 3
    num_from_2_loop = num_insert_on_double_box + num_insert_on_penta_box + num_insert_on_box_bubble
    print(f"2. Graphs from 1-loop self-energy insertions on 2-loop graphs: {num_insert_on_double_box} + {num_insert_on_penta_box} + {num_insert_on_box_bubble} = {num_from_2_loop}")

    # Category 3: 2-loop self-energy insertions on 1-loop graphs
    # The 1-loop graph is the box, which has one type of propagator by symmetry.
    # There are 2 types of 2-loop self-energy graphs: the primitive 'sunrise' diagram, and a non-primitive line with two 1-loop bubbles.
    num_from_1_loop = 2
    print(f"3. Graphs from 2-loop self-energy insertions on the 1-loop box graph: {num_from_1_loop}")
    
    # Final Calculation for Part 1
    total_graphs = num_primitive_3_loop + num_from_2_loop + num_from_1_loop
    print("\nFinal Calculation (Part 1):")
    print(f"Total number of graphs = (primitive) + (insertions on 2-loop) + (insertions on 1-loop)")
    print(f"Total number of graphs = {num_primitive_3_loop} + {num_from_2_loop} + {num_from_1_loop} = {total_graphs}")
    
    print("\n" + "="*50 + "\n")

    # --- Part 2: Power of the leading divergence ---

    print("--- Part 2: Power of the Leading Divergent Term ---")
    print("The task is to find the power of the leading pole in the epsilon expansion of the 3-loop massless on-shell amplitude.")
    print("\nAnalysis:")
    print("For scattering amplitudes of massless particles, infrared (soft and collinear) divergences dominate.")
    print("In dimensional regularization (d = 4 - 2*epsilon), these appear as poles in epsilon.")
    print("A general result of QFT states that the leading divergence for an L-loop amplitude has a pole of order 2L.")
    
    # Calculation for Part 2
    L = 3
    leading_power = 2 * L
    
    print("\nFinal Calculation (Part 2):")
    print("The power of the leading divergent term is given by the formula: 2 * L")
    print(f"For L = {L} loops, the calculation is:")
    print(f"Power = 2 * {L} = {leading_power}")
    print(f"\nThis means the leading divergent term behaves as 1/epsilon^{leading_power}.")
    
    # Storing final numerical answers
    # This part is for the final answer block, not for user output.
    global answer1, answer2
    answer1 = total_graphs
    answer2 = leading_power


# Execute the function to print the explanations and results
solve_feynman_queries()

# The final answer is then formatted as requested.
# For multiple answers, they are presented in a list.
sys.stdout.write(f"\n<<<{[answer1, answer2]}>>>\n")