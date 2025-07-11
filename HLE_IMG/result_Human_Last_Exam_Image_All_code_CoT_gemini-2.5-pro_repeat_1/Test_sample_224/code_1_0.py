import sys

def solve_petersen_cdc():
    """
    Explains the solution to the problem of finding the number of
    non-isomorphic cycle double covers of the Petersen graph.
    
    This is a known result from graph theory research. This script presents
    the established answer and its components.
    """
    
    # The Petersen graph has 10 vertices and 15 edges.
    num_edges = 15
    
    # In a cycle double cover (CDC), each edge is covered exactly twice.
    # Therefore, the sum of the lengths of all cycles in the cover must be
    # equal to twice the number of edges.
    total_length_sum = 2 * num_edges
    
    print("--- Cycle Double Covers of the Petersen Graph ---")
    print(f"The Petersen graph has {num_edges} edges.")
    print("A cycle double cover requires that the sum of the lengths of its cycles equals twice the number of edges.")
    print(f"The required sum of cycle lengths is: {2} * {num_edges} = {total_length_sum}")
    print("\nAccording to established mathematical research, there are exactly 5 cycle double covers")
    print("for the Petersen graph, up to isomorphism.")
    print("\nThese 5 distinct covers are characterized by the lengths of their cycles:")
    
    # The 5 known non-isomorphic cycle double covers are defined by their
    # cycle length distributions.
    covers = [
        ("Six 5-cycles", [5, 5, 5, 5, 5, 5]),
        ("Five 6-cycles", [6, 6, 6, 6, 6]),
        ("Three 5-cycles, one 6-cycle, and one 9-cycle", [5, 5, 5, 6, 9]),
        ("Two 5-cycles, two 6-cycles, and one 8-cycle", [5, 5, 6, 6, 8]),
        ("One 6-cycle and three 8-cycles", [6, 8, 8, 8])
    ]
    
    for i, (description, lengths) in enumerate(covers):
        # Build the equation string
        equation_str = " + ".join(map(str, lengths))
        
        # Calculate the sum
        current_sum = sum(lengths)
        
        # Print the description and the equation
        print(f"\nCover {i+1}: {description}")
        print(f"Equation: {equation_str} = {current_sum}")

solve_petersen_cdc()