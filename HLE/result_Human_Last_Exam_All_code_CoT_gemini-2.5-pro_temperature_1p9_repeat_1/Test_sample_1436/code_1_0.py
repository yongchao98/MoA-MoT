import math

def calculate_feynman_symmetry_factor():
    """
    Calculates the symmetry factor for a Feynman diagram with two vertices
    connected by four internal lines.
    """
    
    # Number of vertices in the diagram
    num_vertices = 2
    
    # Number of identical lines (field contractions) connecting the vertices
    num_lines = 4
    
    # --- Step 1: Calculate symmetry from interchanging vertices ---
    # The two vertices are identical, so they can be swapped without
    # changing the diagram's topology.
    # The number of ways to permute them is num_vertices!
    vertex_symmetry = math.factorial(num_vertices)
    
    # --- Step 2: Calculate symmetry from interchanging lines ---
    # The four lines connecting the two vertices are identical. We can
    # permute them in any order without changing the diagram.
    # The number of ways to permute them is num_lines!
    line_symmetry = math.factorial(num_lines)
    
    # --- Step 3: Calculate the total symmetry factor ---
    # The total symmetry factor S is the product of the symmetries from
    # all independent sources.
    total_symmetry_factor = vertex_symmetry * line_symmetry
    
    print("Calculating the symmetry factor (S) for a Feynman diagram with 2 vertices and 4 connecting lines.")
    print("-" * 70)
    print("Sources of symmetry:")
    print(f"1. Permutations of the {num_vertices} identical vertices: {num_vertices}! = {vertex_symmetry}")
    print(f"2. Permutations of the {num_lines} identical lines: {num_lines}! = {line_symmetry}")
    print("\nThe total symmetry factor is the product of these individual symmetry factors.")
    print("\nFinal Equation:")
    print(f"S = {num_vertices}! * {num_lines}!")
    print(f"S = {vertex_symmetry} * {line_symmetry}")
    print(f"S = {total_symmetry_factor}")
    
calculate_feynman_symmetry_factor()