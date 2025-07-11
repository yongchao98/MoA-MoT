import math

def calculate_symmetry_factor():
    """
    Calculates the symmetry factor for a Feynman diagram with two vertices
    connected by four propagators.
    """
    # Number of vertices in the diagram
    num_vertices = 2

    # Number of identical lines (propagators) connecting the vertices
    num_lines = 4

    # The symmetry factor includes permutations of indistinguishable vertices.
    # For 2 vertices, this is 2!
    vertex_permutations = math.factorial(num_vertices)

    # The symmetry factor also includes permutations of indistinguishable lines.
    # For 4 lines connecting the same two vertices, this is 4!
    line_permutations = math.factorial(num_lines)

    # The total symmetry factor is the product of these independent symmetries.
    total_symmetry_factor = vertex_permutations * line_permutations

    print("Calculating the symmetry factor for a diagram with 2 vertices and 4 connecting lines.")
    print("-" * 70)
    print(f"1. Permutations of the {num_vertices} indistinguishable vertices: {num_vertices}! = {vertex_permutations}")
    print(f"2. Permutations of the {num_lines} indistinguishable lines: {num_lines}! = {line_permutations}")
    print("\nThe total symmetry factor S is the product of these factors.")
    print(f"Final Equation: S = {num_vertices}! * {num_lines}!")
    print(f"Calculation: S = {vertex_permutations} * {line_permutations} = {total_symmetry_factor}")
    print("-" * 70)

calculate_symmetry_factor()