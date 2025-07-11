import math

def calculate_symmetry_factor():
    """
    Calculates the symmetry factor for a Feynman diagram with 2 vertices
    and 4 field contractions (internal lines) between them.
    """
    num_vertices = 2
    num_lines = 4

    # Calculate the number of permutations for identical vertices
    vertex_permutations = math.factorial(num_vertices)

    # Calculate the number of permutations for identical lines connecting the vertices
    line_permutations = math.factorial(num_lines)

    # The total symmetry factor is the product of these factors
    symmetry_factor = vertex_permutations * line_permutations

    # Print the explanation and the final calculation
    print("The symmetry factor S is calculated by multiplying the number of permutations of identical vertices by the number of permutations of identical internal lines.")
    print(f"Number of vertices = {num_vertices}")
    print(f"Number of internal lines = {num_lines}")
    print("\nCalculation:")
    print(f"S = (permutations of vertices) * (permutations of lines)")
    print(f"S = {num_vertices}! * {num_lines}!")
    print(f"S = {vertex_permutations} * {line_permutations}")
    print(f"S = {symmetry_factor}")

calculate_symmetry_factor()