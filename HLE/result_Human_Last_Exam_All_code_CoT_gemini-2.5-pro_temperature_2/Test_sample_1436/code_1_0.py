import math

def calculate_feynman_symmetry_factor():
    """
    Calculates the symmetry factor for a Feynman diagram with two vertices
    connected by four identical propagators.
    """
    # Number of indistinguishable vertices in the diagram
    num_vertices = 2

    # Number of identical propagators (lines) connecting the vertices
    num_propagators = 4

    # Calculate the symmetry factor from permuting the vertices
    # This is n! for n identical vertices.
    vertex_symmetry_factor = math.factorial(num_vertices)

    # Calculate the symmetry factor from permuting the identical lines
    # connecting the same two vertices. This is m! for m identical lines.
    propagator_symmetry_factor = math.factorial(num_propagators)

    # The total symmetry factor is the product of all individual symmetry factors.
    total_symmetry_factor = vertex_symmetry_factor * propagator_symmetry_factor

    print("Calculating the symmetry factor (S) for a diagram with 2 vertices connected by 4 propagators.")
    print("This corresponds to a vacuum bubble diagram with two sources of symmetry:")
    print(f"1. The {num_vertices} vertices are indistinguishable and can be permuted.")
    print(f"2. The {num_propagators} propagators connecting them are identical and can be permuted.")
    print("-" * 30)
    print(f"Vertex permutation factor = {num_vertices}! = {vertex_symmetry_factor}")
    print(f"Propagator permutation factor = {num_propagators}! = {propagator_symmetry_factor}")
    print("-" * 30)
    print("The total symmetry factor is the product of these factors.")
    print(f"S = {num_vertices}! * {num_propagators}! = {vertex_symmetry_factor} * {propagator_symmetry_factor} = {total_symmetry_factor}")

if __name__ == '__main__':
    calculate_feynman_symmetry_factor()