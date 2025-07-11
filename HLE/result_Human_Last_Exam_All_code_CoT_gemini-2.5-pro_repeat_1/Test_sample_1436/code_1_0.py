import math

def calculate_symmetry_factor():
    """
    Calculates the symmetry factor for a Feynman diagram with two vertices
    connected by four identical propagators.
    """
    # Number of identical vertices that can be permuted
    num_vertices = 2

    # Number of identical propagators connecting the vertices that can be permuted
    num_propagators = 4

    # The symmetry contribution from swapping vertices is n!
    vertex_symmetry = math.factorial(num_vertices)

    # The symmetry contribution from swapping identical propagators is m!
    propagator_symmetry = math.factorial(num_propagators)

    # The total symmetry factor is the product of all symmetry contributions
    total_symmetry_factor = vertex_symmetry * propagator_symmetry

    print("Calculating the symmetry factor 'S' for a diagram with 2 vertices and 4 propagators connecting them.")
    print("-" * 80)
    print("The total symmetry is the product of the symmetries of its components:")
    print(f"1. Symmetry from permuting the {num_vertices} identical vertices: {num_vertices}!")
    print(f"2. Symmetry from permuting the {num_propagators} identical propagators: {num_propagators}!")
    print("-" * 80)
    print("The final equation is:")
    print(f"S = ({num_vertices}!) * ({num_propagators}!) = {vertex_symmetry} * {propagator_symmetry} = {total_symmetry_factor}")

calculate_symmetry_factor()