import math

def calculate_symmetry_factor():
    """
    Calculates the symmetry factor for a Feynman diagram with two vertices
    connected by four propagators.
    """

    # Number of identical vertices that can be permuted
    num_vertices = 2
    # Number of identical propagators connecting the two vertices
    num_propagators = 4

    # Calculate the symmetry factor from permuting vertices
    vertex_symmetry = math.factorial(num_vertices)

    # Calculate the symmetry factor from permuting propagators
    propagator_symmetry = math.factorial(num_propagators)

    # The total symmetry factor is the product of these individual factors
    total_symmetry_factor = vertex_symmetry * propagator_symmetry

    print("The symmetry factor (S) is the product of the symmetries from permuting identical vertices and identical propagators.")
    print("\nS = (Symmetry of vertices) * (Symmetry of propagators)")
    print(f"S = {num_vertices}! * {num_propagators}!")
    print(f"S = {vertex_symmetry} * {propagator_symmetry}")
    print(f"\nThe calculated symmetry factor is: {total_symmetry_factor}")

calculate_symmetry_factor()