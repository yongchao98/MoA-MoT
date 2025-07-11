import math

def calculate_symmetry_factor():
    """
    Calculates the symmetry factor for a Feynman diagram with two vertices
    connected by four propagators.
    """
    # Number of identical vertices that can be interchanged
    num_vertices = 2
    # Number of identical propagators connecting the same two vertices
    num_propagators = 4

    # Calculate the contribution from permuting the vertices
    vertex_symmetry = math.factorial(num_vertices)
    
    # Calculate the contribution from permuting the propagators
    propagator_symmetry = math.factorial(num_propagators)

    # The total symmetry factor is the product of these contributions
    total_symmetry_factor = vertex_symmetry * propagator_symmetry

    print("Calculating the symmetry factor for a diagram with 2 vertices and 4 propagators:")
    print("1. Symmetry from interchanging the two identical vertices:")
    print(f"   Factor = {num_vertices}! = {vertex_symmetry}")
    print("2. Symmetry from permuting the four identical propagators connecting the vertices:")
    print(f"   Factor = {num_propagators}! = {propagator_symmetry}")
    print("\nTotal symmetry factor S is the product of these factors.")
    print(f"S = {num_vertices}! * {num_propagators}! = {vertex_symmetry} * {propagator_symmetry} = {total_symmetry_factor}")
    
    # Final answer in the required format
    # print(f"\n<<<{total_symmetry_factor}>>>")

calculate_symmetry_factor()