import math

def calculate_symmetry_factor():
    """
    Calculates the symmetry factor for a Feynman diagram with two vertices
    connected by four identical propagators.
    """
    # Number of identical vertices that can be permuted
    num_vertices = 2
    
    # Number of identical propagators (lines) connecting the two vertices
    num_propagators = 4

    # The symmetry from permuting identical vertices is n!
    vertex_symmetry = math.factorial(num_vertices)

    # The symmetry from permuting identical propagators is k!
    propagator_symmetry = math.factorial(num_propagators)

    # The total symmetry factor is the product of these individual factors.
    total_symmetry_factor = vertex_symmetry * propagator_symmetry

    print("Assuming the diagram has 2 vertices connected by 4 identical propagators ('melon' diagram):")
    print("\nThe symmetry factor S is the product of the symmetries of its components:")
    print(f"1. Symmetry from permuting {num_vertices} identical vertices = {num_vertices}! = {vertex_symmetry}")
    print(f"2. Symmetry from permuting {num_propagators} identical propagators = {num_propagators}! = {propagator_symmetry}")
    print("\nFinal Calculation:")
    print(f"S = (Vertex Symmetry) * (Propagator Symmetry)")
    print(f"S = {vertex_symmetry} * {propagator_symmetry} = {total_symmetry_factor}")

calculate_symmetry_factor()