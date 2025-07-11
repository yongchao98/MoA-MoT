import math

def calculate_symmetry_factor():
    """
    Calculates the symmetry factor for a Feynman diagram with 2 vertices
    connected by 4 propagators.
    """
    
    # Number of identical vertices that can be interchanged
    num_vertices = 2
    
    # Number of identical propagators connecting the vertices
    num_propagators = 4
    
    # Calculate the symmetry factor from permuting the vertices
    vertex_factor = math.factorial(num_vertices)
    
    # Calculate the symmetry factor from permuting the propagators
    propagator_factor = math.factorial(num_propagators)
    
    # The total symmetry factor is the product of these individual factors
    total_symmetry_factor = vertex_factor * propagator_factor
    
    print("This Feynman diagram has two key symmetries:")
    print(f"1. The {num_vertices} vertices can be interchanged, giving a factor of {num_vertices}! = {vertex_factor}")
    print(f"2. The {num_propagators} propagators connecting them can be interchanged, giving a factor of {num_propagators}! = {propagator_factor}")
    print("\nThe total symmetry factor 'S' is the product of these factors.")
    print(f"S = {vertex_factor} * {propagator_factor} = {total_symmetry_factor}")

calculate_symmetry_factor()