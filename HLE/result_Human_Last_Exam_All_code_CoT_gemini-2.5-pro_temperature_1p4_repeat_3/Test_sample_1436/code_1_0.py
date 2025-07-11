import math

def calculate_feynman_symmetry_factor():
    """
    Calculates the symmetry factor for a Feynman diagram with two vertices
    connected by four propagators.
    """
    # Number of interchangeable vertices
    num_vertices = 2
    # Number of identical propagators connecting the vertices
    num_propagators = 4

    # The symmetry factor from permuting identical vertices is n!
    vertex_symmetry = math.factorial(num_vertices)

    # The symmetry factor from permuting identical propagators is m!
    propagator_symmetry = math.factorial(num_propagators)

    # The total symmetry factor is the product of the individual factors
    total_symmetry_factor = vertex_symmetry * propagator_symmetry

    # Print the explanation and step-by-step calculation
    print("Calculating the symmetry factor (S) for a Feynman diagram with 2 vertices and 4 propagators:")
    print("-" * 80)
    
    print(f"1. Vertex Symmetry: There are {num_vertices} identical vertices that can be interchanged.")
    print(f"   Contribution = {num_vertices}! = {vertex_symmetry}")
    
    print(f"\n2. Propagator Symmetry: There are {num_propagators} identical propagators connecting the vertices.")
    print(f"   Contribution = {num_propagators}! = {propagator_symmetry}")
    
    print("\n3. Total Symmetry Factor (S): This is the product of the individual symmetry contributions.")
    print("   S = (Vertex Symmetry) * (Propagator Symmetry)")
    
    # Final output showing the equation with numbers
    print(f"\nFinal Equation:")
    print(f"S = {vertex_symmetry} * {propagator_symmetry} = {total_symmetry_factor}")

if __name__ == '__main__':
    calculate_feynman_symmetry_factor()