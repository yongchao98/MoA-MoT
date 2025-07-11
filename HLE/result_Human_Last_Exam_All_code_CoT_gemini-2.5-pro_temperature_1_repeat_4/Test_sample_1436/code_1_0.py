import math

def calculate_symmetry_factor():
    """
    Calculates the symmetry factor for a Feynman diagram with two vertices
    connected by four propagators.
    """
    
    # Define the parameters of the diagram
    num_vertices = 2
    num_propagators = 4

    # Explain the context and method
    print("The Feynman diagram described has two vertices and four field contractions.")
    print("This corresponds to a vacuum bubble diagram in a phi^4-type theory, where two vertices are connected by four propagators.")
    print("\nThe symmetry factor, S, is the number of permutations of the diagram's internal parts that leave the topology unchanged.")
    print("S can be calculated as the product of the symmetries from its vertices and propagators (lines).")

    # Step 1: Calculate the symmetry from interchanging identical vertices.
    # Since the two vertices are identical, they can be swapped.
    vertex_symmetry = math.factorial(num_vertices)
    print(f"\nStep 1: Calculate Vertex Symmetry")
    print(f"There are {num_vertices} identical vertices which can be interchanged.")
    print(f"This contributes a factor of {num_vertices}! = {vertex_symmetry}.")

    # Step 2: Calculate the symmetry from interchanging identical propagators.
    # The four propagators connecting the vertices are identical and can be permuted.
    propagator_symmetry = math.factorial(num_propagators)
    print(f"\nStep 2: Calculate Propagator Symmetry")
    print(f"There are {num_propagators} identical propagators connecting the two vertices which can be permuted.")
    print(f"This contributes a factor of {num_propagators}! = {propagator_symmetry}.")

    # Step 3: Calculate the total symmetry factor.
    # The total symmetry factor is the product of the individual symmetries.
    total_symmetry_factor = vertex_symmetry * propagator_symmetry
    print(f"\nStep 3: Calculate the Total Symmetry Factor")
    print("The total symmetry factor S is the product of the vertex and propagator symmetries.")
    print(f"S = (Vertex Symmetry) * (Propagator Symmetry)")
    print(f"S = {vertex_symmetry} * {propagator_symmetry} = {total_symmetry_factor}")
    
    print(f"\nThe final symmetry factor associated with the diagram is {total_symmetry_factor}.")

if __name__ == "__main__":
    calculate_symmetry_factor()