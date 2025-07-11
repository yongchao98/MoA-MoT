import math

def calculate_symmetry_factor():
    """
    Calculates the symmetry factor for a Feynman diagram with two vertices
    connected by four propagators.
    """

    # Number of identical vertices that can be permuted
    num_vertices = 2
    # Number of identical propagators connecting the vertices
    num_propagators = 4

    # Calculate the symmetry factor from permuting vertices
    vertex_symmetry = math.factorial(num_vertices)

    # Calculate the symmetry factor from permuting propagators
    propagator_symmetry = math.factorial(num_propagators)

    # The total symmetry factor is the product of these individual symmetries
    total_symmetry_factor = vertex_symmetry * propagator_symmetry

    print("This script calculates the symmetry factor (S) for a Feynman diagram with 2 vertices connected by 4 propagators.")
    print("-" * 50)

    # Explain vertex symmetry
    print("Step 1: Calculate the symmetry from permuting identical vertices.")
    print(f"There are {num_vertices} identical vertices that can be swapped, leading to a factor of {num_vertices}!.")
    print(f"Vertex Symmetry Factor = {vertex_symmetry}")
    print("")

    # Explain propagator symmetry
    print("Step 2: Calculate the symmetry from permuting identical propagators.")
    print(f"There are {num_propagators} identical propagators connecting the vertices, which can be interchanged.")
    print(f"This leads to a factor of {num_propagators}!.")
    print(f"Propagator Symmetry Factor = {propagator_symmetry}")
    print("")

    # Explain and show the final calculation
    print("Step 3: Calculate the total symmetry factor.")
    print("The total symmetry factor S is the product of the vertex and propagator symmetries.")
    print("S = (Vertex Symmetry) * (Propagator Symmetry)")
    print(f"S = {vertex_symmetry} * {propagator_symmetry}")
    print(f"S = {total_symmetry_factor}")
    print("-" * 50)
    print(f"The final symmetry factor for the diagram is {total_symmetry_factor}.")

if __name__ == '__main__':
    calculate_symmetry_factor()