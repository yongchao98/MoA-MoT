import math

def calculate_symmetry_factor():
    """
    Calculates the symmetry factor for a Feynman diagram with two vertices
    connected by four propagators.
    """
    # Define the characteristics of the diagram
    num_interchangeable_vertices = 2
    num_interchangeable_propagators = 4

    # Calculate the symmetry factor from permuting the identical vertices.
    # Since there are 2 vertices, they can be swapped. This corresponds to 2! permutations.
    vertex_symmetry = math.factorial(num_interchangeable_vertices)

    # Calculate the symmetry factor from permuting the identical propagators.
    # The 4 propagators connecting the two vertices are equivalent and can be
    # permuted among themselves in 4! ways.
    propagator_symmetry = math.factorial(num_interchangeable_propagators)

    # The total symmetry factor is the product of the individual symmetry factors.
    total_symmetry_factor = vertex_symmetry * propagator_symmetry

    print("The symmetry factor (S) of a Feynman diagram is the number of permutations of its internal components that leave the diagram unchanged.")
    print("\nFor a diagram with two vertices connected by four propagators:\n")

    print(f"1. Vertex Symmetry: The diagram has {num_interchangeable_vertices} identical vertices.")
    print("   Swapping these vertices leaves the diagram topologically invariant.")
    print(f"   The number of permutations is {num_interchangeable_vertices}! = {vertex_symmetry}.")

    print(f"\n2. Propagator Symmetry: The diagram has {num_interchangeable_propagators} identical propagators connecting the vertices.")
    print("   These propagators can be interchanged without changing the diagram.")
    print(f"   The number of permutations is {num_interchangeable_propagators}! = {propagator_symmetry}.")

    print("\n3. Total Symmetry Factor:")
    print("   The total symmetry factor is the product of the vertex and propagator symmetries.")
    print("   S = (Vertex Permutations) * (Propagator Permutations)")
    print(f"   S = {vertex_symmetry} * {propagator_symmetry} = {total_symmetry_factor}")


if __name__ == '__main__':
    calculate_symmetry_factor()