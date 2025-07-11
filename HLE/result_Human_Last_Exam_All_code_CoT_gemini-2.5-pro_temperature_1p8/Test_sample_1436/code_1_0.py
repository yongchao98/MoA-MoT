import math

def calculate_symmetry_factor():
    """
    Calculates the symmetry factor for a Feynman diagram with two vertices
    connected by four identical propagators.
    """
    # Number of identical vertices
    num_vertices = 2

    # Number of identical propagators connecting the vertices
    num_propagators = 4

    # The symmetry factor from permuting the identical vertices.
    # Swapping the two vertices leaves the diagram unchanged.
    vertex_symmetry = math.factorial(num_vertices)

    # The symmetry factor from permuting the identical propagators.
    # The four lines are indistinguishable and can be interchanged.
    propagator_symmetry = math.factorial(num_propagators)

    # The total symmetry factor is the product of the individual symmetry factors.
    total_symmetry_factor = vertex_symmetry * propagator_symmetry

    print("The symmetry factor (S) of a Feynman diagram is the number of permutations of its components that leave the diagram unchanged.")
    print("\nFor a diagram with 2 vertices connected by 4 identical propagators:")
    print(f"1. Symmetry from permuting the {num_vertices} identical vertices = {num_vertices}! = {vertex_symmetry}")
    print(f"2. Symmetry from permuting the {num_propagators} identical propagators = {num_propagators}! = {propagator_symmetry}")
    print("\nThe total symmetry factor is the product of these values.")
    print("S = (Vertex Symmetry) * (Propagator Symmetry)")

    print("\nFinal Equation:")
    print(f"{vertex_symmetry} * {propagator_symmetry} = {total_symmetry_factor}")


if __name__ == "__main__":
    calculate_symmetry_factor()