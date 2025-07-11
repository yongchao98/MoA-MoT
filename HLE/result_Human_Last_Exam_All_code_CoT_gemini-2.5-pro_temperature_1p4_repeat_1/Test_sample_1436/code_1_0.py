import math

def calculate_symmetry_factor():
    """
    Calculates the symmetry factor for a Feynman diagram with two vertices
    connected by four propagators.
    """

    # 1. Define the components of the diagram
    num_vertices = 2
    num_propagators = 4

    # 2. Explain the physical context
    print("This script calculates the symmetry factor for a specific Feynman diagram.")
    print("The diagram is assumed to be a vacuum bubble in phi^4 theory, featuring:")
    print(f"- {num_vertices} vertices")
    print(f"- {num_propagators} propagators (field contractions) connecting these vertices.")
    print("-" * 30)

    # 3. Calculate the symmetry from permuting vertices
    # Swapping the two vertices leaves the diagram unchanged.
    vertex_symmetry = math.factorial(num_vertices)
    print("Symmetry from permuting identical vertices:")
    print(f"Factor = {num_vertices}! = {vertex_symmetry}")
    print()

    # 4. Calculate the symmetry from permuting propagators
    # The four propagators connecting the vertices are identical and can be interchanged.
    propagator_symmetry = math.factorial(num_propagators)
    print("Symmetry from permuting identical propagators:")
    print(f"Factor = {num_propagators}! = {propagator_symmetry}")
    print()

    # 5. Calculate the total symmetry factor
    total_symmetry = vertex_symmetry * propagator_symmetry

    # 6. Print the final equation and result
    print("-" * 30)
    print("The total symmetry factor S is the product of these factors.")
    print("Final Equation:")
    print(f"S = ({num_vertices}!) * ({num_propagators}!) = {vertex_symmetry} * {propagator_symmetry} = {total_symmetry}")

if __name__ == "__main__":
    calculate_symmetry_factor()