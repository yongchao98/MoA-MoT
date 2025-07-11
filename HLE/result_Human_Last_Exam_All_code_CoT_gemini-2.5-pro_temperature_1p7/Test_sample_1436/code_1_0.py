import math

def calculate_feynman_symmetry_factor():
    """
    Calculates and explains the symmetry factor for a Feynman diagram
    with two vertices connected by four propagators.
    """
    # Define the characteristics of the diagram
    num_vertices = 2
    num_propagators = 4

    print("Calculating the symmetry factor (S) for a Feynman diagram with two vertices connected by four identical propagators.")
    print("-" * 80)

    # Explain the sources of symmetry
    print(f"There are two sources of symmetry in this diagram:")
    print(f"1. The {num_vertices} vertices are identical and can be swapped. This contributes a factor of {num_vertices}!")
    print(f"2. The {num_propagators} propagators are identical and can be permuted. This contributes a factor of {num_propagators}!")
    print("")

    # Calculate the value for each symmetry component
    vertex_symmetry = math.factorial(num_vertices)
    propagator_symmetry = math.factorial(num_propagators)

    # Calculate the total symmetry factor
    total_symmetry_factor = vertex_symmetry * propagator_symmetry

    # Print the final calculation, showing each step
    print("The total symmetry factor S is the product of these individual factors.")
    print(f"S = {num_vertices}! * {num_propagators}!")
    print(f"S = {vertex_symmetry} * {propagator_symmetry}")
    print(f"S = {total_symmetry_factor}")


# Execute the function to get the result
calculate_feynman_symmetry_factor()
