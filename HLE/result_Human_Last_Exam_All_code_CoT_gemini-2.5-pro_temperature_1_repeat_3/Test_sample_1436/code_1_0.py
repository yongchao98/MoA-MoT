import math

def calculate_symmetry_factor():
    """
    Calculates the symmetry factor for a Feynman diagram with two vertices
    connected by four propagators.
    """
    # Define the parameters of the Feynman diagram
    num_vertices = 2
    num_propagators = 4

    print("This script calculates the symmetry factor 'S' for a Feynman diagram.")
    print("The diagram in question has two identical vertices and four identical propagators connecting them.\n")
    print("The formula for the symmetry factor is the product of the permutations of identical components:")
    print("S = (permutations of identical vertices) * (permutations of identical propagators)\n")

    # --- Step 1: Calculate permutations for vertices ---
    # The two vertices are identical and can be swapped.
    vertex_permutations = math.factorial(num_vertices)
    print(f"Step 1: Calculate the vertex permutation factor.")
    print(f"There are {num_vertices} identical vertices, so the factor is {num_vertices}! = {vertex_permutations}\n")

    # --- Step 2: Calculate permutations for propagators ---
    # The four propagators are identical and can be permuted.
    propagator_permutations = math.factorial(num_propagators)
    print(f"Step 2: Calculate the propagator permutation factor.")
    print(f"There are {num_propagators} identical propagators connecting the vertices, so the factor is {num_propagators}! = {propagator_permutations}\n")

    # --- Step 3: Calculate the total symmetry factor ---
    symmetry_factor = vertex_permutations * propagator_permutations
    print("Step 3: Multiply the factors to get the total symmetry factor 'S'.")
    print("\n--- Final Equation ---")
    print(f"S = {num_vertices}! * {num_propagators}!")
    print(f"S = {vertex_permutations} * {propagator_permutations}")
    print(f"S = {symmetry_factor}")

if __name__ == '__main__':
    calculate_symmetry_factor()