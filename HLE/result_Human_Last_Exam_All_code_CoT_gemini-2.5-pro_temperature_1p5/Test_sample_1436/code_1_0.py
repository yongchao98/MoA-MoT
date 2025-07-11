import math

def calculate_symmetry_factor():
    """
    Calculates the symmetry factor for a Feynman diagram with 2 vertices
    and 4 propagators connecting them (a "figure-eight" diagram).
    """
    
    # Define the diagram's components
    num_vertices = 2
    num_propagators = 4

    # --- Introduction ---
    print("Calculating the symmetry factor (S) for a 'figure-eight' Feynman diagram.")
    print("This diagram has 2 vertices and 4 propagators connecting them.\n")
    
    # --- Step 1: Vertex Permutation Symmetry ---
    print(f"Step 1: Calculate the factor from permuting the {num_vertices} identical vertices.")
    # The two vertices are identical and can be swapped. This is a permutation of 2 items.
    vertex_factor = math.factorial(num_vertices)
    print(f"The symmetry factor from vertex permutation is {num_vertices}! = {vertex_factor}.\n")

    # --- Step 2: Propagator Permutation Symmetry ---
    print(f"Step 2: Calculate the factor from permuting the {num_propagators} identical propagators.")
    # The four propagators all connect the same two vertices and are indistinguishable.
    propagator_factor = math.factorial(num_propagators)
    print(f"The symmetry factor from propagator permutation is {num_propagators}! = {propagator_factor}.\n")
    
    # --- Step 3: Total Symmetry Factor Calculation ---
    print("Step 3: Multiply the individual factors to get the total symmetry factor S.")
    
    total_symmetry_factor = vertex_factor * propagator_factor
    
    # The final output must show the numbers in the equation itself.
    print(f"The total symmetry factor S = (vertex factor) * (propagator factor)")
    print(f"S = {vertex_factor} * {propagator_factor} = {total_symmetry_factor}")


if __name__ == "__main__":
    calculate_symmetry_factor()