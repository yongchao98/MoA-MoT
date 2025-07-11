import math

def calculate_symmetry_factor():
    """
    Calculates the symmetry factor for a Feynman diagram with
    two vertices connected by four internal lines.
    """

    # Define the parameters of the Feynman diagram
    num_vertices = 2
    num_internal_lines = 4

    # --- Explanation ---
    print("Calculating the symmetry factor (S) for a Feynman diagram.")
    print(f"The diagram has {num_vertices} vertices and {num_internal_lines} identical lines connecting them.")
    print("The total symmetry factor is the product of the number of permutations")
    print("of identical components that leave the diagram unchanged.")
    print("-" * 50)

    # --- Step 1: Calculate symmetry from permuting vertices ---
    print("Step 1: Symmetry from permuting identical vertices.")
    print(f"Since the {num_vertices} vertices are identical, they can be swapped.")
    sym_vertices = math.factorial(num_vertices)
    print(f"This contributes a factor of {num_vertices}! = {sym_vertices}.")
    print()

    # --- Step 2: Calculate symmetry from permuting internal lines ---
    print("Step 2: Symmetry from permuting identical internal lines.")
    print(f"The {num_internal_lines} lines connecting the vertices are identical and can be permuted.")
    sym_lines = math.factorial(num_internal_lines)
    print(f"This contributes a factor of {num_internal_lines}! = {sym_lines}.")
    print()

    # --- Step 3: Calculate the total symmetry factor ---
    total_symmetry_factor = sym_vertices * sym_lines

    print("-" * 50)
    print("The total symmetry factor S is the product of these factors.")
    print("Final Equation:")
    print(f"S = (Vertex Permutations) * (Line Permutations)")
    print(f"S = {sym_vertices} * {sym_lines} = {total_symmetry_factor}")
    print("-" * 50)

if __name__ == "__main__":
    calculate_symmetry_factor()