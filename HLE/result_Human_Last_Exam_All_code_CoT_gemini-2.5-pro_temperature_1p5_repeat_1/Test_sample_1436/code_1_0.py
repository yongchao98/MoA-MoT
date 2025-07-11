import math

def calculate_symmetry_factor():
    """
    Calculates the symmetry factor for a Feynman diagram with two vertices
    and four lines connecting them.
    """
    num_vertices = 2
    num_lines = 4

    print("This Feynman diagram consists of two vertices connected by four identical lines (propagators).")
    print("The symmetry factor 'S' is calculated by considering the permutations of identical components that leave the diagram unchanged.")
    print("-" * 50)

    # 1. Symmetry from permuting the vertices
    vertex_permutations = math.factorial(num_vertices)
    print(f"First, the two vertices are identical. They can be swapped without changing the diagram.")
    print(f"This contributes a factor of {num_vertices}! = {vertex_permutations}.")
    
    # 2. Symmetry from permuting the lines
    line_permutations = math.factorial(num_lines)
    print(f"\nSecond, the four lines connecting the vertices are identical.")
    print(f"They can be permuted among themselves in any way.")
    print(f"This contributes a factor of {num_lines}! = {line_permutations}.")
    
    # 3. Total symmetry factor calculation
    symmetry_factor = vertex_permutations * line_permutations
    
    print("\nThe total symmetry factor is the product of these contributions.")
    print(f"Final Equation: S = ({num_vertices}!) * ({num_lines}!)")
    print(f"S = {vertex_permutations} * {line_permutations}")
    print(f"The symmetry factor S is: {symmetry_factor}")


if __name__ == "__main__":
    calculate_symmetry_factor()