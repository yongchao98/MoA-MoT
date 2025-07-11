import math

def calculate_symmetry_factor():
    """
    Calculates the symmetry factor for a Feynman diagram with two vertices
    connected by four identical propagators.
    """
    # Number of identical vertices in the diagram
    num_vertices = 2

    # Number of identical propagators (lines) connecting the two vertices
    num_lines = 4

    # The number of ways to permute the identical vertices is n!
    vertex_permutations = math.factorial(num_vertices)

    # The number of ways to permute the identical lines is k!
    line_permutations = math.factorial(num_lines)

    # The total symmetry factor is the product of these individual permutations
    total_symmetry_factor = vertex_permutations * line_permutations

    print("To find the symmetry factor (S) for the described Feynman diagram, we follow these steps:")
    print("1. Identify the number of identical vertices that can be interchanged.")
    print(f"   - Number of identical vertices: {num_vertices}")
    print(f"   - Number of permutations for vertices: {num_vertices}! = {vertex_permutations}\n")
    print("2. Identify the number of identical lines (propagators) that can be interchanged.")
    print(f"   - Number of identical lines: {num_lines}")
    print(f"   - Number of permutations for lines: {num_lines}! = {line_permutations}\n")
    print("3. The total symmetry factor is the product of these permutation counts.")
    print(f"Final Equation: S = ({num_vertices}!) * ({num_lines}!)")
    print(f"Calculation: S = {vertex_permutations} * {line_permutations} = {total_symmetry_factor}")


if __name__ == "__main__":
    calculate_symmetry_factor()