import math

def calculate_symmetry_factor():
    """
    Calculates the symmetry factor for a Feynman diagram with two vertices
    connected by four internal lines.

    The symmetry factor (S) is calculated as:
    S = (permutations of identical vertices) * (permutations of identical lines)
    """

    num_vertices = 2
    num_lines = 4

    # Calculate the factorial for vertex and line permutations
    vertex_permutations = math.factorial(num_vertices)
    line_permutations = math.factorial(num_lines)

    # The total symmetry factor is the product of these values
    symmetry_factor = vertex_permutations * line_permutations

    print("The symmetry factor S is calculated by multiplying the number of ways to permute identical components.")
    print(f"Number of permutations for {num_vertices} identical vertices = {num_vertices}! = {vertex_permutations}")
    print(f"Number of permutations for {num_lines} identical lines = {num_lines}! = {line_permutations}")
    print("\nThe final equation is:")
    print(f"S = {vertex_permutations} * {line_permutations} = {symmetry_factor}")


if __name__ == "__main__":
    calculate_symmetry_factor()