import math

def calculate_symmetry_factor():
    """
    Calculates the symmetry factor for a Feynman diagram with two vertices
    connected by four propagators.
    """
    # The diagram has two identical vertices that can be swapped.
    num_identical_vertices = 2
    vertex_permutations = math.factorial(num_identical_vertices)

    # The diagram has four identical lines connecting the two vertices.
    # These lines can be permuted among themselves.
    num_identical_lines = 4
    line_permutations = math.factorial(num_identical_lines)

    # The total symmetry factor is the product of the permutations of
    # identical vertices and identical lines.
    symmetry_factor = vertex_permutations * line_permutations

    print("This Feynman diagram consists of two vertices connected by four internal lines.")
    print("The symmetry factor 'S' is the product of the number of permutations of its identical components.")
    print(f"1. Permutations of {num_identical_vertices} identical vertices: {num_identical_vertices}! = {vertex_permutations}")
    print(f"2. Permutations of {num_identical_lines} identical lines: {num_identical_lines}! = {line_permutations}")
    print("\nFinal Equation:")
    print(f"S = {vertex_permutations} * {line_permutations} = {symmetry_factor}")

if __name__ == "__main__":
    calculate_symmetry_factor()