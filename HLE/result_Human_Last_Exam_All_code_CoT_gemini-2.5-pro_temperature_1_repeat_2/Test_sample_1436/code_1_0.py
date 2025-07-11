import math

def calculate_symmetry_factor():
    """
    Calculates the symmetry factor for a Feynman diagram with two vertices
    connected by four propagators.
    """
    # Number of identical vertices that can be permuted
    num_vertices = 2
    # Number of identical lines connecting the vertices that can be permuted
    num_lines = 4

    # Calculate the symmetry contributions
    # 1. Permuting the two identical vertices
    vertex_symmetry = math.factorial(num_vertices)
    # 2. Permuting the four identical lines
    line_symmetry = math.factorial(num_lines)

    # The total symmetry factor is the product of the individual factors
    total_symmetry_factor = vertex_symmetry * line_symmetry

    # Output the explanation and calculation
    print("This diagram consists of 2 vertices connected by 4 lines.")
    print("The symmetry factor S is calculated by considering the permutations of identical components.")
    print("-" * 30)

    print(f"Symmetry from permuting {num_vertices} identical vertices = {num_vertices}! = {vertex_symmetry}")
    print(f"Symmetry from permuting {num_lines} identical lines = {num_lines}! = {line_symmetry}")
    print("-" * 30)

    print("The total symmetry factor is the product:")
    print(f"S = (Vertex Symmetry) * (Line Symmetry)")
    print(f"S = {vertex_symmetry} * {line_symmetry}")
    print(f"S = {total_symmetry_factor}")

if __name__ == '__main__':
    calculate_symmetry_factor()