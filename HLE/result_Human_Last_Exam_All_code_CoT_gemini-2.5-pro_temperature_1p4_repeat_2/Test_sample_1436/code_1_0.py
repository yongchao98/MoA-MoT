import math

def calculate_symmetry_factor():
    """
    Calculates the symmetry factor for a Feynman diagram with
    two vertices and four internal lines connecting them.
    """
    # Number of identical vertices in the diagram
    num_vertices = 2

    # Number of identical lines (propagators) connecting the vertices
    num_lines = 4

    print("Calculating the symmetry factor for a Feynman diagram with two vertices and four connecting lines.")
    print("-----------------------------------------------------------------------------------------")

    # 1. Contribution from permuting identical vertices
    # The two vertices are identical and can be swapped.
    vertex_symmetry = math.factorial(num_vertices)
    print(f"There are {num_vertices} identical vertices. The number of permutations is {num_vertices}! = {vertex_symmetry}.")

    # 2. Contribution from permuting identical lines
    # The four lines connecting the vertices are identical and can be permuted.
    line_symmetry = math.factorial(num_lines)
    print(f"There are {num_lines} identical lines connecting them. The number of permutations is {num_lines}! = {line_symmetry}.")

    # 3. The total symmetry factor is the product of these individual factors
    total_symmetry_factor = vertex_symmetry * line_symmetry

    print("\nThe total symmetry factor 'S' is the product of these two numbers.")
    # Final equation with numbers
    print(f"S = (vertex permutations) * (line permutations)")
    print(f"S = {vertex_symmetry} * {line_symmetry}")
    print(f"S = {total_symmetry_factor}")

    print("\nTherefore, the final symmetry factor is:", total_symmetry_factor)

if __name__ == "__main__":
    calculate_symmetry_factor()