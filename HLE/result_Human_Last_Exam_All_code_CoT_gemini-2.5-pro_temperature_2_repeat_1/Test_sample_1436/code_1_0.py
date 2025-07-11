import math

def calculate_symmetry_factor():
    """
    Calculates the symmetry factor for a Feynman diagram with two vertices
    connected by four identical lines.
    """
    
    # The diagram has two identical vertices. The number of ways to permute them is 2!
    num_vertices = 2
    vertex_symmetry = math.factorial(num_vertices)
    
    # The diagram has four identical lines (propagators) connecting the two vertices.
    # The number of ways to permute these lines is 4!
    num_lines = 4
    line_symmetry = math.factorial(num_lines)
    
    # The total symmetry factor is the product of the vertex and line symmetries.
    total_symmetry_factor = vertex_symmetry * line_symmetry
    
    print("This Feynman diagram has two identical vertices and four identical lines connecting them.")
    print("The symmetry factor 'S' is the product of the number of permutations of these identical components.")
    print("\nCalculation steps:")
    print(f"1. Symmetry from interchanging the {num_vertices} vertices = {num_vertices}! = {vertex_symmetry}")
    print(f"2. Symmetry from interchanging the {num_lines} lines = {num_lines}! = {line_symmetry}")
    print("\nFinal Equation:")
    print(f"S = {num_vertices}! * {num_lines}! = {vertex_symmetry} * {line_symmetry} = {total_symmetry_factor}")

calculate_symmetry_factor()
<<<48>>>