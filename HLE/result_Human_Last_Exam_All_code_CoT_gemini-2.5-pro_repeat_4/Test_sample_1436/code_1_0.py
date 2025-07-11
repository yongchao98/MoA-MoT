import math

def calculate_symmetry_factor():
    """
    Calculates the symmetry factor for a Feynman diagram with two identical
    vertices connected by four identical lines.
    """
    # Number of identical vertices in the diagram
    num_identical_vertices = 2
    
    # Number of identical lines (propagators) connecting the vertices
    num_identical_lines = 4

    # Calculate the symmetry contribution from permuting the vertices
    vertex_symmetry = math.factorial(num_identical_vertices)
    
    # Calculate the symmetry contribution from permuting the lines
    line_symmetry = math.factorial(num_identical_lines)
    
    # The total symmetry factor is the product of these contributions
    total_symmetry_factor = vertex_symmetry * line_symmetry

    print("This calculation is for a Feynman diagram with 2 identical vertices and 4 identical lines.")
    print("-" * 50)
    print(f"1. Symmetry from permuting the {num_identical_vertices} identical vertices: {num_identical_vertices}! = {vertex_symmetry}")
    print(f"2. Symmetry from permuting the {num_identical_lines} identical lines: {num_identical_lines}! = {line_symmetry}")
    print("-" * 50)
    print("The total symmetry factor (S) is the product of these values.")
    
    # Print the final equation with each number explicitly
    print(f"S = {num_identical_vertices}! x {num_identical_lines}! = {vertex_symmetry} x {line_symmetry} = {total_symmetry_factor}")

if __name__ == "__main__":
    calculate_symmetry_factor()