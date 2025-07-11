import math

def calculate_symmetry_factor():
    """
    Calculates the symmetry factor for a Feynman diagram with two vertices
    connected by four internal lines.
    """
    # Number of identical vertices that can be interchanged
    num_vertices = 2
    
    # Number of identical lines connecting the vertices
    num_lines = 4
    
    # The symmetry factor arises from permutations of identical components.
    # 1. Permutations of the two identical vertices.
    vertex_permutations = math.factorial(num_vertices)
    
    # 2. Permutations of the four identical lines connecting the vertices.
    line_permutations = math.factorial(num_lines)
    
    # The total symmetry factor is the product of these individual factors.
    symmetry_factor = vertex_permutations * line_permutations
    
    print("Calculating the symmetry factor (S) for a diagram with 2 vertices and 4 internal lines.")
    print("The formula for the symmetry factor is S = (permutations of identical vertices) * (permutations of identical lines).")
    print(f"Number of identical vertices = {num_vertices}")
    print(f"Number of identical lines = {num_lines}")
    print(f"S = {num_vertices}! * {num_lines}!")
    print("The final equation is:")
    print(f"S = {vertex_permutations} * {line_permutations} = {symmetry_factor}")

if __name__ == "__main__":
    calculate_symmetry_factor()
    # The final answer is the result of the calculation.
    # For the given diagram, S = 2! * 4! = 2 * 24 = 48.
    final_answer = 48
    # The following line is for automated result extraction.
    # print(f'<<<{final_answer}>>>')