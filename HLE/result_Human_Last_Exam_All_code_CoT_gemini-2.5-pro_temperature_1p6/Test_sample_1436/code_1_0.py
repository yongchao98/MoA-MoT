import math

def calculate_symmetry_factor():
    """
    Calculates the symmetry factor for a Feynman diagram with two vertices
    connected by four lines.
    """
    num_vertices = 2
    num_lines = 4

    # Calculate the number of permutations for vertices and lines
    vertex_permutations = math.factorial(num_vertices)
    line_permutations = math.factorial(num_lines)

    # The total symmetry factor is the product of these permutations
    symmetry_factor = vertex_permutations * line_permutations

    # Print the explanation and the final calculation
    print("The symmetry factor (S) for a Feynman diagram is the number of ways its components can be interchanged without changing the diagram's structure.")
    print("\nThis diagram consists of:")
    print(f"- {num_vertices} identical vertices")
    print(f"- {num_lines} identical lines (propagators) connecting the vertices")
    print("\nThe calculation involves multiplying the factorials of the counts of these identical components.")
    print("\nFinal Equation:")
    
    # The final equation showing all the numbers involved
    print(f"S = ({num_vertices}!) * ({num_lines}!) = {vertex_permutations} * {line_permutations} = {symmetry_factor}")


if __name__ == "__main__":
    calculate_symmetry_factor()