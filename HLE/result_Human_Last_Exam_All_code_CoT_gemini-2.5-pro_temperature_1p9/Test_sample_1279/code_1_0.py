import math

def calculate_symmetry_factors():
    """
    Calculates and sums the symmetry factors for the two second-order
    vacuum bubble diagrams in phi^4 scalar field theory.
    """
    # Introduction
    print("In phi^4 theory, there are two connected second-order vacuum bubble diagrams.")
    print("-" * 30)

    # Diagram 1: Sunset diagram
    print("1. The 'Sunset' diagram:")
    print("   This diagram has two vertices connected by four propagators.")
    
    # Calculation for S1
    s1_vertex_perms = 2
    s1_line_perms = math.factorial(4)
    s1 = s1_vertex_perms * s1_line_perms
    
    print(f"   Symmetry factor S1 = (vertex permutations) * (line permutations)")
    print(f"   S1 = {s1_vertex_perms} * {s1_line_perms} = {s1}")
    print("-" * 30)
    
    # Diagram 2: Figure-Eight diagram
    print("2. The 'Figure-Eight' diagram:")
    print("   This diagram has two vertices, each with one loop, connected by two propagators.")
    
    # Calculation for S2
    s2_vertex_perms = 2
    s2_line_perms = 2
    s2_loop1_flip = 2
    s2_loop2_flip = 2
    s2 = s2_vertex_perms * s2_line_perms * s2_loop1_flip * s2_loop2_flip
    
    print(f"   Symmetry factor S2 = (vertex perms) * (line perms) * (loop1 flip) * (loop2 flip)")
    print(f"   S2 = {s2_vertex_perms} * {s2_line_perms} * {s2_loop1_flip} * {s2_loop2_flip} = {s2}")
    print("-" * 30)

    # Summing the factors
    total_sum = s1 + s2
    print("The sum of all symmetry factors is S1 + S2.")
    print(f"The final equation is: {s1} + {s2} = {total_sum}")
    print(f"Total Sum = {total_sum}")

if __name__ == "__main__":
    calculate_symmetry_factors()
    # The final answer for the sum is 64.
    # To conform to the requested output format:
    # print("<<<64>>>")