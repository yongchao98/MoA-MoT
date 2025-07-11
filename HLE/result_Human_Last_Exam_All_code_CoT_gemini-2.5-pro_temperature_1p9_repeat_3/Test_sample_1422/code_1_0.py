import math

def solve_min_vertices():
    """
    Calculates the minimum number of vertices for a two-loop Feynman diagram
    in an interacting scalar field theory (phi^4).
    """

    # Parameters for the problem
    num_loops = 2           # L, given as "two-loop"
    lines_per_vertex = 4    # n, for a standard phi^4 interaction

    print("To find the minimum number of vertices (V), we use two standard relations for connected Feynman diagrams:")
    print("1. The loop formula: L = I - V + 1")
    print(f"2. The vertex-line rule for phi^4 theory: n*V = 2*I + E, where n = {lines_per_vertex}")
    print("   L=loops, V=vertices, I=internal lines, E=external lines\n")

    print("We can solve these equations for V in terms of L and E.")
    print("From (1), we get: I = L + V - 1")
    print("Substituting I into (2): 4V = 2*(L + V - 1) + E")
    print("Solving for V gives the general formula: V = L - 1 + E/2\n")

    print(f"For a two-loop diagram, we are given L = {num_loops}.")
    print("Substituting L=2 into the formula for V:")
    # Using float for division in the print statement
    equation_str = f"V = {num_loops} - 1 + E/2  =>  V = 1 + E/2\n"
    print(equation_str)

    print("To minimize V, we must find the minimum valid number of external lines E.")
    print("E must be a non-negative even integer. The smallest possible value is E = 0 (for a vacuum diagram).")
    min_external_lines = 0
    print(f"Using E = {min_external_lines}.\n")

    # Calculate the minimum number of vertices using the final derived formula
    min_vertices = (num_loops - 1) + (min_external_lines / 2)

    print("Plugging the values into the equation:")
    print(f"Minimum V = {num_loops} - 1 + {min_external_lines}/2")
    final_calc_str = f"Minimum V = 1 + {min_external_lines / 2}"
    print(final_calc_str)
    final_result_str = f"Minimum V = {int(min_vertices)}"
    print(final_result_str + "\n")

    print(f"So, the minimum number of vertices in a two-loop diagram is {int(min_vertices)}.")

solve_min_vertices()