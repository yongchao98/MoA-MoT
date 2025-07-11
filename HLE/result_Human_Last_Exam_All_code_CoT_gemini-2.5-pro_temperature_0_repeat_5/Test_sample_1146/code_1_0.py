import math

def get_gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    while b:
        a, b = b, a % b
    return a

def find_ring_graph_n(limit):
    """
    Finds values of n up to a given limit for which AG(Z_n) is a ring graph.
    A graph is a ring graph if every vertex has a degree of 2.
    """
    solutions = []
    for n in range(2, limit + 1):
        # Vertices are non-zero elements of Z_n
        vertices = list(range(1, n))
        if not vertices:
            continue

        # Find the units of Z_n
        units = {u for u in vertices if get_gcd(u, n) == 1}

        is_ring_graph = True
        all_degrees = []

        # Calculate the degree of each vertex
        for v in vertices:
            # The associates of v form its equivalence class
            associates = {(v * u) % n for u in units}
            # The degree is the number of other vertices in its class
            degree = len(associates) - 1
            all_degrees.append(degree)

            # For a ring graph, every vertex must have degree 2
            if degree != 2:
                is_ring_graph = False
                # No need to check other vertices for this n
                break
        
        if is_ring_graph:
            solutions.append(n)

    return solutions

# Set a limit for checking
N_LIMIT = 100
solutions = find_ring_graph_n(N_LIMIT)

# As derived in the explanation, the set of solutions is empty.
# The code will confirm this for the checked range.
# The final output format requires printing the equation.
if not solutions:
    print("n in { }")
else:
    # This part of the code is unlikely to be executed.
    solution_str = ', '.join(map(str, solutions))
    print(f"n in {{ {solution_str} }}")
