import itertools

def solve_loopless_cycle_cover_parity(adj_matrix):
    """
    Solves the LooplessCycleCover problem by brute force and returns the parity.

    A cycle cover is a selection of arcs corresponding to a permutation pi of the
    vertices, where each arc is (i, pi(i)).
    A cycle cover is loopless if the permutation has no 2-cycles.
    The problem asks for the parity of the number of such loopless cycle covers.
    A multigraph is handled by the adjacency matrix entries, which count the
    number of parallel arcs.
    """
    n = len(adj_matrix)
    if n == 0:
        print("Graph has no vertices. Number of loopless cycle covers is 1 (the empty cover).")
        print("Parity: 1")
        return

    vertices = range(n)
    permutations = itertools.permutations(vertices)
    
    loopless_cycle_cover_count = 0

    # Explanation of the method
    print("Finding the parity of the number of loopless cycle covers.")
    print(f"The graph has {n} vertices.")
    print("We will iterate through all n! permutations of the vertices.")
    print("For each permutation pi, we check if it's a valid and loopless cycle cover.")
    print("-" * 20)

    for p_tuple in permutations:
        p = list(p_tuple) # permutation p where p[i] is the vertex that i maps to.
        
        # A self-loop (pi(i) = i) is not a 2-cycle, but problem states graph has no self-loops.
        # We handle this by ensuring adj_matrix[i][i] is 0.
        # However, a permutation can have fixed points, leading to a product of 0 if no self-loops exist.

        # Check for looplessness (no 2-cycles)
        is_loopless = True
        for i in range(n):
            if p[p[i]] == i and p[i] != i:
                is_loopless = False
                break
        
        if not is_loopless:
            continue

        # This permutation is loopless. Now check if it's a valid cycle cover
        # and calculate the number of ways it can be formed.
        # This is the product of the number of arcs (i, p[i]) for all i.
        num_ways = 1
        for i in range(n):
            num_ways *= adj_matrix[i][p[i]]
        
        if num_ways > 0:
            loopless_cycle_cover_count += num_ways

    parity = loopless_cycle_cover_count % 2
    
    print(f"Total number of loopless cycle covers found: {loopless_cycle_cover_count}")
    print(f"The parity is {loopless_cycle_cover_count} mod 2 = {parity}")


# Example from the problem description's context
# A directed multigraph G without self-loops.
# Let's define a simple 4-vertex graph.
# M = [[0, 1, 1, 0],
#      [1, 0, 0, 1],
#      [1, 0, 0, 1],
#      [0, 1, 1, 0]]
# We found two loopless cycle covers: (0->1->3->2->0) and (0->2->3->1->0)
# Let's use a slightly different graph to avoid identical rows.
# G has arcs: (0,1), (1,0), (0,2), (2,3), (3,1)
# And let's add (1,2) and (2,1) to create more possibilities
adj_matrix_example = [
    [0, 1, 1, 0], # Arcs from 0: (0,1), (0,2)
    [1, 0, 1, 0], # Arcs from 1: (1,0), (1,2)
    [0, 1, 0, 1], # Arcs from 2: (2,1), (2,3)
    [0, 1, 0, 0]  # Arcs from 3: (3,1)
]

solve_loopless_cycle_cover_parity(adj_matrix_example)
