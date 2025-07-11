import math
from collections import Counter

def solve_functional_graphs(n=4):
    """
    Calculates the number of non-isomorphic functional graphs on n vertices.
    This corresponds to the number of equivalence classes of endomorphisms on a set of size n.
    """

    # t[k] = number of non-isomorphic rooted trees with k vertices. (OEIS A000081)
    t = {1: 1, 2: 1, 3: 2, 4: 4}

    # c[k] = number of non-isomorphic connected functional graphs with k vertices.
    c = {}

    # c[1]: a single fixed point
    c[1] = t[1]
    
    # c[2]: a 2-cycle or a t_2 rooted tree on a fixed point.
    c[2] = 1 + t[2]
    
    # c[3]: a 3-cycle, a 2-cycle with one vertex attached, or a t_3 tree on a fixed point.
    c[3] = 1 + 1 + t[3]

    # c[4]: a 4-cycle, 3-cycle with 1 vertex, 2-cycle with 2 vertices, or a t_4 tree on a fixed point.
    # For a 2-cycle, the 2 remaining vertices can be partitioned as {2} (one t_2 tree) or {1,1} (two t_1 trees).
    # The two t_1 trees can be placed both on one node or one on each. This gives 1+2=3 ways for the 2-cycle case.
    c[4] = 1 + 1 + 3 + t[4]
    
    print("This problem is equivalent to counting non-isomorphic functional graphs on 4 vertices.")
    print("We sum the possibilities over all integer partitions of 4.\n")
    
    def get_partitions(num):
        # Generates integer partitions of a number
        if num == 0:
            yield []
            return
        for i in range(1, num + 1):
            for p in get_partitions(num - i):
                if not p or i >= p[0]:
                    yield [i] + p
    
    partitions = sorted(list(get_partitions(n)), reverse=True)
    
    total_classes = 0
    final_values = []

    print(f"Number of connected component types for size k=1,2,3,4: c_k = {c[1]}, {c[2]}, {c[3]}, {c[4]}\n")
    print("Calculation for each partition of 4:")

    for p in partitions:
        counts = Counter(p)
        term_ways = 1
        
        # Using combinations with replacement for repeated partition parts
        for k, m in counts.items():
            term_ways *= math.comb(c[k] + m - 1, m)
            
        print(f"Partition {p}: Contributes {term_ways} classes.")
        total_classes += term_ways
        final_values.append(term_ways)
        
    print("\nThe total number of classes is the sum of these contributions.")
    print("Final Equation:")
    print(f"{' + '.join(map(str, final_values))} = {total_classes}")

solve_functional_graphs()
<<<19>>>