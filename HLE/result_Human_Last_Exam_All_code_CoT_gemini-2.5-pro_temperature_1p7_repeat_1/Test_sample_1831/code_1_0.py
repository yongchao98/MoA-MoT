import math
from collections import Counter

def combinations(n, k):
    """Calculates n choose k."""
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

def get_integer_partitions(n, max_parts=None):
    """
    Generates all integer partitions of n.
    If max_parts is specified, it limits the number of parts in the partition.
    """
    partitions = set()
    stack = [(n, [])]
    while stack:
        current_sum, current_partition = stack.pop()
        if current_sum == 0:
            partitions.add(tuple(sorted(current_partition, reverse=True)))
            continue

        start = current_partition[-1] if current_partition else 1
        for i in range(start, current_sum + 1):
            new_partition = current_partition + [i]
            if max_parts is None or len(new_partition) <= max_parts:
                stack.append((current_sum - i, new_partition))
    return partitions

def calculate_c(k, T_values):
    """
    Calculates C(k), the number of non-isomorphic connected functional graphs on k vertices.
    C(k) = sum over cycle lengths c from 1 to k of {ways to attach k-c tree nodes to a c-cycle}.
    """
    if k == 0:
        return 1
    total_c = 0
    # Sum over possible cycle lengths 'c' from 1 to k
    for c in range(1, k + 1):
        num_tree_nodes = k - c
        
        # This sub-problem is to count the ways to form a forest of `num_tree_nodes`
        # total vertices, whose roots are the `c` (indistinguishable) cycle vertices.
        # This is equivalent to partitioning the integer `num_tree_nodes` into at most `c` parts.
        
        partitions = get_integer_partitions(num_tree_nodes, max_parts=c)
        
        ways_for_c = 0
        if num_tree_nodes == 0:
             ways_for_c = 1 # Just the bare c-cycle
        else:
            for p in partitions:
                # For each partition, count the distinct structures.
                # A partition like [s1, s2, s2] means we need to choose 1 tree of size s1+1
                # and 2 trees of size s2+1, from the available non-isomorphic rooted trees.
                counts = Counter(p)
                term = 1
                for size, mult in counts.items():
                    # The number of non-isomorphic rooted trees on `size+1` vertices is T[size+1]
                    # We are choosing `mult` trees from `T[size+1]` types, with replacement.
                    term *= combinations(T_values[size + 1] + mult - 1, mult)
                ways_for_c += term
                
        total_c += ways_for_c
    return total_c

def solve():
    """
    Main function to solve the problem for N=4.
    """
    N = 4
    # T(n): Number of non-isomorphic rooted trees on n vertices.
    # We only need up to N. (OEIS A000081)
    T = {1: 1, 2: 1, 3: 2, 4: 4, 5: 9}
    
    print(f"The problem is to find the number of equivalence classes of endomorphisms of a set of size {N}.")
    print("This is equivalent to counting the number of non-isomorphic functional graphs on 4 vertices, let's call it F(4).\n")
    print("First, we calculate C(k), the number of *connected* non-isomorphic functional graphs on k vertices.")

    C = {}
    for k in range(1, N + 1):
        C[k] = calculate_c(k, T)
        print(f"C({k}) = {C[k]}")

    print("\nNow, we use these values to find F(4). A graph on 4 vertices can be decomposed")
    print("into connected components. We sum the possibilities over all integer partitions of 4.\n")
    
    # Calculate F(4) using partitions of 4
    total_f = 0
    
    # Partition 4
    p1 = C[4]
    total_f += p1
    print(f"For partition 4 (one component of size 4): ways = C(4) = {p1}")

    # Partition 3+1
    p2 = C[3] * C[1]
    total_f += p2
    print(f"For partition 3+1 (one of size 3, one of size 1): ways = C(3) * C(1) = {C[3]} * {C[1]} = {p2}")
    
    # Partition 2+2
    p3 = combinations(C[2] + 2 - 1, 2)
    total_f += p3
    print(f"For partition 2+2 (two of size 2): ways = (C(2) + 2 - 1 choose 2) = ({C[2]} + 1 choose 2) = {p3}")
    
    # Partition 2+1+1
    p4 = C[2] * combinations(C[1] + 2 - 1, 2)
    total_f += p4
    print(f"For partition 2+1+1 (one of size 2, two of size 1): ways = C(2) * (C(1) + 2 - 1 choose 2) = {C[2]} * ({C[1]} + 1 choose 2) = {p4}")
    
    # Partition 1+1+1+1
    p5 = combinations(C[1] + 4 - 1, 4)
    total_f += p5
    print(f"For partition 1+1+1+1 (four of size 1): ways = (C(1) + 4 - 1 choose 4) = ({C[1]} + 3 choose 4) = {p5}")
    
    print(f"\nThe total number of equivalence classes is the sum: {p1} + {p2} + {p3} + {p4} + {p5} = {total_f}")

solve()
<<<19>>>