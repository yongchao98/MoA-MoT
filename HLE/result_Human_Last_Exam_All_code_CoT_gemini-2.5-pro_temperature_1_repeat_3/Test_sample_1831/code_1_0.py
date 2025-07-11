import math

def number_of_rooted_trees(n):
    """
    Returns the number of non-isomorphic rooted trees on n vertices.
    These are pre-calculated values from OEIS A000081.
    r(1)=1, r(2)=1, r(3)=2, r(4)=4.
    """
    if n == 1: return 1
    if n == 2: return 1
    if n == 3: return 2
    if n == 4: return 4
    return 0 # Only need up to n=4

def ways_to_attach_trees(num_vertices, num_cycle_nodes):
    """
    Calculates the number of ways to form trees on `num_vertices`
    and attach them to `num_cycle_nodes` distinct cycle vertices.
    This is a simplified direct enumeration for small numbers.
    """
    if num_vertices == 0:
        return 1
    if num_vertices == 1:
        # One vertex, can be attached to any cycle node. Symmetric.
        return 1
    if num_vertices == 2:
        # The 2 vertices can form a single tree of size 2, or two trees of size 1.
        # A single tree of size 2 (r(2)=1 type) attached to one cycle node.
        ways_one_tree = 1
        # Two trees of size 1. Distribute 2 identical items into num_cycle_nodes bins.
        # For c=2: (2,0), (1,1) -> 2 ways.
        # For c=1: (2) -> 1 way.
        ways_two_trees = num_cycle_nodes
        if num_cycle_nodes == 1:
             ways_two_trees = 1
        if num_cycle_nodes == 2:
             ways_two_trees = 2 # (2 on one node) or (1 on each node)
        return ways_one_tree + ways_two_trees
    return 0

def num_connected_components(m):
    """
    Calculates N(m), the number of non-isomorphic connected functional graphs on m vertices.
    A connected functional graph has exactly one cycle.
    """
    if m == 0: return 0
    total_ways = 0
    # Iterate through all possible cycle lengths c from 1 to m
    for c in range(1, m + 1):
        remaining_vertices = m - c
        if remaining_vertices == 0:
            # Just a simple cycle of length c
            ways = 1
        elif c == 1:
            # A single fixed point with trees attached. This is a rooted tree.
            ways = number_of_rooted_trees(m)
        else:
            # A cycle of length c with trees on `remaining_vertices` attached.
            ways = ways_to_attach_trees(remaining_vertices, c)
        
        # print(f"  For a component of size {m} with cycle length {c}, there are {ways} ways.")
        total_ways += ways
    return total_ways

def main():
    """
    Main function to calculate the number of equivalence classes.
    """
    print("Step 1: Calculate N(m), the number of non-isomorphic connected functional graphs on m vertices.")
    
    n1 = num_connected_components(1)
    print(f"N(1): A single vertex must be a fixed point (cycle of length 1). N(1) = {n1}")
    
    n2 = num_connected_components(2)
    print(f"N(2): A 2-cycle, or a rooted tree on 2 vertices. N(2) = 1 + 1 = {n2}")
    
    n3 = num_connected_components(3)
    print(f"N(3): A 3-cycle (1 way), a 2-cycle with one leaf (1 way), or a rooted tree on 3 vertices (2 ways). N(3) = 1 + 1 + 2 = {n3}")
    
    n4 = num_connected_components(4)
    print(f"N(4): A 4-cycle (1), 3-cycle+leaf (1), 2-cycle+trees (3), or rooted tree (4). N(4) = 1 + 1 + 3 + 4 = {n4}")

    print("\nStep 2: Calculate the total number of classes by considering partitions of 4.")
    
    # Partition 4: One component of size 4
    p4 = n4
    print(f"Partition '4': One component of size 4. This is N(4). Number of ways = {p4}")
    
    # Partition 3+1: One component of size 3 and one of size 1
    p31 = n3 * n1
    print(f"Partition '3+1': One component of size 3 (N(3) ways) and one of size 1 (N(1) ways). Number of ways = {n3} * {n1} = {p31}")
    
    # Partition 2+2: Two components of size 2.
    # N(2) = 2 types. We need to choose a multiset of size 2 from these 2 types.
    # (types + k - 1) choose k = (2 + 2 - 1) choose 2 = 3 choose 2 = 3
    p22 = 3
    print(f"Partition '2+2': Two components of size 2. There are N(2)={n2} types. Choosing a multiset of 2 gives 3 ways.")
    
    # Partition 2+1+1: One component of size 2, two of size 1
    p211 = n2
    print(f"Partition '2+1+1': One component of size 2 (N(2)={n2} ways) and two of size 1 (N(1)=1 way). Number of ways = {p211}")
    
    # Partition 1+1+1+1: Four components of size 1
    p1111 = 1
    print(f"Partition '1+1+1+1': Four components of size 1 (N(1)=1 way). Number of ways = {p1111}")
    
    total_classes = p4 + p31 + p22 + p211 + p1111
    
    print("\nStep 3: Sum the counts for all partitions.")
    print("The final equation is:")
    print(f"{p4} + {p31} + {p22} + {p211} + {p1111} = {total_classes}")
    
    print(f"\nThus, there are {total_classes} elements of E represented by F.")

if __name__ == "__main__":
    main()
