import collections
from itertools import combinations

def get_neighbors(n, r, c):
    """Gets the valid neighbors of a vertex (r, c) in an n x n grid."""
    for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
        nr, nc = r + dr, c + dc
        if 0 <= nr < n and 0 <= nc < n:
            yield (nr, nc)

def is_connected(n, vertices):
    """Checks if a set of vertices induces a connected subgraph in an n x n grid."""
    if not vertices:
        return True
    q = collections.deque([next(iter(vertices))])
    seen = {next(iter(vertices))}
    while q:
        r, c = q.popleft()
        for neighbor in get_neighbors(n, r, c):
            if neighbor in vertices and neighbor not in seen:
                seen.add(neighbor)
                q.append(neighbor)
    return len(seen) == len(vertices)

def solve():
    """
    Calculates P_n for n=3, the number of distinct partitions of the vertices
    of an n x n grid graph into 3 sets of equal size, each inducing
    a connected subgraph.
    """
    n = 3
    if (n * n) % 3 != 0:
        print(f"P_{n} = 0, because n*n is not divisible by 3.")
        return

    partition_size = (n * n) // 3
    all_vertices = frozenset((r, c) for r in range(n) for c in range(n))
    
    ordered_partition_count = 0
    
    # Iterate through all possible sets for the first partition A
    for a_vertices in combinations(all_vertices, partition_size):
        a_set = frozenset(a_vertices)
        if not is_connected(n, a_set):
            continue
            
        remaining_vertices = all_vertices - a_set
        
        # Iterate through all possible sets for the second partition B
        for b_vertices in combinations(remaining_vertices, partition_size):
            b_set = frozenset(b_vertices)
            if not is_connected(n, b_set):
                continue

            # The third partition C is what's left
            c_set = remaining_vertices - b_set
            if not is_connected(n, c_set):
                continue
            
            # Found a valid ordered partition (A, B, C)
            ordered_partition_count += 1
            
    # The number of distinct (unordered) partitions is the ordered count / 3!
    # since the sets A, B, C are indistinguishable.
    if ordered_partition_count % 6 != 0:
         print("Error: Ordered count should be divisible by 6.")
         p_n = -1 # Should not happen
    else:
        p_n = ordered_partition_count // 6

    print(f"For n = {n}:")
    print(f"The number of vertices is {n*n}.")
    print(f"Each of the 3 partitions must have {partition_size} vertices.")
    print(f"The number of ways to form valid ordered partitions (A, B, C) is {ordered_partition_count}.")
    print(f"The number of distinct unordered partitions P_{n} is {ordered_partition_count} / 6 = {p_n}.")

    if p_n % 2 == 0:
        print(f"\nP_{n} = {p_n} is an even number.")
    else:
        print(f"\nP_{n} = {p_n} is an odd number.")

solve()