import collections

def get_canonical(poly):
    """
    Calculates the canonical representation of a polyomino.

    A polyomino is represented as a frozenset of (x, y) coordinates.
    This function generates all 8 transformations (4 rotations, 4 reflections)
    of the polyomino, translates each to the origin, and returns the one
    that is lexicographically smallest.
    """
    # frozenset of tuples is hashable and comparable
    # Python's comparison for sets/frozensets works by first comparing sizes,
    # then elements. For elements (tuples), it's lexicographical.
    # We sort the tuples to ensure consistent comparison.
    
    # Generate 8 symmetries (rotations and reflections)
    symmetries = []
    current_poly = list(poly) # Use a list for mutable operations
    
    for _ in range(2): # Original and reflection
        for _ in range(4): # 4 rotations
            # Translate to origin
            min_x = min(p[0] for p in current_poly)
            min_y = min(p[1] for p in current_poly)
            normalized_poly = frozenset(sorted([(p[0] - min_x, p[1] - min_y) for p in current_poly]))
            symmetries.append(normalized_poly)
            
            # Rotate 90 degrees: (x, y) -> (-y, x)
            current_poly = [(-y, x) for x, y in current_poly]
            
        # Reflect across y-axis: (x, y) -> (-x, y)
        current_poly = [(-x, y) for x, y in current_poly]
        
    return min(symmetries)

def generate_free_polyominoes(n):
    """
    Generates all free polyominoes of a given size n.
    """
    if n == 1:
        return {frozenset([(0, 0)])}

    # Start with smaller polyominoes and add blocks to them
    smaller_polys = generate_free_polyominoes(n - 1)
    new_polys = set()
    
    for poly in smaller_polys:
        # For each polyomino, find all adjacent empty squares
        neighbors = set()
        for x, y in poly:
            for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                neighbor = (x + dx, y + dy)
                if neighbor not in poly:
                    neighbors.add(neighbor)
        
        # Create new polyominoes by adding each neighbor
        for neighbor in neighbors:
            new_poly_raw = poly.union({neighbor})
            canonical_form = get_canonical(new_poly_raw)
            new_polys.add(canonical_form)
            
    return new_polys

def check_domino_tiling(poly):
    """
    Checks if a polyomino can be tiled by 1x2 dominoes.
    This is true if it has an equal number of black and white squares
    on a chessboard. This condition is sufficient for polyominoes of size 8.
    """
    white_squares = 0
    black_squares = 0
    for x, y in poly:
        if (x + y) % 2 == 0:
            white_squares += 1
        else:
            black_squares += 1
    return white_squares == black_squares

def check_hamiltonian_path(poly):
    """
    Checks if a polyomino contains a Hamiltonian path.
    """
    nodes = list(poly)
    num_nodes = len(nodes)
    if num_nodes <= 1:
        return True

    # Build adjacency list for the graph
    adj = collections.defaultdict(list)
    node_map = {node: i for i, node in enumerate(nodes)}
    
    for i in range(num_nodes):
        for j in range(i + 1, num_nodes):
            x1, y1 = nodes[i]
            x2, y2 = nodes[j]
            if abs(x1 - x2) + abs(y1 - y2) == 1:
                adj[nodes[i]].append(nodes[j])
                adj[nodes[j]].append(nodes[i])

    # Backtracking search for a Hamiltonian path
    def find_path_from(start_node):
        path = [start_node]
        visited = {start_node}
        
        def search():
            if len(path) == num_nodes:
                return True
            
            last_node = path[-1]
            for neighbor in adj[last_node]:
                if neighbor not in visited:
                    visited.add(neighbor)
                    path.append(neighbor)
                    if search():
                        return True
                    # Backtrack
                    path.pop()
                    visited.remove(neighbor)
            return False

        return search()

    # Try starting the path from each node
    for start_node in nodes:
        if find_path_from(start_node):
            return True
            
    return False

def solve_t4():
    """
    Calculates T(4) by generating, filtering, and counting polyforms.
    """
    # For T(4), we need polyforms of order 4, which are octominoes (size 8).
    n = 8
    
    print(f"Step 1: Generating all unique free polyominoes of size {n} (octominoes)...")
    all_octominoes = generate_free_polyominoes(n)
    print(f"Found {len(all_octominoes)} unique free octominoes.")
    
    print("\nStep 2: Filtering for octominoes that can be tiled by dominoes...")
    tileable_octominoes = [p for p in all_octominoes if check_domino_tiling(p)]
    print(f"Found {len(tileable_octominoes)} octominoes that can be tiled.")

    print("\nStep 3: Filtering for tileable octominoes that have a Hamiltonian path...")
    count = 0
    final_polyforms = []
    for i, poly in enumerate(tileable_octominoes):
        if check_hamiltonian_path(poly):
            count += 1
            final_polyforms.append(poly)
    
    print(f"\nCalculation complete.")
    print(f"The number of non-equivalent polyforms of order 4, T(4), is: {count}")

if __name__ == '__main__':
    solve_t4()