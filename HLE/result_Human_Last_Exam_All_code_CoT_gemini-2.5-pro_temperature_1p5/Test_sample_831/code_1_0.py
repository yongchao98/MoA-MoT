import collections

def solve():
    """
    Calculates T(4), the number of non-equivalent traversable polyforms of order 4.
    This corresponds to the number of free octominoes that have a Hamiltonian path.
    """
    N = 8  # For T(4), we consider polyominoes of size 2*4=8 (octominoes).

    def get_symmetries(polyomino):
        """Generates all 8 symmetries (rotations and reflections) of a polyomino."""
        symmetries = []
        # A polyomino is represented as a set of (x, y) coordinate tuples.
        current_poly = list(polyomino)
        
        for _ in range(4):
            # Rotation by 90 degrees: (x, y) -> (y, -x)
            current_poly = sorted([(y, -x) for x, y in current_poly])
            symmetries.append(tuple(current_poly))
            
            # Reflection across y-axis: (x, y) -> (-x, y)
            reflected_poly = sorted([(-x, y) for x, y in current_poly])
            symmetries.append(tuple(reflected_poly))
            
        return symmetries

    def get_canonical(polyomino):
        """Finds the canonical representation of a polyomino."""
        symmetries = get_symmetries(polyomino)
        canonical_forms = []
        for sym in symmetries:
            min_x = min(p[0] for p in sym)
            min_y = min(p[1] for p in sym)
            # Translate to origin to make it comparable
            normalized = frozenset(sorted((x - min_x, y - min_y) for x, y in sym))
            canonical_forms.append(normalized)
        
        # The canonical form is the lexicographically smallest one.
        return min(canonical_forms)

    def has_hamiltonian_path(polyomino):
        """Checks if the polyomino's graph has a Hamiltonian path."""
        n = len(polyomino)
        if n <= 1:
            return True
        
        points = list(polyomino)
        point_to_idx = {p: i for i, p in enumerate(points)}
        adj = collections.defaultdict(list)
        for i, (x, y) in enumerate(points):
            for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                neighbor_coord = (x + dx, y + dy)
                if neighbor_coord in point_to_idx:
                    adj[i].append(point_to_idx[neighbor_coord])

        # Backtracking search for a Hamiltonian path
        for start_node in range(n):
            visited = [False] * n
            
            def find_path(u, count):
                visited[u] = True
                if count == n:
                    return True
                
                for v in adj[u]:
                    if not visited[v]:
                        if find_path(v, count + 1):
                            return True
                
                visited[u] = False # Backtrack
                return False

            if find_path(start_node, 1):
                return True
        return False

    # 1. Generate all fixed octominoes (N=8)
    # We generate polyominoes level by level by size.
    polys = {frozenset([(0, 0)])}
    for _ in range(N - 1):
        new_polys = set()
        for poly in polys:
            for x, y in poly:
                for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                    neighbor = (x + dx, y + dy)
                    if neighbor not in poly:
                        new_poly = poly.union({neighbor})
                        new_polys.add(new_poly)
        polys = new_polys
    
    fixed_octominoes = polys
    
    # 2. Filter, Canonicalize, and Count
    unique_traversable_polys = set()
    for poly in fixed_octominoes:
        if has_hamiltonian_path(poly):
            canonical_form = get_canonical(poly)
            unique_traversable_polys.add(canonical_form)
    
    result = len(unique_traversable_polys)
    
    # Final Output
    print(f"The number of non-equivalent traversable polyforms of order 4 is T(4).")
    print(f"This is equivalent to the number of free octominoes with a Hamiltonian path.")
    print(f"Calculation: T(4) = {result}")

solve()