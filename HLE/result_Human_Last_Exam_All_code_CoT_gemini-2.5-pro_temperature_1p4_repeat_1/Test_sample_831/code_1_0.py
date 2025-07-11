import sys

def solve_T4():
    """
    This script calculates T(4), the number of non-equivalent polyforms of order 4
    (4-dominominoes) that have a Hamiltonian path.
    """

    # It's good practice to increase the recursion limit for deep recursive calls,
    # though it may not be strictly necessary for this problem's scale.
    sys.setrecursionlimit(2000)

    def normalize(squares):
        """Translates a polyomino to a canonical position (min_x=0, min_y=0)."""
        if not squares:
            return []
        min_x = min(s[0] for s in squares)
        min_y = min(s[1] for s in squares)
        # Return a sorted list for consistent comparison
        return sorted([(s[0] - min_x, s[1] - min_y) for s in squares])

    def get_symmetries(squares):
        """Generates all 8 symmetries (rotations and reflections) of a polyomino."""
        # A frozenset of coordinates (x, y)
        s_coords = set(squares)
        
        # Transformations
        symmetries_set = [
            s_coords,                                  # Identity
            {(-y, x) for x, y in s_coords},            # Rot 90
            {(-x, -y) for x, y in s_coords},           # Rot 180
            {(y, -x) for x, y in s_coords},            # Rot 270
            {(-x, y) for x, y in s_coords},            # Flip horizontal
            {(y, x) for x, y in s_coords},             # Flip -> Rot 90
            {(x, -y) for x, y in s_coords},            # Flip -> Rot 180
            {(-y, -x) for x, y in s_coords}            # Flip -> Rot 270
        ]
        
        return [normalize(s) for s in symmetries_set]

    def get_canonical_form(squares):
        """Finds the canonical representation, the lexicographically smallest symmetry."""
        return frozenset(min(get_symmetries(squares)))

    def has_hamiltonian_path(squares):
        """Checks if a polyomino's grid graph has a Hamiltonian path using backtracking."""
        n_squares = len(squares)
        if n_squares == 0:
            return True
        
        adj = {sq: [] for sq in squares}
        square_list = list(squares)

        for i in range(n_squares):
            for j in range(i + 1, n_squares):
                s1 = square_list[i]
                s2 = square_list[j]
                if abs(s1[0] - s2[0]) + abs(s1[1] - s2[1]) == 1:
                    adj[s1].append(s2)
                    adj[s2].append(s1)

        # visited must be reset for each starting node attempt
        def find_path_from(u, visited):
            visited.add(u)
            if len(visited) == n_squares:
                return True

            for v in adj[u]:
                if v not in visited:
                    if find_path_from(v, visited):
                        return True
            
            # Backtrack
            visited.remove(u)
            return False

        # Try starting the path from every square in the polyomino
        for start_node in squares:
            if find_path_from(start_node, set()):
                return True
                
        return False

    # --- Main Logic ---

    # Start with the single domino shape (order 1)
    current_polys = {frozenset([(0, 0), (0, 1)])}

    # Generate polyforms of order 2, 3, and finally 4
    for n in range(2, 5):
        next_polys_canonical = set()
        for poly in current_polys:
            # Find all empty squares adjacent to the current polyomino
            adj_empty = set()
            for x, y in poly:
                for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                    neighbor = (x + dx, y + dy)
                    if neighbor not in poly:
                        adj_empty.add(neighbor)
            
            # Try to add a new domino using one of the adjacent empty squares
            for s1 in adj_empty:
                x1, y1 = s1
                # Try to form a domino with one of s1's neighbors
                for dx, dy in [(0, 1), (1, 0)]: # Check only two directions to avoid duplicates
                    s2 = (x1 + dx, y1 + dy)
                    if s2 not in poly:
                        new_domino = {s1, s2}
                        new_poly_squares = poly.union(new_domino)
                        canonical_form = get_canonical_form(new_poly_squares)
                        next_polys_canonical.add(canonical_form)
        
        current_polys = next_polys_canonical

    dominominoes_order_4 = current_polys
    
    hamiltonian_count = 0
    for shape in dominominoes_order_4:
        if has_hamiltonian_path(shape):
            hamiltonian_count += 1
    
    # Per the instructions to "output each number in the final equation"
    print(f"Total unique 4-dominominoes = {len(dominominoes_order_4)}")
    print(f"Number of traversable 4-dominominoes = {hamiltonian_count}")
    print(f"T(4) = {hamiltonian_count}")

solve_T4()