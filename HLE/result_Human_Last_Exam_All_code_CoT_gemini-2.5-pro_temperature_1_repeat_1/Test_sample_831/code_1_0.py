import collections

def solve_t4():
    """
    This script calculates T(4), the number of non-equivalent, traversable
    polyforms of order 4.
    """

    def normalize(squares):
        """
        Normalizes a polyomino to its canonical representation by checking all 8 symmetries.
        """
        forms = []
        original_squares = list(squares)
        
        # Symmetries: identity, 3 rotations, and their reflections
        symmetries = [
            lambda x, y: (x, y),   # identity
            lambda x, y: (-y, x),  # rotate 90
            lambda x, y: (-x, -y), # rotate 180
            lambda x, y: (y, -x),  # rotate 270
            lambda x, y: (-x, y),  # reflect across y-axis
            lambda x, y: (y, x),   # reflect and rotate 90
            lambda x, y: (x, -y),  # reflect and rotate 180
            lambda x, y: (-y, -x)  # reflect and rotate 270
        ]

        for sym in symmetries:
            transformed = [sym(x, y) for x, y in original_squares]
            min_x = min(p[0] for p in transformed)
            min_y = min(p[1] for p in transformed)
            # Translate to origin, sort, and convert to tuple for consistent comparison
            forms.append(tuple(sorted([(x - min_x, y - min_y) for x, y in transformed])))

        # The canonical form is the lexicographically smallest one
        return frozenset(min(forms))

    def has_hamiltonian_path(squares):
        """
        Checks if a polyomino (represented as a set of squares) has a Hamiltonian path.
        """
        n = len(squares)
        if n == 0:
            return False
        
        adj = collections.defaultdict(list)
        square_list = list(squares)
        # Create a mapping from coordinate to index for faster adjacency list creation
        square_to_idx = {sq: i for i, sq in enumerate(square_list)}

        for i, (r, c) in enumerate(square_list):
            for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                neighbor = (r + dr, c + dc)
                if neighbor in square_to_idx:
                    adj[i].append(square_to_idx[neighbor])
        
        # Backtracking function to find a path
        def find_path_recursive(u, path):
            path.add(u)
            if len(path) == n:
                return True
            
            for v in adj[u]:
                if v not in path:
                    if find_path_recursive(v, path):
                        return True
            
            path.remove(u) # Backtrack
            return False

        # Try starting the path from each square (vertex)
        for i in range(n):
            if find_path_recursive(i, set()):
                return True
        return False

    # --- Main Logic ---

    # 1. Generate all unique polyforms of order 4 (octadominoes)
    # Start with a single domino (order 1)
    current_polys = {frozenset({(0, 0), (1, 0)})}

    for order in range(2, 5):
        next_polys = set()
        for poly in current_polys:
            # Find all empty squares adjacent to the current polyomino
            candidate_squares = set()
            for r, c in poly:
                for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                    neighbor = (r + dr, c + dc)
                    if neighbor not in poly:
                        candidate_squares.add(neighbor)
            
            # Try to place a new domino on each candidate square
            for r_cand, c_cand in candidate_squares:
                for dr, dc in [(0, 1), (1, 0)]: # Check only 2 directions to avoid duplicate dominoes
                    neighbor_of_cand = (r_cand + dr, c_cand + dc)
                    if neighbor_of_cand not in poly:
                        new_domino = frozenset({(r_cand, c_cand), neighbor_of_cand})
                        new_poly = poly.union(new_domino)
                        next_polys.add(normalize(new_poly))
        current_polys = next_polys
    
    order_4_polys = current_polys
    total_polyforms = len(order_4_polys)

    # 2. Count how many of them are traversable
    traversable_count = 0
    for poly in order_4_polys:
        if has_hamiltonian_path(poly):
            traversable_count += 1

    # 3. Print the results
    non_traversable_count = total_polyforms - traversable_count
    
    print(f"Total non-equivalent polyforms of order 4: {total_polyforms}")
    print(f"Number of traversable polyforms (T(4)): {traversable_count}")
    print(f"Number of non-traversable polyforms: {non_traversable_count}")
    print(f"The final equation: {traversable_count} (traversable) + {non_traversable_count} (non-traversable) = {total_polyforms} (total)")

solve_t4()