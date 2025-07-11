import sys

def solve_T4():
    """
    This function calculates T(4), the number of non-equivalent polyforms of order 4 (4 dominoes)
    that can be traversed by visiting each square exactly once.
    """
    # It's good practice to increase the recursion limit for deep search problems,
    # although it may not be strictly necessary for this small problem size.
    sys.setrecursionlimit(2000)

    # --- Helper Functions ---

    def get_neighbors(square):
        """Returns the four grid neighbors of a square."""
        x, y = square
        return [(x + 1, y), (x - 1, y), (x, y + 1), (x, y - 1)]

    def normalize(shape_coords):
        """
        Translates a shape to the origin and sorts its coordinates to create a consistent representation.
        Returns a tuple of tuples.
        """
        if not shape_coords:
            return tuple()
        min_x = min(s[0] for s in shape_coords)
        min_y = min(s[1] for s in shape_coords)
        return tuple(sorted([(x - min_x, y - min_y) for x, y in shape_coords]))

    def get_canonical(shape):
        """
        Finds the canonical representation of a shape by checking all 8 symmetries.
        The canonical form is the lexicographically smallest of all normalized symmetric forms.
        """
        symmetries = set()
        s_coords = list(shape)
        for _ in range(4):  # Apply 4 rotations
            # Add normalized form of the current orientation
            symmetries.add(normalize(s_coords))
            # Add normalized form of its reflection
            symmetries.add(normalize([(x, -y) for x, y in s_coords]))
            # Rotate the coordinates by 90 degrees
            s_coords = [(-y, x) for x, y in s_coords]
        return min(symmetries)

    def has_hamiltonian_path(shape_coords):
        """
        Checks if a shape's grid graph has a Hamiltonian path using a backtracking search.
        """
        num_squares = len(shape_coords)
        if num_squares == 0:
            return False
        
        # Build the adjacency list for the graph
        adj = {sq: [] for sq in shape_coords}
        for sq in shape_coords:
            for neighbor in get_neighbors(sq):
                if neighbor in shape_coords:
                    adj[sq].append(neighbor)

        def do_search(u, path_set):
            """Recursive search function."""
            path_set.add(u)
            if len(path_set) == num_squares:
                return True
            for v in adj[u]:
                if v not in path_set:
                    if do_search(v, path_set):
                        return True
            path_set.remove(u)  # Backtrack
            return False

        # A Hamiltonian path must start somewhere. We try starting from every square.
        for start_node in shape_coords:
            # For each starting attempt, we use a new, empty set for the path.
            if do_search(start_node, set()):
                return True
        return False

    # --- Main Logic ---

    n = 4
    num_squares = 2 * n

    # Step 1: Generate all "fixed" polyforms of order 4
    # We use a breadth-first approach, starting with one domino.
    # `q` holds the shapes to be expanded at the current level.
    q = {frozenset([(0, 0), (0, 1)])}
    
    for _ in range(n - 1):
        next_q = set()
        for shape in q:
            # Find perimeter squares to add the next domino
            perimeter = set()
            for square in shape:
                for neighbor in get_neighbors(square):
                    if neighbor not in shape:
                        perimeter.add(neighbor)
            
            # Try adding a new domino at each possible perimeter position
            for p_square in perimeter:
                for neighbor in get_neighbors(p_square):
                    if neighbor not in shape:
                        new_domino = frozenset([p_square, neighbor])
                        new_shape = shape.union(new_domino)
                        next_q.add(new_shape)
        q = next_q

    fixed_polyforms = q
    
    # Step 2 & 3: Check traversability and find unique canonical forms
    traversable_canonical_shapes = set()
    for shape in fixed_polyforms:
        if has_hamiltonian_path(shape):
            canonical_form = get_canonical(shape)
            traversable_canonical_shapes.add(canonical_form)
            
    # Step 4: Final Count
    result = len(traversable_canonical_shapes)
    
    # As requested, outputting the numbers in the final equation
    print(f"The number of squares per polyform is 2 * 4 = {num_squares}.")
    print(f"The number of non-equivalent, traversable polyforms of order 4 is T(4).")
    print(f"T(4) = {result}")

solve_T4()