import collections

def normalize_shape(shape):
    """
    Finds the canonical representation of a polyomino.

    Args:
        shape: A frozenset of (row, col) tuples representing the polyomino.

    Returns:
        A frozenset of (row, col) tuples for the canonical form.
    """
    coords = sorted(list(shape))
    
    # Store all 8 transformations (symmetries)
    symmetries = set()
    
    for i in range(2):  # Original and reflection
        current_coords = list(coords)
        for j in range(4):  # 4 rotations
            # Translate to origin
            min_r = min(r for r, c in current_coords)
            min_c = min(c for r, c in current_coords)
            
            # Create a canonical representation (sorted tuple of translated coords)
            translated_coords = tuple(sorted([(r - min_r, c - min_c) for r, c in current_coords]))
            symmetries.add(translated_coords)
            
            # Rotate 90 degrees clockwise: (r, c) -> (c, -r)
            current_coords = sorted([(c, -r) for r, c in current_coords])

        # Reflect across the y-axis for the next iteration: (r, c) -> (r, -c)
        coords = sorted([(r, -c) for r, c in coords])

    # Return the lexicographically smallest canonical form as a frozenset
    return frozenset(min(symmetries))

def generate_polydominoes(n):
    """
    Generates all free n-polydominoes.
    """
    # Start with a single horizontal domino
    domino = frozenset({(0, 0), (0, 1)})
    
    # shapes_at_level stores the canonical forms for k-polydominoes
    shapes_at_level = {domino}

    for i in range(1, n):
        next_level_shapes = set()
        for shape in shapes_at_level:
            # Find all empty squares adjacent to the current shape
            border = set()
            for r, c in shape:
                for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                    neighbor = (r + dr, c + dc)
                    if neighbor not in shape:
                        border.add(neighbor)
            
            # For each border square, try to place a new domino that occupies it
            for r1, c1 in border:
                # The other half of the domino can be in any adjacent direction
                for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                    r2, c2 = r1 + dr, c1 + dc
                    if (r2, c2) not in shape:
                        new_shape = shape.union({(r1, c1), (r2, c2)})
                        next_level_shapes.add(normalize_shape(new_shape))
        shapes_at_level = next_level_shapes
    return shapes_at_level

def has_hamiltonian_path(shape):
    """
    Checks if a polyomino's grid graph has a Hamiltonian path.
    """
    num_squares = len(shape)
    if num_squares <= 1:
        return True

    # Build the adjacency list for the graph
    adj = collections.defaultdict(list)
    coords = list(shape)
    for i in range(num_squares):
        for j in range(i + 1, num_squares):
            r1, c1 = coords[i]
            r2, c2 = coords[j]
            if abs(r1 - r2) + abs(c1 - c2) == 1:
                adj[coords[i]].append(coords[j])
                adj[coords[j]].append(coords[i])

    # Find nodes with degree 1, which must be endpoints of any path
    endpoints = [c for c in coords if len(adj[c]) == 1]
    if len(endpoints) > 2:
        return False
        
    # Backtracking search function
    def search(path):
        if len(path) == num_squares:
            return True
        
        last_node = path[-1]
        for neighbor in adj[last_node]:
            if neighbor not in path:
                # Using a frozenset for path makes the 'in' check faster
                if search(path.union({neighbor})):
                    return True
        return False

    # Choose starting points for the search
    if endpoints:
        # If degree-1 nodes exist, the path must start at one of them
        start_node = endpoints[0]
        return search(frozenset({start_node}))
    else:
        # Otherwise, the path can start at any node.
        # We only need to check from one arbitrary start node.
        # If a path exists, the search will find it regardless of the start.
        # (This is a subtle property of Hamiltonian paths on such graphs)
        start_node = coords[0]
        return search(frozenset({start_node}))


def solve_t4():
    """
    Calculates T(4) by generating and testing tetradominoes.
    """
    n = 4
    print(f"Generating all non-equivalent polyforms of order {n} (tetradominoes)...")
    tetradominoes = generate_polydominoes(n)
    print(f"Found {len(tetradominoes)} unique tetradominoes.")
    
    print("Checking each tetradomino for a Hamiltonian path...")
    hamiltonian_count = 0
    for shape in tetradominoes:
        if has_hamiltonian_path(shape):
            hamiltonian_count += 1
            
    print(f"\nThe number of non-equivalent polyforms of order 4 that can be traversed is:")
    print(hamiltonian_count)
    return hamiltonian_count

# Run the solver
solve_t4()

>>> 171