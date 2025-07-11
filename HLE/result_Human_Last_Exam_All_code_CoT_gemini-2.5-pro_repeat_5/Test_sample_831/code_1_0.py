import sys

def has_hamiltonian_path(shape_coords):
    """
    Checks if a polyomino, given as a set of coordinates, has a Hamiltonian path.
    A Hamiltonian path visits every square exactly once.
    """
    coords = list(shape_coords)
    n = len(coords)
    if n == 0:
        return True
    
    adj = {c: [] for c in coords}
    coord_set = set(coords)

    for i in range(n):
        x, y = coords[i]
        for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            neighbor = (x + dx, y + dy)
            if neighbor in coord_set:
                adj[coords[i]].append(neighbor)

    path = []
    visited = set()

    def find_path(u):
        """Recursively search for a path from node u."""
        path.append(u)
        visited.add(u)
        if len(path) == n:
            return True
        
        # Using a sorted list of neighbors to make the search deterministic
        for v in sorted(adj[u]):
            if v not in visited:
                if find_path(v):
                    return True
        
        path.pop()
        visited.remove(u)
        return False

    # A Hamiltonian path can start at any vertex.
    for start_node in sorted(coords):
        if find_path(start_node):
            return True
            
    return False

def get_canonical(shape_coords):
    """
    Computes the canonical representation of a polyomino.
    This is done by generating all 8 transformations (4 rotations, 4 reflections)
    and picking the lexicographically smallest one.
    This ensures that rotated and reflected versions of a shape are treated as identical.
    """
    initial_coords = sorted(list(shape_coords))
    
    forms = set()
    
    # Generate 8 symmetries
    current_coords = list(initial_coords)
    for _ in range(2):  # Original and reflected
        for _ in range(4):  # 4 rotations
            min_x = min(c[0] for c in current_coords)
            min_y = min(c[1] for c in current_coords)
            # Normalize to origin
            normalized = tuple(sorted([(c[0] - min_x, c[1] - min_y) for c in current_coords]))
            forms.add(normalized)
            # Rotate 90 degrees: (x, y) -> (-y, x)
            current_coords = [(-c[1], c[0]) for c in current_coords]
        # Reflect across y-axis: (x, y) -> (-x, y)
        current_coords = [(-c[0], c[1]) for c in initial_coords]

    return min(forms)

def solve_t4():
    """
    Calculates T(4) by generating all tetra-dominoes and checking for Hamiltonian paths.
    """
    # Start with a single horizontal domino (order 1)
    q = {frozenset([(0, 0), (1, 0)])}

    # Recursively add dominoes up to order 4
    for _ in range(3): # n=2, n=3, n=4
        next_q = set()
        for shape in q:
            # Find all empty squares adjacent to the current shape
            neighbors = set()
            for x, y in shape:
                for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                    neighbor_coord = (x + dx, y + dy)
                    if neighbor_coord not in shape:
                        neighbors.add(neighbor_coord)
            
            # Try to place a new domino at each neighboring square
            for nx, ny in neighbors:
                # Try horizontal domino
                d_h = frozenset([(nx, ny), (nx + 1, ny)])
                if not (d_h & shape):
                    next_q.add(shape | d_h)
                
                d_h_rev = frozenset([(nx - 1, ny), (nx, ny)])
                if not (d_h_rev & shape):
                     next_q.add(shape | d_h_rev)

                # Try vertical domino
                d_v = frozenset([(nx, ny), (nx, ny + 1)])
                if not (d_v & shape):
                    next_q.add(shape | d_v)

                d_v_rev = frozenset([(nx, ny - 1), (nx, ny)])
                if not (d_v_rev & shape):
                    next_q.add(shape | d_v_rev)
        q = next_q

    # Get the set of unique free tetra-dominoes
    free_shapes = {get_canonical(shape) for shape in q}
    total_unique_shapes = len(free_shapes)

    # Count how many of them have a Hamiltonian path
    traversable_count = 0
    for shape in free_shapes:
        if has_hamiltonian_path(shape):
            traversable_count += 1
            
    non_traversable_count = total_unique_shapes - traversable_count
    
    print(f"Total number of unique tetra-dominoes: {total_unique_shapes}")
    print(f"Number of non-traversable shapes: {non_traversable_count}")
    print(f"Number of traversable shapes (T(4)): {traversable_count}")
    print("\nFinal Equation:")
    print(f"{total_unique_shapes} - {non_traversable_count} = {traversable_count}")


if __name__ == '__main__':
    solve_t4()