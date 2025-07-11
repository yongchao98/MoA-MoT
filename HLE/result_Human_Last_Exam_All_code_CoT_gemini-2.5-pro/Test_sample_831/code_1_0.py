import collections

def get_shape_from_polybar(polybar_set):
    """Converts a set of horizontal bar positions to a set of cell coordinates."""
    shape_cells = set()
    for x, y in polybar_set:
        shape_cells.add((x, y))
        shape_cells.add((x + 1, y))
    return frozenset(shape_cells)

def get_canonical(shape_cells):
    """Finds the canonical representation of a shape among its 8 symmetries."""
    symmetries = []
    points = sorted(list(shape_cells))
    
    for _ in range(2):  # Normal and reflected
        for _ in range(4):  # 4 rotations
            min_x = min(p[0] for p in points)
            min_y = min(p[1] for p in points)
            symmetries.append(frozenset((p[0] - min_x, p[1] - min_y) for p in points))
            points = sorted([(-y, x) for x, y in points])
        points = sorted([(x, -y) for x, y in points])
        
    return min(symmetries)

def generate_free_polybar_shapes(n):
    """Generates all unique free polybar shapes of order n."""
    # Start with a single horizontal domino at the origin
    fixed_polybars = {frozenset({(0, 0)})}
    
    for _ in range(n - 1):
        next_fixed_polybars = set()
        for polybar in fixed_polybars:
            # For each bar in the current polybar, find adjacent spots for a new bar
            for x, y in polybar:
                # Potential neighbor bar positions
                neighbors = {
                    (x - 2, y), (x + 2, y),  # Horizontal neighbors
                    (x - 1, y - 1), (x, y - 1), (x + 1, y - 1),  # Neighbors below
                    (x - 1, y + 1), (x, y + 1), (x + 1, y + 1)   # Neighbors above
                }
                for neighbor_pos in neighbors:
                    if neighbor_pos not in polybar:
                        next_fixed_polybars.add(polybar.union({neighbor_pos}))
        fixed_polybars = next_fixed_polybars
        
    free_shapes = set()
    for polybar in fixed_polybars:
        shape = get_shape_from_polybar(polybar)
        free_shapes.add(get_canonical(shape))
        
    return free_shapes

def has_hamiltonian_path(shape_cells):
    """Checks if a shape (set of cells) has a Hamiltonian path."""
    if not shape_cells:
        return False
    
    n = len(shape_cells)
    nodes = list(shape_cells)
    adj = {node: [] for node in nodes}
    for x, y in nodes:
        for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            neighbor = (x + dx, y + dy)
            if neighbor in adj:
                adj[(x, y)].append(neighbor)

    path = []
    
    def solve(current_node):
        path.append(current_node)
        if len(path) == n:
            return True
        for neighbor in adj[current_node]:
            if neighbor not in path:
                if solve(neighbor):
                    return True
        path.pop()
        return False

    # Try starting the path from each node
    for start_node in nodes:
        if solve(start_node):
            return True
    
    return False

def main():
    """
    Calculates T(4) by generating polyforms and checking for Hamiltonian paths.
    """
    n = 4
    
    # Generate the unique shapes
    polybar_shapes = generate_free_polybar_shapes(n)
    total_shapes = len(polybar_shapes)
    
    # Count how many have a Hamiltonian path
    hp_count = 0
    for shape in polybar_shapes:
        if has_hamiltonian_path(shape):
            hp_count += 1
            
    print(f"For n = {n}:")
    print(f"Total number of unique polyforms with fixed domino orientation = {total_shapes}")
    print(f"Number of these polyforms with a Hamiltonian path = {hp_count}")
    print(f"Therefore, T(4) = {hp_count}")

if __name__ == '__main__':
    main()