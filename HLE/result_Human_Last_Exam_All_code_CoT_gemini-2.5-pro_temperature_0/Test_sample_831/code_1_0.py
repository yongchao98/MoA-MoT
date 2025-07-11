import collections

def get_transformations(poly):
    """
    Generates all 8 symmetries (rotations and reflections) of a polyomino.
    The polyomino is represented as a set of (x, y) coordinate tuples.
    """
    symmetries = []
    current_poly = set(poly)
    
    # Generate 4 rotations, then reflect and generate 4 more rotations.
    for _ in range(2):
        for _ in range(4):
            symmetries.append(frozenset(current_poly))
            current_poly = {(-y, x) for x, y in current_poly} # 90-degree rotation
        current_poly = {(-x, y) for x, y in current_poly} # Reflection across y-axis
        
    return symmetries

def normalize(poly):
    """
    Translates a polyomino to the origin and returns it as a sorted tuple of coordinates.
    This sorted tuple form is the "canonical representation".
    """
    if not poly:
        return tuple()
    min_x = min(p[0] for p in poly)
    min_y = min(p[1] for p in poly)
    return tuple(sorted([(x - min_x, y - min_y) for x, y in poly]))

def canonical(poly):
    """
    Finds the canonical form of a polyomino by taking the lexicographically
    smallest of its 8 normalized symmetries.
    """
    return min(normalize(p) for p in get_transformations(poly))

def get_neighbors(cell):
    """Returns the four neighbors of a grid cell."""
    x, y = cell
    return [(x + 1, y), (x - 1, y), (x, y + 1), (x, y - 1)]

def generate_polydominoes(n):
    """
    Generates all non-equivalent polydominoes of order n.
    """
    if n == 1:
        # The base case is the canonical form of a single domino.
        return {canonical(frozenset([(0, 0), (0, 1)]))}

    smaller_polys = generate_polydominoes(n - 1)
    new_polys = set()

    for poly_tuple in smaller_polys:
        poly_set = frozenset(poly_tuple)
        perimeter = set()
        for cell in poly_set:
            for neighbor in get_neighbors(cell):
                if neighbor not in poly_set:
                    perimeter.add(neighbor)
        
        for p_cell in perimeter:
            x, y = p_cell
            # Try to add a horizontal domino
            cell2_h = (x + 1, y)
            if cell2_h not in poly_set:
                new_poly = poly_set.union([p_cell, cell2_h])
                new_polys.add(canonical(new_poly))

            # Try to add a vertical domino
            cell2_v = (x, y + 1)
            if cell2_v not in poly_set:
                new_poly = poly_set.union([p_cell, cell2_v])
                new_polys.add(canonical(new_poly))
                
    return new_polys

def has_hamiltonian_path(poly_tuple):
    """
    Checks if a polyomino has a Hamiltonian path using a backtracking search.
    """
    squares = list(poly_tuple)
    num_squares = len(squares)
    if num_squares <= 1:
        return True
    
    adj = collections.defaultdict(list)
    for i in range(num_squares):
        for j in range(i + 1, num_squares):
            s1 = squares[i]
            s2 = squares[j]
            if abs(s1[0] - s2[0]) + abs(s1[1] - s2[1]) == 1:
                adj[s1].append(s2)
                adj[s2].append(s1)

    # A graph with more than 2 nodes of degree 1 cannot have a Hamiltonian path.
    if [len(v) for v in adj.values()].count(1) > 2:
        return False

    path = []
    visited = set()

    def backtrack(u):
        path.append(u)
        visited.add(u)

        if len(path) == num_squares:
            return True

        for v in adj[u]:
            if v not in visited:
                if backtrack(v):
                    return True
        
        path.pop()
        visited.remove(u)
        return False

    # Start the search from any node.
    return backtrack(squares[0])

def solve_T(n):
    """
    Calculates T(n) by generating all polydominoes of order n and
    counting how many have a Hamiltonian path.
    """
    polydominoes = generate_polydominoes(n)
    
    traversable_count = 0
    for poly in polydominoes:
        if has_hamiltonian_path(poly):
            traversable_count += 1
            
    return traversable_count

# Main execution
if __name__ == "__main__":
    n = 4
    result = solve_T(n)
    print(f"T({n}) = {result}")
