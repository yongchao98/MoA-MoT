import collections

def canonical(polyomino):
    """
    Generates the canonical representation of a polyomino among its 8 symmetries.
    A polyomino is represented as a frozenset of (x, y) coordinates.
    The canonical form is the lexicographically smallest tuple of sorted coordinates,
    normalized to have min x and y coordinates as 0.
    """
    forms = []
    points = list(polyomino)
    
    # Generate 8 symmetries: 4 rotations, and their reflections
    for _ in range(2): # Original and reflection
        # Reflection over y-axis for the second iteration: (x, y) -> (-x, y)
        points = [(-x, y) for x, y in points]
        
        for _ in range(4): # 4 rotations
            # Rotation by 90 degrees clockwise: (x, y) -> (y, -x)
            # This is equivalent to counter-clockwise (-y, x) used in planning
            # as long as we are consistent.
            points = [(y, -x) for x, y in points]
            
            # Normalize the transformed polyomino to the origin
            min_x = min(p[0] for p in points)
            min_y = min(p[1] for p in points)
            normalized_points = tuple(sorted([(x - min_x, y - min_y) for x, y in points]))
            
            forms.append(normalized_points)
            
    return min(forms)

def get_neighbors(p):
    """Returns the 4 orthogonal neighbors of a point."""
    x, y = p
    return [(x + 1, y), (x - 1, y), (x, y + 1), (x, y - 1)]

def generate_polyominoes(n):
    """
    Generates all free polyominoes of size n using an iterative approach.
    """
    # Start with the set of (n-1)-ominoes
    if n == 1:
        # Base case: the monomino
        return {frozenset([(0, 0)])}
        
    smaller_polyominoes = generate_polyominoes(n - 1)
    
    # Generate n-ominoes by adding one square to the border of each (n-1)-omino
    new_polyominoes = set()
    for p in smaller_polyominoes:
        border = set()
        for square in p:
            for neighbor in get_neighbors(square):
                if neighbor not in p:
                    border.add(neighbor)
        
        for new_square in border:
            new_poly = p.union({new_square})
            new_polyominoes.add(frozenset(new_poly))

    # Reduce to canonical forms to count only unique "free" polyominoes
    canonical_forms = {canonical(p) for p in new_polyominoes}
    
    # Return polyominoes as frozensets for the next iteration
    return {frozenset(p_tuple) for p_tuple in canonical_forms}


def has_hamiltonian_path(polyomino):
    """
    Checks if a polyomino's graph has a Hamiltonian path using backtracking.
    """
    points = list(polyomino)
    n = len(points)
    adj = collections.defaultdict(list)
    point_to_idx = {p: i for i, p in enumerate(points)}

    for i, p1 in enumerate(points):
        for p2 in get_neighbors(p1):
            if p2 in point_to_idx:
                adj[i].append(point_to_idx[p2])

    def backtrack(path, visited):
        if len(path) == n:
            return True
        
        last_node_idx = path[-1]
        for neighbor_idx in adj[last_node_idx]:
            if not visited[neighbor_idx]:
                visited[neighbor_idx] = True
                path.append(neighbor_idx)
                if backtrack(path, visited):
                    return True
                path.pop()
                visited[neighbor_idx] = False
        return False

    # Try starting the path from every square (node)
    for start_node_idx in range(n):
        path = [start_node_idx]
        visited = [False] * n
        visited[start_node_idx] = True
        if backtrack(path, visited):
            return True
            
    return False

def is_tileable(polyomino):
    """
    Checks if a polyomino can be tiled by 1x2 dominoes.
    This is equivalent to finding a perfect matching in its bipartite graph.
    """
    white_squares = {p for p in polyomino if (p[0] + p[1]) % 2 == 0}
    black_squares = {p for p in polyomino if (p[0] + p[1]) % 2 != 0}

    if len(white_squares) != len(black_squares):
        return False

    # Adjacency list for the bipartite graph (from white to black squares)
    adj = {w: [n for n in get_neighbors(w) if n in black_squares] for w in white_squares}

    # Bipartite matching algorithm (using augmenting paths via DFS)
    match = {}  # Stores matching: black_square -> white_square

    def dfs(w_node, visited_b):
        for b_node in adj[w_node]:
            if b_node not in visited_b:
                visited_b.add(b_node)
                if b_node not in match or dfs(match[b_node], visited_b):
                    match[b_node] = w_node
                    return True
        return False

    match_count = 0
    for w_node in white_squares:
        # For each attempt to match a white node, use a fresh set of visited black nodes
        if dfs(w_node, set()):
            match_count += 1
            
    return match_count == len(white_squares)


def solve_t4():
    """
    Calculates T(4) by generating all free octominoes and checking them
    for the two required properties: having a Hamiltonian path and being
    tileable by dominoes.
    """
    # Polyforms of order n=4 correspond to polyominoes of size 2*n = 8
    polyomino_size = 8
    
    # Generate all unique free octominoes
    octominoes = generate_polyominoes(polyomino_size)
    
    valid_polyforms_count = 0
    
    for p in octominoes:
        # A shape is a solution if it has a Hamiltonian path AND is tileable
        if has_hamiltonian_path(p):
            if is_tileable(p):
                valid_polyforms_count += 1
    
    return valid_polyforms_count

# Calculate and print the final answer for T(4)
t4_value = solve_t4()
print(f"T(4) = {t4_value}")
<<<106>>>