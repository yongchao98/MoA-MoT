import collections

def normalize_polyomino(polyomino):
    """
    Finds the canonical representation of a polyomino by checking all 8 orientations
    (4 rotations and 4 reflections) and picking the lexicographically smallest one.
    A polyomino is a frozenset of (x, y) tuples.
    """
    representations = []
    # Use a list to allow for transformations
    points = list(polyomino)

    for _ in range(2):  # Original and reflected
        # 4 rotations
        for _ in range(4):
            # Translate to origin
            min_x = min(p[0] for p in points)
            min_y = min(p[1] for p in points)
            translated = frozenset((x - min_x, y - min_y) for x, y in points)
            
            # Sort to create a canonical tuple for comparison
            representations.append(tuple(sorted(list(translated))))
            
            # Rotate 90 degrees: (x, y) -> (-y, x)
            points = [(-y, x) for x, y in points]
            
        # Reflect across y-axis: (x, y) -> (-x, y)
        points = [(-x, y) for x, y in points]
        
    return min(representations)

def generate_free_polyominoes(n):
    """
    Generates all free polyominoes of a given size n.
    """
    if n == 1:
        return {frozenset({(0, 0)})}
        
    smaller_polys = generate_free_polyominoes(n - 1)
    new_polys = set()
    
    for poly in smaller_polys:
        neighbors = set()
        for x, y in poly:
            neighbors.add((x + 1, y))
            neighbors.add((x - 1, y))
            neighbors.add((x, y + 1))
            neighbors.add((x, y - 1))
        
        for neighbor in neighbors:
            if neighbor not in poly:
                new_poly = poly.union({neighbor})
                # Add the canonical form to the set to ensure uniqueness
                new_polys.add(normalize_polyomino(new_poly))
                
    return new_polys

def is_tileable_by_dominoes(polyomino):
    """
    Checks if a polyomino can be tiled by 1x2 dominoes using a backtracking algorithm.
    """
    # Condition 1: Must have an even number of squares.
    if len(polyomino) % 2 != 0:
        return False
        
    # Condition 2: Chessboard coloring must be balanced.
    num_white_squares = sum(1 for x, y in polyomino if (x + y) % 2 != 0)
    if num_white_squares != len(polyomino) / 2:
        return False

    # Backtracking solver for the exact cover problem.
    memo = {}
    def can_tile(uncovered_squares):
        if not uncovered_squares:
            return True
        
        # Memoization: use a frozenset of the tuple as key
        state = frozenset(uncovered_squares)
        if state in memo:
            return memo[state]

        # Pick a square to cover (the lexicographically smallest one for consistency)
        p_x, p_y = min(uncovered_squares)
        
        # Try placing a domino horizontally
        neighbor_h = (p_x + 1, p_y)
        if neighbor_h in uncovered_squares:
            remaining_h = tuple(sorted(list(uncovered_squares - {(p_x, p_y), neighbor_h})))
            if can_tile(remaining_h):
                memo[state] = True
                return True

        # Try placing a domino vertically
        neighbor_v = (p_x, p_y + 1)
        if neighbor_v in uncovered_squares:
            remaining_v = tuple(sorted(list(uncovered_squares - {(p_x, p_y), neighbor_v})))
            if can_tile(remaining_v):
                memo[state] = True
                return True
        
        memo[state] = False
        return False

    return can_tile(tuple(sorted(list(polyomino))))

def has_hamiltonian_path(polyomino):
    """
    Checks if the graph of a polyomino has a Hamiltonian path.
    """
    n = len(polyomino)
    adj = {p: [] for p in polyomino}
    
    # Build adjacency list for the graph
    for x, y in polyomino:
        for dx, dy in [(0, 1), (1, 0)]: # Check only two directions to avoid duplicates
            neighbor = (x + dx, y + dy)
            if neighbor in polyomino:
                adj[(x, y)].append(neighbor)
                adj[neighbor].append((x, y))

    def find_path(current_node, visited):
        visited.add(current_node)
        if len(visited) == n:
            return True
        
        for neighbor in adj[current_node]:
            if neighbor not in visited:
                if find_path(neighbor, visited):
                    return True
        
        visited.remove(current_node)  # Backtrack
        return False

    # Try starting the path from each square
    for start_node in polyomino:
        if find_path(start_node, set()):
            return True
    return False

def solve_t4():
    """
    Calculates T(4) by generating, filtering, and counting octominoes.
    """
    # T(4) corresponds to polyominoes of size 8 (octominoes).
    N_SQUARES = 8
    N_DOMINOES = 4
    
    # Step 1: Generate all 369 free octominoes.
    free_octominoes = generate_free_polyominoes(N_SQUARES)
    
    count = 0
    # Steps 2 & 3: Filter and count
    for p_tuple in free_octominoes:
        p_frozenset = frozenset(p_tuple)
        
        # Check for domino tileability
        if is_tileable_by_dominoes(p_frozenset):
            # Check for Hamiltonian path
            if has_hamiltonian_path(p_frozenset):
                count += 1
    
    # Step 4: Output the final result
    print(f"The number of non-equivalent polyforms of order {N_DOMINOES} (octominoes) that are both domino-tileable and have a Hamiltonian path is:")
    print(f"T({N_DOMINOES}) = {count}")

if __name__ == '__main__':
    solve_t4()