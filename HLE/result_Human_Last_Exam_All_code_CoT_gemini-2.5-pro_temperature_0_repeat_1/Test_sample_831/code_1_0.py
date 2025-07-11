import sys

# Increase recursion limit for the search algorithms
sys.setrecursionlimit(2000)

# Memoization caches to store results for previously seen shapes
memo_tile = {}
memo_hamiltonian = {}

def get_canonical(points):
    """
    Computes the canonical form of a polyomino.
    The canonical form is the lexicographically smallest representation among all 8 symmetries.
    """
    # points is a frozenset of (x, y) tuples
    symmetries = []
    p_list = list(points)

    # Generate 8 symmetries (4 rotations of original + 4 rotations of reflection)
    for _ in range(2):
        # Reflect across y-axis: (x, y) -> (-x, y)
        p_list = [(-x, y) for x, y in p_list]
        
        for _ in range(4):
            # Rotate 90 degrees clockwise: (x, y) -> (y, -x)
            p_list = [(y, -x) for x, y in p_list]
            
            # Normalize by translating to origin (min_x=0, min_y=0)
            min_x = min(p[0] for p in p_list)
            min_y = min(p[1] for p in p_list)
            
            normalized_poly = tuple(sorted([(x - min_x, y - min_y) for x, y in p_list]))
            symmetries.append(normalized_poly)
            
    return min(symmetries)

def generate_free_polyominoes(n, memo_gen={}):
    """
    Generates all free polyominoes of size n using a recursive approach.
    """
    if n in memo_gen:
        return memo_gen[n]
    if n == 1:
        return {frozenset([(0, 0)])}
    
    smaller_polys = generate_free_polyominoes(n - 1, memo_gen)
    new_fixed_polys = set()
    
    # Add a new square to each smaller polyomino
    for poly in smaller_polys:
        for x, y in poly:
            for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                neighbor = (x + dx, y + dy)
                if neighbor not in poly:
                    new_poly = poly.union({neighbor})
                    new_fixed_polys.add(new_poly)
    
    # Convert to canonical forms to count only free polyominoes
    canonical_forms = {get_canonical(poly) for poly in new_fixed_polys}
    memo_gen[n] = {frozenset(p) for p in canonical_forms}
    return memo_gen[n]

def can_be_tiled(points):
    """
    Checks if a polyomino can be tiled by 1x2 dominoes.
    """
    points_tuple = tuple(sorted(list(points)))
    if points_tuple in memo_tile:
        return memo_tile[points_tuple]

    # A polyomino must have an even number of squares to be tiled.
    if len(points) % 2 != 0:
        memo_tile[points_tuple] = False
        return False

    # It must have an equal number of black and white squares.
    white_squares = sum(1 for x, y in points if (x + y) % 2 == 0)
    black_squares = len(points) - white_squares
    if white_squares != black_squares:
        memo_tile[points_tuple] = False
        return False

    # Backtracking search for a perfect matching (tiling)
    def solve(current_points):
        if not current_points:
            return True
        
        p1 = min(current_points)  # Pick a canonical point to start
        x1, y1 = p1
        
        # Try to tile it with a neighbor
        for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
            p2 = (x1 + dx, y1 + dy)
            if p2 in current_points:
                remaining = current_points - {p1, p2}
                if solve(remaining):
                    return True
        return False

    result = solve(points)
    memo_tile[points_tuple] = result
    return result

def has_hamiltonian_path(points):
    """
    Checks if a polyomino's graph has a Hamiltonian path.
    """
    points_tuple = tuple(sorted(list(points)))
    if points_tuple in memo_hamiltonian:
        return memo_hamiltonian[points_tuple]

    n = len(points)
    adj = {p: [] for p in points}
    for x1, y1 in points:
        p1 = (x1, y1)
        for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            p2 = (x1 + dx, y1 + dy)
            if p2 in points:
                adj[p1].append(p2)

    # Backtracking search for a path visiting all n vertices
    def find_path(path, visited_set):
        if len(path) == n:
            return True
        
        last_point = path[-1]
        for neighbor in adj[last_point]:
            if neighbor not in visited_set:
                if find_path(path + [neighbor], visited_set | {neighbor}):
                    return True
        return False

    # Try starting the path from each point in the polyomino
    for start_point in points:
        if find_path([start_point], {start_point}):
            memo_hamiltonian[points_tuple] = True
            return True
    
    memo_hamiltonian[points_tuple] = False
    return False

def solve_t4():
    """
    Calculates T(4) by generating, filtering, and counting polyominoes.
    """
    # For T(4), we consider polyominoes of size 2*4 = 8 (octominoes).
    print("Generating all 369 free octominoes...")
    octominoes = generate_free_polyominoes(8)
    print(f"Generated {len(octominoes)} octominoes. Now checking properties...")

    count = 0
    for i, poly in enumerate(octominoes):
        # A polyomino for T(n) must be tilable by dominoes
        if can_be_tiled(poly):
            # It must also have a Hamiltonian path
            if has_hamiltonian_path(poly):
                count += 1
    
    return count

if __name__ == '__main__':
    t4_value = solve_t4()
    print("\nThe number of non-equivalent polyforms of order 4, T(4), is:")
    print(t4_value)