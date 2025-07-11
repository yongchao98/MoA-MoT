import sys

# Increase recursion limit for the path generation
sys.setrecursionlimit(2000)

def canonical(poly):
    """
    Calculates the canonical form of a polyomino.
    The canonical form is the lexicographically smallest representation
    among all 8 rotations and reflections.
    """
    forms = set()
    
    # Start with the given polyomino (as a list of sorted tuples)
    current_poly = sorted(list(poly))

    # Generate all 8 transformations (4 rotations + 4 reflected rotations)
    for _ in range(2):  # Original and reflected
        poly_to_transform = list(current_poly)
        for _ in range(4):  # 4 rotations
            # Normalize the coordinates to start at (0,0)
            min_x = min(p[0] for p in poly_to_transform)
            min_y = min(p[1] for p in poly_to_transform)
            normalized_form = tuple(sorted([(x - min_x, y - min_y) for x, y in poly_to_transform]))
            forms.add(normalized_form)
            
            # Rotate 90 degrees counter-clockwise: (x, y) -> (-y, x)
            poly_to_transform = [(-y, x) for x, y in poly_to_transform]
        
        # Reflect across y-axis for the next iteration: (x, y) -> (-x, y)
        current_poly = [(-x, y) for x, y in current_poly]
        
    return min(forms)

def is_tileable(poly):
    """
    Checks if a polyomino can be tiled by 1x2 dominoes using a backtracking algorithm.
    """
    # A polyomino must have an even number of squares to be tileable.
    if len(poly) % 2 != 0:
        return False

    # A quick check for checkerboard coloring. Necessary but not always sufficient.
    white_squares = sum(1 for x, y in poly if (x + y) % 2 == 0)
    black_squares = len(poly) - white_squares
    if white_squares != black_squares:
        return False

    # A recursive function to find a perfect matching (tiling)
    def can_tile_recursive(remaining_poly):
        if not remaining_poly:
            return True
        
        # Pick the lexicographically smallest square to tile.
        p1 = min(remaining_poly)
        x1, y1 = p1
        
        # Define its horizontal and vertical neighbors
        neighbors = [(x1 + 1, y1), (x1, y1 + 1)]
        
        for p2 in neighbors:
            if p2 in remaining_poly:
                # If a neighbor is part of the polyomino, try to place a domino
                new_remaining = remaining_poly - {p1, p2}
                if can_tile_recursive(new_remaining):
                    return True
        
        # If no tiling worked from this square, backtrack.
        return False

    return can_tile_recursive(frozenset(poly))

def solve():
    """
    Finds all non-equivalent octominoes that have a Hamiltonian path
    and can be tiled by dominoes.
    """
    solutions = set()
    # Cache to avoid re-checking tiling for already processed shapes
    checked_canonical_forms = {} 

    def extend_path(path, visited_squares):
        if len(path) == 8:
            poly = frozenset(path)
            canon_poly = canonical(poly)
            
            if canon_poly not in checked_canonical_forms:
                if is_tileable(canon_poly):
                    checked_canonical_forms[canon_poly] = True
                    solutions.add(canon_poly)
                    print(f"Found solution #{len(solutions)}:")
                    
                    # Pretty print the shape
                    coords = sorted(list(canon_poly))
                    max_x = max(c[0] for c in coords)
                    max_y = max(c[1] for c in coords)
                    grid = [['.' for _ in range(max_x + 1)] for _ in range(max_y + 1)]
                    for x, y in coords:
                        grid[y][x] = '#'
                    for row in reversed(grid):
                        print(" ".join(row))
                    print("-" * 20)

                else:
                    checked_canonical_forms[canon_poly] = False
            elif checked_canonical_forms[canon_poly]:
                 solutions.add(canon_poly)

            return

        last_x, last_y = path[-1]
        for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            next_pos = (last_x + dx, last_y + dy)
            if next_pos not in visited_squares:
                extend_path(path + [next_pos], visited_squares | {next_pos})

    print("Searching for traversable and tileable octominoes...")
    # Start the search from a single square (0,0)
    extend_path([(0, 0)], {(0, 0)})
    
    count = len(solutions)
    print(f"The number of non-equivalent polyforms of order 4, T(4), is {count}.")
    print(f"\nFinal Equation: T(4) = {count}")


solve()

<<<8>>>