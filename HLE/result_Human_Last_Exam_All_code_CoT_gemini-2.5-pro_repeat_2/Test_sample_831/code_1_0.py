import collections

def solve_t4():
    """
    This function calculates T(4) by generating all free octominoes and
    checking them for domino tileability and the existence of a Hamiltonian path.
    """

    def get_canonical(poly):
        """
        Computes the canonical representation of a polyomino.
        A polyomino is a frozenset of (x, y) coordinates.
        The canonical form is the lexicographically smallest representation
        among all 8 symmetries (rotations and reflections).
        """
        symmetries = []
        points = sorted(list(poly))

        for _ in range(4):  # 4 rotations
            # Normalize by translating min x, min y to 0,0
            min_x = min(p[0] for p in points)
            min_y = min(p[1] for p in points)
            symmetries.append(frozenset((p[0] - min_x, p[1] - min_y) for p in points))

            # Also add reflection across the y-axis
            reflected_points = [(-p[0], p[1]) for p in points]
            min_x = min(p[0] for p in reflected_points)
            min_y = min(p[1] for p in reflected_points)
            symmetries.append(frozenset((p[0] - min_x, p[1] - min_y) for p in reflected_points))

            # Rotate points 90 degrees: (x, y) -> (y, -x)
            points = sorted([(p[1], -p[0]) for p in points])

        return min(symmetries, key=lambda p: sorted(list(p)))

    def check_tileability(poly):
        """
        Checks if a polyomino can be tiled by dominoes using backtracking.
        """
        squares = tuple(sorted(list(poly)))
        
        # Quick check for black/white square balance
        if sum(1 for x, y in squares if (x + y) % 2 != 0) != len(squares) / 2:
            return False

        adj = collections.defaultdict(list)
        square_set = set(squares)
        for x, y in squares:
            for dx, dy in [(0, 1), (1, 0)]: # Check only two directions to avoid duplicate pairs
                neighbor = (x + dx, y + dy)
                if neighbor in square_set:
                    adj[(x, y)].append(neighbor)
        
        memo = {}
        def can_tile(current_squares_tuple):
            if not current_squares_tuple:
                return True
            if current_squares_tuple in memo:
                return memo[current_squares_tuple]

            s1 = current_squares_tuple[0]
            current_squares_set = set(current_squares_tuple)
            
            for s2 in adj[s1]:
                if s2 in current_squares_set:
                    remaining_squares = list(current_squares_set - {s1, s2})
                    remaining_squares.sort()
                    if can_tile(tuple(remaining_squares)):
                        memo[current_squares_tuple] = True
                        return True
            
            memo[current_squares_tuple] = False
            return False

        return can_tile(squares)

    def check_hamiltonian_path(poly):
        """
        Checks if the graph of a polyomino has a Hamiltonian path.
        """
        squares = list(poly)
        n = len(squares)
        adj = collections.defaultdict(list)
        square_to_idx = {sq: i for i, sq in enumerate(squares)}

        for i, (x, y) in enumerate(squares):
            for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                neighbor = (x + dx, y + dy)
                if neighbor in square_to_idx:
                    adj[i].append(square_to_idx[neighbor])

        memo = {}
        def find_path(last_node, visited_mask):
            if visited_mask == (1 << n) - 1:
                return True
            
            state = (last_node, visited_mask)
            if state in memo:
                return memo[state]

            for neighbor in adj[last_node]:
                if not ((visited_mask >> neighbor) & 1):
                    if find_path(neighbor, visited_mask | (1 << neighbor)):
                        memo[state] = True
                        return True

            memo[state] = False
            return False

        for start_node in range(n):
            if find_path(start_node, 1 << start_node):
                return True
        return False

    # Step 1: Generate all free octominoes
    polys = {frozenset([(0, 0)])}
    for _ in range(7):  # Expand 7 times to get from size 1 to 8
        new_polys = set()
        for poly in polys:
            for x, y in poly:
                for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                    neighbor = (x + dx, y + dy)
                    if neighbor not in poly:
                        new_polys.add(poly.union({neighbor}))
        
        canonical_polys = set()
        for poly in new_polys:
            canonical_polys.add(get_canonical(poly))
        polys = canonical_polys
    
    free_octominoes = polys
    
    # Step 2 & 3: Filter and find solutions
    solutions = []
    for p_canonical in free_octominoes:
        if check_tileability(p_canonical):
            if check_hamiltonian_path(p_canonical):
                solutions.append(p_canonical)
    
    # Step 4: Group solutions by bounding box and print the result
    grouped_solutions = collections.defaultdict(list)
    for poly in solutions:
        points = list(poly)
        min_x = min(p[0] for p in points)
        max_x = max(p[0] for p in points)
        min_y = min(p[1] for p in points)
        max_y = max(p[1] for p in points)
        bbox = tuple(sorted((max_x - min_x + 1, max_y - min_y + 1)))
        grouped_solutions[bbox].append(poly)

    counts = [len(v) for k, v in sorted(grouped_solutions.items())]
    total = sum(counts)
    equation_parts = [str(c) for c in counts]
    print(" + ".join(equation_parts) + f" = {total}")

solve_t4()
<<<8>>>