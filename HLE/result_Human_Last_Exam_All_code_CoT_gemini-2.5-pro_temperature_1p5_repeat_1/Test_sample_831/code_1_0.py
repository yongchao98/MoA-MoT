import collections

def canonicalize_fixed(dominoes):
    """
    Translates a list of dominoes to a canonical form for a 'fixed' arrangement.
    The form is a sorted tuple of dominoes, translated so that the
    lexicographically smallest square is at the origin (0,0).
    """
    if not dominoes:
        return tuple()
    
    all_squares = []
    for d_part1, d_part2 in dominoes:
        all_squares.append(d_part1)
        all_squares.append(d_part2)

    min_x = min(s[0] for s in all_squares)
    min_y = min(s[1] for s in all_squares)

    new_dominoes = []
    for d_part1, d_part2 in dominoes:
        new_d1 = (d_part1[0] - min_x, d_part1[1] - min_y)
        new_d2 = (d_part2[0] - min_x, d_part2[1] - min_y)
        # Sort squares within a domino for a consistent representation
        if new_d1 > new_d2:
            new_d1, new_d2 = new_d2, new_d1
        new_dominoes.append((new_d1, new_d2))

    # Sort the dominoes themselves
    new_dominoes.sort()
    return tuple(new_dominoes)

def generate_fixed(n, _cache={}):
    """
    Generates all 'fixed' arrangements of n dominoes.
    This is done recursively, adding one domino at a time.
    """
    if n in _cache:
        return _cache[n]
    if n == 0:
        return {tuple()}
    if n == 1:
        # Start with a horizontal and a vertical domino
        return {(((0, 0), (0, 1)),), (((0, 0), (1, 0)),)}

    prev_arrangements = generate_fixed(n - 1)
    new_arrangements = set()

    for arrangement in prev_arrangements:
        all_squares = {sq for d in arrangement for sq in d}
        
        # Find all empty squares adjacent to the current shape
        candidate_squares = set()
        for x, y in all_squares:
            for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                neighbor = (x + dx, y + dy)
                if neighbor not in all_squares:
                    candidate_squares.add(neighbor)
        
        # Try to place a new domino on each candidate square
        for x1, y1 in candidate_squares:
            # Horizontal domino
            x2, y2 = x1 + 1, y1
            if (x2, y2) not in all_squares:
                new_domino = ((x1, y1), (x2, y2))
                new_arrangements.add(canonicalize_fixed(list(arrangement) + [new_domino]))
            
            # Vertical domino
            x2, y2 = x1, y1 + 1
            if (x2, y2) not in all_squares:
                new_domino = ((x1, y1), (x2, y2))
                new_arrangements.add(canonicalize_fixed(list(arrangement) + [new_domino]))
    
    _cache[n] = new_arrangements
    return new_arrangements

def canonicalize_free(arrangement):
    """
    Finds the canonical representation for a 'free' polyform by checking
    all 8 symmetries (rotations and reflections).
    """
    symmetries = set()
    current_arrangement = arrangement

    for _ in range(4): # 4 rotations
        # Reflection
        reflected_arrangement = [((-p1[0], p1[1]), (-p2[0], p2[1])) for p1, p2 in current_arrangement]
        symmetries.add(canonicalize_fixed(reflected_arrangement))
        
        symmetries.add(canonicalize_fixed(list(current_arrangement)))
        
        # Rotation (90 degrees counter-clockwise)
        current_arrangement = [((-p1[1], p1[0]), (-p2[1], p2[0])) for p1, p2 in current_arrangement]
        
    return min(symmetries)

def has_hamiltonian_path(squares):
    """
    Checks if a polyomino shape (represented by a set of squares)
    has a Hamiltonian path using a backtracking algorithm.
    """
    n = len(squares)
    if n <= 1:
        return True
    
    adj = {s: [] for s in squares}
    for x, y in squares:
        for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            neighbor = (x + dx, y + dy)
            if neighbor in squares:
                adj[x, y].append(neighbor)

    # Any node that is a single point of connection between parts of the graph
    # (a cut vertex) makes a Hamiltonian path impossible for n > 2.
    # This check helps prune some shapes early.
    for s in squares:
        # Temporarily remove s and check if the graph splits
        remaining_squares = squares - {s}
        if not remaining_squares: continue
        
        q = collections.deque([next(iter(remaining_squares))])
        visited = {next(iter(remaining_squares))}
        count = 0
        while q:
            count += 1
            u = q.popleft()
            for v_neighbor in adj[u]:
                if v_neighbor in remaining_squares and v_neighbor not in visited:
                    visited.add(v_neighbor)
                    q.append(v_neighbor)
        if count < n - 1:
            return False


    # Backtracking search for the path
    path = []
    path_set = set()
    def search():
        if len(path) == n:
            return True
        last_node = path[-1]
        for neighbor in adj[last_node]:
            if neighbor not in path_set:
                path.append(neighbor)
                path_set.add(neighbor)
                if search():
                    return True
                path_set.remove(neighbor) # backtrack
                path.pop() # backtrack
        return False

    # Try starting a path from every square
    for start_node in sorted(list(squares)):
        path = [start_node]
        path_set = {start_node}
        if search():
            return True
            
    return False

def solve_t4():
    """
    Main function to calculate T(4).
    """
    N = 4
    
    # 1. Generate fixed arrangements
    fixed_arrangements = generate_fixed(N)
    
    # 2. Find unique free arrangements
    free_arrangements = set()
    for arr in fixed_arrangements:
        canon_free_form = canonicalize_free(arr)
        free_arrangements.add(canon_free_form)
        
    # 3. Check for Hamiltonian path and count
    hp_count = 0
    traversable_polyforms = []
    for free_arr in sorted(list(free_arrangements)):
        squares = {sq for d in free_arr for sq in d}
        if has_hamiltonian_path(squares):
            hp_count += 1
            traversable_polyforms.append(free_arr)
            
    # 4. Output the result
    sum_str = " + ".join(["1"] * hp_count)
    if not sum_str: sum_str = "0"
    print(f"Found {len(free_arrangements)} non-equivalent polyforms of order {N}.")
    print(f"Found {hp_count} of them can be traversed.")
    print(f"T({N}) = {sum_str} = {hp_count}")


if __name__ == '__main__':
    solve_t4()
<<<26>>>