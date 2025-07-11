def solve_puzzle(input_str):
    """
    Solves the puzzle by applying the deduced transformation rules.
    """
    # 1. Parse the input string into a 2D grid
    rows = input_str.split(',')
    grid = [list(row) for row in rows]
    height = len(grid)
    width = len(grid[0])

    # 2. Find the coordinates of the '2'
    r2, c2 = -1, -1
    try:
        for r_idx, row in enumerate(grid):
            c_idx = row.index('2')
            r2, c2 = r_idx, c_idx
            break
    except ValueError:
        # No '2' found, which is an invalid state based on examples
        print("MISSING")
        return

    # 3. Find candidate '0's adjacent to the '2'
    deltas = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]
    candidate_coords = []
    for dr, dc in deltas:
        nr, nc = r2 + dr, c2 + dc
        if 0 <= nr < height and 0 <= nc < width and grid[nr][nc] == '0':
            candidate_coords.append((nr, nc))

    # If there are no candidate '0's, the result is MISSING
    if not candidate_coords:
        print("MISSING")
        return

    # 4. Score each candidate '0' by counting its '0' neighbors
    scored_candidates = []
    for r_cand, c_cand in candidate_coords:
        score = 0
        for dr, dc in deltas:
            nr, nc = r_cand + dr, c_cand + dc
            if 0 <= nr < height and 0 <= nc < width and grid[nr][nc] == '0':
                score += 1
        scored_candidates.append({'r': r_cand, 'c': c_cand, 'score': score})

    # 5. Determine the best candidate based on the rules
    min_score = min(c['score'] for c in scored_candidates)
    best_candidates = [c for c in scored_candidates if c['score'] == min_score]

    tie_count = len(best_candidates)

    if tie_count > 2:
        # Ambiguous case with more than 2 tied candidates
        print("MISSING")
    else:
        # If tie_count is 1 or 2, sort by reading order to find the target
        best_candidates.sort(key=lambda c: (c['r'], c['c']))
        target = best_candidates[0]
        rt, ct = target['r'], target['c']

        # Perform the swap on the grid
        grid[r2][c2] = '0'
        grid[rt][ct] = '2'

        # Format the result grid back into a string
        output_rows = ["".join(row) for row in grid]
        result = ",".join(output_rows)
        print(result)

# The input for which we need to find the "MISSING" value
input_c = '000000,011120,111111'

# Solve the puzzle for Input_C
solve_puzzle(input_c)