def solve_puzzle(input_str):
    """
    Solves the puzzle based on the derived transformation rule.
    The rule involves swapping the '2' with one of its neighboring '0's.
    The choice of which '0' neighbor depends on a scoring system and the
    parity of the sum of the '2's coordinates.
    """
    grid = [list(map(int, row)) for row in input_str.split(',')]
    rows = len(grid)
    cols = len(grid[0])
    
    # Find the location of '2'
    r2, c2 = -1, -1
    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == 2:
                r2, c2 = r, c
                break
        if r2 != -1:
            break

    # If '2' is not found, something is wrong with the input
    if r2 == -1:
        return "ERROR: '2' not found"

    # Find all neighboring '0' candidates
    candidates = []
    # Relative coordinates of 8 neighbors in reading order
    neighbor_deltas = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]
    for dr, dc in neighbor_deltas:
        nr, nc = r2 + dr, c2 + dc
        if 0 <= nr < rows and 0 <= nc < cols and grid[nr][nc] == 0:
            candidates.append((nr, nc))

    # If no '0' neighbors, the transformation is not possible
    if not candidates:
        return "MISSING"

    # Score each candidate by counting its '1' neighbors
    scores = []
    for r_cand, c_cand in candidates:
        one_neighbors_count = 0
        for dr, dc in neighbor_deltas:
            nr, nc = r_cand + dr, c_cand + dc
            if 0 <= nr < rows and 0 <= nc < cols and grid[nr][nc] == 1:
                one_neighbors_count += 1
        scores.append(one_neighbors_count)

    # Determine which rule to apply based on parity of r2+c2
    if (r2 + c2) % 2 == 0:  # Even -> MIN rule
        min_score = min(scores)
        best_candidates = [cand for i, cand in enumerate(candidates) if scores[i] == min_score]
    else:  # Odd -> MAX rule
        max_score = max(scores)
        best_candidates = [cand for i, cand in enumerate(candidates) if scores[i] == max_score]
    
    # The first candidate in reading order is the target (candidates are already sorted by reading order)
    r0, c0 = best_candidates[0]
    
    # Perform the swap
    grid[r2][c2], grid[r0][c0] = grid[r0][c0], grid[r2][c2]
    
    # Format the grid back to a string
    output_rows = ["".join(map(str, row)) for row in grid]
    return ",".join(output_rows)

# Input for the final question
input_c = '000000,011120,111111'

# Calculate and print the result
missing_value = solve_puzzle(input_c)
print(missing_value)