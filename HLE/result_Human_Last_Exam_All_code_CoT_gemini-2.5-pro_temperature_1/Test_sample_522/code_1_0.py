import math

def solve_puzzle(input_str):
    """
    Solves the puzzle by applying the derived transformation rule.
    """
    rows = input_str.split(',')
    grid = [list(r) for r in rows]
    height = len(grid)
    width = len(grid[0])

    # 1. Find the '2'
    r2, c2 = -1, -1
    for r in range(height):
        for c in range(width):
            if grid[r][c] == '2':
                r2, c2 = r, c
                break
        if r2 != -1:
            break

    if r2 == -1:
        print("No '2' found in the grid.")
        return

    # 2. Determine diagonal type and find candidate '0's
    candidate_zeros = []
    if c2 % 2 != 0:  # Odd column -> main diagonal
        diag_val = r2 - c2
        for r in range(height):
            for c in range(width):
                if r - c == diag_val and grid[r][c] == '0':
                    candidate_zeros.append((r, c))
    else:  # Even column -> anti-diagonal
        diag_val = r2 + c2
        for r in range(height):
            for c in range(width):
                if r + c == diag_val and grid[r][c] == '0':
                    candidate_zeros.append((r, c))

    # 3. Filter candidates: must have a '1' neighbor
    valid_targets = []
    for r0, c0 in candidate_zeros:
        has_one_neighbor = False
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue
                nr, nc = r0 + dr, c0 + dc
                if 0 <= nr < height and 0 <= nc < width and grid[nr][nc] == '1':
                    has_one_neighbor = True
                    break
            if has_one_neighbor:
                break
        if has_one_neighbor:
            valid_targets.append((r0, c0))

    # 4. If no valid targets, output is MISSING
    if not valid_targets:
        print("MISSING")
        return

    # 5. Find the closest valid target
    best_target = None
    min_dist = float('inf')

    for rt, ct in valid_targets:
        dist = abs(rt - r2) + abs(ct - c2)
        if dist < min_dist:
            min_dist = dist
            best_target = (rt, ct)
        elif dist == min_dist:
            # Tie-breaking: smallest row, then smallest column
            if rt < best_target[0] or (rt == best_target[0] and ct < best_target[1]):
                best_target = (rt, ct)

    # 6. Perform the move and print the result
    r0, c0 = best_target
    grid[r0][c0] = '2'
    grid[r2][c2] = '0'

    output_rows = ["".join(r) for r in grid]
    print(",".join(output_rows))


# Solve for the MISSING value using Input_C
input_c = '000000,011120,111111'
solve_puzzle(input_c)