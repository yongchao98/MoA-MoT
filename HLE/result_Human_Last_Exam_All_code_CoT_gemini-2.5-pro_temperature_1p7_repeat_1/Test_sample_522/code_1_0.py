def solve_puzzle():
    """
    Solves the puzzle by finding a '2', calculating its neighbor sum S2,
    then finding the first '0' (in reading order) whose neighbor sum S0
    satisfies the rule S2 + S0 = 9. If found, the values are swapped.
    If not found, the result is 'MISSING'.
    """
    input_str = '000000,011120,111111'

    rows = input_str.split(',')
    try:
        grid = [[int(char) for char in row] for row in rows]
        height = len(grid)
        width = len(grid[0])
    except (ValueError, IndexError):
        print("Invalid grid format.")
        return

    # Find coordinates of the '2'
    r2, c2 = -1, -1
    for r in range(height):
        for c in range(width):
            if grid[r][c] == 2:
                r2, c2 = r, c
                break
        if r2 != -1:
            break
    
    if r2 == -1:
        print("MISSING")
        return

    def get_neighbor_sum(g, r_center, c_center):
        s = 0
        h = len(g)
        w = len(g[0])
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue
                r, c = r_center + dr, c_center + dc
                if 0 <= r < h and 0 <= c < w:
                    s += g[r][c]
        return s

    s2 = get_neighbor_sum(grid, r2, c2)
    s_target_req = 9 - s2

    # Find the first '0' that matches the required neighbor sum
    target_coords = None
    for r in range(height):
        for c in range(width):
            if grid[r][c] == 0:
                s0 = get_neighbor_sum(grid, r, c)
                if s0 == s_target_req:
                    target_coords = (r, c)
                    break
        if target_coords:
            break

    if target_coords is None:
        print("MISSING")
    else:
        rt, ct = target_coords
        grid[r2][c2], grid[rt][ct] = 0, 2
        
        output_rows = ["".join(map(str, row)) for row in grid]
        print(",".join(output_rows))

solve_puzzle()