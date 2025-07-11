def solve_puzzle(input_str):
    """
    Solves the puzzle by applying the discovered transformation rule.

    The rule is:
    1. The '2' moves to one of its '0' neighbors.
    2. The choice of neighbor is the one in the highest row.
    3. If there's a tie (multiple '0' neighbors in the highest row), the choice depends on the
       column of the original '2':
       - If the '2's column is odd, pick the leftmost neighbor.
       - If the '2's column is even, pick the rightmost neighbor.
    """
    grid = [list(row) for row in input_str.split(',')]
    rows = len(grid)
    cols = len(grid[0]) if rows > 0 else 0

    # 1. Find the '2'
    r2, c2 = -1, -1
    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == '2':
                r2, c2 = r, c
                break
        if r2 != -1:
            break
    
    if r2 == -1:
        print("No '2' found in the grid.")
        return

    print(f"Found '2' at row {r2}, column {c2}.")

    # 2. Find all '0' neighbors
    zero_neighbors = []
    for dr in [-1, 0, 1]:
        for dc in [-1, 0, 1]:
            if dr == 0 and dc == 0:
                continue
            rn, cn = r2 + dr, c2 + dc
            if 0 <= rn < rows and 0 <= cn < cols and grid[rn][cn] == '0':
                zero_neighbors.append((rn, cn))
    
    if not zero_neighbors:
        print("The '2' has no '0' neighbors to move to.")
        return

    print(f"'0' neighbors are: {zero_neighbors}")
    
    # 3. Filter for the top-most row
    min_r = min(r for r, c in zero_neighbors)
    top_zero_neighbors = [(r, c) for r, c in zero_neighbors if r == min_r]
    
    print(f"Top-most '0' neighbors are at row {min_r}: {top_zero_neighbors}")

    # 4. Select the target neighbor based on the column of the '2'
    target_neighbor = None
    if c2 % 2 != 0: # Odd column
        print(f"Column of '2' ({c2}) is odd, selecting left-most neighbor.")
        min_c = min(c for r, c in top_zero_neighbors)
        target_neighbor = (min_r, min_c)
    else: # Even column
        print(f"Column of '2' ({c2}) is even, selecting right-most neighbor.")
        max_c = max(c for r, c in top_zero_neighbors)
        target_neighbor = (min_r, max_c)
        
    print(f"Selected neighbor for new '2': {target_neighbor}")
    
    # 5. Apply the transformation
    rt, ct = target_neighbor
    grid[r2][c2] = '0'
    grid[rt][ct] = '2'

    # 6. Generate the final output string
    output_rows = ["".join(row) for row in grid]
    final_string = ",".join(output_rows)

    print("\nFinal grid as a string:")
    print(final_string)
    return final_string


# The input for which we need to find the MISSING value
Input_C = '000000,011120,111111'
missing_value = solve_puzzle(Input_C)