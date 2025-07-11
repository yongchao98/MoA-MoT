def solve_puzzle(input_str):
    """
    Solves the puzzle by moving the number 2 based on a specific rule.

    The rule is:
    1. Find the coordinates of the number 2.
    2. Find all adjacent neighbors (including diagonals) that have a value of 0.
    3. Sort these "zero-neighbors" by their coordinates in reading order (top-to-bottom, left-to-right).
    4. The number 2 moves to the neighbor at index `2 % (number of zero-neighbors)`.
    5. The original position of the 2 becomes 0.
    """
    # 1. Parse the input string into a 2D grid of integers.
    grid_str_rows = input_str.split(',')
    grid = [list(map(int, list(row))) for row in grid_str_rows]
    height = len(grid)
    width = len(grid[0])

    # 2. Find the coordinates (r2, c2) of the number 2.
    r2, c2 = -1, -1
    for r in range(height):
        for c in range(width):
            if grid[r][c] == 2:
                r2, c2 = r, c
                break
        if r2 != -1:
            break

    if r2 == -1:
        print("Error: No '2' found in the input grid.")
        return

    # 3. Find all neighboring '0's in reading order.
    zero_neighbors = []
    # Iterate through potential neighbors in reading order
    for dr in [-1, 0, 1]:
        for dc in [-1, 0, 1]:
            if dr == 0 and dc == 0:
                continue
            
            nr, nc = r2 + dr, c2 + dc
            
            # Check if neighbor is within grid bounds
            if 0 <= nr < height and 0 <= nc < width:
                if grid[nr][nc] == 0:
                    zero_neighbors.append((nr, nc))

    # 4. If no zero-neighbors are found, the move is not possible.
    if not zero_neighbors:
        # This corresponds to the 'MISSING' case if it meant an impossible move.
        print("MISSING")
        return
        
    # 5. Determine the destination based on the modulo rule.
    num_zero_neighbors = len(zero_neighbors)
    target_index = 2 % num_zero_neighbors
    r0, c0 = zero_neighbors[target_index]
    
    # 6. Apply the transformation to the grid.
    grid[r2][c2] = 0
    grid[r0][c0] = 2

    # 7. Convert the modified grid back to a string and print it.
    output_rows = ["".join(map(str, row)) for row in grid]
    result = ",".join(output_rows)
    print(result)

# The input for which we need to find the "MISSING" value.
input_C = '000000,011120,111111'
solve_puzzle(input_C)