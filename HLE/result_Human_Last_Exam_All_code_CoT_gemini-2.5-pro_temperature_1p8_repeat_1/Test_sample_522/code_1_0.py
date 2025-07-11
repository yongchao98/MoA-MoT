def solve_puzzle(input_str):
    """
    Solves the puzzle by moving the '2' based on a derived rule.

    The rule is as follows:
    1. Find the position of the '2'.
    2. Count the total number of '0's in the grid.
    3. Find all neighboring cells of the '2' that contain a '0'.
    4. The destination for the '2' is determined by the formula:
       index = (total number of '0's) % (number of '0'-neighbors).
       The '0'-neighbors are sorted in reading order (top-to-bottom, left-to-right).
    5. The '2' moves to the selected '0'-neighbor, and its original position becomes '0'.
    """
    grid = [list(row) for row in input_str.split(',')]
    height = len(grid)
    width = len(grid[0])

    two_pos = None
    total_zeros = 0

    # Find the '2' and count all '0's
    for r in range(height):
        for c in range(width):
            if grid[r][c] == '2':
                two_pos = (r, c)
            if grid[r][c] == '0':
                total_zeros += 1

    if two_pos is None:
        print("Error: No '2' found in the grid.")
        return

    # Find all '0'-neighbors in reading order
    zero_neighbors = []
    # Relative positions of 8 neighbors in reading order
    directions = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]

    for dr, dc in directions:
        nr, nc = two_pos[0] + dr, two_pos[1] + dc
        if 0 <= nr < height and 0 <= nc < width and grid[nr][nc] == '0':
            zero_neighbors.append((nr, nc))

    if not zero_neighbors:
        print("Error: The '2' has no '0'-neighbors.")
        return
        
    # Calculate the destination index
    dest_index = total_zeros % len(zero_neighbors)
    dest_pos = zero_neighbors[dest_index]

    # Update the grid
    grid[two_pos[0]][two_pos[1]] = '0'
    grid[dest_pos[0]][dest_pos[1]] = '2'

    # Format and print the output
    output_rows = ["".join(row) for row in grid]
    result = ",".join(output_rows)
    print(result)

# Input for the puzzle to be solved
input_c = '000000,011120,111111'
solve_puzzle(input_c)