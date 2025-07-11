def solve_puzzle(input_str):
    """
    Solves the puzzle by applying the derived transformation rule.
    """
    # 1. Parse the input string into a grid of characters
    grid_rows = input_str.split(',')
    grid = [list(row) for row in grid_rows]

    # 2. Find the coordinates (r2, c2) of the '2'
    r2, c2 = -1, -1
    for r, row in enumerate(grid):
        try:
            c = row.index('2')
            r2, c2 = r, c
            break
        except ValueError:
            continue

    if r2 == -1:
        return "Error: No '2' found in the input."
    if r2 == 0:
        return "Error: '2' is in the top row, cannot move up."

    # 3. Count '0's to the left and right of the '2' in its own row
    row_of_2 = grid[r2]
    zeros_left = row_of_2[0:c2].count('0')
    zeros_right = row_of_2[c2 + 1:].count('0')

    # 4. Find all available '0's in the row above
    row_above = grid[r2 - 1]
    zero_indices = [i for i, char in enumerate(row_above) if char == '0']

    if not zero_indices:
        return "Error: No '0' to move to in the row above."

    # 5. Select the target column based on the rule
    if zeros_left > zeros_right:
        # Pick the LAST '0'
        c0 = zero_indices[-1]
    else:
        # Pick the FIRST '0'
        c0 = zero_indices[0]
    
    r0 = r2 - 1

    # 6. Perform the swap
    grid[r2][c2] = '0'
    grid[r0][c0] = '2'

    # 7. Convert the modified grid back to a string and return
    output_rows = ["".join(row) for row in grid]
    return ",".join(output_rows)

# Input for the puzzle
input_c = '000000,011120,111111'

# Calculate the missing value
missing_value = solve_puzzle(input_c)

# Print the final result
print(missing_value)