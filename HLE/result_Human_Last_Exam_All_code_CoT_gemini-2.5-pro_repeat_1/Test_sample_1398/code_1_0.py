import numpy as np

def solve_grid(grid):
    """
    Applies the transformation rule to solve the grid puzzle.
    The rule is: the output cell is 1 if the number of its cardinal neighbors
    is greater than or equal to the number of its diagonal neighbors.
    Otherwise, the output is 0.
    """
    rows, cols = grid.shape
    output_grid = np.zeros_like(grid)

    for r in range(rows):
        for c in range(cols):
            diagonal_sum = 0
            cardinal_sum = 0

            # Iterate through all 8 neighbors
            for dr in [-1, 0, 1]:
                for dc in [-1, 0, 1]:
                    if dr == 0 and dc == 0:
                        continue

                    nr, nc = r + dr, c + dc

                    if 0 <= nr < rows and 0 <= nc < cols:
                        # Check if it's a diagonal neighbor
                        if abs(dr) == 1 and abs(dc) == 1:
                            diagonal_sum += grid[nr, nc]
                        # Otherwise it's a cardinal neighbor
                        else:
                            cardinal_sum += grid[nr, nc]

            # The final equation for each cell's output
            # output = 1 if cardinal_sum >= diagonal_sum else 0
            # For the printout requirement:
            print(f"Cell ({r},{c}): Cardinal Sum = {cardinal_sum}, Diagonal Sum = {diagonal_sum}. Result: {1 if cardinal_sum >= diagonal_sum else 0}")


            if cardinal_sum >= diagonal_sum:
                output_grid[r, c] = 1
            else:
                output_grid[r, c] = 0

    return output_grid

# The 6th test input grid
test_input_grid = np.array([
    [0, 1, 1, 0, 0],
    [1, 0, 0, 1, 0],
    [1, 0, 1, 0, 1],
    [0, 1, 0, 0, 1],
    [1, 0, 1, 1, 0]
])

# Solve the grid and get the result
# The printout of equations will be generated here
final_output_grid = solve_grid(test_input_grid)

# Flatten the 2D array into a 1D list of strings
result_list = final_output_grid.flatten().astype(str).tolist()

# Join the list into a single string
final_answer_string = "".join(result_list)
print("\nFinal Output Grid as a binary string:")
print(final_answer_string)