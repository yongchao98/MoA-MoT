def is_valid_selection(selected_positions, grid):
    # Check if the product condition is satisfied for each row and column
    for i in range(5):
        row_neg_count = sum(1 for j in range(5) if (i, j) in selected_positions and grid[i][j] < 0)
        if row_neg_count % 2 != 0:
            return False
    for j in range(5):
        col_neg_count = sum(1 for i in range(5) if (i, j) in selected_positions and grid[i][j] < 0)
        if col_neg_count % 2 != 0:
            return False
    return True

def find_selection(grid):
    selected_positions = [(0, 4)]
    selected_numbers = [-9]

    def backtrack():
        if len(selected_positions) == 10:
            if is_valid_selection(selected_positions, grid):
                return True
            return False

        for i in range(5):
            for j in range(5):
                if (i, j) not in selected_positions:
                    selected_positions.append((i, j))
                    selected_numbers.append(grid[i][j])
                    if backtrack():
                        return True
                    selected_positions.pop()
                    selected_numbers.pop()
        return False

    backtrack()
    return selected_positions

grid = [
    [-1, 2, -3, 8, -9],
    [2, -8, -5, 1, -3],
    [-10, 6, 1, 1, 10],
    [-2, -2, 6, 9, -6],
    [-3, -10, 2, 6, -9]
]

selected_positions = find_selection(grid)
print(selected_positions)