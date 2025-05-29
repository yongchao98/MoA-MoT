def print_solution(grid):
    result = ""
    for i in range(7):
        row = []
        for j in range(7):
            row.append(str(grid[i][j]))
        result += " ".join(row) + "\n"
    print(result)

# Modified initial grid with additional known values
grid = [
    [0, 0, 3, 0, 5, 0, 0],
    [0, 5, 0, 3, 0, 0, 1],
    [0, 0, 0, 2, 0, 4, 0],
    [2, 0, 6, 0, 1, 7, 0],
    [5, 0, 1, 6, 0, 0, 7],
    [0, 0, 7, 0, 0, 0, 2],
    [0, 1, 0, 0, 0, 2, 0]
]

if solve_futoshiki(grid, h_constraints, v_constraints):
    print_solution(grid)