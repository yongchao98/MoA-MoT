def print_grid(grid):
    for row in grid:
        print(row)

def is_valid(grid, row, col, num):
    # Check row
    for x in range(8):
        if grid[row][x] == num:
            return False
    
    # Check column
    for x in range(8):
        if grid[x][col] == num:
            return False
    
    # Check horizontal constraints (<)
    h_less = {
        (1,6): True,   # row 1, col 6: num < right number
        (3,6): True,   # row 3, col 6: num < right number
        (4,7): True,   # row 4, col 7: num < right number
        (7,4): True    # row 7, col 4: num < right number
    }
    
    # Check vertical constraints (∨=down, ∧=up)
    v_constraints = {
        (0,1): 'down',  # row 0, col 1: top < bottom
        (1,1): 'up',    # row 1, col 1: top > bottom
        (1,7): 'up',    # row 1, col 7: top > bottom
        (2,7): 'down',  # row 2, col 7: top < bottom
        (4,4): 'up',    # row 4, col 4: top > bottom
        (5,3): 'down'   # row 5, col 3: top < bottom
    }
    
    # Check horizontal constraints
    if (row, col) in h_less:
        if grid[row][col+1] != 0 and num >= grid[row][col+1]:
            return False
    
    # Check reverse horizontal constraints
    for (r, c), _ in h_less.items():
        if row == r and col == c + 1 and grid[r][c] != 0:
            if grid[r][c] >= num:
                return False
    
    # Check vertical constraints
    for (r, c), direction in v_constraints.items():
        if col == c:
            if row == r:
                if direction == 'down' and grid[r+1][c] != 0:
                    if num >= grid[r+1][c]:
                        return False
                if direction == 'up' and grid[r+1][c] != 0:
                    if num <= grid[r+1][c]:
                        return False
            elif row == r + 1:
                if direction == 'down' and grid[r][c] != 0:
                    if grid[r][c] >= num:
                        return False
                if direction == 'up' and grid[r][c] != 0:
                    if grid[r][c] <= num:
                        return False
    
    return True

def find_empty(grid):
    for i in range(8):
        for j in range(8):
            if grid[i][j] == 0:
                return (i, j)
    return None

def solve(grid, depth=0):
    pos = find_empty(grid)
    if not pos:
        return True
    
    row, col = pos
    for num in range(1, 9):
        if is_valid(grid, row, col, num):
            grid[row][col] = num
            if solve(grid, depth + 1):
                return True
            grid[row][col] = 0
    
    return False

# Initialize grid
grid = [
    [0,2,1,0,5,7,0,0],
    [0,1,0,0,0,3,0,0],
    [4,0,5,0,0,0,0,0],
    [0,7,4,6,0,0,0,5],
    [0,3,7,0,0,0,1,0],
    [0,0,0,2,0,0,5,0],
    [6,0,0,0,2,4,0,0],
    [0,0,0,0,4,1,7,0]
]

# Verify initial setup
def verify_initial_setup():
    for i in range(8):
        for j in range(8):
            if grid[i][j] != 0:
                temp = grid[i][j]
                grid[i][j] = 0
                if not is_valid(grid, i, j, temp):
                    print(f"Contradiction found at position ({i},{j}) with value {temp}")
                    return False
                grid[i][j] = temp
    return True

print("Verifying initial setup...")
if verify_initial_setup():
    print("Initial setup is valid, solving puzzle...")
    if solve(grid):
        print("\nSolution found:")
        print("<<<")
        for i in range(8):
            row = ""
            for j in range(7):
                row += str(grid[i][j])
                if (i == 1 and j == 6) or (i == 3 and j == 6) or (i == 4 and j == 7) or (i == 7 and j == 4):
                    row += " < "
                else:
                    row += "   "
            row += str(grid[i][7])
            print(row)
            if i < 7:
                v_row = ["   "] * 8
                if i == 0: v_row[1] = " ∨ "
                if i == 1: v_row[1] = v_row[7] = " ∧ "
                if i == 2: v_row[7] = " ∨ "
                if i == 4: v_row[4] = " ∧ "
                if i == 5: v_row[3] = " ∨ "
                print("".join(v_row))
        print(">>>")
    else:
        print("No solution exists")
else:
    print("Initial setup contains contradictions")