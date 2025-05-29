def print_solution(grid):
    result = ""
    for i in range(8):
        for j in range(8):
            result += str(grid[i][j])
            result += " " if j < 7 else "\n"
    print(result)

def is_safe(grid, row, col, num, h_constraints, v_constraints):
    # Check row
    if num in grid[row]:
        return False
    
    # Check column
    if num in [grid[i][col] for i in range(8)]:
        return False
    
    # Check horizontal constraints
    for r, c1, c2 in h_constraints:
        if r == row:
            if col == c1 and grid[r][c2] != 0:  # Left side of '<'
                if num >= grid[r][c2]:
                    return False
            if col == c2 and grid[r][c1] != 0:  # Right side of '<'
                if num <= grid[r][c1]:
                    return False

    # Check vertical constraints
    for r1, r2, c in v_constraints:
        if col == c:
            if row == r1 and grid[r2][c] != 0:  # Upper side of '∨'
                if num >= grid[r2][c]:
                    return False
            if row == r2 and grid[r1][c] != 0:  # Lower side of '∨'
                if num <= grid[r1][c]:
                    return False

    return True

def find_empty(grid):
    min_possibilities = 9
    best_pos = None
    
    for i in range(8):
        for j in range(8):
            if grid[i][j] == 0:
                count = sum(1 for num in range(1, 9) if is_safe(grid, i, j, num, h_constraints, v_constraints))
                if count < min_possibilities:
                    min_possibilities = count
                    best_pos = (i, j)
                if count == 1:  # If we find a cell with only one possibility, return it immediately
                    return (i, j)
    
    return best_pos

def solve(grid, h_constraints, v_constraints):
    pos = find_empty(grid)
    if not pos:
        return True
    
    row, col = pos
    
    # Try numbers in an order that's more likely to succeed
    numbers = list(range(1, 9))
    # If there are constraints, try to order numbers more intelligently
    if any(r == row and (c1 == col or c2 == col) for r, c1, c2 in h_constraints):
        numbers.sort()  # For '<' constraints, try smaller numbers first
    elif any(c == col and (r1 == row or r2 == row) for r1, r2, c in v_constraints):
        numbers.sort()  # For '∨' constraints, try smaller numbers first
    
    for num in numbers:
        if is_safe(grid, row, col, num, h_constraints, v_constraints):
            grid[row][col] = num
            if solve(grid, h_constraints, v_constraints):
                return True
            grid[row][col] = 0
    
    return False

# Initial grid
grid = [
    [0, 1, 5, 0, 7, 0, 2, 8],
    [5, 7, 0, 1, 0, 0, 0, 3],
    [0, 0, 0, 0, 0, 8, 3, 0],
    [2, 0, 6, 0, 0, 4, 0, 0],
    [1, 0, 0, 6, 0, 0, 0, 0],
    [0, 3, 0, 0, 0, 0, 4, 1],
    [8, 2, 0, 0, 0, 5, 0, 6],
    [4, 8, 0, 0, 0, 0, 0, 5]
]

# Define constraints
h_constraints = [(6, 2, 3)]  # Row 7, col 3 < col 4
v_constraints = [(3, 4, 7), (6, 7, 3)]  # (Row 4, col 8) < (Row 5, col 8) and (Row 7, col 4) < (Row 8, col 4)

if solve(grid, h_constraints, v_constraints):
    print_solution(grid)
    
    # Format the answer in the required format
    answer = ""
    for i in range(8):
        for j in range(8):
            answer += str(grid[i][j])
            if j < 7:
                if (i, j, j+1) in h_constraints:
                    answer += "<"
                else:
                    answer += " "
        answer += "\n"
        if i < 7:
            for j in range(8):
                if (i, i+1, j) in v_constraints:
                    answer += "∨ "
                else:
                    answer += "  "
            answer += "\n"
    
    print("\n<<<" + answer.rstrip() + ">>>")
else:
    print("No solution exists")