from itertools import permutations

def is_valid_row(row, fixed_positions_in_row):
    # Check if row matches fixed positions
    for pos, letter in fixed_positions_in_row:
        if row[pos] != letter:
            return False
    return True

def get_valid_rows(fixed_positions_in_row, diagonal_pos=None):
    valid_rows = []
    for p in permutations('abcdefg'):
        row = list(p)
        if diagonal_pos is not None:
            if row[diagonal_pos] != 'g':
                continue
        if is_valid_row(row, fixed_positions_in_row):
            valid_rows.append(row)
    return valid_rows

# Initial grid
initial_grid = [
    ['', '', 'a', 'c', '', '', 'g'],
    ['', 'a', 'c', 'b', '', 'g', ''],
    ['a', '', '', '', 'g', '', ''],
    ['c', 'b', 'f', '', 'e', 'd', 'a'],
    ['b', 'f', '', '', '', 'a', ''],
    ['', '', 'e', '', '', 'c', ''],
    ['', '', 'd', 'a', '', 'b', 'f']
]

# Get fixed positions for each row
fixed_positions = []
for i in range(7):
    row_fixed = []
    for j in range(7):
        if initial_grid[i][j] != '':
            row_fixed.append((j, initial_grid[i][j]))
    fixed_positions.append(row_fixed)

# Generate valid rows for each position
valid_rows = []
for i in range(7):
    diagonal_pos = 6 - i  # Position of diagonal in this row
    rows = get_valid_rows(fixed_positions[i], diagonal_pos)
    valid_rows.append(rows)

def check_columns(solution):
    # Check if each column has all different letters
    for j in range(7):
        column = [solution[i][j] for i in range(7)]
        if len(set(column)) != 7:
            return False
    return True

def solve(valid_rows, row=0, current_solution=[]):
    if row == 7:
        if check_columns(current_solution):
            return current_solution
        return None
    
    for valid_row in valid_rows[row]:
        current_solution.append(valid_row)
        result = solve(valid_rows, row + 1, current_solution)
        if result is not None:
            return result
        current_solution.pop()
    
    return None

solution = solve(valid_rows)

if solution:
    for row in solution:
        print(','.join(row))
else:
    print("No solution exists")