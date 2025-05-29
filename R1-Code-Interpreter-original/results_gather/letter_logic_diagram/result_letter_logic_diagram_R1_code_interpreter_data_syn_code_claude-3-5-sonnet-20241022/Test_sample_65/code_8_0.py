def print_grid(grid):
    result = ""
    for row in grid:
        result += ",".join(row) + "\n"
    print("<<<")
    print(result.strip())
    print(">>>")

def is_valid(grid, row, col, letter, initial_grid):
    # Check if this position was pre-filled
    if initial_grid[row][col] != '' and initial_grid[row][col] != letter:
        return False

    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False

    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False

    return True

def solve():
    initial_grid = [
        ['a','','d','','b','g','c'],
        ['f','d','','','','','a'],
        ['d','','','','','',''],
        ['','','','','','f',''],
        ['','g','','a','f','','e'],
        ['','c','','f','','e',''],
        ['','a','','','e','','']
    ]

    grid = [[''] * 7 for _ in range(7)]
    
    # First, copy initial values
    for i in range(7):
        for j in range(7):
            if initial_grid[i][j] != '':
                grid[i][j] = initial_grid[i][j]

    def find_minor_diagonal_options():
        # Find possible letters for minor diagonal based on initial constraints
        possible_letters = set('abcdefg')
        for i in range(7):
            if initial_grid[i][6-i] != '':
                possible_letters = {initial_grid[i][6-i]}
                break
        return list(possible_letters)

    def try_with_minor_diagonal(minor_letter):
        # Fill minor diagonal
        for i in range(7):
            if grid[i][6-i] == '':
                grid[i][6-i] = minor_letter
            elif grid[i][6-i] != minor_letter:
                return False

        def backtrack(pos=0):
            if pos == 49:
                return True

            row = pos // 7
            col = pos % 7

            # Skip pre-filled cells and minor diagonal
            if grid[row][col] != '':
                return backtrack(pos + 1)

            for letter in 'abcdefg':
                if is_valid(grid, row, col, letter, initial_grid):
                    grid[row][col] = letter
                    if backtrack(pos + 1):
                        return True
                    grid[row][col] = ''

            return False

        return backtrack(0)

    # Try each possible minor diagonal letter
    minor_diag_options = find_minor_diagonal_options()
    for minor_letter in minor_diag_options:
        # Create a fresh grid with initial values
        grid = [[''] * 7 for _ in range(7)]
        for i in range(7):
            for j in range(7):
                if initial_grid[i][j] != '':
                    grid[i][j] = initial_grid[i][j]
        
        if try_with_minor_diagonal(minor_letter):
            print_grid(grid)
            return True

    print("No solution found")
    return False

solve()