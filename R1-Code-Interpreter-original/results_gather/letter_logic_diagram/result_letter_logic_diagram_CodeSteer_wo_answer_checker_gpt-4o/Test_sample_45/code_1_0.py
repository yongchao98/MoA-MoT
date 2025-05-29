import heapq

def solve_puzzle(grid):
    def is_valid(grid, row, col, letter):
        # Check if the letter can be placed at grid[row][col]
        for i in range(7):
            if grid[row][i] == letter or grid[i][col] == letter:
                return False
        return True

    def update_possible_values(grid, possible_values, row, col, letter):
        # Update possible values after placing a letter
        for i in range(7):
            possible_values[row][i].discard(letter)
            possible_values[i][col].discard(letter)

    def solve(grid, possible_values):
        # Use a priority queue to prioritize cells with fewer possible values
        pq = []
        for i in range(7):
            for j in range(7):
                if grid[i][j] == '':
                    heapq.heappush(pq, (len(possible_values[i][j]), i, j))

        if not pq:
            return True

        _, row, col = heapq.heappop(pq)

        for letter in possible_values[row][col]:
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                new_possible_values = [row[:] for row in possible_values]
                update_possible_values(grid, new_possible_values, row, col, letter)
                if solve(grid, new_possible_values):
                    return True
                grid[row][col] = ''

        return False

    # Initialize possible values for each cell
    possible_values = [[set('abcdefg') for _ in range(7)] for _ in range(7)]
    for i in range(7):
        for j in range(7):
            if grid[i][j] != '':
                possible_values[i][j] = set(grid[i][j])
                for k in range(7):
                    possible_values[i][k].discard(grid[i][j])
                    possible_values[k][j].discard(grid[i][j])

    # Choose a letter for the minor diagonal that doesn't conflict
    minor_diagonal_letter = None
    for letter in 'abcdefg':
        if all(grid[i][6-i] in ('', letter) for i in range(7)):
            minor_diagonal_letter = letter
            break

    if minor_diagonal_letter is None:
        print("No solution found")
        return

    for i in range(7):
        grid[i][6-i] = minor_diagonal_letter
        possible_values[i][6-i] = set(minor_diagonal_letter)

    if solve(grid, possible_values):
        for row in grid:
            print(','.join(row))
    else:
        print("No solution found")

# Initial grid with empty cells as ''
grid = [
    ['c', '', 'a', 'e', '', 'b', ''],
    ['', 'a', '', 'f', '', '', ''],
    ['', '', 'f', 'b', 'g', '', ''],
    ['', '', '', '', '', '', ''],
    ['', '', 'g', '', '', '', ''],
    ['b', '', '', 'd', 'a', 'e', 'f'],
    ['', '', 'd', 'a', '', 'f', 'b']
]

solve_puzzle(grid)