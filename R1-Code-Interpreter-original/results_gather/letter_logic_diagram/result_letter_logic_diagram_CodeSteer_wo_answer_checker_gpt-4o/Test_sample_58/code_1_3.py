def solve_puzzle(grid):
    from collections import defaultdict, deque

    def initialize_possibilities(grid):
        possibilities = defaultdict(lambda: set('abcdefg'))
        for r in range(7):
            for c in range(7):
                if grid[r][c] != '':
                    possibilities[(r, c)] = {grid[r][c]}
        return possibilities

    def ac3(possibilities):
        queue = deque([(r, c) for r in range(7) for c in range(7)])
        while queue:
            r, c = queue.popleft()
            if len(possibilities[(r, c)]) == 1:
                letter = next(iter(possibilities[(r, c)]))
                for i in range(7):
                    if i != c and letter in possibilities[(r, i)]:
                        possibilities[(r, i)].remove(letter)
                        if len(possibilities[(r, i)]) == 0:
                            return False
                        if len(possibilities[(r, i)]) == 1:
                            queue.append((r, i))
                    if i != r and letter in possibilities[(i, c)]:
                        possibilities[(i, c)].remove(letter)
                        if len(possibilities[(i, c)]) == 0:
                            return False
                        if len(possibilities[(i, c)]) == 1:
                            queue.append((i, c))
        return True

    def fill_grid(grid, possibilities, row=0, col=0):
        if row == 7:
            return True
        if col == 7:
            return fill_grid(grid, possibilities, row + 1, 0)
        if grid[row][col] != '':
            return fill_grid(grid, possibilities, row, col + 1)

        for letter in possibilities[(row, col)]:
            grid[row][col] = letter
            new_possibilities = defaultdict(set, {k: v.copy() for k, v in possibilities.items()})
            new_possibilities[(row, col)] = {letter}
            if ac3(new_possibilities) and fill_grid(grid, new_possibilities, row, col + 1):
                return True
            grid[row][col] = ''
        return False

    # Initialize possibilities
    possibilities = initialize_possibilities(grid)
    if not ac3(possibilities):
        return None

    # Fill the minor diagonal with the same letter
    for letter in 'abcdefg':
        temp_grid = [row[:] for row in grid]
        temp_possibilities = defaultdict(set, {k: v.copy() for k, v in possibilities.items()})
        for i in range(7):
            temp_grid[i][6-i] = letter
            temp_possibilities[(i, 6-i)] = {letter}
        if ac3(temp_possibilities) and fill_grid(temp_grid, temp_possibilities):
            return temp_grid

    return None

# Initial grid setup
grid = [
    ['d', '', '', '', '', 'f', ''],
    ['', 'b', '', '', 'f', 'a', 'd'],
    ['b', '', '', 'f', '', '', 'd'],
    ['', '', 'f', '', '', '', ''],
    ['c', '', '', '', 'e', '', 'g'],
    ['', 'a', 'd', '', '', '', 'c'],
    ['', '', 'e', '', '', 'c', 'f']
]

solution = solve_puzzle(grid)
if solution:
    print("<<<")
    for row in solution:
        print(','.join(row))
    print(">>>")
else:
    print("No solution found.")