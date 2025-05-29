from itertools import permutations

def is_valid(grid, constraints):
    # Check all rows and columns for unique numbers
    for i in range(7):
        if len(set(grid[i])) != 7 or len(set(row[i] for row in grid)) != 7:
            return False
    
    # Check horizontal constraints
    for (r, c1, c2, op) in constraints['horizontal']:
        if op == '<' and not (grid[r][c1] < grid[r][c2]):
            return False
        if op == '>' and not (grid[r][c1] > grid[r][c2]):
            return False
    
    # Check vertical constraints
    for (c, r1, r2, op) in constraints['vertical']:
        if op == '∧' and not (grid[r1][c] < grid[r2][c]):
            return False
        if op == '∨' and not (grid[r1][c] > grid[r2][c]):
            return False
    
    return True

def solve_futoshiki():
    # Initial grid with known values
    grid = [
        [0, 5, 0, 0, 0, 0, 0],
        [4, 0, 0, 3, 0, 0, 0],
        [0, 3, 0, 5, 0, 1, 7],
        [0, 0, 0, 0, 0, 3, 0],
        [1, 0, 0, 0, 5, 0, 0],
        [0, 0, 0, 0, 0, 0, 4],
        [0, 0, 6, 2, 0, 0, 0]
    ]
    
    # Constraints
    constraints = {
        'horizontal': [
            (1, 0, 1, '<'), (1, 2, 3, '>'),
            (2, 3, 4, '<'), (2, 4, 5, '>'),
            (3, 3, 4, '<'),
            (4, 1, 2, '>'), (4, 4, 5, '>'),
            (5, 0, 1, '>'),
            (6, 1, 2, '<')
        ],
        'vertical': [
            (1, 0, 1, '∧'),
            (3, 1, 2, '∧'),
            (4, 2, 3, '∨'),
            (5, 3, 4, '∨'),
            (6, 4, 5, '∧')
        ]
    }
    
    # Try all permutations for each row
    for perm1 in permutations(range(1, 8)):
        grid[0] = list(perm1)
        for perm2 in permutations(range(1, 8)):
            if perm2[0] == 4 and perm2[3] == 3:
                grid[1] = list(perm2)
                for perm3 in permutations(range(1, 8)):
                    if perm3[1] == 3 and perm3[3] == 5 and perm3[5] == 1 and perm3[6] == 7:
                        grid[2] = list(perm3)
                        for perm4 in permutations(range(1, 8)):
                            if perm4[5] == 3:
                                grid[3] = list(perm4)
                                for perm5 in permutations(range(1, 8)):
                                    if perm5[0] == 1 and perm5[4] == 5:
                                        grid[4] = list(perm5)
                                        for perm6 in permutations(range(1, 8)):
                                            if perm6[6] == 4:
                                                grid[5] = list(perm6)
                                                for perm7 in permutations(range(1, 8)):
                                                    if perm7[2] == 6 and perm7[3] == 2:
                                                        grid[6] = list(perm7)
                                                        if is_valid(grid, constraints):
                                                            for row in grid:
                                                                print(row)
                                                            return

solve_futoshiki()