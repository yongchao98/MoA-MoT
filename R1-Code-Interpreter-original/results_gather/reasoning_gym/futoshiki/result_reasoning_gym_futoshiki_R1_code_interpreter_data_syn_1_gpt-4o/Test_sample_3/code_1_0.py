from itertools import permutations

def is_valid(grid, constraints):
    # Check all rows and columns for unique numbers
    for i in range(8):
        if len(set(grid[i])) != 8 or len(set(row[i] for row in grid)) != 8:
            return False
    
    # Check all constraints
    for (x1, y1, x2, y2, op) in constraints:
        if op == '<' and not (grid[x1][y1] < grid[x2][y2]):
            return False
        if op == '>' and not (grid[x1][y1] > grid[x2][y2]):
            return False
        if op == '∧' and not (grid[x1][y1] < grid[x2][y2]):
            return False
        if op == '∨' and not (grid[x1][y1] > grid[x2][y2]):
            return False
    return True

def solve_futoshiki():
    # Initial grid with given numbers
    grid = [
        [1, 4, 5, 0, 0, 0, 0, 3],
        [2, 0, 0, 0, 7, 0, 0, 0],
        [4, 0, 0, 0, 0, 6, 7, 0],
        [0, 0, 0, 4, 0, 1, 0, 8],
        [0, 0, 6, 2, 0, 3, 0, 0],
        [0, 0, 0, 0, 3, 0, 0, 0],
        [0, 0, 0, 0, 8, 0, 0, 5],
        [3, 6, 0, 0, 4, 7, 0, 0]
    ]
    
    # Constraints in the form (x1, y1, x2, y2, op)
    constraints = [
        (0, 2, 0, 3, '<'), (1, 2, 1, 3, '>'), (1, 4, 1, 5, '>'),
        (2, 0, 1, 0, '∧'), (2, 4, 2, 5, '∧'), (2, 5, 2, 6, '∨'),
        (3, 1, 2, 1, '∧'), (3, 2, 3, 3, '>'), (3, 5, 3, 6, '<'),
        (4, 2, 3, 2, '∧'), (5, 3, 4, 3, '∨'), (5, 6, 5, 7, '>'),
        (6, 1, 5, 1, '∧'), (6, 7, 6, 6, '>'), (7, 1, 7, 2, '<')
    ]
    
    # Try all permutations of numbers 1 to 8 for each row
    for perm in permutations(range(1, 9)):
        grid[0] = [1, 4, 5] + list(perm[:5])
        grid[1] = [2] + list(perm[5:8]) + [7] + list(perm[8:])
        grid[2] = [4] + list(perm[8:11]) + [6, 7] + list(perm[11:])
        grid[3] = list(perm[11:14]) + [4] + list(perm[14:16]) + [1, 8]
        grid[4] = list(perm[16:18]) + [6, 2] + list(perm[18:20]) + [3] + list(perm[20:])
        grid[5] = list(perm[20:23]) + [3] + list(perm[23:])
        grid[6] = list(perm[23:27]) + [8] + list(perm[27:29]) + [5]
        grid[7] = [3, 6] + list(perm[29:31]) + [4, 7] + list(perm[31:])
        
        if is_valid(grid, constraints):
            for row in grid:
                print(' '.join(map(str, row)))
            return

solve_futoshiki()