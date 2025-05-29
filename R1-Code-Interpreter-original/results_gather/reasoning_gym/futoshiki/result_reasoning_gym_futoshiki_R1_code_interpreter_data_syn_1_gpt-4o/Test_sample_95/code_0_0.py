from itertools import permutations

def is_valid_solution(grid):
    # Check rows and columns for duplicates
    for i in range(4):
        if len(set(grid[i])) != 4 or len(set(row[i] for row in grid)) != 4:
            return False
    return True

def solve_futoshiki():
    # Initial grid with given numbers
    grid = [
        [3, 0, 0, 1],
        [0, 0, 2, 0],
        [0, 0, 0, 3],
        [0, 0, 0, 0]
    ]
    
    # Possible numbers to fill in
    numbers = [1, 2, 3, 4]
    
    # Try all permutations of numbers for each row
    for perm1 in permutations(numbers):
        grid[0] = [3 if x == 0 else x for x in perm1]
        for perm2 in permutations(numbers):
            grid[1] = [2 if x == 0 else x for x in perm2]
            for perm3 in permutations(numbers):
                grid[2] = [3 if x == 0 else x for x in perm3]
                for perm4 in permutations(numbers):
                    grid[3] = perm4
                    if is_valid_solution(grid):
                        return grid

solution = solve_futoshiki()
for row in solution:
    print(' '.join(map(str, row)))