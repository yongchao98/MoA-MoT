from itertools import permutations

def is_valid_solution(grid):
    # Check rows and columns for unique numbers
    for i in range(5):
        if len(set(grid[i])) != 5:
            return False
        if len(set(row[i] for row in grid)) != 5:
            return False
    
    # Check the inequality constraint
    if not (grid[4][4] > grid[4][3]):
        return False
    
    return True

def solve_futoshiki():
    # Initial grid with known numbers
    grid = [
        [0, 1, 0, 0, 3],
        [0, 0, 0, 3, 0],
        [0, 2, 0, 0, 5],
        [3, 0, 0, 1, 0],
        [0, 0, 4, 5, 0]
    ]
    
    # Possible numbers to fill in
    numbers = [1, 2, 3, 4, 5]
    
    # Try all permutations for each row
    for perm1 in permutations(numbers):
        if perm1[1] == 1 and perm1[4] == 3:
            grid[0] = list(perm1)
            for perm2 in permutations(numbers):
                if perm2[3] == 3:
                    grid[1] = list(perm2)
                    for perm3 in permutations(numbers):
                        if perm3[1] == 2 and perm3[4] == 5:
                            grid[2] = list(perm3)
                            for perm4 in permutations(numbers):
                                if perm4[0] == 3 and perm4[3] == 1:
                                    grid[3] = list(perm4)
                                    for perm5 in permutations(numbers):
                                        if perm5[2] == 4 and perm5[3] == 5 and perm5[4] > perm5[3]:
                                            grid[4] = list(perm5)
                                            if is_valid_solution(grid):
                                                return grid

solution = solve_futoshiki()
if solution:
    for row in solution:
        print(row)
else:
    print("No solution found")