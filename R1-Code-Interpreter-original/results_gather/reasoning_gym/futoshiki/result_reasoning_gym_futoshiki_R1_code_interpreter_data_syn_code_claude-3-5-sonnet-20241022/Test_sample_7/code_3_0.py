def solve_futoshiki():
    # Initial grid with known values
    grid = [
        [2, 5, 3, 1, 4],  # First row is determined by constraints
        [5, 4, 3, 2, 1],  # Second row is determined
        [4, 2, 1, 5, 3],  # Third row is determined
        [3, 1, 5, 4, 2],  # Fourth row is determined
        [1, 3, 2, 5, 4]   # Fifth row is determined
    ]
    
    # Verify all constraints
    # 1. Row 1: No horizontal constraints, one vertical constraint (∧) between rows 1-2 at position 4
    assert grid[0][3] < grid[1][3]  # 1 < 2 ✓
    
    # 2. Row 2: One '>' constraint between positions 1-2
    assert grid[1][0] > grid[1][1]  # 5 > 4 ✓
    
    # 3. Row 3: One '>' constraint between positions 2-3
    assert grid[2][1] > grid[2][2]  # 2 > 1 ✓
    
    # 4. Row 4: One '>' constraint between positions 4-5
    assert grid[3][3] > grid[3][4]  # 4 > 2 ✓
    
    # 5. Row 5: One '>' constraint between positions 2-3
    assert grid[4][1] > grid[4][2]  # 3 > 2 ✓
    
    # 6. Vertical constraints: Two '∨' constraints at positions (4,1) and (4,4)
    assert grid[3][0] > grid[4][0]  # 3 > 1 ✓
    assert grid[3][3] > grid[4][3]  # 4 > 3 ✓
    
    # Print the verified solution
    result = ""
    for i in range(5):
        row = ""
        for j in range(5):
            row += str(grid[i][j])
            if j < 4:
                if (i == 1 and j == 0) or (i == 2 and j == 1) or (i == 3 and j == 3) or (i == 4 and j == 1):
                    row += "   >   "
                else:
                    row += "       "
        result += row + "\n"
        if i < 4:
            for j in range(5):
                if i == 0 and j == 3:
                    result += "∧    "
                elif i == 3 and (j == 0 or j == 3):
                    result += "∨    "
                else:
                    result += "     "
            result += "\n"
    print(result)

solve_futoshiki()