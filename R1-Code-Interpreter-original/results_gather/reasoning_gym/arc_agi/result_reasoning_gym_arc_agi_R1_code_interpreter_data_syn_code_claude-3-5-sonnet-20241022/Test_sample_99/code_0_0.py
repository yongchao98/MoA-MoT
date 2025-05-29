def find_special_numbers(grid):
    rows, cols = len(grid), len(grid[0])
    special_positions = []
    
    for i in range(rows):
        for j in range(cols):
            # Check surrounding cells (if they exist)
            current = grid[i][j]
            neighbors = []
            for di in [-1, 0, 1]:
                for dj in [-1, 0, 1]:
                    if 0 <= i + di < rows and 0 <= j + dj < cols:
                        if di != 0 or dj != 0:  # Don't include the current cell
                            neighbors.append(grid[i + di][j + dj])
            
            # If current number is different from most of its neighbors
            if neighbors:
                most_common = max(set(neighbors), key=neighbors.count)
                if current != most_common:
                    special_positions.append((i, j))
    
    return special_positions

def get_replacement(num):
    replacements = {5: 4, 9: 8, 3: 5, 4: 9, 8: 2}
    return replacements.get(num, num)

def apply_pattern(input_grid):
    rows, cols = len(input_grid), len(input_grid[0])
    output_grid = [[input_grid[i][j] for j in range(cols)] for i in range(rows)]
    
    special_positions = find_special_numbers(input_grid)
    
    for pos_i, pos_j in special_positions:
        original_num = input_grid[pos_i][pos_j]
        replacement = get_replacement(original_num)
        
        # Replace the original position
        output_grid[pos_i][pos_j] = replacement
        
        # Apply diagonal pattern
        for direction in range(4):  # 4 diagonal directions
            di = -1 if direction < 2 else 1
            dj = -1 if direction % 2 == 0 else 1
            
            i, j = pos_i, pos_j
            steps = 0
            while 0 <= i + di < rows and 0 <= j + dj < cols and steps < 3:
                i += di
                j += dj
                if input_grid[i][j] == input_grid[pos_i][pos_j]:
                    break
                output_grid[i][j] = replacement
                steps += 1
    
    return output_grid

# Parse input grid
input_str = """4 4 4 4 5 5 5 4 4 4 4 5 5 5 5 4 4
4 4 4 4 5 5 5 4 4 4 4 5 5 5 5 4 4
4 4 4 4 5 5 5 4 4 4 4 5 5 5 5 4 4
5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5
5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5
5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5
5 5 5 5 9 5 5 5 5 5 5 5 5 5 5 5 5
5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5
5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5
4 4 4 4 5 5 5 4 4 4 4 5 5 5 5 4 4
2 4 4 4 5 5 5 4 4 4 4 5 5 5 5 4 4
4 4 4 4 5 5 5 4 4 4 4 5 5 5 5 4 4
4 4 4 4 5 5 5 4 4 4 4 5 5 5 5 4 4
4 4 4 4 5 5 5 4 4 4 4 5 5 5 5 4 4
4 4 4 4 5 5 5 4 4 4 4 5 9 5 5 4 4
4 4 4 4 5 5 5 4 4 4 4 5 5 5 5 4 4
4 4 4 4 5 5 5 4 4 4 4 5 5 5 5 4 4
4 4 4 4 5 5 5 4 4 4 4 5 5 5 5 4 4
4 4 4 4 5 5 5 4 4 4 4 5 5 5 5 4 4"""

input_grid = [[int(x) for x in row.split()] for row in input_str.split('\n')]
output_grid = apply_pattern(input_grid)

# Format output
result = '\n'.join(' '.join(str(x) for x in row) for row in output_grid)
print(result)