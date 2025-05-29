def count_rectangles(grid):
    if not grid or not grid[0]:
        return 0
    
    rows, cols = len(grid), len(grid[0])
    visited = set()
    count = 0
    
    def is_border(i, j):
        return grid[i][j] in ['#', '█']
    
    def find_rectangle(i, j):
        if (i, j) in visited or not is_border(i, j):
            return False
            
        # Find right edge
        right = j
        while right < cols and is_border(i, right):
            right += 1
        right -= 1
        
        # Find bottom edge
        bottom = i
        while bottom < rows and is_border(bottom, j):
            bottom += 1
        bottom -= 1
        
        # Verify rectangle
        for r in range(i, bottom + 1):
            if not is_border(r, j) or not is_border(r, right):
                return False
        for c in range(j, right + 1):
            if not is_border(i, c) or not is_border(bottom, c):
                return False
        
        # Mark as visited
        for r in range(i, bottom + 1):
            for c in range(j, right + 1):
                if is_border(r, c):
                    visited.add((r, c))
        
        return True
    
    # Scan for top-left corners
    for i in range(rows):
        for j in range(cols):
            if is_border(i, j):
                # Check if it's a potential top-left corner
                if (i == 0 or not is_border(i-1, j)) and (j == 0 or not is_border(i, j-1)):
                    if find_rectangle(i, j):
                        count += 1
    
    return count

# Read the grid
grid = []
try:
    while True:
        line = input()
        if not line.strip():  # Skip empty lines
            continue
        grid.append(line)
except EOFError:
    pass

print(count_rectangles(grid))