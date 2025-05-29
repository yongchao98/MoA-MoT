def find_rectangles(grid):
    rows, cols = len(grid), len(grid[0])
    visited = set()
    
    def is_corner(i, j):
        if grid[i][j] not in '#█':
            return False
        # Check if it's a top-left corner
        left = j == 0 or grid[i][j-1] not in '#█'
        top = i == 0 or grid[i-1][j] not in '#█'
        return left and top
    
    def validate_rectangle(top, left):
        if (top, left) in visited:
            return False
            
        # Find right boundary
        right = left
        while right < cols and grid[top][right] in '#█':
            right += 1
        right -= 1
        
        if right == left:
            return False
            
        # Find bottom boundary
        bottom = top
        while bottom < rows and grid[bottom][left] in '#█':
            bottom += 1
        bottom -= 1
        
        if bottom == top:
            return False
            
        # Quick validation
        for i in range(top, bottom + 1):
            if grid[i][left] not in '#█' or grid[i][right] not in '#█':
                return False
        for j in range(left, right + 1):
            if grid[top][j] not in '#█' or grid[bottom][j] not in '#█':
                return False
                
        # Mark as visited
        for i in range(top, bottom + 1):
            for j in range(left, right + 1):
                if grid[i][j] in '#█':
                    visited.add((i, j))
                    
        return True
    
    count = 0
    # Only check positions that could be corners
    for i in range(rows):
        for j in range(cols):
            if is_corner(i, j) and validate_rectangle(i, j):
                count += 1
                
    return count

# Read input efficiently
grid = []
while True:
    try:
        line = input()
        if not line and grid:
            break
        grid.append(line)
    except EOFError:
        break

print(f"<<<{find_rectangles(grid)}>>>")