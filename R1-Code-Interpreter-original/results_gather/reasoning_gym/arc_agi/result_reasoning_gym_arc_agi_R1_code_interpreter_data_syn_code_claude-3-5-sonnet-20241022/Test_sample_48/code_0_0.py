import numpy as np
from scipy.ndimage import label

def find_patterns(grid):
    # Convert input string to numpy array
    grid = np.array([[int(x) for x in row.split()] for row in grid.strip().split('\n')])
    
    # Create a binary mask for value 1
    mask = (grid == 1)
    
    # Find connected components
    labeled_array, num_features = label(mask)
    
    # For each component, check if it forms a valid pattern
    result = grid.copy()
    for i in range(1, num_features + 1):
        component = (labeled_array == i)
        if is_valid_pattern(component):
            result[component] = 5
    
    return result

def is_valid_pattern(component):
    # Get the bounding box of the component
    rows = np.any(component, axis=1)
    cols = np.any(component, axis=0)
    rmin, rmax = np.where(rows)[0][[0, -1]]
    cmin, cmax = np.where(cols)[0][[0, -1]]
    
    # Extract the component's bounding box
    box = component[rmin:rmax+1, cmin:cmax+1]
    
    # Count the number of '1's in the component
    num_ones = np.sum(component)
    
    # Check if it forms a valid pattern (connected and forms a recognizable shape)
    # Valid patterns typically have at least 4 cells and form a connected shape
    return num_ones >= 4 and is_connected(box)

def is_connected(pattern):
    # Check if all '1's are connected
    h, w = pattern.shape
    visited = np.zeros_like(pattern, dtype=bool)
    
    def dfs(i, j):
        if i < 0 or i >= h or j < 0 or j >= w or visited[i,j] or not pattern[i,j]:
            return 0
        visited[i,j] = True
        count = 1
        for di, dj in [(0,1), (1,0), (0,-1), (-1,0)]:
            count += dfs(i+di, j+dj)
        return count
    
    # Find first '1'
    for i in range(h):
        for j in range(w):
            if pattern[i,j]:
                return dfs(i,j) == np.sum(pattern)
    return False

# Test input
test_input = """6 6 6 6 6 6 6 6 6 6 6 6 6 1 1
6 6 6 1 6 6 6 6 6 6 6 6 6 6 6
6 6 6 1 6 6 6 1 1 1 1 1 6 6 6
6 6 6 1 6 6 6 1 6 6 1 6 6 6 6
6 6 6 1 6 6 6 1 6 6 1 6 6 6 6
1 6 6 1 6 6 6 1 1 1 1 6 6 6 6
1 1 1 1 6 6 6 1 6 6 6 6 6 6 6
1 6 6 1 6 6 6 1 6 6 6 1 1 1 1
1 6 6 6 6 6 6 6 6 6 6 1 6 6 1
1 1 6 6 6 6 6 6 6 6 6 1 6 6 1
6 6 6 6 6 6 6 6 6 6 6 6 6 1 1
6 6 6 6 6 6 6 6 6 6 6 6 6 6 6
6 1 1 6 6 6 1 1 1 1 1 1 6 6 6
6 1 6 6 6 6 1 6 6 6 6 1 6 6 6
6 1 6 6 6 6 1 6 6 6 6 1 6 6 6
6 1 6 6 6 1 1 1 1 1 1 1 6 6 6"""

result = find_patterns(test_input)
output = '\n'.join(' '.join(str(x) for x in row) for row in result)
print(output)