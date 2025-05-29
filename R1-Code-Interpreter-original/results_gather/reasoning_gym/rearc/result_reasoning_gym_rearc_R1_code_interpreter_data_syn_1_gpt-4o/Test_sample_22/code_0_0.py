import numpy as np
from scipy.ndimage import label, find_objects

def find_subgrid(input_grid):
    # Convert input to numpy array for easier manipulation
    grid = np.array(input_grid)
    
    # Find the most common value
    unique, counts = np.unique(grid, return_counts=True)
    most_common_value = unique[np.argmax(counts)]
    
    # Create a mask for values different from the most common value
    mask = grid != most_common_value
    
    # Label connected components in the mask
    labeled_array, num_features = label(mask)
    
    # Find the largest connected component
    max_area = 0
    max_slice = None
    for i in range(1, num_features + 1):
        component_slice = find_objects(labeled_array == i)[0]
        area = (component_slice[0].stop - component_slice[0].start) * (component_slice[1].stop - component_slice[1].start)
        if area > max_area:
            max_area = area
            max_slice = component_slice
    
    # Extract the 4x4 subgrid centered around the largest component
    if max_slice:
        r_start, r_stop = max_slice[0].start, max_slice[0].stop
        c_start, c_stop = max_slice[1].start, max_slice[1].stop
        
        # Calculate the center of the component
        r_center = (r_start + r_stop) // 2
        c_center = (c_start + c_stop) // 2
        
        # Extract a 4x4 subgrid centered around the component
        r_start = max(0, r_center - 2)
        c_start = max(0, c_center - 2)
        r_stop = r_start + 4
        c_stop = c_start + 4
        
        subgrid = grid[r_start:r_stop, c_start:c_stop]
        return subgrid.tolist()
    else:
        return []

# Test input grid
input_grid = [
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 9, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 9, 7, 7, 8, 8, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 0, 7, 7, 2, 7, 7, 7, 7],
    [7, 7, 7, 7, 0, 7, 7, 2, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 4, 4, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7]
]

# Find and print the output grid
output_grid = find_subgrid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))