import numpy as np
from collections import Counter

def find_pattern(input_grid):
    # Convert input grid to numpy array for easier manipulation
    grid = np.array(input_grid)
    
    # Get the shape of the grid
    rows, cols = grid.shape
    
    # Extract the diagonal elements
    diagonal = np.diagonal(grid)
    
    # Count the frequency of each element in the diagonal
    counter = Counter(diagonal)
    
    # Sort elements by frequency and value
    most_common = sorted(counter.items(), key=lambda x: (-x[1], x[0]))
    
    # Create the output grid based on the most common elements
    output_grid = [x[0] for x in most_common]
    
    return output_grid

# Test input grid
input_grid = [
    [3, 3, 8],
    [3, 3, 3],
    [8, 3, 8]
]

# Find the pattern and print the output grid
output_grid = find_pattern(input_grid)
print(output_grid)