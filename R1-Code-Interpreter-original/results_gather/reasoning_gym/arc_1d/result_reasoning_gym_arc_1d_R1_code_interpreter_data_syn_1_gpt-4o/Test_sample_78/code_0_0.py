from collections import Counter

def transform_grid(input_grid):
    # Count the frequency of each number in the input grid
    counter = Counter(input_grid)
    # Remove zero from the counter as we are interested in non-zero numbers
    if 0 in counter:
        del counter[0]
    # Find the most common non-zero number
    most_common_number, _ = counter.most_common(1)[0]
    # Transform the grid
    output_grid = [most_common_number if x != 0 else 0 for x in input_grid]
    return output_grid

# Test input
input_grid = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 9, 2, 2, 4, 1, 1, 1, 1, 4, 1, 1, 9, 1, 1, 1, 1, 0, 0]
output_grid = transform_grid(input_grid)
print(output_grid)