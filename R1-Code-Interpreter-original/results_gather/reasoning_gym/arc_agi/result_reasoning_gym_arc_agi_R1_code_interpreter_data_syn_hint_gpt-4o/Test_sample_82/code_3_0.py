def find_clusters(grid):
    # This function will identify clusters of numbers in the grid
    from collections import defaultdict

    clusters = defaultdict(int)
    for row in grid:
        for num in row:
            if num != 1:  # Assuming 1 is the background number
                clusters[num] += 1
    return clusters

def select_numbers(clusters):
    # This function will select numbers based on the observed pattern
    # For simplicity, let's assume we select the most frequent numbers
    sorted_clusters = sorted(clusters.items(), key=lambda x: -x[1])
    selected_numbers = [num for num, _ in sorted_clusters[:3]]  # Select top 3 frequent numbers
    return sorted(selected_numbers)

# Test input grid
test_input = [
    [1, 8, 1, 1, 1, 1, 1, 1, 1, 1],
    [8, 8, 8, 1, 1, 1, 1, 6, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 6, 6, 1],
    [1, 1, 1, 1, 1, 1, 1, 6, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 3, 3, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 0],
    [1, 1, 1, 1, 1, 9, 9, 1, 1, 0],
    [1, 1, 1, 1, 1, 1, 9, 1, 1, 0],
    [1, 1, 1, 1, 1, 1, 9, 1, 1, 1]
]

clusters = find_clusters(test_input)
selected_numbers = select_numbers(clusters)

# Create the output grid
output_grid = [selected_numbers for _ in range(len(selected_numbers))]

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))