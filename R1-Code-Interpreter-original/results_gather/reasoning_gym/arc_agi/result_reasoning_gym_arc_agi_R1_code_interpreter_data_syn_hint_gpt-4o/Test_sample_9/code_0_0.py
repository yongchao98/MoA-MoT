# Define the input grid
input_grid = [
    [1, 9, 9, 1, 9, 9, 1, 7, 7, 1, 9, 9, 1],
    [1, 9, 9, 1, 7, 7, 1, 9, 9, 1, 2, 2, 1],
    [1, 4, 4, 1, 9, 9, 1, 9, 9, 1, 4, 4, 1],
    [1, 9, 9, 1, 2, 2, 1, 9, 9, 1, 9, 9, 1],
    [1, 9, 9, 1, 9, 9, 1, 9, 9, 1, 7, 7, 1],
    [1, 7, 7, 1, 9, 9, 1, 2, 2, 1, 9, 9, 1],
    [1, 9, 9, 1, 7, 7, 1, 9, 9, 1, 9, 9, 1]
]

# Flatten the grid and find unique numbers excluding 1
unique_numbers = set()
for row in input_grid:
    for num in row:
        if num != 1:
            unique_numbers.add(num)

# Convert the set to a sorted list
output_grid = sorted(unique_numbers)
print(output_grid)