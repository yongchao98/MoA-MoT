def find_most_frequent_number(grid, start_row, start_col):
    # Extract the 3x3 block
    block = [grid[i][j] for i in range(start_row, start_row + 3) for j in range(start_col, start_col + 3)]
    # Count the frequency of each number
    frequency = {}
    for number in block:
        if number != 1:  # Ignore the number 1
            if number in frequency:
                frequency[number] += 1
            else:
                frequency[number] = 1
    # Find the most frequent number, choose the smallest in case of a tie
    most_frequent = min(frequency, key=lambda x: (-frequency[x], x))
    return most_frequent

# Define the input grid
input_grid = [
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 9, 9, 1, 9, 9, 1, 7, 7, 1, 9, 9, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 9, 9, 1, 7, 7, 1, 9, 9, 1, 2, 2, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 4, 4, 1, 9, 9, 1, 9, 9, 1, 4, 4, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 9, 9, 1, 2, 2, 1, 9, 9, 1, 9, 9, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 9, 9, 1, 9, 9, 1, 9, 9, 1, 7, 7, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 7, 7, 1, 9, 9, 1, 2, 2, 1, 9, 9, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 9, 9, 1, 7, 7, 1, 9, 9, 1, 9, 9, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
]

# Calculate the output
output = [
    find_most_frequent_number(input_grid, 1, 1),
    find_most_frequent_number(input_grid, 5, 5),
    find_most_frequent_number(input_grid, 9, 9)
]

# Print the output
print("<<<", *output, ">>>")