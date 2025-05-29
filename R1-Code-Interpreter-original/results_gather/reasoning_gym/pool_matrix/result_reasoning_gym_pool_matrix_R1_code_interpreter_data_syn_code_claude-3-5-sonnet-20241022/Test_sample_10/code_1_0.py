# Input matrix as a list of lists
matrix = [
    [6, 3],
    [7, 4],
    [6, 9]
]

# Find maximum value in the entire matrix
max_value = max(max(row) for row in matrix)

# Format to 2 decimal places
result = "{:.2f}".format(max_value)
print(result)