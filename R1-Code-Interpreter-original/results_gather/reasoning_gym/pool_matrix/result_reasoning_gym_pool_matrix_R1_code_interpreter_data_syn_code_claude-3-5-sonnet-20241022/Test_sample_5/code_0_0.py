# Input matrix
input_str = """5 6 3 3 5 8 5 7 8 8
3 2 8 1 8 3 8 4 1 3
2 3 3 2 9 3 7 5 6 1"""

# Convert input string to matrix
matrix = [[float(num) for num in row.split()] for row in input_str.split('\n')]

# For kernel size 1, max pooling just returns the same matrix
# Format each number to 2 decimal places
result = [[f"{num:.2f}" for num in row] for row in matrix]

# Print result in the required format
for row in result:
    print(" ".join(row))