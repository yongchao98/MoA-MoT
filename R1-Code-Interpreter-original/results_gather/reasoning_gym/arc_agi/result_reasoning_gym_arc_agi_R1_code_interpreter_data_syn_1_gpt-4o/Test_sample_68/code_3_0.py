def reflect_and_rotate(grid):
    # Reflect along the main diagonal
    reflected = [list(row) for row in zip(*grid)]
    # Rotate 90 degrees clockwise
    return [list(reversed(col)) for col in zip(*reflected)]

# Test Input
test_input = [
    [1, 3, 5],
    [3, 5, 3],
    [8, 8, 3]
]
test_output = reflect_and_rotate(test_input)

print("Test Output:", test_output)