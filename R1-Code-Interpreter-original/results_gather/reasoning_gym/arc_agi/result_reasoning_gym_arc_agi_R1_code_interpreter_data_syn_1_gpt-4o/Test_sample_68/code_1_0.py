def rotate_90_counterclockwise(grid):
    return [list(col) for col in zip(*grid[::-1])]

# Example 1
input1 = [
    [1, 1, 1, 5, 8, 1],
    [3, 8, 5, 5, 1, 8],
    [1, 8, 3, 8, 1, 8],
    [8, 3, 5, 1, 8, 1],
    [1, 5, 8, 1, 5, 8],
    [5, 5, 8, 3, 1, 8]
]
output1 = rotate_90_counterclockwise(input1)

# Example 2
input2 = [
    [8, 3, 3, 3],
    [8, 3, 1, 3],
    [8, 5, 8, 5],
    [8, 8, 5, 1]
]
output2 = rotate_90_counterclockwise(input2)

# Example 3
input3 = [
    [8, 3, 3, 3, 5, 5, 5],
    [1, 5, 5, 1, 1, 5, 3],
    [3, 3, 8, 1, 8, 5, 5],
    [8, 8, 5, 5, 5, 8, 8],
    [8, 5, 8, 1, 5, 8, 5],
    [1, 1, 8, 8, 3, 3, 3],
    [3, 3, 3, 3, 3, 8, 3]
]
output3 = rotate_90_counterclockwise(input3)

print("Output 1:", output1)
print("Output 2:", output2)
print("Output 3:", output3)