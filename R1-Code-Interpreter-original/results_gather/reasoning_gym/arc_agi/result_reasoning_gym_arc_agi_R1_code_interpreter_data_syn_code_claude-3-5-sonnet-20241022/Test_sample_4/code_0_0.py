def analyze_changes(input_grid, output_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    
    print("Changes from input to output:")
    for i in range(rows):
        for j in range(cols):
            if input_grid[i][j] != output_grid[i][j]:
                print(f"Position ({i},{j}): {input_grid[i][j]} -> {output_grid[i][j]}")

# Example 1
input1 = [
    [3,0,1,1,1],
    [5,5,9,5,8],
    [1,5,5,5,3],
    [1,1,8,6,9],
    [9,5,2,4,6]
]

output1 = [
    [3,0,5,5,5],
    [1,1,9,1,8],
    [5,1,1,1,3],
    [5,5,8,6,9],
    [9,1,2,4,6]
]

print("Example 1 analysis:")
analyze_changes(input1, output1)

# Example 2
input2 = [
    [9,5,6],
    [5,0,5],
    [1,1,3]
]

output2 = [
    [9,1,6],
    [1,0,1],
    [5,5,3]
]

print("\nExample 2 analysis:")
analyze_changes(input2, output2)