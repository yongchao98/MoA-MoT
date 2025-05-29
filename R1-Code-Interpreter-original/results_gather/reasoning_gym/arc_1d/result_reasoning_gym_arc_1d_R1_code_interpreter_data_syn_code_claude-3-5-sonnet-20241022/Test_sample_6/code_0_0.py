def print_grid(grid, width=5):
    for i in range(0, len(grid), width):
        print(grid[i:i+width])

def find_pattern(input_grid, output_grid):
    print("Input pattern:")
    print_grid(input_grid)
    print("\nOutput pattern:")
    print_grid(output_grid)

# Example 1
input1 = [2,2,2,2,2,0,0,0,2,0,0,0,0,0,0,2,0,0,0,0,0,0,0]
output1 = [2,2,2,2,2,0,2,2,2,2,2,0,0,2,2,2,2,2,0,0,0,0,0]

print("Example 1:")
find_pattern(input1, output1)