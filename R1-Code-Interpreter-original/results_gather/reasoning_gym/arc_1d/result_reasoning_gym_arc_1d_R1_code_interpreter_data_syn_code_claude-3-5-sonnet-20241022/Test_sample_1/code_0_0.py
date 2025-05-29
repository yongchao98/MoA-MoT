def transform_grid(input_grid):
    # Convert string to list of integers
    grid = [int(x) for x in input_grid.split()]
    
    # Get last two numbers
    last_two = grid[-2:]
    
    # Create output by combining last two numbers with rest (shifted right)
    output = last_two + grid[:-2]
    
    # Convert to string format
    return ' '.join(map(str, output))

# Test cases
test1 = "3 9 6 9 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 8 9 9"
test2 = "0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 2 5 7 5 9 8 1 7 5 6 9 4 2 1 6"
test3 = "2 2 3 2 2 6 1 1 0 4 2 6 2 6 4 9 5 6 1 6 9 5 9 7 7 3 2 1 6 9"
test_input = "0 0 0 0 0 0 0 0 6 4 4 7 8 8 6 9 8 5 8 2 1 8 0 0 0 0 0 0 0 0"

# Verify pattern with test cases
print("Test case 1 output:", transform_grid(test1))
print("Test case 2 output:", transform_grid(test2))
print("Test case 3 output:", transform_grid(test3))
print("\nFinal answer for test input:", transform_grid(test_input))