def transform_grid(input_grid):
    # Convert string to list of integers
    grid = [int(x) for x in input_grid.split()]
    
    # Last two numbers go to front
    last_two = grid[-2:]
    
    # Third and fourth from end go to end
    third_fourth_from_end = grid[-4:-2]
    
    # Middle part is everything except last 4 numbers
    middle = grid[:-4]
    
    # Construct output
    output = last_two + middle + third_fourth_from_end
    
    return ' '.join(map(str, output))

# Test input
test_input = "0 0 0 0 0 0 0 0 6 4 4 7 8 8 6 9 8 5 8 2 1 8 0 0 0 0 0 0 0 0"
print(transform_grid(test_input))