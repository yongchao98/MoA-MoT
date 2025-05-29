def transform_grid(input_grid):
    # Convert string to list of integers
    grid = [int(x) for x in input_grid.split()]
    
    # Get relevant parts
    last_two = grid[-2:]  # Last two numbers
    main_part = grid[:-4]  # Everything except last 4 numbers
    third_fourth_from_end = grid[-4:-2]  # Numbers that were 3rd and 4th from end
    
    # Construct output:
    # 1. Start with last two numbers
    # 2. Add main part
    # 3. Add the numbers that were 3rd and 4th from end
    output = last_two + main_part + third_fourth_from_end
    
    # Ensure exactly 30 numbers
    output = output[:30]
    
    return ' '.join(map(str, output))

# Test input
test_input = "0 0 0 0 0 0 0 0 6 4 4 7 8 8 6 9 8 5 8 2 1 8 0 0 0 0 0 0 0 0"

# Print result
print(transform_grid(test_input))