def process_grid(input_grid):
    n = len(input_grid)
    output = input_grid.copy()
    
    # Find all non-zero numbers and their positions
    numbers = {}
    for i, num in enumerate(input_grid):
        if num != 0:
            if num not in numbers:
                numbers[num] = []
            numbers[num].append(i)
    
    # Process each unique number
    for num, positions in numbers.items():
        if len(positions) == 1:  # Single number case
            # Fill all positions to the right
            for i in range(positions[0], n):
                output[i] = num
        elif len(positions) == 2:  # Two numbers case
            # Fill positions between the two numbers
            for i in range(positions[0], positions[1] + 1):
                output[i] = num
    
    return output

# Test input
test_input = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9, 0, 0, 0, 9]

# Process and print result
result = process_grid(test_input)
print(' '.join(map(str, result)))