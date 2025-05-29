def process_grid(input_grid):
    rows = len(input_grid)
    output = [[7]*3 for _ in range(rows)]
    
    # Copy third column
    for i in range(rows):
        output[i][2] = input_grid[i][2]
    
    # Process first two columns
    for col in range(2):
        # Find all significant numbers and their positions
        numbers = []
        for i in range(rows):
            if input_grid[i][col] not in [7, 9]:
                numbers.append((i, input_grid[i][col]))
        
        # Process each number
        for i in range(len(numbers)):
            pos, num = numbers[i]
            
            # Determine the range this number should fill
            start = pos
            # If there's a next number, fill until that position
            if i < len(numbers) - 1:
                end = numbers[i + 1][0]
            else:
                # If it's the last number, fill until the end
                end = rows
            
            # Fill all positions from start to end
            for j in range(start, end):
                output[j][col] = num
                
            # Special case: if this is not the last number and the next number
            # is different, don't fill all the way to the next number's position
            if i < len(numbers) - 1 and numbers[i + 1][1] != num:
                # Find the last non-9 before the next number
                for j in range(numbers[i + 1][0] - 1, pos, -1):
                    if input_grid[j][col] != 9:
                        end = j + 1
                        break
                # Refill with the correct end position
                for j in range(start, end):
                    output[j][col] = num
    
    # Print the result
    for row in output:
        print(' '.join(map(str, row)))

# Test input
test_input = [
    [7, 0, 0],
    [7, 7, 9],
    [7, 7, 7],
    [7, 5, 9],
    [7, 9, 7],
    [7, 7, 7],
    [9, 6, 6],
    [7, 7, 9],
    [7, 7, 7],
    [7, 4, 9],
    [7, 4, 7]
]

process_grid(test_input)