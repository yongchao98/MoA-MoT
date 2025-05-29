def analyze_and_solve(input_grid):
    # Count occurrences of each number
    counts = {0: 0, 1: 0, 2: 0}
    for num in input_grid:
        counts[num] += 1
    
    # Generate output by arranging in descending order
    output = []
    # First add all 2s
    output.extend([2] * counts[2])
    # Then add all 1s
    output.extend([1] * counts[1])
    # Finally add all 0s
    output.extend([0] * counts[0])
    
    # Print counts for verification
    print("Number counts:", counts)
    print("Generated output:", output)
    
    return output

# Test input
test_input = [2, 1, 0, 1, 1, 0, 0, 0, 0, 0, 2, 0, 0, 2, 1, 1]
result = analyze_and_solve(test_input)