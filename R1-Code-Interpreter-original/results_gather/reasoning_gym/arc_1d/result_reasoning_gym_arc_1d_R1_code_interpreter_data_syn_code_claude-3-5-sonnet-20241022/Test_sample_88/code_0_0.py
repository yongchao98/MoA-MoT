def analyze_pattern(inputs, outputs):
    # Convert string inputs to lists
    examples = [
        ([4,4,4,0,4,4,0,4,0,4,4,4,4,0,0,0,0], [4,0,4,4,0,4,4,4,0,4,4,4,4,0,0,0,0]),
        ([4,0,4,4,4,0,4,4,4,4,0,4,4,0,0,0,0], [4,0,4,4,0,4,4,4,0,4,4,4,4,0,0,0,0]),
        ([6,6,6,0,6,6,6,6,0,6,0,6,6,0,0,0,0], [6,0,6,6,0,6,6,6,0,6,6,6,6,0,0,0,0])
    ]
    
    # Analyze the pattern
    test_input = [5,0,5,5,5,0,5,5,5,5,0,5,5,0,0,0,0]
    
    # The pattern appears to be:
    # X 0 X X 0 X X X 0 X X X X 0 0 0 0
    # where X is the non-zero number
    
    # Generate output based on pattern
    result = []
    non_zero = test_input[0]  # Get the non-zero number (5 in this case)
    pattern = [1,0,1,1,0,1,1,1,0,1,1,1,1,0,0,0,0]
    
    for p in pattern:
        if p == 1:
            result.append(non_zero)
        else:
            result.append(0)
    
    print(' '.join(map(str, result)))

analyze_pattern(None, None)