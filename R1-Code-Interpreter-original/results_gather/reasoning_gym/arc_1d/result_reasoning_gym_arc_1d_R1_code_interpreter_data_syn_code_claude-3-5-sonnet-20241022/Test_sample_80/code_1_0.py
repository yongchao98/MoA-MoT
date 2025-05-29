def analyze_patterns(inputs, outputs):
    # Convert string inputs to lists of numbers
    input1 = [0,0,0,0,4,0,4,0,6,3,0,0,7,2,4,2,0,0,0,4,9,1,0,7,0,0,7,2]
    output1 = [0,0,0,0,5,0,5,0,5,5,0,0,5,5,5,5,0,0,0,5,5,5,0,5,0,0,5,5]
    input2 = [5,3,4,0,2,2,0,0,8,9,0,1,0,0,0,0,7,9,5,0,3,0,0,2,7,0,0,6]
    output2 = [5,5,5,0,5,5,0,0,5,5,0,5,0,0,0,0,5,5,5,0,5,0,0,5,5,0,0,5]
    input3 = [3,9,8,9,0,1,0,5,0,1,0,0,0,0,0,7,0,3,4,0,9,5,0,3,6,8,5,0]
    output3 = [5,5,5,5,0,5,0,5,0,5,0,0,0,0,0,5,0,5,5,0,5,5,0,5,5,5,5,0]
    
    # Analyze patterns
    rules = []
    for i in range(28):
        if input1[i] != 0 and output1[i] == 5:
            rules.append("non-zero to 5")
        elif input1[i] == 0 and output1[i] == 0:
            rules.append("zero stays zero")
        else:
            rules.append("special case")
    
    # Test input
    test_input = [0,2,9,6,0,4,7,0,0,8,8,0,0,0,4,2,0,2,0,7,9,0,6,6,4,0,9,0]
    
    # Apply the discovered rule
    result = []
    for num in test_input:
        if num == 0:
            result.append(0)
        else:
            result.append(5)
    
    print(' '.join(map(str, result)))