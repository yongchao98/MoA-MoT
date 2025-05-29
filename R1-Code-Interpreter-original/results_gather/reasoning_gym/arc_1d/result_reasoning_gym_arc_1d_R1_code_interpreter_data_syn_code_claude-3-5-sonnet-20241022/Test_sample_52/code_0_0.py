def analyze_patterns():
    # Example inputs and outputs
    inputs = [
        [2,0,2,0,0,0,2,2,0,2,1],
        [2,0,1,1,0,1,1,0,1,2,1],
        [0,0,0,1,2,0,2,1,0,2,0],
    ]
    outputs = [
        [2,2,2,2,2,1,0,0,0,0,0],
        [2,2,1,1,1,1,1,1,0,0,0],
        [2,2,2,1,1,0,0,0,0,0,0],
    ]
    
    for i in range(len(inputs)):
        input_counts = {0: inputs[i].count(0), 1: inputs[i].count(1), 2: inputs[i].count(2)}
        output_counts = {0: outputs[i].count(0), 1: outputs[i].count(1), 2: outputs[i].count(2)}
        print(f"Example {i+1}:")
        print(f"Input counts: {input_counts}")
        print(f"Output counts: {output_counts}")
        print("---")

    # Test input
    test = [0,0,0,2,2,2,0,2,0,0,0]
    test_counts = {0: test.count(0), 1: test.count(1), 2: test.count(2)}
    print("Test input counts:", test_counts)

analyze_patterns()