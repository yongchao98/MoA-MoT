def analyze_transformations():
    def find_isolated_numbers(arr):
        sequences = []
        current = arr[0]
        count = 1
        start = 0
        
        for i in range(1, len(arr)):
            if arr[i] != current:
                if count >= 2:
                    sequences.append((start, current, count))
                current = arr[i]
                count = 1
                start = i
            else:
                count += 1
        if count >= 2:
            sequences.append((start, current, count))
        return sequences
    
    # Test cases
    examples = [
        ([8,8,0,0,0,0,0,0,0,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8],
         [0,0,0,0,0,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,0,0]),
        ([6,6,6,6,6,6,6,6,0,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6],
         [6,6,6,6,0,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6])
    ]
    
    for i, (inp, out) in enumerate(examples):
        print(f"\nExample {i+1}:")
        print("Input sequences:", find_isolated_numbers(inp))
        print("Output sequences:", find_isolated_numbers(out))
    
    test = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6,6,0,0,0,0,0]
    print("\nTest input sequences:", find_isolated_numbers(test))

analyze_transformations()