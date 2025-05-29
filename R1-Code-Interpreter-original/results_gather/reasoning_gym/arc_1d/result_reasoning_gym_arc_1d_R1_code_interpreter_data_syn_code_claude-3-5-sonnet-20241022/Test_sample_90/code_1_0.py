def analyze_grid_patterns():
    examples = [
        ([5,5,5,5,0,5,5,0,5,5,5,0,5,0,0,0], [5,0,5,5,0,5,5,5,0,5,5,5,5,0,0,0]),
        ([4,4,0,4,0,4,4,4,0,4,4,4,4,0,0,0], [4,0,4,4,0,4,4,4,0,4,4,4,4,0,0,0]),
        ([2,2,2,2,0,2,0,2,2,0,2,2,2,0,0,0], [2,0,2,2,0,2,2,2,0,2,2,2,2,0,0,0])
    ]
    
    for i, (inp, out) in enumerate(examples):
        print(f"\nExample {i+1}:")
        print("Input grid:")
        for j in range(0, 16, 4):
            print(inp[j:j+4])
        print("Output grid:")
        for j in range(0, 16, 4):
            print(out[j:j+4])

    test = [9,9,9,0,9,0,9,9,0,9,9,9,9,0,0,0]
    print("\nTest input grid:")
    for j in range(0, 16, 4):
        print(test[j:j+4])