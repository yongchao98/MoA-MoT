def analyze_patterns():
    # Example pairs
    examples = [
        ([5,5,5,5,0,5,5,0,5,5,5,0,5,0,0,0], [5,0,5,5,0,5,5,5,0,5,5,5,5,0,0,0]),
        ([4,4,0,4,0,4,4,4,0,4,4,4,4,0,0,0], [4,0,4,4,0,4,4,4,0,4,4,4,4,0,0,0]),
        ([2,2,2,2,0,2,0,2,2,0,2,2,2,0,0,0], [2,0,2,2,0,2,2,2,0,2,2,2,2,0,0,0])
    ]
    
    for i, (inp, out) in enumerate(examples):
        print(f"\nExample {i+1}:")
        print("Input:  ", inp)
        print("Output: ", out)
        print("Changes at positions:", [i for i in range(16) if inp[i] != out[i]])

    # Test input
    test = [9,9,9,0,9,0,9,9,0,9,9,9,9,0,0,0]
    print("\nTest input:", test)

analyze_patterns()