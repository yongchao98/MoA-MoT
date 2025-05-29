def analyze_pattern():
    # Example pairs
    examples = [
        ("2 2 5 9 5 9 6 3 1 8 3 0 0 0 0 0 3 8 2 6 5 1",
         "9 5 9 6 3 1 8 3 0 0 0 0 0 3 8 2 6 5 1 2 2 5"),
        ("9 6 8 6 8 6 0 0 0 0 0 0 0 0 0 3 1 4 5 5 7 6",
         "6 8 6 0 0 0 0 0 0 0 0 0 3 1 4 5 5 7 6 9 6 8"),
        ("9 8 7 1 8 3 7 1 6 7 6 7 9 8 5 8 2 2 0 0 0 5",
         "1 8 3 7 1 6 7 6 7 9 8 5 8 2 2 0 0 0 5 9 8 7")
    ]
    
    for i, (inp, out) in enumerate(examples, 1):
        inp = inp.split()
        out = out.split()
        print(f"\nExample {i}:")
        print("Input:  ", inp)
        print("Output: ", out)
        
        # Find the shift
        n = len(inp)
        for shift in range(n):
            if inp[shift:] + inp[:shift] == out:
                print(f"Shift found: {shift} positions to the right")

print("Pattern Analysis:")
analyze_pattern()