def analyze_patterns():
    # Example pairs of input and output grids
    examples = [
        (
            [[4,4,4], [4,8,8], [4,8,8]],
            [[4,4,4], [8,8,4], [8,8,4]]
        ),
        (
            [[5,5,5], [5,8,8], [8,8,5]],
            [[8,5,5], [8,8,5], [5,8,5]]
        ),
        (
            [[4,4,8], [4,4,8], [4,4,4]],
            [[4,4,4], [4,4,4], [4,8,8]]
        ),
        (
            [[8,8,7], [7,8,7], [7,8,7]],
            [[7,7,8], [8,8,8], [7,7,7]]
        )
    ]
    
    for idx, (inp, out) in enumerate(examples):
        print(f"\nExample {idx + 1}:")
        print("Input:")
        for row in inp:
            print(row)
        print("Output:")
        for row in out:
            print(row)
        
        # Analyze number movements
        numbers = set(sum(inp, []))
        for num in numbers:
            print(f"\nNumber {num} positions:")
            inp_pos = [(i,j) for i in range(3) for j in range(3) if inp[i][j] == num]
            out_pos = [(i,j) for i in range(3) for j in range(3) if out[i][j] == num]
            print(f"Input positions: {inp_pos}")
            print(f"Output positions: {out_pos}")

analyze_patterns()