def analyze_vertical_lines():
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
        # Check vertical lines
        for col in range(3):
            print(f"Column {col}:")
            print(f"Input: [{inp[0][col]}, {inp[1][col]}, {inp[2][col]}]")
            print(f"Output: [{out[0][col]}, {out[1][col]}, {out[2][col]}]")

analyze_vertical_lines()