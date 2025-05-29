def analyze_transformations():
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
        # Check if numbers move clockwise or counterclockwise
        print("First row comparison:")
        print(f"Input first row: {inp[0]}")
        print(f"Output first row: {out[0]}")
        print("Last column comparison:")
        print(f"Input last column: [{inp[0][2]}, {inp[1][2]}, {inp[2][2]}]")
        print(f"Output last column: [{out[0][2]}, {out[1][2]}, {out[2][2]}]")

analyze_transformations()