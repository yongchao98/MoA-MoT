def analyze_shifts():
    examples = [
        ([0,0,0,0,0,4,4,4,4,4,4,4,4,4,4,4,4,4], [0,0,4,4,4,4,4,4,4,4,4,4,4,4,4,0,0,0]),
        ([1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1], [1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1]),
        ([3,3,3,3,3,3,3,0,0,0,0,0,0,0,0,0,3,3], [3,3,3,3,0,0,0,0,0,0,0,0,0,3,3,3,3,3])
    ]
    
    test_input = [0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0]
    
    for i, (inp, out) in enumerate(examples, 1):
        non_zero = next(x for x in inp if x != 0)
        inp_positions = [i for i, x in enumerate(inp) if x != 0]
        out_positions = [i for i, x in enumerate(out) if x != 0]
        
        print(f"\nExample {i}:")
        print(f"Input non-zero positions: {inp_positions}")
        print(f"Output non-zero positions: {out_positions}")
        print(f"Input center: {sum(inp_positions)/len(inp_positions):.1f}")
        print(f"Output center: {sum(out_positions)/len(out_positions):.1f}")

    # Analyze test input
    print("\nTest input:")
    test_positions = [i for i, x in enumerate(test_input) if x != 0]
    print(f"Non-zero position: {test_positions}")
    print(f"Center position in grid: {len(test_input)/2 - 0.5}")

analyze_shifts()