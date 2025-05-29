def analyze_all_examples():
    # All examples
    examples = [
        ([0,0,0,0,0,4,4,4,4,4,4,4,4,4,4,4,4,4], [0,0,4,4,4,4,4,4,4,4,4,4,4,4,4,0,0,0]),
        ([1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1], [1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1]),
        ([3,3,3,3,3,3,3,0,0,0,0,0,0,0,0,0,3,3], [3,3,3,3,0,0,0,0,0,0,0,0,0,3,3,3,3,3])
    ]
    
    for i, (inp, out) in enumerate(examples, 1):
        print(f"\nExample {i}:")
        # Find the non-zero number and its count
        non_zero = next(x for x in inp if x != 0)
        count_nonzero = inp.count(non_zero)
        count_zeros = inp.count(0)
        
        # Analyze the distribution
        first_nonzero_in = inp.index(non_zero)
        last_nonzero_in = len(inp) - inp[::-1].index(non_zero) - 1
        first_nonzero_out = out.index(non_zero)
        last_nonzero_out = len(out) - out[::-1].index(non_zero) - 1
        
        print(f"Number {non_zero} count: {count_nonzero}")
        print(f"Zeros count: {count_zeros}")
        print(f"Input span: {first_nonzero_in} to {last_nonzero_in}")
        print(f"Output span: {first_nonzero_out} to {last_nonzero_out}")
        print(f"Input leading zeros: {first_nonzero_in}")
        print(f"Output leading zeros: {first_nonzero_out}")

analyze_all_examples()