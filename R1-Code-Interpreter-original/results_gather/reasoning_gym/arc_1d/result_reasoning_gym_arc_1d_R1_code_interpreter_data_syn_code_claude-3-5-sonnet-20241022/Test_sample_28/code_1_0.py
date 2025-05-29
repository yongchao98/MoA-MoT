def verify_test_case():
    test_input = [1,0,0,1,2,1,0,0,0,0,1,1,2,0,2,0,0,2,0]
    
    # Count each digit
    ones = test_input.count(1)
    zeros = test_input.count(0)
    twos = test_input.count(2)
    
    print(f"Test input analysis:")
    print(f"Number of 1s: {ones}")
    print(f"Number of 0s: {zeros}")
    print(f"Number of 2s: {twos}")
    print(f"Total length: {len(test_input)}")
    
    # Generate output
    output = [1]*ones + [0]*zeros + [2]*twos
    print("\nGenerated output:", " ".join(map(str, output)))

verify_test_case()