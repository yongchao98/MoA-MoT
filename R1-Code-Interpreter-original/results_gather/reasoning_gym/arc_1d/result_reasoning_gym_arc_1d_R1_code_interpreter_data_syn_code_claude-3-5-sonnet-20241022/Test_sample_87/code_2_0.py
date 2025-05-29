def analyze_sequences():
    examples = [
        ("0 0 0 0 0 2 0 1 2 1 2 0 0 2 2", "2 2 2 2 2 1 1 0 0 0 0 0 0 0 0"),
        ("2 0 0 1 0 1 0 1 0 0 0 0 2 1 0", "2 2 1 1 1 1 0 0 0 0 0 0 0 0 0"),
        ("0 1 1 2 0 2 0 0 0 0 0 1 2 0 1", "2 2 2 1 1 1 1 0 0 0 0 0 0 0 0")
    ]
    
    test = "0 1 0 0 2 1 0 0 0 1 0 1 1 2 0"
    test_list = [int(x) for x in test.split()]
    
    print("Test input analysis:")
    print(f"2's count: {test_list.count(2)} (positions in input: {[i for i, x in enumerate(test_list) if x == 2]})")
    print(f"1's count: {test_list.count(1)} (positions in input: {[i for i, x in enumerate(test_list) if x == 1]})")
    print(f"0's count: {test_list.count(0)} (positions in input: {[i for i, x in enumerate(test_list) if x == 0]})")

analyze_sequences()