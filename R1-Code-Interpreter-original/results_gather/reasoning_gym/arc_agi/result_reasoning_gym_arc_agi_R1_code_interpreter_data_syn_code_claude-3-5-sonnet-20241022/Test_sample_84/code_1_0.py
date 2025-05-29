def analyze_pattern():
    examples = [
        {"regions": 2, "output_length": 3},
        {"regions": 4, "output_length": 5},
        {"regions": 1, "output_length": 2}
    ]
    
    for ex in examples:
        print(f"Regions: {ex['regions']}, Output length: {ex['output_length']}")
        print(f"Difference: {ex['output_length'] - ex['regions']}")
        
    print("\nTest input regions:", 4)
    print("Expected output length:", 5)

analyze_pattern()