def analyze_patterns():
    examples = [
        ("1 1 1 1 0 0 0 0 2 0 0 2 2 1 0 2 2 1 2 2 0 0 0",
         "2 2 2 2 2 2 2 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1"),
        ("0 0 0 0 1 2 0 1 2 2 2 2 0 0 1 0 0 0 0 1 2 1 2",
         "2 2 2 2 2 2 2 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1"),
        ("0 1 1 0 0 0 1 0 0 2 1 0 0 1 0 0 2 0 0 0 2 0 1",
         "2 2 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1")
    ]
    
    for i, (inp, out) in enumerate(examples):
        inp = [int(x) for x in inp.split()]
        out = [int(x) for x in out.split()]
        print(f"\nExample {i+1}:")
        print(f"Input counts: 0s={inp.count(0)}, 1s={inp.count(1)}, 2s={inp.count(2)}")
        print(f"Output counts: 0s={out.count(0)}, 1s={out.count(1)}, 2s={out.count(2)}")

    # Analyze test input
    test = "2 2 1 0 0 0 0 2 0 0 0 0 0 0 1 0 1 0 2 0 0 0 0"
    test = [int(x) for x in test.split()]
    print("\nTest input counts:")
    print(f"0s={test.count(0)}, 1s={test.count(1)}, 2s={test.count(2)}")

analyze_patterns()