def analyze_grids():
    examples = [
        ("0 0 0 0 0 2 0 1 2 1 2 0 0 2 2", "2 2 2 2 2 1 1 0 0 0 0 0 0 0 0"),
        ("2 0 0 1 0 1 0 1 0 0 0 0 2 1 0", "2 2 1 1 1 1 0 0 0 0 0 0 0 0 0"),
        ("0 1 1 2 0 2 0 0 0 0 0 1 2 0 1", "2 2 2 1 1 1 1 0 0 0 0 0 0 0 0")
    ]
    
    for i, (inp, out) in enumerate(examples):
        inp_list = [int(x) for x in inp.split()]
        out_list = [int(x) for x in out.split()]
        
        print(f"\nExample {i+1}:")
        print(f"Input counts:  0s={inp_list.count(0)}, 1s={inp_list.count(1)}, 2s={inp_list.count(2)}")
        print(f"Output counts: 0s={out_list.count(0)}, 1s={out_list.count(1)}, 2s={out_list.count(2)}")

    # Analyze test input
    test = "0 1 0 0 2 1 0 0 0 1 0 1 1 2 0"
    test_list = [int(x) for x in test.split()]
    print(f"\nTest input counts: 0s={test_list.count(0)}, 1s={test_list.count(1)}, 2s={test_list.count(2)}")

analyze_grids()