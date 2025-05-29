def analyze_grids():
    examples = [
        ([1,0,1,1,2,0,0,0,0,1,1,1,2,0,1,0,2,2,0], [1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,2,2,2,2]),
        ([0,2,0,0,2,0,0,2,1,0,0,1,0,0,0,0,1,0,2], [1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,2,2,2,2]),
        ([1,0,2,1,0,0,0,0,1,0,2,0,1,1,1,2,0,0,0], [1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,2,2,2])
    ]
    
    for i, (inp, out) in enumerate(examples):
        print(f"\nExample {i+1}:")
        print(f"Input counts: 1s={inp.count(1)}, 0s={inp.count(0)}, 2s={inp.count(2)}")
        print(f"Output counts: 1s={out.count(1)}, 0s={out.count(0)}, 2s={out.count(2)}")

analyze_grids()