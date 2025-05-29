def analyze_pattern():
    examples = [
        ("0 0 0 0 0 2 0 1 2 1 2 0 0 2 2", "2 2 2 2 2 1 1 0 0 0 0 0 0 0 0"),
        ("2 0 0 1 0 1 0 1 0 0 0 0 2 1 0", "2 2 1 1 1 1 0 0 0 0 0 0 0 0 0"),
        ("0 1 1 2 0 2 0 0 0 0 0 1 2 0 1", "2 2 2 1 1 1 1 0 0 0 0 0 0 0 0")
    ]
    
    for i, (inp, out) in enumerate(examples):
        inp_list = [int(x) for x in inp.split()]
        out_list = [int(x) for x in out.split()]
        
        print(f"\nExample {i+1}:")
        print(f"Input:  {inp}")
        print(f"Output: {out}")
        print(f"Pattern: {out_list[:out_list.count(2)]} | {out_list[out_list.count(2):out_list.count(2)+out_list.count(1)]} | {out_list[-out_list.count(0):]}")

analyze_pattern()