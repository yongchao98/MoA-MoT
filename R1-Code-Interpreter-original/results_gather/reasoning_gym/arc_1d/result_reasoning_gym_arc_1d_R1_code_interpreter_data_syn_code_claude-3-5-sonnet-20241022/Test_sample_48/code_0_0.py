def analyze_patterns(examples):
    for i, (inp, out) in enumerate(examples, 1):
        print(f"\nExample {i}:")
        print(f"Input counts:  0s={inp.count('0')}, 1s={inp.count('1')}, 2s={inp.count('2')}")
        print(f"Output counts: 0s={out.count('0')}, 1s={out.count('1')}, 2s={out.count('2')}")
        
example1_in = "0 0 0 2 2 0 0 1 1 1 0 2 0 0 0 1 1 0 2 2 0 2"
example1_out = "2 2 2 2 2 2 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1"
example2_in = "0 2 0 0 1 2 0 1 0 1 1 1 0 1 0 2 0 0 0 0 1 2"
example2_out = "2 2 2 2 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1"
example3_in = "0 0 0 1 0 0 0 2 1 0 0 0 2 1 2 1 0 0 1 0 2 0"
example3_out = "2 2 2 2 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1"
test_in = "0 0 0 1 1 0 0 0 1 0 0 0 1 0 1 0 2 1 2 2 0 0"

examples = [
    (example1_in.split(), example1_out.split()),
    (example2_in.split(), example2_out.split()),
    (example3_in.split(), example3_out.split())
]

analyze_patterns(examples)
print("\nTest input counts:")
test_input = test_in.split()
print(f"0s={test_input.count('0')}, 1s={test_input.count('1')}, 2s={test_input.count('2')}")