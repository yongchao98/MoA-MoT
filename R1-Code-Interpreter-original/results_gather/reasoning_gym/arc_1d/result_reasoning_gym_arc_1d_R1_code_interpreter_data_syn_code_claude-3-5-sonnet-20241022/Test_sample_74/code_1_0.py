def verify_rule(input_str, output_str):
    input_nums = [int(x) for x in input_str.split()]
    output_nums = [int(x) for x in output_str.split()]
    
    print("Rule verification:")
    all_follow_rule = True
    for i, (inp, out) in enumerate(zip(input_nums, output_nums)):
        follows_rule = (inp == 0 and out == 0) or (inp != 0 and out == 5)
        print(f"Position {i}: {inp} -> {out} : {'Follows rule' if follows_rule else 'Breaks rule'}")
        all_follow_rule = all_follow_rule and follows_rule
    
    print(f"\nOverall rule compliance: {all_follow_rule}")

# Test all examples
examples = [
    ("0 0 5 0 8 8 0 4 0 8", "0 0 5 0 5 5 0 5 0 5"),
    ("0 6 2 6 0 7 2 0 3 0", "0 5 5 5 0 5 5 0 5 0"),
    ("9 5 7 3 3 6 0 0 0 3", "5 5 5 5 5 5 0 0 0 5")
]

for i, (inp, out) in enumerate(examples):
    print(f"\nExample {i+1}:")
    verify_rule(inp, out)