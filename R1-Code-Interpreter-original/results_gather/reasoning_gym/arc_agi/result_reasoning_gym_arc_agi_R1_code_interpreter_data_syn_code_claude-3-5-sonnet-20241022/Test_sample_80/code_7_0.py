def analyze_all_patterns():
    examples = [
        # Example 1 - single 9 at (2,2)
        {"distribution": "bottom edge twos", "output": [(2,2)]},
        # Example 2 - single 9 at (2,2)
        {"distribution": "bottom edge twos", "output": [(2,2)]},
        # Example 3 - full diagonal
        {"distribution": "symmetric distribution", "output": [(0,0), (1,1), (2,2)]},
        # Example 4 - two 9s
        {"distribution": "bottom heavy, balanced sides", "output": [(1,1), (2,2)]},
        # Example 5 - two 9s
        {"distribution": "bottom heavy, balanced sides", "output": [(1,1), (2,2)]}
    ]
    
    patterns = {}
    for i, example in enumerate(examples, 1):
        key = f"bottom heavy + balanced sides"
        if key in patterns:
            patterns[key].append(example["output"])
        else:
            patterns[key] = [example["output"]]
    
    return patterns

# Analyze test input characteristics
test_characteristics = {
    "left_right_balance": "5:4 (nearly balanced)",
    "top_bottom_balance": "1:8 (bottom heavy)",
    "cluster_pattern": "two large clusters + one small",
    "similar_examples": "4 and 5"
}

print("Test input characteristics:", test_characteristics)