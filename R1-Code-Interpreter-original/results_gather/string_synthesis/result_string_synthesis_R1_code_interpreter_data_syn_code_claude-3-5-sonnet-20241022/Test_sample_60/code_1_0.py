def simulate_synthesis():
    blocks = {
        '[A]': 1, '[B]': 1, '[C]': 3,
        '{A}': 0, '{B}': 0, '{C}': 0,
        '(A)': 0, '(B)': 0, '(C)': 0
    }
    
    def print_state():
        print("Current state:", {k: v for k, v in blocks.items() if v > 0})
    
    # Step 1: Method 4 ([C] + [C] -> {C})
    blocks['[C]'] -= 2
    blocks['{C}'] += 1
    print_state()
    
    # Step 2: Method 2 ([A] + [B] -> {C})
    blocks['[A]'] -= 1
    blocks['[B]'] -= 1
    blocks['{C}'] += 1
    print_state()
    
    # Step 3: Method 4 ([C] + remaining [C] -> {C})
    blocks['[C]'] -= 1
    blocks['{C}'] += 1
    print_state()
    
    # Step 4: Method 5 ({A} + {C} -> (A) + (B))
    blocks['{C}'] -= 2
    blocks['(B)'] += 1
    print_state()
    
    print("\nFinal result:", "1 (B)")

simulate_synthesis()