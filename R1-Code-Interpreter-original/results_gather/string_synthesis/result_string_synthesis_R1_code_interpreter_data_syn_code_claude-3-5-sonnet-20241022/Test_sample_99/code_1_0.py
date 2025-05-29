def simulate_synthesis(blocks):
    def try_method(blocks):
        # Method 2: [A] + [B] -> {C}
        if blocks['[A]'] >= 1 and blocks['[B]'] >= 1:
            blocks['[A]'] -= 1
            blocks['[B]'] -= 1
            blocks['{C}'] += 1
            return True
            
        # Method 5: {A} + {C} -> (A) + (B)
        if blocks['{A}'] >= 1 and blocks['{C}'] >= 1:
            blocks['{A}'] -= 1
            blocks['{C}'] -= 1
            blocks['(A)'] += 1
            blocks['(B)'] += 1
            return True
            
        return False

    changed = True
    while changed:
        changed = try_method(blocks)
    
    result = ""
    for block, count in blocks.items():
        if count > 0:
            result += f"{count} {block}, "
    print(result.rstrip(", "))

initial = {
    '[A]': 5, '[B]': 7, '[C]': 1,
    '{A}': 0, '{B}': 0, '{C}': 0,
    '(A)': 0, '(B)': 0, '(C)': 0
}

simulate_synthesis(initial)