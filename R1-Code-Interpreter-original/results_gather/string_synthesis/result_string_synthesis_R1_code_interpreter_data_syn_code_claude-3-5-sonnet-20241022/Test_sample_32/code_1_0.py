def simulate():
    blocks = {
        '[A]': 5, '[B]': 5, '[C]': 3,
        '{A}': 0, '{B}': 0, '{C}': 0,
        '(A)': 0, '(B)': 0, '(C)': 0
    }
    
    # Method 1: [A] + [B] + [C] -> {A}
    blocks['[A]'] -= 1
    blocks['[B]'] -= 1
    blocks['[C]'] -= 1
    blocks['{A}'] += 1
    
    # Method 2: [A] + [B] -> {C} (twice)
    blocks['[A]'] -= 2
    blocks['[B]'] -= 2
    blocks['{C}'] += 2
    
    # Method 4: [C] + [C] -> {C}
    blocks['[C]'] -= 2
    blocks['{C}'] += 1
    
    # Method 5: {A} + {C} -> (A) + (B) (twice)
    blocks['{A}'] -= 1
    blocks['{C}'] -= 3
    blocks['(A)'] += 2
    blocks['(B)'] += 2
    
    result = ""
    for block, count in sorted(blocks.items()):
        if count > 0:
            result += f"{count} {block}, "
    print(result[:-2])

simulate()