def synthesize():
    blocks = {
        '[A]': 2, '[B]': 5, '[C]': 1,
        '{A}': 0, '{B}': 0, '{C}': 0,
        '(A)': 0, '(B)': 0, '(C)': 0
    }
    
    while True:
        initial_state = blocks.copy()
        
        # Try each method in sequence
        # Method 1: [A] + [B] + [C] -> {A}
        if blocks['[A]'] >= 1 and blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
            blocks['[A]'] -= 1
            blocks['[B]'] -= 1
            blocks['[C]'] -= 1
            blocks['{A}'] += 1
            continue

        # Method 2: [A] + [B] -> {C}
        if blocks['[A]'] >= 1 and blocks['[B]'] >= 1:
            blocks['[A]'] -= 1
            blocks['[B]'] -= 1
            blocks['{C}'] += 1
            continue

        # Method 3: [B] + [C] -> {B}
        if blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
            blocks['[B]'] -= 1
            blocks['[C]'] -= 1
            blocks['{B}'] += 1
            continue

        # Method 4: [C] + [C] -> {C}
        if blocks['[C]'] >= 2:
            blocks['[C]'] -= 2
            blocks['{C}'] += 1
            continue

        # Method 5: {A} + {C} -> (A) + (B)
        if blocks['{A}'] >= 1 and blocks['{C}'] >= 1:
            blocks['{A}'] -= 1
            blocks['{C}'] -= 1
            blocks['(A)'] += 1
            blocks['(B)'] += 1
            continue

        # Method 6: {B} + {B} -> (C)
        if blocks['{B}'] >= 2:
            blocks['{B}'] -= 2
            blocks['(C)'] += 1
            continue

        if blocks == initial_state:
            break
    
    result = ""
    for block_type, count in sorted(blocks.items()):
        if count > 0:
            result += f"{count} {block_type}, "
    print(result[:-2])

synthesize()