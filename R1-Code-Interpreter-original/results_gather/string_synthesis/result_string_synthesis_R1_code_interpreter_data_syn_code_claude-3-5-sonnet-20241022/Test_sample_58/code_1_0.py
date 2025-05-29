def simulate_synthesis(blocks):
    changes = True
    while changes:
        changes = False
        # Store initial state for comparison
        initial_state = str(blocks)
        
        # Method 1: [A] + [B] + [C] -> {A}
        if blocks['[A]'] >= 1 and blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
            blocks['[A]'] -= 1
            blocks['[B]'] -= 1
            blocks['[C]'] -= 1
            blocks['{A}'] += 1
            changes = True
            continue
            
        # Method 2: [A] + [B] -> {C}
        if blocks['[A]'] >= 1 and blocks['[B]'] >= 1:
            blocks['[A]'] -= 1
            blocks['[B]'] -= 1
            blocks['{C}'] += 1
            changes = True
            continue
            
        # Method 3: [B] + [C] -> {B}
        if blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
            blocks['[B]'] -= 1
            blocks['[C]'] -= 1
            blocks['{B}'] += 1
            changes = True
            continue
            
        # Method 4: [C] + [C] -> {C}
        if blocks['[C]'] >= 2:
            blocks['[C]'] -= 2
            blocks['{C}'] += 1
            changes = True
            continue
            
        # Method 5: {A} + {C} -> (A) + (B)
        if blocks['{A}'] >= 1 and blocks['{C}'] >= 1:
            blocks['{A}'] -= 1
            blocks['{C}'] -= 1
            blocks['(A)'] += 1
            blocks['(B)'] += 1
            changes = True
            continue
            
        # Method 6: {B} + {B} -> (C)
        if blocks['{B}'] >= 2:
            blocks['{B}'] -= 2
            blocks['(C)'] += 1
            changes = True
            continue

    return blocks

initial_blocks = {
    '[A]': 1, '[B]': 7, '[C]': 3,
    '{A}': 0, '{B}': 0, '{C}': 0,
    '(A)': 0, '(B)': 0, '(C)': 0
}

result = simulate_synthesis(initial_blocks)
output = ', '.join(f"{v} {k}" for k, v in result.items() if v > 0)
print(output)