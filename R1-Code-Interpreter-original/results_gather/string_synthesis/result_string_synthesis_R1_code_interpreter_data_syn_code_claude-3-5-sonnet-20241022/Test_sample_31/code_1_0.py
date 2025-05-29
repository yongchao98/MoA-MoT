def simulate_synthesis():
    # Initial blocks
    blocks = {
        '[A]': 2, '[B]': 5, '[C]': 4,
        '{A}': 0, '{B}': 0, '{C}': 0,
        '(A)': 0, '(B)': 0, '(C)': 0
    }
    
    changed = True
    while changed:
        changed = False
        # Method 1: [A] + [B] + [C] -> {A}
        while blocks['[A]'] >= 1 and blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
            blocks['[A]'] -= 1
            blocks['[B]'] -= 1
            blocks['[C]'] -= 1
            blocks['{A}'] += 1
            changed = True
            
        # Method 4: [C] + [C] -> {C}
        while blocks['[C]'] >= 2:
            blocks['[C]'] -= 2
            blocks['{C}'] += 1
            changed = True
            
        # Method 5: {A} + {C} -> (A) + (B)
        while blocks['{A}'] >= 1 and blocks['{C}'] >= 1:
            blocks['{A}'] -= 1
            blocks['{C}'] -= 1
            blocks['(A)'] += 1
            blocks['(B)'] += 1
            changed = True
            
        # Method 6: {B} + {B} -> (C)
        while blocks['{B}'] >= 2:
            blocks['{B}'] -= 2
            blocks['(C)'] += 1
            changed = True
            
        # Method 2: [A] + [B] -> {C}
        while blocks['[A]'] >= 1 and blocks['[B]'] >= 1:
            blocks['[A]'] -= 1
            blocks['[B]'] -= 1
            blocks['{C}'] += 1
            changed = True
            
        # Method 3: [B] + [C] -> {B}
        while blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
            blocks['[B]'] -= 1
            blocks['[C]'] -= 1
            blocks['{B}'] += 1
            changed = True
    
    result = ""
    for block_type in ['[A]', '[B]', '[C]', '{A}', '{B}', '{C}', '(A)', '(B)', '(C)']:
        if blocks[block_type] > 0:
            result += f"{blocks[block_type]} {block_type}, "
    print(result[:-2])

simulate_synthesis()