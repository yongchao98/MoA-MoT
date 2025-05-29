def simulate_final():
    blocks = {
        '[A]': 2, '[B]': 5, '[C]': 4,
        '{A}': 0, '{B}': 0, '{C}': 0,
        '(A)': 0, '(B)': 0, '(C)': 0
    }
    
    # Step 1: Apply Method 1 twice
    for _ in range(2):
        blocks['[A]'] -= 1
        blocks['[B]'] -= 1
        blocks['[C]'] -= 1
        blocks['{A}'] += 1
    
    # Step 2: Apply Method 4 once
    blocks['[C]'] -= 2
    blocks['{C}'] += 1
    
    # Step 3: Apply Method 5 once
    blocks['{A}'] -= 1
    blocks['{C}'] -= 1
    blocks['(A)'] += 1
    blocks['(B)'] += 1
    
    result = ""
    for block_type in ['[A]', '[B]', '[C]', '{A}', '{B}', '{C}', '(A)', '(B)', '(C)']:
        if blocks[block_type] > 0:
            result += f"{blocks[block_type]} {block_type}, "
    print(result[:-2])

simulate_final()