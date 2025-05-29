def synthesize():
    blocks = {
        '[A]': 5, '[B]': 7, '[C]': 5,
        '{A}': 0, '{B}': 0, '{C}': 0,
        '(A)': 0, '(B)': 0, '(C)': 0
    }
    
    # First, try to make as many {A} as possible using Rule 1
    while blocks['[A]'] >= 1 and blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
        blocks['[A]'] -= 1
        blocks['[B]'] -= 1
        blocks['[C]'] -= 1
        blocks['{A}'] += 1
    
    result = ""
    for block in ['[A]', '[B]', '[C]', '{A}', '{B}', '{C}', '(A)', '(B)', '(C)']:
        if blocks[block] > 0:
            result += f"{blocks[block]} {block}, "
    print(result[:-2])

synthesize()