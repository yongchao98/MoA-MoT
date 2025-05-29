def synthesize_from_state():
    blocks = {
        '[A]': 0, '[B]': 3, '[C]': 0,
        '{A}': 3, '{B}': 0, '{C}': 0,
        '(A)': 0, '(B)': 0, '(C)': 0
    }
    
    while True:
        any_change = False
        methods = [apply_method_1, apply_method_2, apply_method_3, 
                  apply_method_4, apply_method_5, apply_method_6]
        
        for method in methods:
            if method(blocks):
                any_change = True
                break
        
        if not any_change:
            break
    
    result = ""
    for block_type in ['[A]', '[B]', '[C]', '{A}', '{B}', '{C}', '(A)', '(B)', '(C)']:
        if blocks[block_type] > 0:
            result += f"{blocks[block_type]}{block_type} "
    
    print(result.strip())

synthesize_from_state()