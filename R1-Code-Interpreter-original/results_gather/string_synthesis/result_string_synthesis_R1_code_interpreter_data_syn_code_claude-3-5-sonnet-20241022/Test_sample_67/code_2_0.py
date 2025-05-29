def synthesize():
    blocks = {
        '[A]': 3, '[B]': 5, '[C]': 7,
        '{A}': 0, '{B}': 0, '{C}': 0,
        '(A)': 0, '(B)': 0, '(C)': 0
    }
    
    print("Initial:", blocks)
    
    # Rule application sequence with state tracking
    while True:
        old_state = blocks.copy()
        
        # Rule 1: [A] + [B] + [C] -> {A}
        while blocks['[A]'] >= 1 and blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
            blocks['[A]'] -= 1
            blocks['[B]'] -= 1
            blocks['[C]'] -= 1
            blocks['{A}'] += 1
            print("After Rule 1:", blocks)
            
        # Rule 2: [A] + [B] -> {C}
        while blocks['[A]'] >= 1 and blocks['[B]'] >= 1:
            blocks['[A]'] -= 1
            blocks['[B]'] -= 1
            blocks['{C}'] += 1
            print("After Rule 2:", blocks)
            
        # Rule 3: [B] + [C] -> {B}
        while blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
            blocks['[B]'] -= 1
            blocks['[C]'] -= 1
            blocks['{B}'] += 1
            print("After Rule 3:", blocks)
            
        # Rule 4: [C] + [C] -> {C}
        while blocks['[C]'] >= 2:
            blocks['[C]'] -= 2
            blocks['{C}'] += 1
            print("After Rule 4:", blocks)
            
        # Rule 5: {A} + {C} -> (A) + (B)
        while blocks['{A}'] >= 1 and blocks['{C}'] >= 1:
            blocks['{A}'] -= 1
            blocks['{C}'] -= 1
            blocks['(A)'] += 1
            blocks['(B)'] += 1
            print("After Rule 5:", blocks)
            
        # Rule 6: {B} + {B} -> (C)
        while blocks['{B}'] >= 2:
            blocks['{B}'] -= 2
            blocks['(C)'] += 1
            print("After Rule 6:", blocks)
            
        if blocks == old_state:
            break
    
    print("\nFinal state:", blocks)

synthesize()