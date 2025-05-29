def synthesize():
    blocks = {
        '[A]': 1, '[B]': 5, '[C]': 7,
        '{A}': 0, '{B}': 0, '{C}': 0,
        '(A)': 0, '(B)': 0, '(C)': 0
    }
    
    print("Initial:", dict(blocks))
    
    # Rule 1: [A]+[B]+[C] -> {A}
    blocks['[A]'] -= 1
    blocks['[B]'] -= 1
    blocks['[C]'] -= 1
    blocks['{A}'] += 1
    print("After Rule 1:", dict(blocks))
    
    # Rule 3 (4 times): [B]+[C] -> {B}
    for _ in range(4):
        blocks['[B]'] -= 1
        blocks['[C]'] -= 1
        blocks['{B}'] += 1
    print("After Rule 3 x4:", dict(blocks))
    
    # Rule 4: [C]+[C] -> {C}
    blocks['[C]'] -= 2
    blocks['{C}'] += 1
    print("After Rule 4:", dict(blocks))
    
    # Rule 5: {A}+{C} -> (A)+(B)
    blocks['{A}'] -= 1
    blocks['{C}'] -= 1
    blocks['(A)'] += 1
    blocks['(B)'] += 1
    print("After Rule 5:", dict(blocks))
    
    # Rule 6 (twice): {B}+{B} -> (C)
    blocks['{B}'] -= 4
    blocks['(C)'] += 2
    print("Final state:", dict(blocks))

synthesize()