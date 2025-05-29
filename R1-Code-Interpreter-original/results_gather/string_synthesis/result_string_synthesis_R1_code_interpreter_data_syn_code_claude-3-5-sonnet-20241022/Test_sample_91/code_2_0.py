def apply_method_1(blocks):
    # [A] + [B] + [C] -> {A}
    if blocks['[A]'] >= 1 and blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
        blocks['[A]'] -= 1
        blocks['[B]'] -= 1
        blocks['[C]'] -= 1
        blocks['{A}'] += 1
        return True
    return False

def apply_method_2(blocks):
    # [A] + [B] -> {C}
    if blocks['[A]'] >= 1 and blocks['[B]'] >= 1:
        blocks['[A]'] -= 1
        blocks['[B]'] -= 1
        blocks['{C}'] += 1
        return True
    return False

def apply_method_3(blocks):
    # [B] + [C] -> {B}
    if blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
        blocks['[B]'] -= 1
        blocks['[C]'] -= 1
        blocks['{B}'] += 1
        return True
    return False

def apply_method_4(blocks):
    # [C] + [C] -> {C}
    if blocks['[C]'] >= 2:
        blocks['[C]'] -= 2
        blocks['{C}'] += 1
        return True
    return False

def apply_method_5(blocks):
    # {A} + {C} -> (A) + (B)
    if blocks['{A}'] >= 1 and blocks['{C}'] >= 1:
        blocks['{A}'] -= 1
        blocks['{C}'] -= 1
        blocks['(A)'] += 1
        blocks['(B)'] += 1
        return True
    return False

def apply_method_6(blocks):
    # {B} + {B} -> (C)
    if blocks['{B}'] >= 2:
        blocks['{B}'] -= 2
        blocks['(C)'] += 1
        return True
    return False

def simulate():
    blocks = {
        '[A]': 1, '[B]': 6, '[C]': 6,
        '{A}': 0, '{B}': 0, '{C}': 0,
        '(A)': 0, '(B)': 0, '(C)': 0
    }
    
    while True:
        old_state = blocks.copy()
        
        # Apply methods in specific order to maximize combinations
        # First, create as many {C} as possible using method 4
        while apply_method_4(blocks): pass
        
        # Then try to create {A} using method 1
        while apply_method_1(blocks): pass
        
        # Create more {C} using method 2
        while apply_method_2(blocks): pass
        
        # Create {B} using method 3
        while apply_method_3(blocks): pass
        
        # Create (A) and (B) using method 5
        while apply_method_5(blocks): pass
        
        # Finally, create (C) using method 6
        while apply_method_6(blocks): pass
        
        # If no changes were made in this iteration, we're done
        if blocks == old_state:
            break
    
    result = ""
    for block, count in sorted(blocks.items()):
        if count > 0:
            result += f"{count} {block}, "
    print(result[:-2])

simulate()