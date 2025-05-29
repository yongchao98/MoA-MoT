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

def synthesize(initial_blocks):
    blocks = initial_blocks.copy()
    while True:
        any_change = False
        # Try each method in order
        if apply_method_1(blocks): any_change = True
        elif apply_method_2(blocks): any_change = True
        elif apply_method_3(blocks): any_change = True
        elif apply_method_4(blocks): any_change = True
        elif apply_method_5(blocks): any_change = True
        elif apply_method_6(blocks): any_change = True
        
        if not any_change:
            break
    
    return blocks

# Initial blocks
initial_blocks = {
    '[A]': 4, '[B]': 1, '[C]': 5,
    '{A}': 0, '{B}': 0, '{C}': 0,
    '(A)': 0, '(B)': 0, '(C)': 0
}

final_blocks = synthesize(initial_blocks)
result = ""
for block_type in ['[A]', '[B]', '[C]', '{A}', '{B}', '{C}', '(A)', '(B)', '(C)']:
    if final_blocks[block_type] > 0:
        result += f"{final_blocks[block_type]} {block_type}, "
result = result.rstrip(", ")
print(result)