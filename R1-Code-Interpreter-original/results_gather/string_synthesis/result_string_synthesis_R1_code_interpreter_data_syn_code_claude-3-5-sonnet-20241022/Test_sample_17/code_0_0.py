def apply_rule_1(blocks):
    # [A] + [B] + [C] -> {A}
    if blocks['[A]'] >= 1 and blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
        blocks['[A]'] -= 1
        blocks['[B]'] -= 1
        blocks['[C]'] -= 1
        blocks['{A}'] += 1
        return True
    return False

def apply_rule_2(blocks):
    # [A] + [B] -> {C}
    if blocks['[A]'] >= 1 and blocks['[B]'] >= 1:
        blocks['[A]'] -= 1
        blocks['[B]'] -= 1
        blocks['{C}'] += 1
        return True
    return False

def apply_rule_3(blocks):
    # [B] + [C] -> {B}
    if blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
        blocks['[B]'] -= 1
        blocks['[C]'] -= 1
        blocks['{B}'] += 1
        return True
    return False

def apply_rule_4(blocks):
    # [C] + [C] -> {C}
    if blocks['[C]'] >= 2:
        blocks['[C]'] -= 2
        blocks['{C}'] += 1
        return True
    return False

def apply_rule_5(blocks):
    # {A} + {C} -> (A) + (B)
    if blocks['{A}'] >= 1 and blocks['{C}'] >= 1:
        blocks['{A}'] -= 1
        blocks['{C}'] -= 1
        blocks['(A)'] += 1
        blocks['(B)'] += 1
        return True
    return False

def apply_rule_6(blocks):
    # {B} + {B} -> (C)
    if blocks['{B}'] >= 2:
        blocks['{B}'] -= 2
        blocks['(C)'] += 1
        return True
    return False

def synthesize():
    blocks = {
        '[A]': 2, '[B]': 3, '[C]': 3,
        '{A}': 0, '{B}': 0, '{C}': 0,
        '(A)': 0, '(B)': 0, '(C)': 0
    }
    
    while True:
        changes_made = False
        for rule in [apply_rule_1, apply_rule_2, apply_rule_3, apply_rule_4, apply_rule_5, apply_rule_6]:
            if rule(blocks):
                changes_made = True
                break
        if not changes_made:
            break
    
    result = ""
    for block_type in ['[A]', '[B]', '[C]', '{A}', '{B}', '{C}', '(A)', '(B)', '(C)']:
        if blocks[block_type] > 0:
            result += f"{blocks[block_type]} {block_type}, "
    print(result[:-2])

synthesize()