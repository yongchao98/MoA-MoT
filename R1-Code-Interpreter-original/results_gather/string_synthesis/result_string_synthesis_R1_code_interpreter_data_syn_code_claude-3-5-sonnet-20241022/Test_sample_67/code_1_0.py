def apply_rule1(blocks):
    # [A] + [B] + [C] -> {A}
    if blocks['[A]'] >= 1 and blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
        blocks['[A]'] -= 1
        blocks['[B]'] -= 1
        blocks['[C]'] -= 1
        blocks['{A}'] += 1
        return True, "Rule 1: [A] + [B] + [C] -> {A}"
    return False, ""

def apply_rule2(blocks):
    # [A] + [B] -> {C}
    if blocks['[A]'] >= 1 and blocks['[B]'] >= 1:
        blocks['[A]'] -= 1
        blocks['[B]'] -= 1
        blocks['{C}'] += 1
        return True, "Rule 2: [A] + [B] -> {C}"
    return False, ""

def apply_rule3(blocks):
    # [B] + [C] -> {B}
    if blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
        blocks['[B]'] -= 1
        blocks['[C]'] -= 1
        blocks['{B}'] += 1
        return True, "Rule 3: [B] + [C] -> {B}"
    return False, ""

def apply_rule4(blocks):
    # [C] + [C] -> {C}
    if blocks['[C]'] >= 2:
        blocks['[C]'] -= 2
        blocks['{C}'] += 1
        return True, "Rule 4: [C] + [C] -> {C}"
    return False, ""

def apply_rule5(blocks):
    # {A} + {C} -> (A) + (B)
    if blocks['{A}'] >= 1 and blocks['{C}'] >= 1:
        blocks['{A}'] -= 1
        blocks['{C}'] -= 1
        blocks['(A)'] += 1
        blocks['(B)'] += 1
        return True, "Rule 5: {A} + {C} -> (A) + (B)"
    return False, ""

def apply_rule6(blocks):
    # {B} + {B} -> (C)
    if blocks['{B}'] >= 2:
        blocks['{B}'] -= 2
        blocks['(C)'] += 1
        return True, "Rule 6: {B} + {B} -> (C)"
    return False, ""

def synthesize():
    blocks = {
        '[A]': 3, '[B]': 5, '[C]': 7,
        '{A}': 0, '{B}': 0, '{C}': 0,
        '(A)': 0, '(B)': 0, '(C)': 0
    }
    
    step = 1
    print(f"Initial: {blocks}")
    
    while True:
        changes = False
        for rule_func in [apply_rule1, apply_rule2, apply_rule3, apply_rule4, apply_rule5, apply_rule6]:
            while True:
                changed, message = rule_func(blocks)
                if changed:
                    changes = True
                    print(f"Step {step}: {message}")
                    print(f"State: {blocks}")
                    step += 1
                else:
                    break
        if not changes:
            break
    
    result = ""
    for block_type in ['[A]', '[B]', '[C]', '{A}', '{B}', '{C}', '(A)', '(B)', '(C)']:
        if blocks[block_type] > 0:
            result += f"{blocks[block_type]} {block_type}, "
    print(f"\nFinal result: {result[:-2]}")

synthesize()