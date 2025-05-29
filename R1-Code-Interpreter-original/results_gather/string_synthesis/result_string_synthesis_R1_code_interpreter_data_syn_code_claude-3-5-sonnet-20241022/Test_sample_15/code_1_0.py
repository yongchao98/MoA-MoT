def apply_rule1(blocks):
    # [A] + [B] + [C] -> {A}
    if blocks['[A]'] >= 1 and blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
        blocks['[A]'] -= 1
        blocks['[B]'] -= 1
        blocks['[C]'] -= 1
        blocks['{A}'] += 1
        print("Rule 1: [A] + [B] + [C] -> {A}")
        return True
    return False

def apply_rule2(blocks):
    # [A] + [B] -> {C}
    if blocks['[A]'] >= 1 and blocks['[B]'] >= 1:
        blocks['[A]'] -= 1
        blocks['[B]'] -= 1
        blocks['{C}'] += 1
        print("Rule 2: [A] + [B] -> {C}")
        return True
    return False

def apply_rule3(blocks):
    # [B] + [C] -> {B}
    if blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
        blocks['[B]'] -= 1
        blocks['[C]'] -= 1
        blocks['{B}'] += 1
        print("Rule 3: [B] + [C] -> {B}")
        return True
    return False

def apply_rule4(blocks):
    # [C] + [C] -> {C}
    if blocks['[C]'] >= 2:
        blocks['[C]'] -= 2
        blocks['{C}'] += 1
        print("Rule 4: [C] + [C] -> {C}")
        return True
    return False

def apply_rule5(blocks):
    # {A} + {C} -> (A) + (B)
    if blocks['{A}'] >= 1 and blocks['{C}'] >= 1:
        blocks['{A}'] -= 1
        blocks['{C}'] -= 1
        blocks['(A)'] += 1
        blocks['(B)'] += 1
        print("Rule 5: {A} + {C} -> (A) + (B)")
        return True
    return False

def apply_rule6(blocks):
    # {B} + {B} -> (C)
    if blocks['{B}'] >= 2:
        blocks['{B}'] -= 2
        blocks['(C)'] += 1
        print("Rule 6: {B} + {B} -> (C)")
        return True
    return False

def synthesize():
    blocks = {
        '[A]': 1, '[B]': 5, '[C]': 7,
        '{A}': 0, '{B}': 0, '{C}': 0,
        '(A)': 0, '(B)': 0, '(C)': 0
    }
    
    print("Initial:", blocks)
    cycle = 1
    while True:
        changes_made = False
        print(f"\nCycle {cycle}:")
        for rule in [apply_rule1, apply_rule2, apply_rule3, apply_rule4, apply_rule5, apply_rule6]:
            while rule(blocks):
                changes_made = True
                print("Current state:", blocks)
        
        if not changes_made:
            break
        cycle += 1
    
    print("\nFinal state:", blocks)
    result = ""
    for block in ['[A]', '[B]', '[C]', '{A}', '{B}', '{C}', '(A)', '(B)', '(C)']:
        if blocks[block] > 0:
            result += f"{blocks[block]} {block}, "
    print("\nFinal result:", result[:-2])

synthesize()