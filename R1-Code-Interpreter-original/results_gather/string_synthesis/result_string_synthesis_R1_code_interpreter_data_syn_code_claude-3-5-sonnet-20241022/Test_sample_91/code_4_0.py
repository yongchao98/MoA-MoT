def simulate():
    # Initialize blocks
    blocks = {
        '[A]': 1, '[B]': 6, '[C]': 6,
        '{A}': 0, '{B}': 0, '{C}': 0,
        '(A)': 0, '(B)': 0, '(C)': 0
    }
    
    def apply_rule(rule_num):
        nonlocal blocks
        if rule_num == 1:  # [A] + [B] + [C] -> {A}
            if blocks['[A]'] >= 1 and blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
                blocks['[A]'] -= 1
                blocks['[B]'] -= 1
                blocks['[C]'] -= 1
                blocks['{A}'] += 1
                return True
        elif rule_num == 2:  # [A] + [B] -> {C}
            if blocks['[A]'] >= 1 and blocks['[B]'] >= 1:
                blocks['[A]'] -= 1
                blocks['[B]'] -= 1
                blocks['{C}'] += 1
                return True
        elif rule_num == 3:  # [B] + [C] -> {B}
            if blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
                blocks['[B]'] -= 1
                blocks['[C]'] -= 1
                blocks['{B}'] += 1
                return True
        elif rule_num == 4:  # [C] + [C] -> {C}
            if blocks['[C]'] >= 2:
                blocks['[C]'] -= 2
                blocks['{C}'] += 1
                return True
        elif rule_num == 5:  # {A} + {C} -> (A) + (B)
            if blocks['{A}'] >= 1 and blocks['{C}'] >= 1:
                blocks['{A}'] -= 1
                blocks['{C}'] -= 1
                blocks['(A)'] += 1
                blocks['(B)'] += 1
                return True
        elif rule_num == 6:  # {B} + {B} -> (C)
            if blocks['{B}'] >= 2:
                blocks['{B}'] -= 2
                blocks['(C)'] += 1
                return True
        return False

    # First phase: Convert [C] pairs to {C}
    while apply_rule(4): pass
    
    # Second phase: Convert [A] and [B] to {C}
    while apply_rule(2): pass
    
    # Third phase: Convert [B] and remaining [C] to {B}
    while apply_rule(3): pass
    
    # Fourth phase: Convert {B} pairs to (C)
    while apply_rule(6): pass
    
    # Fifth phase: Convert {A} and {C} to (A) and (B)
    while apply_rule(5): pass
    
    result = ""
    for block, count in sorted(blocks.items()):
        if count > 0:
            result += f"{count} {block}, "
    print(result[:-2])

simulate()