def apply_rules(blocks):
    def try_rule(rule_num):
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

    changed = True
    while changed:
        changed = False
        # First, try to combine [C]s (Rule 4)
        while try_rule(4):
            changed = True
        
        # Then try Rule 1 to create {A}
        while try_rule(1):
            changed = True
            
        # Then Rule 2 to create {C}
        while try_rule(2):
            changed = True
            
        # Then Rule 3 to create {B}
        while try_rule(3):
            changed = True
            
        # Then Rule 5 to create (A) and (B)
        while try_rule(5):
            changed = True
            
        # Finally Rule 6 to create (C)
        while try_rule(6):
            changed = True

blocks = {
    '[A]': 1, '[B]': 6, '[C]': 6,
    '{A}': 0, '{B}': 0, '{C}': 0,
    '(A)': 0, '(B)': 0, '(C)': 0
}

apply_rules(blocks)

result = ""
for block, count in sorted(blocks.items()):
    if count > 0:
        result += f"{count} {block}, "
print(result[:-2])