def apply_rules(blocks):
    def try_rule1():
        # [A] + [B] + [C] -> {A}
        if blocks['[A]'] >= 1 and blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
            blocks['[A]'] -= 1
            blocks['[B]'] -= 1
            blocks['[C]'] -= 1
            blocks['{A}'] += 1
            return True
        return False

    def try_rule2():
        # [A] + [B] -> {C}
        if blocks['[A]'] >= 1 and blocks['[B]'] >= 1:
            blocks['[A]'] -= 1
            blocks['[B]'] -= 1
            blocks['{C}'] += 1
            return True
        return False

    def try_rule3():
        # [B] + [C] -> {B}
        if blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
            blocks['[B]'] -= 1
            blocks['[C]'] -= 1
            blocks['{B}'] += 1
            return True
        return False

    def try_rule4():
        # [C] + [C] -> {C}
        if blocks['[C]'] >= 2:
            blocks['[C]'] -= 2
            blocks['{C}'] += 1
            return True
        return False

    def try_rule5():
        # {A} + {C} -> (A) + (B)
        if blocks['{A}'] >= 1 and blocks['{C}'] >= 1:
            blocks['{A}'] -= 1
            blocks['{C}'] -= 1
            blocks['(A)'] += 1
            blocks['(B)'] += 1
            return True
        return False

    def try_rule6():
        # {B} + {B} -> (C)
        if blocks['{B}'] >= 2:
            blocks['{B}'] -= 2
            blocks['(C)'] += 1
            return True
        return False

    changed = True
    while changed:
        changed = False
        rules = [try_rule1, try_rule2, try_rule3, try_rule4, try_rule5, try_rule6]
        for rule in rules:
            if rule():
                changed = True
                break

blocks = {
    '[A]': 3, '[B]': 6, '[C]': 7,
    '{A}': 0, '{B}': 0, '{C}': 0,
    '(A)': 0, '(B)': 0, '(C)': 0
}

apply_rules(blocks)
result = ""
for block_type in ['[A]', '[B]', '[C]', '{A}', '{B}', '{C}', '(A)', '(B)', '(C)']:
    if blocks[block_type] > 0:
        result += f"{blocks[block_type]} {block_type}, "
print(result[:-2])