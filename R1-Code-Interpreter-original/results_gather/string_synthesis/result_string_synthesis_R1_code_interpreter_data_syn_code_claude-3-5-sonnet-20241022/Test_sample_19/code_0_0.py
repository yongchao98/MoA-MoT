def can_apply_rule(blocks, rule_num):
    if rule_num == 1:
        return blocks['[A]'] >= 1 and blocks['[B]'] >= 1 and blocks['[C]'] >= 1
    elif rule_num == 2:
        return blocks['[A]'] >= 1 and blocks['[B]'] >= 1
    elif rule_num == 3:
        return blocks['[B]'] >= 1 and blocks['[C]'] >= 1
    elif rule_num == 4:
        return blocks['[C]'] >= 2
    elif rule_num == 5:
        return blocks['{A}'] >= 1 and blocks['{C}'] >= 1
    elif rule_num == 6:
        return blocks['{B}'] >= 2
    return False

def apply_rule(blocks, rule_num):
    if rule_num == 1:
        blocks['[A]'] -= 1
        blocks['[B]'] -= 1
        blocks['[C]'] -= 1
        blocks['{A}'] += 1
    elif rule_num == 2:
        blocks['[A]'] -= 1
        blocks['[B]'] -= 1
        blocks['{C}'] += 1
    elif rule_num == 3:
        blocks['[B]'] -= 1
        blocks['[C]'] -= 1
        blocks['{B}'] += 1
    elif rule_num == 4:
        blocks['[C]'] -= 2
        blocks['{C}'] += 1
    elif rule_num == 5:
        blocks['{A}'] -= 1
        blocks['{C}'] -= 1
        blocks['(A)'] += 1
        blocks['(B)'] += 1
    elif rule_num == 6:
        blocks['{B}'] -= 2
        blocks['(C)'] += 1

def synthesize(initial_blocks):
    blocks = initial_blocks.copy()
    while True:
        rule_applied = False
        for rule in range(1, 7):
            while can_apply_rule(blocks, rule):
                apply_rule(blocks, rule)
                rule_applied = True
        if not rule_applied:
            break
    return blocks

# Initial blocks
initial_blocks = {
    '[A]': 3, '[B]': 7, '[C]': 2,
    '{A}': 0, '{B}': 0, '{C}': 0,
    '(A)': 0, '(B)': 0, '(C)': 0
}

final_blocks = synthesize(initial_blocks)
result = ""
for block, count in final_blocks.items():
    if count > 0:
        result += f"{count} {block}, "
result = result.rstrip(", ")
print(result)