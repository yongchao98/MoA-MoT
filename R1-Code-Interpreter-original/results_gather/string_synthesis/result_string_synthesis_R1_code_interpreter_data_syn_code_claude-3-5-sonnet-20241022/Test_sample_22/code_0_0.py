def can_apply_method(blocks, method):
    if method == 1:
        return blocks['[A]'] >= 1 and blocks['[B]'] >= 1 and blocks['[C]'] >= 1
    elif method == 2:
        return blocks['[A]'] >= 1 and blocks['[B]'] >= 1
    elif method == 3:
        return blocks['[B]'] >= 1 and blocks['[C]'] >= 1
    elif method == 4:
        return blocks['[C]'] >= 2
    elif method == 5:
        return blocks['{A}'] >= 1 and blocks['{C}'] >= 1
    elif method == 6:
        return blocks['{B}'] >= 2
    return False

def apply_method(blocks, method):
    if method == 1:
        blocks['[A]'] -= 1
        blocks['[B]'] -= 1
        blocks['[C]'] -= 1
        blocks['{A}'] += 1
    elif method == 2:
        blocks['[A]'] -= 1
        blocks['[B]'] -= 1
        blocks['{C}'] += 1
    elif method == 3:
        blocks['[B]'] -= 1
        blocks['[C]'] -= 1
        blocks['{B}'] += 1
    elif method == 4:
        blocks['[C]'] -= 2
        blocks['{C}'] += 1
    elif method == 5:
        blocks['{A}'] -= 1
        blocks['{C}'] -= 1
        blocks['(A)'] += 1
        blocks['(B)'] += 1
    elif method == 6:
        blocks['{B}'] -= 2
        blocks['(C)'] += 1

def synthesize(initial_blocks):
    blocks = initial_blocks.copy()
    changed = True
    while changed:
        changed = False
        for method in range(1, 7):
            while can_apply_method(blocks, method):
                apply_method(blocks, method)
                changed = True
    return blocks

# Initial blocks
initial = {
    '[A]': 4, '[B]': 1, '[C]': 3,
    '{A}': 0, '{B}': 0, '{C}': 0,
    '(A)': 0, '(B)': 0, '(C)': 0
}

final = synthesize(initial)
result = ""
for block in sorted(final.keys()):
    if final[block] > 0:
        result += f"{final[block]} {block}, "
result = result.rstrip(", ")
print(result)