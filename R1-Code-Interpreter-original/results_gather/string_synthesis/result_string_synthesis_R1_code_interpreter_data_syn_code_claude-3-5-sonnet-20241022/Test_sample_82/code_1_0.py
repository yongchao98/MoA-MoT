def apply_method_1(blocks):
    if blocks['[A]'] >= 1 and blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
        blocks['[A]'] -= 1
        blocks['[B]'] -= 1
        blocks['[C]'] -= 1
        blocks['{A}'] += 1
        return True
    return False

def apply_method_2(blocks):
    if blocks['[A]'] >= 1 and blocks['[B]'] >= 1:
        blocks['[A]'] -= 1
        blocks['[B]'] -= 1
        blocks['{C}'] += 1
        return True
    return False

def apply_method_3(blocks):
    if blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
        blocks['[B]'] -= 1
        blocks['[C]'] -= 1
        blocks['{B}'] += 1
        return True
    return False

def apply_method_4(blocks):
    if blocks['[C]'] >= 2:
        blocks['[C]'] -= 2
        blocks['{C}'] += 1
        return True
    return False

def apply_method_5(blocks):
    if blocks['{A}'] >= 1 and blocks['{C}'] >= 1:
        blocks['{A}'] -= 1
        blocks['{C}'] -= 1
        blocks['(A)'] += 1
        blocks['(B)'] += 1
        return True
    return False

def apply_method_6(blocks):
    if blocks['{B}'] >= 2:
        blocks['{B}'] -= 2
        blocks['(C)'] += 1
        return True
    return False

def synthesize(blocks):
    changed = True
    while changed:
        changed = False
        for method in [apply_method_1, apply_method_2, apply_method_3, 
                      apply_method_4, apply_method_5, apply_method_6]:
            if method(blocks):
                changed = True
                break

initial = {
    '[A]': 1, '[B]': 7, '[C]': 4,
    '{A}': 0, '{B}': 0, '{C}': 0,
    '(A)': 0, '(B)': 0, '(C)': 0
}

synthesize(initial)
result = ", ".join(f"{v} {k}" for k, v in initial.items() if v > 0)
print(result)