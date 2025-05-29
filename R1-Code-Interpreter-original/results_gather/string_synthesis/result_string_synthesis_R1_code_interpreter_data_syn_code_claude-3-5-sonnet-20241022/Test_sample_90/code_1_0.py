def apply_rules(blocks):
    changed = True
    iterations = 0
    while changed and iterations < 100:  # Added iteration limit for safety
        changed = False
        iterations += 1
        
        # Rule 1: [A] + [B] + [C] -> {A}
        if blocks['[A]'] >= 1 and blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
            blocks['[A]'] -= 1
            blocks['[B]'] -= 1
            blocks['[C]'] -= 1
            blocks['{A}'] += 1
            changed = True
            continue
        
        # Rule 2: [A] + [B] -> {C}
        if blocks['[A]'] >= 1 and blocks['[B]'] >= 1:
            blocks['[A]'] -= 1
            blocks['[B]'] -= 1
            blocks['{C}'] += 1
            changed = True
            continue
        
        # Rule 3: [B] + [C] -> {B}
        if blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
            blocks['[B]'] -= 1
            blocks['[C]'] -= 1
            blocks['{B}'] += 1
            changed = True
            continue
        
        # Rule 4: [C] + [C] -> {C}
        if blocks['[C]'] >= 2:
            blocks['[C]'] -= 2
            blocks['{C}'] += 1
            changed = True
            continue
        
        # Rule 5: {A} + {C} -> (A) + (B)
        if blocks['{A}'] >= 1 and blocks['{C}'] >= 1:
            blocks['{A}'] -= 1
            blocks['{C}'] -= 1
            blocks['(A)'] += 1
            blocks['(B)'] += 1
            changed = True
            continue
        
        # Rule 6: {B} + {B} -> (C)
        if blocks['{B}'] >= 2:
            blocks['{B}'] -= 2
            blocks['(C)'] += 1
            changed = True
            continue

blocks = {
    '[A]': 3, '[B]': 3, '[C]': 6,
    '{A}': 0, '{B}': 0, '{C}': 0,
    '(A)': 0, '(B)': 0, '(C)': 0
}

apply_rules(blocks)

result = ""
for block, count in sorted(blocks.items()):
    if count > 0:
        result += f"{count} {block}, "
result = result.rstrip(", ")
print(result)