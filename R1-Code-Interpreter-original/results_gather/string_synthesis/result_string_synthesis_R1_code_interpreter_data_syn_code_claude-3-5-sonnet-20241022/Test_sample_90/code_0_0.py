def apply_rules(blocks):
    changed = True
    while changed:
        changed = False
        # Store initial state to check if any changes occurred
        initial_state = blocks.copy()
        
        # Rule 1: [A] + [B] + [C] -> {A}
        while blocks['[A]'] >= 1 and blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
            blocks['[A]'] -= 1
            blocks['[B]'] -= 1
            blocks['[C]'] -= 1
            blocks['{A}'] += 1
            changed = True
        
        # Rule 2: [A] + [B] -> {C}
        while blocks['[A]'] >= 1 and blocks['[B]'] >= 1:
            blocks['[A]'] -= 1
            blocks['[B]'] -= 1
            blocks['{C}'] += 1
            changed = True
        
        # Rule 3: [B] + [C] -> {B}
        while blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
            blocks['[B]'] -= 1
            blocks['[C]'] -= 1
            blocks['{B}'] += 1
            changed = True
        
        # Rule 4: [C] + [C] -> {C}
        while blocks['[C]'] >= 2:
            blocks['[C]'] -= 2
            blocks['{C}'] += 1
            changed = True
        
        # Rule 5: {A} + {C} -> (A) + (B)
        while blocks['{A}'] >= 1 and blocks['{C}'] >= 1:
            blocks['{A}'] -= 1
            blocks['{C}'] -= 1
            blocks['(A)'] += 1
            blocks['(B)'] += 1
            changed = True
        
        # Rule 6: {B} + {B} -> (C)
        while blocks['{B}'] >= 2:
            blocks['{B}'] -= 2
            blocks['(C)'] += 1
            changed = True

# Initial blocks
blocks = {
    '[A]': 3, '[B]': 3, '[C]': 6,
    '{A}': 0, '{B}': 0, '{C}': 0,
    '(A)': 0, '(B)': 0, '(C)': 0
}

# Apply rules until no more changes
apply_rules(blocks)

# Print final result
result = ""
for block, count in blocks.items():
    if count > 0:
        result += f"{count} {block}, "
result = result.rstrip(", ")
print(result)