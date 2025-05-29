# Initialize the counts of each block type
blocks = {
    '[A]': 1,
    '[B]': 5,
    '[C]': 1,
    '{A}': 0,
    '{B}': 0,
    '{C}': 0,
    '(A)': 0,
    '(B)': 0,
    '(C)': 0
}

# Function to apply the rules
def apply_rules(blocks):
    changes = True
    while changes:
        changes = False
        
        # Rule 1: [A] + [B] + [C] -> {A}
        while blocks['[A]'] >= 1 and blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
            blocks['[A]'] -= 1
            blocks['[B]'] -= 1
            blocks['[C]'] -= 1
            blocks['{A}'] += 1
            changes = True
        
        # Rule 2: [A] + [B] -> {C}
        while blocks['[A]'] >= 1 and blocks['[B]'] >= 1:
            blocks['[A]'] -= 1
            blocks['[B]'] -= 1
            blocks['{C}'] += 1
            changes = True
        
        # Rule 3: [B] + [C] -> {B}
        while blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
            blocks['[B]'] -= 1
            blocks['[C]'] -= 1
            blocks['{B}'] += 1
            changes = True
        
        # Rule 4: 2 * [C] -> {C}
        while blocks['[C]'] >= 2:
            blocks['[C]'] -= 2
            blocks['{C}'] += 1
            changes = True
        
        # Rule 5: {A} + {C} -> (A) + (B)
        while blocks['{A}'] >= 1 and blocks['{C}'] >= 1:
            blocks['{A}'] -= 1
            blocks['{C}'] -= 1
            blocks['(A)'] += 1
            blocks['(B)'] += 1
            changes = True
        
        # Rule 6: 2 * {B} -> (C)
        while blocks['{B}'] >= 2:
            blocks['{B}'] -= 2
            blocks['(C)'] += 1
            changes = True

# Apply the rules
apply_rules(blocks)

# Format the output
output = ''.join(f'{count}{block}' for block, count in blocks.items() if count > 0)
print(f'Your answer: {output}')