# Initial counts of blocks
blocks = {
    '[A]': 4,
    '[B]': 1,
    '[C]': 5,
    '{A}': 0,
    '{B}': 0,
    '{C}': 0,
    '(A)': 0,
    '(B)': 0,
    '(C)': 0
}

# Function to apply rules
def apply_rules(blocks):
    changes = True
    while changes:
        changes = False
        # Rule 1
        while blocks['[A]'] >= 1 and blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
            blocks['[A]'] -= 1
            blocks['[B]'] -= 1
            blocks['[C]'] -= 1
            blocks['{A}'] += 1
            changes = True
        # Rule 2
        while blocks['[A]'] >= 1 and blocks['[B]'] >= 1:
            blocks['[A]'] -= 1
            blocks['[B]'] -= 1
            blocks['{C}'] += 1
            changes = True
        # Rule 3
        while blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
            blocks['[B]'] -= 1
            blocks['[C]'] -= 1
            blocks['{B}'] += 1
            changes = True
        # Rule 4
        while blocks['[C]'] >= 2:
            blocks['[C]'] -= 2
            blocks['{C}'] += 1
            changes = True
        # Rule 5
        while blocks['{A}'] >= 1 and blocks['{C}'] >= 1:
            blocks['{A}'] -= 1
            blocks['{C}'] -= 1
            blocks['(A)'] += 1
            blocks['(B)'] += 1
            changes = True
        # Rule 6
        while blocks['{B}'] >= 2:
            blocks['{B}'] -= 2
            blocks['(C)'] += 1
            changes = True

apply_rules(blocks)

# Output the final configuration of blocks
print(blocks)