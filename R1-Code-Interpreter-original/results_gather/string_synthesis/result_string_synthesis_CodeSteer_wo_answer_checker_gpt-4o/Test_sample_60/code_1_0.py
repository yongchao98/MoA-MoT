# Initialize the counts of each block type
blocks = {
    '[A]': 1,
    '[B]': 1,
    '[C]': 3,
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
        # Rule 1
        if blocks['[A]'] >= 1 and blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
            blocks['[A]'] -= 1
            blocks['[B]'] -= 1
            blocks['[C]'] -= 1
            blocks['{A}'] += 1
            changes = True
        # Rule 2
        if blocks['[A]'] >= 1 and blocks['[B]'] >= 1:
            blocks['[A]'] -= 1
            blocks['[B]'] -= 1
            blocks['{C}'] += 1
            changes = True
        # Rule 3
        if blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
            blocks['[B]'] -= 1
            blocks['[C]'] -= 1
            blocks['{B}'] += 1
            changes = True
        # Rule 4
        if blocks['[C]'] >= 2:
            blocks['[C]'] -= 2
            blocks['{C}'] += 1
            changes = True
        # Rule 5
        if blocks['{A}'] >= 1 and blocks['{C}'] >= 1:
            blocks['{A}'] -= 1
            blocks['{C}'] -= 1
            blocks['(A)'] += 1
            blocks['(B)'] += 1
            changes = True
        # Rule 6
        if blocks['{B}'] >= 2:
            blocks['{B}'] -= 2
            blocks['(C)'] += 1
            changes = True

# Apply the rules
apply_rules(blocks)

# Prepare the final answer
final_answer = ''.join(f'{count}{block}' for block, count in blocks.items() if count > 0)

# Output the final answer
print(f'Your answer: {final_answer}')