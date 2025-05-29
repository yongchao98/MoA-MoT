# Initialize the counts of each block type
blocks = {
    '[A]': 3,
    '[B]': 6,
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
    changes = False
    
    # Rule 1: One [A], one [B], and one [C] can be combined to form one {A}
    if blocks['[A]'] >= 1 and blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
        blocks['[A]'] -= 1
        blocks['[B]'] -= 1
        blocks['[C]'] -= 1
        blocks['{A}'] += 1
        changes = True
    
    # Rule 2: One [A] and one [B] can be combined to form one {C}
    if blocks['[A]'] >= 1 and blocks['[B]'] >= 1:
        blocks['[A]'] -= 1
        blocks['[B]'] -= 1
        blocks['{C}'] += 1
        changes = True
    
    # Rule 3: One [B] and one [C] can be combined to form one {B}
    if blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
        blocks['[B]'] -= 1
        blocks['[C]'] -= 1
        blocks['{B}'] += 1
        changes = True
    
    # Rule 4: Two [C] can be combined to form one {C}
    if blocks['[C]'] >= 2:
        blocks['[C]'] -= 2
        blocks['{C}'] += 1
        changes = True
    
    # Rule 5: One {A} and one {C} can be combined to form one (A) and one (B)
    if blocks['{A}'] >= 1 and blocks['{C}'] >= 1:
        blocks['{A}'] -= 1
        blocks['{C}'] -= 1
        blocks['(A)'] += 1
        blocks['(B)'] += 1
        changes = True
    
    # Rule 6: Two {B} can be combined to form one (C)
    if blocks['{B}'] >= 2:
        blocks['{B}'] -= 2
        blocks['(C)'] += 1
        changes = True
    
    return changes

# Apply the rules until no more changes can be made
while apply_rules(blocks):
    pass

# Print the final counts of each block type
result = ''.join(f'{count}{block}' for block, count in blocks.items() if count > 0)
print(f'Your answer: {result}')