# Initialize the blocks
blocks = {
    '[A]': 5,
    '[B]': 7,
    '[C]': 2,
    '{A}': 0,
    '{B}': 0,
    '{C}': 0,
    '(A)': 0,
    '(B)': 0,
    '(C)': 0
}

# Define a flag to check if any changes were made
changes_made = True

# Repeat the process until no more changes can be made
while changes_made:
    changes_made = False
    
    # Method 1: [A] + [B] + [C] -> {A}
    while blocks['[A]'] >= 1 and blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
        blocks['[A]'] -= 1
        blocks['[B]'] -= 1
        blocks['[C]'] -= 1
        blocks['{A}'] += 1
        changes_made = True
    
    # Method 2: [A] + [B] -> {C}
    while blocks['[A]'] >= 1 and blocks['[B]'] >= 1:
        blocks['[A]'] -= 1
        blocks['[B]'] -= 1
        blocks['{C}'] += 1
        changes_made = True
    
    # Method 3: [B] + [C] -> {B}
    while blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
        blocks['[B]'] -= 1
        blocks['[C]'] -= 1
        blocks['{B}'] += 1
        changes_made = True
    
    # Method 4: [C] + [C] -> {C}
    while blocks['[C]'] >= 2:
        blocks['[C]'] -= 2
        blocks['{C}'] += 1
        changes_made = True
    
    # Method 5: {A} + {C} -> (A) + (B)
    while blocks['{A}'] >= 1 and blocks['{C}'] >= 1:
        blocks['{A}'] -= 1
        blocks['{C}'] -= 1
        blocks['(A)'] += 1
        blocks['(B)'] += 1
        changes_made = True
    
    # Method 6: {B} + {B} -> (C)
    while blocks['{B}'] >= 2:
        blocks['{B}'] -= 2
        blocks['(C)'] += 1
        changes_made = True

# Output the final configuration of blocks
print(blocks)