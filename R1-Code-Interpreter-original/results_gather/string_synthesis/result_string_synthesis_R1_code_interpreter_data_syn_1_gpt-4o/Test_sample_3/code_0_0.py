# Initialize the blocks
blocks = {
    '[A]': 5,
    '[B]': 1,
    '[C]': 1,
    '{A}': 0,
    '{B}': 0,
    '{C}': 0,
    '(A)': 0,
    '(B)': 0,
    '(C)': 0
}

# Define a flag to check if any changes were made
changes_made = True

# Cycle through the methods until no more changes can be made
while changes_made:
    changes_made = False
    
    # Method 1: One [A], one [B], and one [C] can be combined to form one {A}
    while blocks['[A]'] >= 1 and blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
        blocks['[A]'] -= 1
        blocks['[B]'] -= 1
        blocks['[C]'] -= 1
        blocks['{A}'] += 1
        changes_made = True
    
    # Method 2: One [A] and one [B] can be combined to form one {C}
    while blocks['[A]'] >= 1 and blocks['[B]'] >= 1:
        blocks['[A]'] -= 1
        blocks['[B]'] -= 1
        blocks['{C}'] += 1
        changes_made = True
    
    # Method 3: One [B] and one [C] can be combined to form one {B}
    while blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
        blocks['[B]'] -= 1
        blocks['[C]'] -= 1
        blocks['{B}'] += 1
        changes_made = True
    
    # Method 4: Two [C] can be combined to form one {C}
    while blocks['[C]'] >= 2:
        blocks['[C]'] -= 2
        blocks['{C}'] += 1
        changes_made = True
    
    # Method 5: One {A} and one {C} can be combined to form one (A) and one (B)
    while blocks['{A}'] >= 1 and blocks['{C}'] >= 1:
        blocks['{A}'] -= 1
        blocks['{C}'] -= 1
        blocks['(A)'] += 1
        blocks['(B)'] += 1
        changes_made = True
    
    # Method 6: Two {B} can be combined to form one (C)
    while blocks['{B}'] >= 2:
        blocks['{B}'] -= 2
        blocks['(C)'] += 1
        changes_made = True

# Print the final count of each block
print(blocks)