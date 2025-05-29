# Initialize the blocks
blocks = {
    '[A]': 3,
    '[B]': 6,
    '[C]': 2,
    '{A}': 0,
    '{B}': 0,
    '{C}': 0,
    '(A)': 0,
    '(B)': 0,
    '(C)': 0
}

# Define a function to apply the rules
def apply_rules(blocks):
    changes = True
    while changes:
        changes = False
        # Rule 1: One [A], one [B], and one [C] can be combined to form one {A}
        while blocks['[A]'] >= 1 and blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
            blocks['[A]'] -= 1
            blocks['[B]'] -= 1
            blocks['[C]'] -= 1
            blocks['{A}'] += 1
            changes = True
        
        # Rule 2: One [A] and one [B] can be combined to form one {C}
        while blocks['[A]'] >= 1 and blocks['[B]'] >= 1:
            blocks['[A]'] -= 1
            blocks['[B]'] -= 1
            blocks['{C}'] += 1
            changes = True
        
        # Rule 3: One [B] and one [C] can be combined to form one {B}
        while blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
            blocks['[B]'] -= 1
            blocks['[C]'] -= 1
            blocks['{B}'] += 1
            changes = True
        
        # Rule 4: Two [C] can be combined to form one {C}
        while blocks['[C]'] >= 2:
            blocks['[C]'] -= 2
            blocks['{C}'] += 1
            changes = True
        
        # Rule 5: One {A} and one {C} can be combined to form one (A) and one (B)
        while blocks['{A}'] >= 1 and blocks['{C}'] >= 1:
            blocks['{A}'] -= 1
            blocks['{C}'] -= 1
            blocks['(A)'] += 1
            blocks['(B)'] += 1
            changes = True
        
        # Rule 6: Two {B} can be combined to form one (C)
        while blocks['{B}'] >= 2:
            blocks['{B}'] -= 2
            blocks['(C)'] += 1
            changes = True

# Apply the rules
apply_rules(blocks)

# Print the final count of each block
print(blocks)