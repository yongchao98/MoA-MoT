# Initialize the blocks
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

# Define a function to apply the methods
def apply_methods(blocks):
    changes = True
    while changes:
        changes = False
        # Method 1: One [A], one [B], and one [C] can be combined to form one {A}
        while blocks['[A]'] >= 1 and blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
            blocks['[A]'] -= 1
            blocks['[B]'] -= 1
            blocks['[C]'] -= 1
            blocks['{A}'] += 1
            changes = True
        
        # Method 2: One [A] and one [B] can be combined to form one {C}
        while blocks['[A]'] >= 1 and blocks['[B]'] >= 1:
            blocks['[A]'] -= 1
            blocks['[B]'] -= 1
            blocks['{C}'] += 1
            changes = True
        
        # Method 3: One [B] and one [C] can be combined to form one {B}
        while blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
            blocks['[B]'] -= 1
            blocks['[C]'] -= 1
            blocks['{B}'] += 1
            changes = True
        
        # Method 4: Two [C] can be combined to form one {C}
        while blocks['[C]'] >= 2:
            blocks['[C]'] -= 2
            blocks['{C}'] += 1
            changes = True
        
        # Method 5: One {A} and one {C} can be combined to form one (A) and one (B)
        while blocks['{A}'] >= 1 and blocks['{C}'] >= 1:
            blocks['{A}'] -= 1
            blocks['{C}'] -= 1
            blocks['(A)'] += 1
            blocks['(B)'] += 1
            changes = True
        
        # Method 6: Two {B} can be combined to form one (C)
        while blocks['{B}'] >= 2:
            blocks['{B}'] -= 2
            blocks['(C)'] += 1
            changes = True

    return blocks

# Apply the methods
final_blocks = apply_methods(blocks)

# Prepare the final answer
answer = f"{final_blocks['[A]']} [A], {final_blocks['[B]']} [B], {final_blocks['[C]']} [C], {final_blocks['{A}']} {{A}}, {final_blocks['{B}']} {{B}}, {final_blocks['{C}']} {{C}}, {final_blocks['(A)']} (A), {final_blocks['(B)']} (B), {final_blocks['(C)']} (C)"
print(answer)