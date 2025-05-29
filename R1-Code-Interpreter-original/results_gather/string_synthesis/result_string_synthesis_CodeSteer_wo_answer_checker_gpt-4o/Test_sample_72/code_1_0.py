# Initialize the counts of each block type
blocks = {
    '[A]': 3,
    '[B]': 6,
    '[C]': 7,
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
        # Rule 1: One [A], one [B], and one [C] -> one {A}
        while blocks['[A]'] >= 1 and blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
            blocks['[A]'] -= 1
            blocks['[B]'] -= 1
            blocks['[C]'] -= 1
            blocks['{A}'] += 1
            changes = True

        # Rule 2: One [A] and one [B] -> one {C}
        while blocks['[A]'] >= 1 and blocks['[B]'] >= 1:
            blocks['[A]'] -= 1
            blocks['[B]'] -= 1
            blocks['{C}'] += 1
            changes = True

        # Rule 3: One [B] and one [C] -> one {B}
        while blocks['[B]'] >= 1 and blocks['[C]'] >= 1:
            blocks['[B]'] -= 1
            blocks['[C]'] -= 1
            blocks['{B}'] += 1
            changes = True

        # Rule 4: Two [C] -> one {C}
        while blocks['[C]'] >= 2:
            blocks['[C]'] -= 2
            blocks['{C}'] += 1
            changes = True

        # Rule 5: One {A} and one {C} -> one (A) and one (B)
        while blocks['{A}'] >= 1 and blocks['{C}'] >= 1:
            blocks['{A}'] -= 1
            blocks['{C}'] -= 1
            blocks['(A)'] += 1
            blocks['(B)'] += 1
            changes = True

        # Rule 6: Two {B} -> one (C)
        while blocks['{B}'] >= 2:
            blocks['{B}'] -= 2
            blocks['(C)'] += 1
            changes = True

# Apply the rules
apply_rules(blocks)

# Format the output
output = f"{blocks['[A]']}{blocks['[B]']}{blocks['[C]']}{blocks['{A}']}{blocks['{B}']}{blocks['{C}']}{blocks['(A)']}{blocks['(B)']}{blocks['(C)']}"
print(f"Your answer: {output}")