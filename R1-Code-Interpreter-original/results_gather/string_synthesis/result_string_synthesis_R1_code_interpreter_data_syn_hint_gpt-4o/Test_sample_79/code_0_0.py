# Initial counts of blocks
A = 5
B = 5
C = 7
curly_A = 0
curly_B = 0
curly_C = 0
paren_A = 0
paren_B = 0
paren_C = 0

# Function to apply the methods
def apply_methods():
    global A, B, C, curly_A, curly_B, curly_C, paren_A, paren_B, paren_C
    changes = True
    while changes:
        changes = False
        # Method 1
        while A >= 1 and B >= 1 and C >= 1:
            A -= 1
            B -= 1
            C -= 1
            curly_A += 1
            changes = True
        # Method 2
        while A >= 1 and B >= 1:
            A -= 1
            B -= 1
            curly_C += 1
            changes = True
        # Method 3
        while B >= 1 and C >= 1:
            B -= 1
            C -= 1
            curly_B += 1
            changes = True
        # Method 4
        while C >= 2:
            C -= 2
            curly_C += 1
            changes = True
        # Method 5
        while curly_A >= 1 and curly_C >= 1:
            curly_A -= 1
            curly_C -= 1
            paren_A += 1
            paren_B += 1
            changes = True
        # Method 6
        while curly_B >= 2:
            curly_B -= 2
            paren_C += 1
            changes = True

apply_methods()

# Output the final counts of each block type
print(f"[A]: {A}, [B]: {B}, [C]: {C}, {{A}}: {curly_A}, {{B}}: {curly_B}, {{C}}: {curly_C}, (A): {paren_A}, (B): {paren_B}, (C): {paren_C}")