# Initial counts of blocks
A = 1
B = 5
C = 7
A_curly = 0
B_curly = 0
C_curly = 0
A_paren = 0
B_paren = 0
C_paren = 0

# Function to apply the methods
def apply_methods():
    global A, B, C, A_curly, B_curly, C_curly, A_paren, B_paren, C_paren
    changes = True
    while changes:
        changes = False
        # Method 1
        while A >= 1 and B >= 1 and C >= 1:
            A -= 1
            B -= 1
            C -= 1
            A_curly += 1
            changes = True
        # Method 2
        while A >= 1 and B >= 1:
            A -= 1
            B -= 1
            C_curly += 1
            changes = True
        # Method 3
        while B >= 1 and C >= 1:
            B -= 1
            C -= 1
            B_curly += 1
            changes = True
        # Method 4
        while C >= 2:
            C -= 2
            C_curly += 1
            changes = True
        # Method 5
        while A_curly >= 1 and C_curly >= 1:
            A_curly -= 1
            C_curly -= 1
            A_paren += 1
            B_paren += 1
            changes = True
        # Method 6
        while B_curly >= 2:
            B_curly -= 2
            C_paren += 1
            changes = True

apply_methods()

# Output the final state
print(f"[A]: {A}, [B]: {B}, [C]: {C}, {{A}}: {A_curly}, {{B}}: {B_curly}, {{C}}: {C_curly}, (A): {A_paren}, (B): {B_paren}, (C): {C_paren}")