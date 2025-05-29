# Initial counts of blocks
A = 2
B = 4
C = 1
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
        if A >= 1 and B >= 1 and C >= 1:
            A -= 1
            B -= 1
            C -= 1
            A_curly += 1
            changes = True
        # Method 2
        if A >= 1 and B >= 1:
            A -= 1
            B -= 1
            C_curly += 1
            changes = True
        # Method 3
        if B >= 1 and C >= 1:
            B -= 1
            C -= 1
            B_curly += 1
            changes = True
        # Method 4
        if C >= 2:
            C -= 2
            C_curly += 1
            changes = True
        # Method 5
        if A_curly >= 1 and C_curly >= 1:
            A_curly -= 1
            C_curly -= 1
            A_paren += 1
            B_paren += 1
            changes = True
        # Method 6
        if B_curly >= 2:
            B_curly -= 2
            C_paren += 1
            changes = True

apply_methods()

# Output the final configuration
print(f"({A_paren})({B_paren})({C_paren})")