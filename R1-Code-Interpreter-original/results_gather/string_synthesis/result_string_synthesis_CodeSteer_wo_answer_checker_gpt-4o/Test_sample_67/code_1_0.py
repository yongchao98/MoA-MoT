# Initial counts of blocks
A = 3
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

# Apply the methods
apply_methods()

# Print the final counts
print(f"Your answer: {paren_A}{paren_B}{paren_C}")