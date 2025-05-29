# Initial counts
A = 4
B = 7
C = 2
curly_A = 0
curly_B = 0
curly_C = 0
paren_A = 0
paren_B = 0
paren_C = 0

# Apply methods until no more changes can be made
while True:
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
    
    # If no changes were made, break the loop
    if not changes:
        break

# Output the final counts
print(f"[A]: {A}, [B]: {B}, [C]: {C}, {{A}}: {curly_A}, {{B}}: {curly_B}, {{C}}: {curly_C}, (A): {paren_A}, (B): {paren_B}, (C): {paren_C}")