# Initial counts of blocks
A = 1
B = 6
C = 6
A_curly = 0
B_curly = 0
C_curly = 0
A_paren = 0
B_paren = 0
C_paren = 0

# Apply methods until no more changes can be made
while True:
    changes = False
    
    # Method 1: [A] + [B] + [C] -> {A}
    if A >= 1 and B >= 1 and C >= 1:
        A -= 1
        B -= 1
        C -= 1
        A_curly += 1
        changes = True
    
    # Method 2: [A] + [B] -> {C}
    if A >= 1 and B >= 1:
        A -= 1
        B -= 1
        C_curly += 1
        changes = True
    
    # Method 3: [B] + [C] -> {B}
    while B >= 1 and C >= 1:
        B -= 1
        C -= 1
        B_curly += 1
        changes = True
    
    # Method 4: [C] + [C] -> {C}
    while C >= 2:
        C -= 2
        C_curly += 1
        changes = True
    
    # Method 5: {A} + {C} -> (A) + (B)
    while A_curly >= 1 and C_curly >= 1:
        A_curly -= 1
        C_curly -= 1
        A_paren += 1
        B_paren += 1
        changes = True
    
    # Method 6: {B} + {B} -> (C)
    while B_curly >= 2:
        B_curly -= 2
        C_paren += 1
        changes = True
    
    # If no changes were made, break the loop
    if not changes:
        break

# Output the final state of the blocks
print(f"{A} [A], {B} [B], {C} [C], {A_curly} {{A}}, {B_curly} {{B}}, {C_curly} {{C}}, {A_paren} (A), {B_paren} (B), {C_paren} (C)")