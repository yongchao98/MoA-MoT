# Initial counts
A = 5
B = 7
C = 1
A_curly = 0
B_curly = 0
C_curly = 0
A_paren = 0
B_paren = 0
C_paren = 0

# Apply the rules in sequence until no more changes can be made
while True:
    # Track changes
    changes = False
    
    # Rule 1: [A] + [B] + [C] -> {A}
    if A > 0 and B > 0 and C > 0:
        A -= 1
        B -= 1
        C -= 1
        A_curly += 1
        changes = True
    
    # Rule 2: [A] + [B] -> {C}
    while A > 0 and B > 0:
        A -= 1
        B -= 1
        C_curly += 1
        changes = True
    
    # Rule 3: [B] + [C] -> {B}
    while B > 0 and C > 0:
        B -= 1
        C -= 1
        B_curly += 1
        changes = True
    
    # Rule 4: [C] + [C] -> {C}
    while C > 1:
        C -= 2
        C_curly += 1
        changes = True
    
    # Rule 5: {A} + {C} -> (A) + (B)
    while A_curly > 0 and C_curly > 0:
        A_curly -= 1
        C_curly -= 1
        A_paren += 1
        B_paren += 1
        changes = True
    
    # Rule 6: {B} + {B} -> (C)
    while B_curly > 1:
        B_curly -= 2
        C_paren += 1
        changes = True
    
    # If no changes, break the loop
    if not changes:
        break

# Output the final state
print(f"({A_paren})({B_paren})({C_paren})")