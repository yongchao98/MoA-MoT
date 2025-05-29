# Initial counts of blocks
A = 1
B = 1
C = 7
curly_A = 0
curly_B = 0
curly_C = 0
paren_A = 0
paren_B = 0
paren_C = 0

# Function to apply rules
def apply_rules():
    global A, B, C, curly_A, curly_B, curly_C, paren_A, paren_B, paren_C
    changes = True
    while changes:
        changes = False
        # Rule 4: Two [C] can be combined to form one {C}.
        while C >= 2:
            C -= 2
            curly_C += 1
            changes = True
        # Rule 2: One [A] and one [B] can be combined to form one {C}.
        if A >= 1 and B >= 1:
            A -= 1
            B -= 1
            curly_C += 1
            changes = True
        # Rule 3: One [B] and one [C] can be combined to form one {B}.
        if B >= 1 and C >= 1:
            B -= 1
            C -= 1
            curly_B += 1
            changes = True
        # Rule 1: One [A], one [B], and one [C] can be combined to form one {A}.
        if A >= 1 and B >= 1 and C >= 1:
            A -= 1
            B -= 1
            C -= 1
            curly_A += 1
            changes = True
        # Rule 5: One {A} and one {C} can be combined to form one (A) and one (B).
        if curly_A >= 1 and curly_C >= 1:
            curly_A -= 1
            curly_C -= 1
            paren_A += 1
            paren_B += 1
            changes = True
        # Rule 6: Two {B} can be combined to form one (C).
        while curly_B >= 2:
            curly_B -= 2
            paren_C += 1
            changes = True

apply_rules()

# Output the final counts of blocks
print(f"[A]: {A}, [B]: {B}, [C]: {C}, {{A}}: {curly_A}, {{B}}: {curly_B}, {{C}}: {curly_C}, (A): {paren_A}, (B): {paren_B}, (C): {paren_C}")