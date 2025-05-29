# Initial counts of blocks
count_A = 5
count_B = 6
count_C = 6
count_curly_A = 0
count_curly_B = 0
count_curly_C = 0
count_paren_A = 0
count_paren_B = 0
count_paren_C = 0

# Function to apply the rules
def apply_rules():
    global count_A, count_B, count_C
    global count_curly_A, count_curly_B, count_curly_C
    global count_paren_A, count_paren_B, count_paren_C
    
    changes = True
    while changes:
        changes = False
        
        # Rule 1: [A] + [B] + [C] -> {A}
        while count_A >= 1 and count_B >= 1 and count_C >= 1:
            count_A -= 1
            count_B -= 1
            count_C -= 1
            count_curly_A += 1
            changes = True
        
        # Rule 2: [A] + [B] -> {C}
        while count_A >= 1 and count_B >= 1:
            count_A -= 1
            count_B -= 1
            count_curly_C += 1
            changes = True
        
        # Rule 3: [B] + [C] -> {B}
        while count_B >= 1 and count_C >= 1:
            count_B -= 1
            count_C -= 1
            count_curly_B += 1
            changes = True
        
        # Rule 4: 2 [C] -> {C}
        while count_C >= 2:
            count_C -= 2
            count_curly_C += 1
            changes = True
        
        # Rule 5: {A} + {C} -> (A) + (B)
        while count_curly_A >= 1 and count_curly_C >= 1:
            count_curly_A -= 1
            count_curly_C -= 1
            count_paren_A += 1
            count_paren_B += 1
            changes = True
        
        # Rule 6: 2 {B} -> (C)
        while count_curly_B >= 2:
            count_curly_B -= 2
            count_paren_C += 1
            changes = True

# Apply the rules
apply_rules()

# Output the final counts
result = f"{count_A} [A], {count_B} [B], {count_C} [C], {count_curly_A} {{A}}, {count_curly_B} {{B}}, {count_curly_C} {{C}}, {count_paren_A} (A), {count_paren_B} (B), {count_paren_C} (C)"
print(result)