# Initial counts
machines_A = 3
machines_B = 6
machines_C = 6
parts_X = 0
parts_Y = 0
parts_Z = 0

# Function to apply methods
def apply_methods():
    global machines_A, machines_B, machines_C, parts_X, parts_Y, parts_Z
    changes = True
    while changes:
        changes = False
        # Method 1: Dismantle machine A
        if machines_A >= 1:
            machines_A -= 1
            parts_X += 2
            parts_Y += 1
            changes = True
        # Method 2: Dismantle machines B
        elif machines_B >= 2:
            machines_B -= 2
            parts_X += 1
            changes = True
        # Method 3: Dismantle machines C
        elif machines_C >= 2:
            machines_C -= 2
            parts_Y += 1
            changes = True
        # Method 4: Combine B and C into A
        elif machines_B >= 1 and machines_C >= 1:
            machines_B -= 1
            machines_C -= 1
            machines_A += 1
            changes = True
        # Method 5: Combine X and Y into Z
        elif parts_X >= 1 and parts_Y >= 1:
            parts_X -= 1
            parts_Y -= 1
            parts_Z += 1
            changes = True

# Apply the methods
apply_methods()

# Output the final counts
result = [str(machines_A), str(machines_B), str(machines_C), str(parts_X), str(parts_Y), str(parts_Z)]
print(result)