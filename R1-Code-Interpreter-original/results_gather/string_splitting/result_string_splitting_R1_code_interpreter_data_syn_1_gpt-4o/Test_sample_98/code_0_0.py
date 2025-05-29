# Initial counts
machines_A = 4
machines_B = 4
machines_C = 3
parts_X = 0
parts_Y = 0
parts_Z = 0

# Function to apply methods
def apply_methods():
    global machines_A, machines_B, machines_C, parts_X, parts_Y, parts_Z
    changes = True
    while changes:
        changes = False
        # Method 1: Dismantle one machine A
        if machines_A >= 1:
            machines_A -= 1
            parts_X += 2
            parts_Y += 1
            changes = True
        # Method 2: Dismantle two machines B
        elif machines_B >= 2:
            machines_B -= 2
            parts_X += 1
            changes = True
        # Method 3: Dismantle two machines C
        elif machines_C >= 2:
            machines_C -= 2
            parts_Y += 1
            changes = True
        # Method 4: Combine one machine B and one machine C to form one machine A
        elif machines_B >= 1 and machines_C >= 1:
            machines_B -= 1
            machines_C -= 1
            machines_A += 1
            changes = True
        # Method 5: Combine one part X and one part Y to form one part Z
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