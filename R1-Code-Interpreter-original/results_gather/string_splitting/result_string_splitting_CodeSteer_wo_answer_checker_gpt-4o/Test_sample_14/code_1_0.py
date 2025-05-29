# Initial quantities
machines_A = 3
machines_B = 1
machines_C = 2
parts_X = 0
parts_Y = 0
parts_Z = 0

# Function to apply each method
def apply_methods():
    global machines_A, machines_B, machines_C, parts_X, parts_Y, parts_Z
    changes = False

    # Method 1: Dismantle machine A
    if machines_A >= 1:
        machines_A -= 1
        parts_X += 2
        parts_Y += 1
        changes = True

    # Method 2: Dismantle machines B
    if machines_B >= 2:
        machines_B -= 2
        parts_X += 1
        changes = True

    # Method 3: Dismantle machines C
    if machines_C >= 2:
        machines_C -= 2
        parts_Y += 1
        changes = True

    # Method 4: Combine machine B and C into machine A
    if machines_B >= 1 and machines_C >= 1:
        machines_B -= 1
        machines_C -= 1
        machines_A += 1
        changes = True

    # Method 5: Combine parts X and Y into part Z
    if parts_X >= 1 and parts_Y >= 1:
        parts_X -= 1
        parts_Y -= 1
        parts_Z += 1
        changes = True

    return changes

# Apply methods until no changes occur
while apply_methods():
    pass

# Output the final quantities
result = [str(machines_A), str(machines_B), str(machines_C), str(parts_X), str(parts_Y), str(parts_Z)]
print(result)