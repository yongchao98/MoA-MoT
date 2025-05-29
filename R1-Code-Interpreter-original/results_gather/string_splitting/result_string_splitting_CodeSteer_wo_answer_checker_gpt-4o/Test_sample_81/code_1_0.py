# Initial quantities
machines_A = 6
machines_B = 4
machines_C = 5
parts_X = 0
parts_Y = 0
parts_Z = 0

# Function to apply methods in sequence
def apply_methods():
    global machines_A, machines_B, machines_C, parts_X, parts_Y, parts_Z
    while True:
        # Method 1: Dismantle machine A
        if machines_A >= 1:
            machines_A -= 1
            parts_X += 2
            parts_Y += 1
        # Method 2: Dismantle machine B
        elif machines_B >= 2:
            machines_B -= 2
            parts_X += 1
        # Method 3: Dismantle machine C
        elif machines_C >= 2:
            machines_C -= 2
            parts_Y += 1
        # Method 4: Combine machine B and C into machine A
        elif machines_B >= 1 and machines_C >= 1:
            machines_B -= 1
            machines_C -= 1
            machines_A += 1
        # Method 5: Combine parts X and Y into part Z
        elif parts_X >= 1 and parts_Y >= 1:
            parts_X -= 1
            parts_Y -= 1
            parts_Z += 1
        else:
            # If no method can be executed, break the loop
            break

# Apply the methods
apply_methods()

# Output the final quantities
result = [str(machines_A), str(machines_B), str(machines_C), str(parts_X), str(parts_Y), str(parts_Z)]
print(result)