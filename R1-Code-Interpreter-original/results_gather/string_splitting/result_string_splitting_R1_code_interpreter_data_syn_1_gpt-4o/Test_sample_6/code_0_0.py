# Initial counts
machines_A = 2
machines_B = 3
machines_C = 3
parts_X = 0
parts_Y = 0
parts_Z = 0

# Function to apply the methods
def apply_methods():
    global machines_A, machines_B, machines_C, parts_X, parts_Y, parts_Z
    while True:
        # Method 1: Dismantle a single machine A
        if machines_A >= 1:
            machines_A -= 1
            parts_X += 2
            parts_Y += 1
        # Method 2: Dismantle two machine B units
        elif machines_B >= 2:
            machines_B -= 2
            parts_X += 1
        # Method 3: Dismantle two machine C units
        elif machines_C >= 2:
            machines_C -= 2
            parts_Y += 1
        # Method 4: Combine one unit of machine B and one unit of machine C into one unit of machine A
        elif machines_B >= 1 and machines_C >= 1:
            machines_B -= 1
            machines_C -= 1
            machines_A += 1
        # Method 5: Combine one unit of part X and one unit of part Y into one unit of part Z
        elif parts_X >= 1 and parts_Y >= 1:
            parts_X -= 1
            parts_Y -= 1
            parts_Z += 1
        else:
            # No more methods can be applied
            break

# Apply the methods
apply_methods()

# Output the final counts
result = [str(machines_A), str(machines_B), str(machines_C), str(parts_X), str(parts_Y), str(parts_Z)]
print(result)