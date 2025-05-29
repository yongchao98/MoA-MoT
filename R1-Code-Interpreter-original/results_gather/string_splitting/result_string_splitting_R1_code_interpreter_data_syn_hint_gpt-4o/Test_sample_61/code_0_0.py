# Initial counts
machines_A = 4
machines_B = 2
machines_C = 2
parts_X = 0
parts_Y = 0
parts_Z = 0

# Function to apply methods in order
def apply_methods():
    global machines_A, machines_B, machines_C, parts_X, parts_Y, parts_Z
    while True:
        # Method 1: One unit of machine A can be dismantled into two units of part X and one unit of part Y
        if machines_A >= 1:
            machines_A -= 1
            parts_X += 2
            parts_Y += 1
        # Method 2: Two units of machine B can be dismantled into one unit of part X
        elif machines_B >= 2:
            machines_B -= 2
            parts_X += 1
        # Method 3: Each round, two machines C break down into one part Y
        elif machines_C >= 2:
            machines_C -= 2
            parts_Y += 1
        # Method 4: Combining one machine B and one machine C results in one machine A
        elif machines_B >= 1 and machines_C >= 1:
            machines_B -= 1
            machines_C -= 1
            machines_A += 1
        # Method 5: Each round, one part X and one part Y merge to form one part Z
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