# Initialize the counts of machines and parts
machines_A = 2
machines_B = 3
machines_C = 2
parts_X = 0
parts_Y = 0
parts_Z = 0

# Function to apply the methods in sequence
def apply_methods():
    global machines_A, machines_B, machines_C, parts_X, parts_Y, parts_Z
    while True:
        # Method 1: Each round, one machine A breaks down into two part X and one part Y.
        if machines_A >= 1:
            machines_A -= 1
            parts_X += 2
            parts_Y += 1
        # Method 2: Dismantling two machine B units produces one part X each cycle.
        elif machines_B >= 2:
            machines_B -= 2
            parts_X += 1
        # Method 3: Each round, two machines C break down into one part Y.
        elif machines_C >= 2:
            machines_C -= 2
            parts_Y += 1
        # Method 4: Combining one machine B and one machine C results in one machine A each cycle.
        elif machines_B >= 1 and machines_C >= 1:
            machines_B -= 1
            machines_C -= 1
            machines_A += 1
        # Method 5: Combining one part X with one part Y produces one part Z each cycle.
        elif parts_X >= 1 and parts_Y >= 1:
            parts_X -= 1
            parts_Y -= 1
            parts_Z += 1
        else:
            # If no method can be executed, break the loop
            break

# Apply the methods
apply_methods()

# Output the final counts in the specified format
result = [str(machines_A), str(machines_B), str(machines_C), str(parts_X), str(parts_Y), str(parts_Z)]
print(result)