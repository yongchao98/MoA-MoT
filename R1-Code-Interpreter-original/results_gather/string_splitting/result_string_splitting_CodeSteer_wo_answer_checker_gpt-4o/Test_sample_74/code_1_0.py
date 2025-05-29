# Initialize the counts of machines and parts
machines_A = 4
machines_B = 4
machines_C = 3
parts_X = 0
parts_Y = 0
parts_Z = 0

# Function to apply the methods in sequence
def apply_methods():
    global machines_A, machines_B, machines_C, parts_X, parts_Y, parts_Z
    changes_made = False

    # Method 1: Each round, one machine A breaks down into two part X and one part Y.
    if machines_A >= 1:
        machines_A -= 1
        parts_X += 2
        parts_Y += 1
        changes_made = True

    # Method 2: Dismantling two machine B units produces one part X each cycle.
    if machines_B >= 2:
        machines_B -= 2
        parts_X += 1
        changes_made = True

    # Method 3: Each round, two machines C break down into one part Y.
    if machines_C >= 2:
        machines_C -= 2
        parts_Y += 1
        changes_made = True

    # Method 4: Combining one machine B and one machine C results in one machine A each cycle.
    if machines_B >= 1 and machines_C >= 1:
        machines_B -= 1
        machines_C -= 1
        machines_A += 1
        changes_made = True

    # Method 5: One unit of part X and one unit of part Y can be combined into one unit of part Z per round.
    if parts_X >= 1 and parts_Y >= 1:
        parts_X -= 1
        parts_Y -= 1
        parts_Z += 1
        changes_made = True

    return changes_made

# Continuously apply methods until no changes can be made
while apply_methods():
    pass

# Output the final counts in the specified format
result = [str(machines_A), str(machines_B), str(machines_C), str(parts_X), str(parts_Y), str(parts_Z)]
print(result)