# Define the positions of the books
positions = ["D", "B", "A", "E", "G", "F", "C"]

# Check constraints
def check_constraints(positions):
    # Constraint 1: G must be adjacent to F
    if abs(positions.index("G") - positions.index("F")) != 1:
        return False
    # Constraint 2: B must be to the right of D
    if positions.index("B") <= positions.index("D"):
        return False
    # Constraint 3: G must be in position 5
    if positions.index("G") != 4:
        return False
    # Constraint 4: B must be adjacent to D
    if abs(positions.index("B") - positions.index("D")) != 1:
        return False
    # Constraint 5: A cannot be in position 7
    if positions.index("A") == 6:
        return False
    # Constraint 6: E must be to the right of A
    if positions.index("E") <= positions.index("A"):
        return False
    # Constraint 7: E must be adjacent to A
    if abs(positions.index("E") - positions.index("A")) != 1:
        return False
    # Constraint 8: D must be in position 1
    if positions.index("D") != 0:
        return False
    # Constraint 9: B cannot be in position 7
    if positions.index("B") == 6:
        return False
    return True

# Check if the current arrangement satisfies all constraints
if check_constraints(positions):
    print(positions)
else:
    print("No valid arrangement found.")