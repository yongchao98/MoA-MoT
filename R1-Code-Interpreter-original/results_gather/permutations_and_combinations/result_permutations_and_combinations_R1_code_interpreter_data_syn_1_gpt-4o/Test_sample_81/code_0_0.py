# Define the positions of the books
positions = ["C", "D", "E", "A", "G", "F", "B"]

# Check the constraints
def check_constraints(positions):
    # Constraint 1: D must be to the left of G
    if positions.index("D") >= positions.index("G"):
        return False
    # Constraint 2: B must be to the right of A
    if positions.index("B") <= positions.index("A"):
        return False
    # Constraint 3: E must be adjacent to A
    if abs(positions.index("E") - positions.index("A")) != 1:
        return False
    # Constraint 4: F must be in position 6
    if positions[5] != "F":
        return False
    # Constraint 5: D must be adjacent to E
    if abs(positions.index("D") - positions.index("E")) != 1:
        return False
    # Constraint 6: F must be adjacent to G
    if abs(positions.index("F") - positions.index("G")) != 1:
        return False
    # Constraint 7: C must be to the right of A
    if positions.index("C") <= positions.index("A"):
        return False
    # Constraint 8: A must be adjacent to E
    if abs(positions.index("A") - positions.index("E")) != 1:
        return False
    return True

# Check if the arrangement satisfies all constraints
is_valid = check_constraints(positions)
print(is_valid)