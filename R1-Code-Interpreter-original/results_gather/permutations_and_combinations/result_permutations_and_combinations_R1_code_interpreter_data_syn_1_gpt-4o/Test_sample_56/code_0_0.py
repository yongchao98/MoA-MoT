# Define the positions of the books
positions = ["C", "G", "A", "B", "D", "F", "E"]

# Check the constraints
def check_constraints(positions):
    # Constraint 1: Book D must be adjacent to book B
    if abs(positions.index("D") - positions.index("B")) != 1:
        return False
    # Constraint 2: Book F must be to the right of book A
    if positions.index("F") <= positions.index("A"):
        return False
    # Constraint 3: Book B cannot be placed in position 1
    if positions.index("B") == 0:
        return False
    # Constraint 4: Book D must be placed in position 5
    if positions.index("D") != 4:
        return False
    # Constraint 5: Book E must be to the right of book A
    if positions.index("E") <= positions.index("A"):
        return False
    # Constraint 6: Book C must be adjacent to book G
    if abs(positions.index("C") - positions.index("G")) != 1:
        return False
    # Constraint 7: Book G must be to the left of book A
    if positions.index("G") >= positions.index("A"):
        return False
    # Constraint 8: Book C must be placed in position 1
    if positions.index("C") != 0:
        return False
    return True

# Check if the arrangement satisfies all constraints
if check_constraints(positions):
    print(positions)
else:
    print("No valid arrangement found.")