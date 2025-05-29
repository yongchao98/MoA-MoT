# Define the positions of the books
positions = ["B", "C", "A", "E", "D", "F"]

# Check constraints
def check_constraints(positions):
    # Constraint 1: Book F cannot be in position 5
    if positions[4] == "F":
        return False
    # Constraint 2: Book A cannot be in position 2
    if positions[1] == "A":
        return False
    # Constraint 3: Book E must be adjacent to book D
    if abs(positions.index("E") - positions.index("D")) != 1:
        return False
    # Constraint 4: Book C must be to the left of book F
    if positions.index("C") > positions.index("F"):
        return False
    # Constraint 5: Book A must be in position 3
    if positions[2] != "A":
        return False
    # Constraint 6: Book C cannot be in position 5
    if positions[4] == "C":
        return False
    # Constraint 7: Book C must be adjacent to book A
    if abs(positions.index("C") - positions.index("A")) != 1:
        return False
    return True

# Check if the arrangement satisfies all constraints
if check_constraints(positions):
    print(positions)
else:
    print("No valid arrangement found.")