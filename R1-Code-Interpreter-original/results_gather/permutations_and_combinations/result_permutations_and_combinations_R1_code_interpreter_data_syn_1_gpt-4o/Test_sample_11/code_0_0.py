# Define the positions of the books
books = ["F", "D", "C", "A", "E", "B"]

# Check the constraints
def check_constraints(books):
    # Constraint 1: Book C must be in position 3
    if books[2] != "C":
        return False
    # Constraint 2: Book E must be to the right of book F
    if books.index("E") <= books.index("F"):
        return False
    # Constraint 3: Book A cannot be in position 3
    if books[2] == "A":
        return False
    # Constraint 4: Book C must be adjacent to book D
    if abs(books.index("C") - books.index("D")) != 1:
        return False
    # Constraint 5: Book C cannot be in position 1
    if books[0] == "C":
        return False
    # Constraint 6: Book E must be adjacent to book A
    if abs(books.index("E") - books.index("A")) != 1:
        return False
    # Constraint 7: Book F must be to the left of book C
    if books.index("F") >= books.index("C"):
        return False
    return True

# Check if the arrangement satisfies all constraints
if check_constraints(books):
    print(books)
else:
    print("No valid arrangement found.")