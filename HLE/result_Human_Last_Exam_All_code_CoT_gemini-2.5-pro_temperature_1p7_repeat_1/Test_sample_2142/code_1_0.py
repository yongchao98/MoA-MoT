def solve_puzzle():
    """
    Solves the puzzle by mapping X1...X10 to digits 0...9 and
    constructing the numbers Y and Z based on the given formulas.
    """
    
    # Step 1: Establish the mapping for X1...X10 to digits 0...9.
    # Xn maps to the digit n-1. We use a list where the index corresponds to n-1.
    X = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

    # Step 2: Define the formulas for Y and Z using 1-based indices.
    y_indices = [1, 2, 3, 4, 2, 5]
    z_indices = [3, 6, 1, 7, 8, 9, 10, 2, 6]

    # Step 3: Construct Y and Z by substituting the digits.
    # We convert the 1-based indices to 0-based for list access (i-1).
    y_digits = [str(X[i - 1]) for i in y_indices]
    z_digits = [str(X[i - 1]) for i in z_indices]
    
    # The final prompt asks to "write Y Z". We will concatenate the digits for each
    # and print them separated by a space.
    y_final = "".join(y_digits)
    z_final = "".join(z_digits)

    # To satisfy the "output each number in the final equation" requirement,
    # we first show the composition of Y and Z from the individual numbers.
    print(f"Based on the mapping Xn = n-1:")
    print(f"Y = X1 X2 X3 X4 X2 X5 = {' '.join(y_digits)}")
    print(f"Z = X3 X6 X1 X7 X8 X9 X10 X2 X6 = {' '.join(z_digits)}")
    
    # Finally, print the solution "Y Z" as requested.
    print("\nFinal Answer:")
    print(f"{y_final} {z_final}")

solve_puzzle()