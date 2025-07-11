def solve_puzzle():
    """
    Solves the puzzle by deciphering Xn as n and constructing Y and Z.
    """
    # Step 1: Decipher the characters X1 through X10 by assigning n to Xn.
    # We use a dictionary to map the character names to their string values.
    # Note that X10 is mapped to '10'.
    x_values = {f"X{i}": str(i) for i in range(1, 11)}

    # Step 2: Define the formulas for Y and Z as lists of character names.
    y_formula = ["X1", "X2", "X3", "X4", "X2", "X5"]
    z_formula = ["X3", "X6", "X1", "X7", "X8", "X9", "X10", "X2", "X6"]

    # Step 3: Get the numerical substitutions for each part of the equations.
    y_equation_parts = [x_values[token] for token in y_formula]
    z_equation_parts = [x_values[token] for token in z_formula]

    # Step 4: As requested, print each number in the final equations.
    print("The deciphered equations are:")
    print(f"Y = X1 X2 X3 X4 X2 X5 = {' '.join(y_equation_parts)}")
    print(f"Z = X3 X6 X1 X7 X8 X9 X10 X2 X6 = {' '.join(z_equation_parts)}")

    # Step 5: Construct the final Y and Z strings and print the result.
    y_result = "".join(y_equation_parts)
    z_result = "".join(z_equation_parts)
    
    print("\nThe final result Y Z is:")
    print(f"{y_result} {z_result}")

# Run the solver function
solve_puzzle()