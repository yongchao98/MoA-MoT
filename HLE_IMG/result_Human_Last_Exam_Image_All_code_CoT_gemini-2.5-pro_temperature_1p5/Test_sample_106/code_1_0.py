def solve_substituents():
    """
    Identifies the substituents at the specified positions in the product molecule
    and prints the result in the required format.
    """
    substituents = {
        1: "CH3",
        2: "CH3",
        3: "H",
        4: "CH3",
        5: "H"
    }

    # Format the output string as "1 = V, 2 = W, 3 = X, 4 = Y, 5 = Z"
    output_string = ", ".join([f"{position} = {group}" for position, group in sorted(substituents.items())])
    
    print(output_string)

solve_substituents()