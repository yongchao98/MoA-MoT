def solve_rearrangement():
    """
    This function determines the substituents at the specified positions in the product molecule.
    The product structure is identified as an olean-12-ene skeleton. The substituents
    at the numbered positions are determined by analyzing this standard structure.
    """
    substituents = {
        1: "CH3",
        2: "CH3",
        3: "H",
        4: "CH3",
        5: "CH3"
    }

    # Print the result in the specified format
    print(f"1 = {substituents[1]}, 2 = {substituents[2]}, 3 = {substituents[3]}, 4 = {substituents[4]}, 5 = {substituents[5]}")

solve_rearrangement()