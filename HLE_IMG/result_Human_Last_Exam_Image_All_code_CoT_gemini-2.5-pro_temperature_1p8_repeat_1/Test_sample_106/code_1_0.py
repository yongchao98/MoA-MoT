def solve_rearrangement():
    """
    This function determines and prints the substituents at the five specified
    positions in the product of the chemical rearrangement.
    """
    # The substituents are identified by visually inspecting the product molecule.
    # In chemical drawings, a line without a specific atomic label at its end
    # represents a methyl group (CH3). A bond shown without any label represents a hydrogen atom (H).
    substituents = {
        1: "CH3",
        2: "CH3",
        3: "H",
        4: "CH3",
        5: "H"
    }

    # Print the answer in the required format.
    for position, substituent in substituents.items():
        print(f"{position} = {substituent}")

solve_rearrangement()