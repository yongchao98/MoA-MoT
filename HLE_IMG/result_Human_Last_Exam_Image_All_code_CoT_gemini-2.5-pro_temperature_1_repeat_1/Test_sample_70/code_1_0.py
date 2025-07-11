def analyze_electrocyclic_reaction(pi_electrons, reaction_type, condition="thermal"):
    """
    Analyzes an electrocyclic reaction using Woodward-Hoffmann rules.

    Args:
        pi_electrons (int): The number of pi electrons involved (e.g., 4, 6).
        reaction_type (str): The type of reaction ('opening' or 'closure').
        condition (str): The reaction condition ('thermal' or 'photochemical').

    Returns:
        str: A description of the reaction.
    """
    mode = ""
    # Woodward-Hoffmann rules for thermal reactions
    if condition == "thermal":
        if pi_electrons % 4 == 0:  # 4n system
            mode = "conrotatory"
        elif pi_electrons % 4 == 2:  # 4n+2 system
            mode = "disrotatory"
        else:
            return "Invalid number of pi electrons."
    else:
        # Rules for photochemical reactions (not needed here but included for completeness)
        if pi_electrons % 4 == 0:  # 4n system
            mode = "disrotatory"
        elif pi_electrons % 4 == 2:  # 4n+2 system
            mode = "conrotatory"
        else:
            return "Invalid number of pi electrons."

    return f"a {pi_electrons}π electron {mode} electrocyclic ring {reaction_type}"

def main():
    """
    Identifies the two pericyclic reactions in the given transformation.
    """
    print("The thermal transformation proceeds via two sequential pericyclic reactions:")

    # Step 1: Electrocyclic ring opening of the cyclobutene ring.
    # This involves the 2 electrons from the pi bond and 2 from the sigma bond.
    pi_electrons_step1 = 4
    reaction1_description = analyze_electrocyclic_reaction(pi_electrons_step1, "opening")
    print(f"1. The first step is {reaction1_description}.")

    # Step 2: Electrocyclic ring closure of the hexatriene part of the intermediate.
    # This involves 3 pi bonds from the intermediate.
    pi_electrons_step2 = 6
    reaction2_description = analyze_electrocyclic_reaction(pi_electrons_step2, "closure")
    print(f"2. The intermediate then undergoes {reaction2_description} to form the final product.")

if __name__ == "__main__":
    main()