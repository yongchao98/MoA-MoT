def analyze_covalency():
    """
    Analyzes the relative covalency of CeF6(2-) and CeCl6(2-) based on orbital overlap.
    """
    # Step 1: Define the core chemical principle.
    principle = "Covalency, which is the degree of electron sharing in a bond, is directly proportional to the extent of orbital overlap between the bonding atoms."
    print("Chemical Principle:")
    print(principle)
    print("-" * 30)

    # Step 2: State the information given in the problem.
    # We are comparing two chemical species: CeF6(2-) and CeCl6(2-).
    # The numbers in these formulas are 6 (for the ligands) and 2 (for the charge).
    compound_1 = "CeF6(2-)"
    compound_2 = "CeCl6(2-)"
    
    print("Information Provided:")
    print(f"In {compound_1}, the overlap is between Cerium 4f and Fluorine 2p orbitals.")
    print(f"In {compound_2}, the overlap is between Cerium 4f and Chlorine 3p orbitals.")
    print(f"It is given that the orbital overlap in {compound_1} is greater than in {compound_2}.")
    print("-" * 30)

    # Step 3: Apply the principle to the given information and draw a conclusion.
    print("Conclusion:")
    print(f"Because {compound_1} has greater orbital overlap than {compound_2}, it means there is more effective electron sharing between the central Cerium atom and the Fluorine ligands.")
    print(f"Therefore, {compound_1} displays stronger covalency compared to {compound_2}.")

# Run the analysis
analyze_covalency()