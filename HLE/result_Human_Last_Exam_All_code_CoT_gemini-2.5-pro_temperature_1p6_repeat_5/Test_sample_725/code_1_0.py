def explain_covalency():
    """
    Explains which compound has stronger covalency based on orbital overlap.
    """
    
    # Principle of Covalency and Orbital Overlap
    print("Principle: The strength of covalency in a chemical bond is directly related to the degree of orbital overlap.")
    print("A greater overlap between atomic orbitals results in a stronger covalent character.\n")

    # Given information from the problem
    compound1 = "CeF6(2-)"
    overlap1_orbitals = "Ce(IV) 4f and F 2p"
    compound2 = "CeCl6(2-)"
    overlap2_orbitals = "Ce(IV) 4f and Cl 3p"
    
    print("Given Information:")
    print(f"1. In {compound1}, the overlap is between {overlap1_orbitals} orbitals.")
    print(f"2. In {compound2}, the overlap is between {overlap2_orbitals} orbitals.")
    print(f"3. It is observed that the orbital overlap in {compound1} is greater than in {compound2}.\n")

    # Conclusion based on the principle and given information
    print("Conclusion:")
    print(f"Because {compound1} has greater orbital overlap than {compound2}, it will exhibit stronger covalent bonding.")
    print(f"Therefore, {compound1} displays stronger covalency compared to {compound2}.")

# Execute the explanation
explain_covalency()