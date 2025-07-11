def solve_covalency_comparison():
    """
    Analyzes the relative covalency of CeF6(2-) and CeCl6(2-) based on given orbital overlap data.
    """
    compound1 = "CeF6(2-)"
    compound2 = "CeCl6(2-)"
    
    # The numbers from the problem description are:
    # Ce(IV) state, 4f orbitals
    # Ligand 1: F, 2p orbitals, forming CeF6(2-)
    # Ligand 2: Cl, 3p orbitals, forming CeCl6(2-)

    # Principle: Covalency is a measure of electron sharing, which is directly related to orbital overlap.
    # Greater orbital overlap results in stronger covalency.
    
    # Given information from the problem:
    # CeF6(2-) has greater 4f-2p orbital overlap than CeCl6(2-) has with its 4f-3p overlap.
    
    conclusion = "stronger"
    
    print("In chemical bonding, the degree of covalency is directly proportional to the extent of orbital overlap.")
    print("Greater orbital overlap leads to more electron sharing and thus stronger covalency.")
    print(f"Given that the 4f orbital overlap with the ligand is greater in {compound1} (with F 2p orbitals) than in {compound2} (with Cl 3p orbitals),")
    print(f"it can be concluded that {compound1} displays {conclusion} covalency.")

solve_covalency_comparison()