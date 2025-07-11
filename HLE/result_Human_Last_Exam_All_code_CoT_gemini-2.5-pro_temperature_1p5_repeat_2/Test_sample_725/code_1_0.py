def explain_covalency():
    """
    Explains the relative covalency of CeF6(2-) and CeCl6(2-)
    based on the provided information about orbital overlap.
    """

    # Step 1: State the fundamental chemical principle.
    print("Step 1: Understand the relationship between orbital overlap and covalency.")
    print("Principle: In chemical bonding, a greater degree of orbital overlap between atoms leads to more effective electron sharing. This results in a stronger covalent interaction and a higher degree of covalency.")
    print("-" * 50)

    # Step 2: State the information given in the problem.
    # The numbers in the chemical formulas and orbitals are included as requested.
    compound_1 = "CeF6(2-)"
    compound_2 = "CeCl6(2-)"
    ce_orbital = "4f"
    f_orbital = "2p"
    cl_orbital = "3p"

    print("Step 2: Analyze the given information.")
    print(f"It is stated that the orbital overlap in {compound_1} (between Cerium's {ce_orbital} orbital and Fluorine's {f_orbital} orbitals) is greater than the overlap in {compound_2} (between Cerium's {ce_orbital} orbital and Chlorine's {cl_orbital} orbitals).")
    print("-" * 50)

    # Step 3: Apply the principle to the information to draw a conclusion.
    conclusion = "stronger"
    print("Step 3: Draw a conclusion.")
    print(f"Applying the principle from Step 1 to the information in Step 2:")
    print(f"Because the orbital overlap is greater in {compound_1}, it displays {conclusion} covalency compared to {compound_2}.")
    print("\n--- Final Answer Derivation ---")
    print(f"Greater Overlap ({compound_1}) ==> {conclusion.capitalize()} Covalency")


explain_covalency()