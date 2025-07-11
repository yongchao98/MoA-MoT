def analyze_covalency():
    """
    Analyzes and explains the relative covalency of CeF6(2-) and CeCl6(2-).
    This function will print the step-by-step reasoning based on the principles
    of chemical bonding.
    """
    
    # Principle relating covalency and orbital overlap
    print("Core Principle: Covalency, the sharing of electrons in a bond, is directly proportional to the extent of orbital overlap. Greater overlap results in a stronger covalent character.")
    print("-" * 80)
    
    # Information provided in the problem
    print("Given Information:")
    print("1. In CeF6(2-), the relevant orbitals are Cerium 4f and Fluorine 2p.")
    print("2. In CeCl6(2-), the relevant orbitals are Cerium 4f and Chlorine 3p.")
    print("3. The orbital overlap between Ce 4f and F 2p is greater than the overlap between Ce 4f and Cl 3p.")
    print("-" * 80)

    # Drawing the conclusion
    print("Conclusion:")
    print("Because the orbital overlap in CeF6(2-) is greater than in CeCl6(2-), the sharing of electrons between Cerium and the ligand is more effective in CeF6(2-).")
    
    # Final Answer
    final_answer = "Stronger"
    print(f"\nTherefore, CeF6(2-) displays {final_answer} covalency compared to CeCl6(2-).")

# Execute the analysis function
analyze_covalency()