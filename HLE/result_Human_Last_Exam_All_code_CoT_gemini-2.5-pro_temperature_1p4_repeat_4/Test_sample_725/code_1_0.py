import sys

def analyze_covalency():
    """
    Analyzes the relative covalency of two cerium compounds based on orbital overlap.
    """
    
    # Principle of covalency
    print("Principle: Covalency is the measure of electron sharing in a chemical bond.")
    print("This sharing is physically represented by the overlap of atomic orbitals.")
    print("Therefore, greater orbital overlap leads to a stronger covalent character.\n")

    # Given information from the problem
    print("Given Information:")
    fact1 = "The 4f orbital overlap with ligand orbitals is greater in CeF6(2-) than in CeCl6(2-)."
    print(f"- {fact1}\n")
    
    # Logical Deduction
    print("Reasoning:")
    print("1. We start with the premise that: (Greater Overlap) ==> (Stronger Covalency).")
    print("2. We are given that the orbital overlap in CeF6(2-) is greater than in CeCl6(2-).")
    
    # Conclusion
    print("\nConclusion:")
    conclusion = "CeF6(2-) displays stronger covalency compared to CeCl6(2-)."
    print(f"Based on the reasoning, it can be concluded that {conclusion}")

if __name__ == '__main__':
    analyze_covalency()