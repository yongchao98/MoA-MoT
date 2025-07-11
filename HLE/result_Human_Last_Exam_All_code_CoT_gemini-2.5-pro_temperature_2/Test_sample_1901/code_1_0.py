def design_liquid_crystal_candidate(n):
    """
    Generates information about a single-ring liquid crystal candidate
    based on the general structure CnH2n+1-Ph-CN.

    Args:
        n (int): The number of carbon atoms in the flexible alkyl chain.
    """
    if not isinstance(n, int) or n < 1:
        print("Error: Chain length 'n' must be a positive integer.")
        return

    # Dictionary for alkyl chain names
    alkyl_names = {
        1: "methyl", 2: "ethyl", 3: "propyl", 4: "butyl",
        5: "pentyl", 6: "hexyl", 7: "heptyl", 8: "octyl",
        9: "nonyl", 10: "decyl"
    }
    alkyl_name = alkyl_names.get(n, f"C{n}")

    # --- Molecular Information ---
    print(f"--- Analyzing Candidate for n={n} ---")
    print("Design based on key requirements:")
    print("- Single benzene ring core")
    print("- Para (1,4) substitution")
    print("- General structure: C(n)H(2n+1)-[Benzene]-CN\n")

    # --- Molecular Formula Calculation ---
    # As requested, showing each number in the final equation.
    print("Calculating Molecular Formula:")
    
    # Part 1: Alkyl chain (CnH2n+1)
    h_alkyl = 2 * n + 1
    print(f"Alkyl Group (CnH2n+1): C({n})H({h_alkyl})")

    # Part 2: Benzene Core (para-substituted, so C6H4)
    c_benzene = 6
    h_benzene = 4
    print(f"Benzene Core (-C6H4-): C({c_benzene})H({h_benzene})")

    # Part 3: Cyano Group (-CN)
    c_cyano = 1
    n_cyano = 1
    print(f"Cyano Group (-CN): C({c_cyano})N({n_cyano})\n")

    # Final Equation and Result
    c_total = n + c_benzene + c_cyano
    h_total = h_alkyl + h_benzene
    print("Final Equation:")
    print(f"C({n} + {c_benzene} + {c_cyano}) H({h_alkyl} + {h_benzene}) N({n_cyano})")
    print(f"Result: C({c_total})H({h_total})N({n_cyano})\n")
    
    # --- Summary ---
    molecule_name = f"4-{alkyl_name}cyanobenzene"
    molecule_formula = f"C{c_total}H{h_total}N"

    print("--- Suggested Candidate Molecule ---")
    print(f"Name: {molecule_name}")
    print(f"Formula: {molecule_formula}\n")

    # --- Design Tuning Advice ---
    print("Tuning Advice (from Rule F):")
    print("- This molecule with n=5 is a great starting point.")
    print("- If experimental transition temperatures are too high, try increasing chain length 'n'.")
    print("- If too low, try decreasing chain length 'n'.")
    print("- Fine-tune by adding lateral substituents (e.g., -F) to the benzene ring.")

# Run the design script for the primary example n=5 (pentyl chain)
if __name__ == "__main__":
    design_liquid_crystal_candidate(5)