import sys

def design_liquid_crystal():
    """
    Designs a single-ring liquid crystal molecule based on a set of rules.
    """
    # A. Key requirements: Single benzene ring, nematic/smectic, near room temp.
    # B. Core structure: Single benzene ring, para (1,4) substitution.
    core_structure = "Benzene Ring (Phenyl group, Ph)"

    # C. Recommended molecular features:
    flexible_tail = "Alkyl Chain (CnH2n+1)"
    polar_head = "Cyano Group (-CN)"

    # E. The general structure is defined as: CnH2n+1-Ph-CN

    # F. To achieve room temperature transitions, start with a pentyl chain.
    n = 5  # Start with a pentyl chain (n=5)

    # Calculate the specific alkyl chain formula for n=5
    carbon_atoms = n
    hydrogen_atoms = 2 * n + 1
    specific_alkyl_chain = f"C{carbon_atoms}H{hydrogen_atoms}"

    # D. Define the specific example structure based on the rules
    # This is a single-ring analogue of 5CB, known as 4-pentylbenzonitrile.
    final_molecule_name = "4-pentylbenzonitrile"
    final_molecular_formula = f"{specific_alkyl_chain}-Ph-{polar_head.split(' ')[1]}"

    # --- Output the Design ---
    print("--- Liquid Crystal Design Steps ---")
    print(f"1. Core: {core_structure}")
    print(f"2. Tail Group: {flexible_tail}")
    print(f"3. Head Group: {polar_head}")
    print("\n--- Tuning for Room Temperature ---")
    print("Strategy: Start with a specific alkyl chain length and adjust if necessary.")
    print(f"Initial choice for chain length 'n': {n}")
    
    print("\n--- Final Proposed Molecular Equation ---")
    # As requested, printing each number in the final equation/formula.
    print(f"The general formula is: C(n)H(2n+1)-Ph-CN")
    print(f"For n = {n}:")
    print(f"Number of Carbon atoms in alkyl chain: {carbon_atoms}")
    print(f"Number of Hydrogen atoms in alkyl chain (2*{n}+1): {hydrogen_atoms}")
    print("-----------------------------------------")
    print(f"Final Structure: C{carbon_atoms}H{hydrogen_atoms}-Ph-CN")
    print(f"Chemical Name: {final_molecule_name}")
    print("-----------------------------------------")

# Execute the design function
design_liquid_crystal()