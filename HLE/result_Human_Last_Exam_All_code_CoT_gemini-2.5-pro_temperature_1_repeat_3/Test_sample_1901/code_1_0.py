def generate_lc_formula(n, core="Ph", polar_group="CN"):
    """
    Generates the chemical formula for a single-ring liquid crystal molecule.

    Args:
        n (int): The number of carbon atoms in the alkyl chain.
        core (str): The representation for the rigid core.
        polar_group (str): The polar terminal group.

    Returns:
        str: The formatted chemical formula.
    """
    if n <= 0:
        return "Invalid chain length"
    
    # Formula for an alkyl group: CnH2n+1
    alkyl_chain = f"C{n}H{2*n + 1}"
    
    return f"{alkyl_chain}-{core}-{polar_group}"

def design_liquid_crystal():
    """
    Guides the user through the design process of a single-ring liquid crystal
    and prints the results.
    """
    print("--- 1. Liquid Crystal Design Requirements ---")
    print("- Core: Single benzene ring")
    print("- Phase: Nematic or Smectic")
    print("- Transition Temperature: Near room temperature (~20-25°C)")
    print("\n" + "="*50 + "\n")

    print("--- 2. Proposed Molecular Design Strategy ---")
    print("The general structure combines a flexible chain, a rigid core, and a polar group.")
    print("The design 'equation' or general formula is: C(n)H(2*n+1)-Ph-CN")
    print("Where:")
    print("  - C(n)H(2*n+1) is a flexible alkyl chain.")
    print("  - Ph represents the rigid benzene ring core.")
    print("  - CN is a polar cyano group, which promotes liquid crystallinity.")
    print("\n" + "="*50 + "\n")

    print("--- 3. Tuning Strategy for Room Temperature Transitions ---")
    print("The primary method for tuning the transition temperature is by adjusting the alkyl chain length (n).")
    print("- A good starting point for room temperature liquid crystals is n=5.")
    print("- To decrease transition temperatures, generally increase chain length (n).")
    print("- To increase transition temperatures, generally decrease chain length (n).")
    print("\n" + "="*50 + "\n")
    
    print("--- 4. Example Molecular Formulas ---")
    
    # Starting point example
    n_start = 5
    print(f"A) Starting Point Example (n = {n_start}):")
    formula_start = generate_lc_formula(n_start)
    print(f"   The general formula is: C({n_start})H(2*{n_start}+1)-Ph-CN")
    print(f"   The final structure is: {formula_start}")
    print("   (This corresponds to the molecule 4-pentylcyanobenzene)\n")
    
    # Longer chain example
    n_longer = 8
    print(f"B) Longer Chain Example (n = {n_longer}):")
    formula_longer = generate_lc_formula(n_longer)
    print(f"   The general formula is: C({n_longer})H(2*{n_longer}+1)-Ph-CN")
    print(f"   The final structure is: {formula_longer}")
    print("   (This corresponds to the molecule 4-octylcyanobenzene, known as 8CB)\n")

    # Shorter chain example
    n_shorter = 3
    print(f"C) Shorter Chain Example (n = {n_shorter}):")
    formula_shorter = generate_lc_formula(n_shorter)
    print(f"   The general formula is: C({n_shorter})H(2*{n_shorter}+1)-Ph-CN")
    print(f"   The final structure is: {formula_shorter}")
    print("   (This corresponds to the molecule 4-propylcyanobenzene)")


if __name__ == '__main__':
    design_liquid_crystal()