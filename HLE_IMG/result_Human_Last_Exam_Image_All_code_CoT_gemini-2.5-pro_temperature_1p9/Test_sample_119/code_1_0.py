def solve_structure_puzzle():
    """
    This script explains the step-by-step reasoning for identifying the compound
    and verifies the molecular formula.
    """
    
    print("--- Structural Analysis Summary ---")

    print("\nStep 1: Molecular Formula Determination from Mass Spectrum and IR")
    mw = 135  # From (M+) m/z 135
    print(f"The molecular weight is {mw} g/mol. Odd MW suggests an odd number of Nitrogen atoms.")
    print("Based on all spectral data, the molecular formula is determined to be C9H13N.")

    # Calculate and verify the molecular weight of the proposed formula C9H13N
    mass_C = 12.011
    mass_H = 1.008
    mass_N = 14.007
    num_C = 9
    num_H = 13
    num_N = 1
    
    calculated_mw = (num_C * mass_C) + (num_H * mass_H) + (num_N * mass_N)
    
    print("\nStep 2: Verifying the Molecular Weight Calculation")
    print(f"The final proposed structure is 1-phenylpropan-2-amine with the formula C{num_C}H{num_H}N{num_N}.")
    print("Let's calculate its molecular weight:")
    print(f"({num_C} * {mass_C}) + ({num_H} * {mass_H}) + ({num_N} * {mass_N}) = {calculated_mw:.2f}")
    print(f"The calculated molecular weight ({calculated_mw:.2f}) matches the experimental value from the mass spectrum ({mw}).")

    print("\nStep 3: Structure Confirmation with NMR")
    print("1H, 13C, DEPT, and HSQC data are fully consistent with the structure of 1-phenylpropan-2-amine:")
    print("- C6H5- group (aromatic signals)")
    print("- -CH2- group (benzylic, δC=43.5, δH=~2.8)")
    print("- -CH(NH2)- group (δC=49.6, δH=~3.1)")
    print("- -CH3 group (δC=19.2, δH=1.2)")

    print("\n--- Final Answer ---")
    iupac_name = "1-phenylpropan-2-amine"
    print(f"The IUPAC name of the compound is: {iupac_name}")

if __name__ == '__main__':
    solve_structure_puzzle()
