def solve_structure():
    """
    This function performs the calculations to determine the empirical and molecular formulas
    for substances B and X, as described in the problem.
    """

    print("--- Part 1: Analysis of Substance B ---")
    # Given mass fractions for B
    mass_frac_C_B = 0.5
    mass_frac_H_B = 0.1
    mass_frac_O_B = 0.4

    # Molar masses of elements
    M_C = 12.01
    M_H = 1.008
    M_O = 16.00

    # Calculate moles assuming 100g of substance B
    moles_C_B = mass_frac_C_B * 100 / M_C
    moles_H_B = mass_frac_H_B * 100 / M_H
    moles_O_B = mass_frac_O_B * 100 / M_O

    print(f"Molar amounts in B (relative): C = {moles_C_B:.4f}, H = {moles_H_B:.4f}, O = {moles_O_B:.4f}")

    # Find the simplest ratio
    min_moles_B = min(moles_C_B, moles_H_B, moles_O_B)
    ratio_C_B = moles_C_B / min_moles_B
    ratio_H_B = moles_H_B / min_moles_B
    ratio_O_B = moles_O_B / min_moles_B

    print(f"Molar ratio in B: C = {ratio_C_B:.2f}, H = {ratio_H_B:.2f}, O = {ratio_O_B:.2f}")
    # To get integers, we multiply by 3
    final_C_B = round(ratio_C_B * 3)
    final_H_B = round(ratio_H_B * 3)
    final_O_B = round(ratio_O_B * 3)
    print(f"Integer ratio -> Empirical formula of B: C{final_C_B}H{final_H_B}O{final_O_B}")
    print("\n")

    print("--- Part 2: Analysis of Substance X ---")
    # Given combustion data for X
    mass_CO2 = 0.7472  # g
    mass_H2O = 0.1834  # g

    # Molar masses of compounds
    M_CO2 = 44.01
    M_H2O = 18.016

    # Calculate moles of C and H in the sample of X
    moles_CO2 = mass_CO2 / M_CO2
    moles_C_X = moles_CO2
    moles_H2O = mass_H2O / M_H2O
    moles_H_X = 2 * moles_H2O

    print(f"Molar amounts in X sample: C = {moles_C_X:.5f}, H = {moles_H_X:.5f}")

    # Calculate the H/C ratio
    H_C_ratio = moles_H_X / moles_C_X
    print(f"Molar ratio H/C in X: {H_C_ratio:.3f}")
    # 1.2 is the ratio 12/10
    print("This ratio suggests a molecular formula with a H:C ratio of 12:10.")

    # Check potential formulas against the molar mass (M ~ 150, range 135-165)
    # Candidate formula based on ratio: C10H12Oz
    C = 10
    H = 12
    for O_count in range(1, 4):
        molar_mass_X = C * M_C + H * M_H + O_count * M_O
        print(f"Testing formula C{C}H{H}O{O_count}: Molar Mass = {molar_mass_X:.2f} g/mol")
        if 135 <= molar_mass_X <= 165:
            print(f"  -> C{C}H{H}O{O_count} is a strong candidate for the molecular formula of X.")

    print("\n--- Final Conclusion ---")
    print("Based on the calculations and chemical reasoning:")
    print("Formula of B is C5H12O3.")
    print("Formula of X is C10H12O2.")
    print("Structure of X is Cyclodeca-1,6-diene-3,8-dione.")
    print("The final ozonolysis equation is:")
    print("C10H12O2 (X) + 2 O3 --(Zn/CH3COOH)--> 2 C5H6O3 (A)")
    print("where each number in the equation is:")
    print("1 molecule of C10H12O2 reacts with")
    print("2 molecules of O3 to produce")
    print("2 molecules of C5H6O3")

solve_structure()