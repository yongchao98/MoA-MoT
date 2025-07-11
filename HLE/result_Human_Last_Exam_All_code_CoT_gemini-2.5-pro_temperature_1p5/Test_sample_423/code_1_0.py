def solve_organic_puzzle():
    """
    Solves the organic chemistry puzzle by calculating molecular formulas
    and proposing a structure based on the provided clues.
    """
    
    # --- Step 1: Determine Molecular Formula of X ---
    mass_co2 = 0.7472  # g
    mass_h2o = 0.1834  # g
    
    # Molar masses
    M_co2 = 12.01 + 2 * 16.00
    M_h2o = 2 * 1.008 + 16.00
    M_c = 12.01
    M_h = 1.008
    M_o = 16.00
    
    # Moles of C and H
    moles_co2 = mass_co2 / M_co2
    moles_c_in_x = moles_co2
    moles_h2o = mass_h2o / M_h2o
    moles_h_in_x = 2 * moles_h2o
    
    # Ratio of C:H
    h_c_ratio = moles_h_in_x / moles_c_in_x  # Should be ~1.2
    
    # We found the empirical ratio is approx C5H6
    
    # Molar Mass of X is ~150 (range 135 to 165)
    # Trying multiples of empirical formula C5H6 + Oxygen
    n = 2  # to get into the molar mass range
    c_count_x = 5 * n
    h_count_x = 6 * n
    
    # Let's test possible molecular formulas
    # X = C10H12O (M = 148.21) -> Fits range and C:H ratio
    # X = C10H12O2 (M = 164.21) -> Fits range and C:H ratio
    # The problem describes reactions consistent with two oxygen atoms (e.g., phenol/enol + another group)
    # So we'll proceed with C10H12O2
    mf_x_c = 10
    mf_x_h = 12
    mf_x_o = 2
    molar_mass_x = mf_x_c * M_c + mf_x_h * M_h + mf_x_o * M_o
    
    print("Step 1: Analysis of substance X")
    print(f"Moles of C in sample = {moles_c_in_x:.4f} mol")
    print(f"Moles of H in sample = {moles_h_in_x:.4f} mol")
    print(f"H:C ratio = {h_c_ratio:.2f}:1, which simplifies to a 6:5 or 12:10 ratio.")
    print(f"Based on M ~ 150, the most plausible molecular formula for X is C{mf_x_c}H{mf_x_h}O{mf_x_o}")
    print(f"Calculated Molar Mass of C10H12O2 = {molar_mass_x:.2f} g/mol\n")
    
    # --- Step 2: Determine Molecular Formula and Structure of B ---
    # Mass fractions in B: C=0.5, H=0.1, O=0.4
    # Relative moles in 100g sample
    moles_c_in_b = 50 / M_c
    moles_h_in_b = 10 / M_h
    moles_o_in_b = 40 / M_o
    
    # Simplest ratio
    min_moles_b = min(moles_c_in_b, moles_h_in_b, moles_o_in_b)
    ef_b_c = round(moles_c_in_b / min_moles_b * 3) # Multiplied by 3 to get integers
    ef_b_h = round(moles_h_in_b / min_moles_b * 3)
    ef_b_o = round(moles_o_in_b / min_moles_b * 3)
    
    # B reduces to n-pentane, so it has 5 carbons. The EF is the MF.
    mf_b_c = 5
    mf_b_h = 12
    mf_b_o = 3
    molar_mass_b = mf_b_c * M_c + mf_b_h * M_h + mf_b_o * M_o

    print("Step 2: Analysis of substance B")
    print(f"The empirical formula of B based on mass fractions is C{ef_b_c}H{ef_b_h}O{ef_b_o}.")
    print("Since B is reduced to n-pentane, it must have 5 carbons.")
    print(f"Therefore, the molecular formula of B is C{mf_b_c}H{mf_b_h}O{mf_b_o}.")
    print(f"Calculated Molar Mass of C5H12O3 = {molar_mass_b:.2f} g/mol\n")

    print("Step 3: Deduction of Structures")
    print("Structure of B: The properties (reduction to n-pentane, 2 bromo-derivatives, optical activity) point to meso-pentane-2,3,4-triol.")
    print("Structure of A: A precursor to B formed via ozonolysis. Logically, A is 3-hydroxypentane-2,4-dione.")
    print("Structure of X: The problem contains contradictory information regarding the atom balance in the X->A->B reaction sequence.")
    print("However, based on all the chemical properties (C10H12O2, symmetric, enol giving red color with FeCl3, reduces Tollens, has a cleavable C=C bond), a plausible structure for X, despite formula inconsistencies in its reaction pathway, is (3Z,7Z)-4,7-dihydroxydeca-3,7-diene-2,9-dione. This structure is a symmetric enol-containing dimer.")
    print("\nFinal Proposed Structure for X:")
    print("Due to contradictions in the problem data, a definitive structure cannot be logically derived without assumptions. A highly plausible candidate based on the properties of X being a symmetric enol dimer is (3Z,7Z)-4,7-dihydroxydeca-3,7-diene-2,9-dione.")


solve_organic_puzzle()