def calculate_glycan_mass():
    """
    Calculates the m/z of a specific derivatized N-glycan.

    The glycan is a biantennary (A2G2S2) structure that has been
    amidated and permethylated, and is detected as a singly-sodiated ion.
    """
    
    # Monoisotopic atomic masses (amu)
    atomic_masses = {
        'C': 12.000000,
        'H': 1.007825,
        'N': 14.003074,
        'O': 15.994915,
        'Na': 22.989770
    }

    # The chemical formula is derived as follows:
    # 1. Base glycan (3 Man, 2 Gal, 2 GlcNAc, 2 Neu5Ac) after polymerization (8 water molecules lost): C68 H112 N4 O52
    # 2. Amidation of 2 sialic acid carboxyls (-COOH -> -CONH2): adds 2*H, 2*N, removes 2*O. Formula becomes C68 H114 N6 O50
    # 3. Permethylation (replacing H with CH3) at 33 sites (25 -OH, 8 amide -NH/-NH2): adds 33*C, 33*2*H. Formula becomes C101 H180 N6 O50
    
    final_formula = {
        'C': 101,
        'H': 180,
        'N': 6,
        'O': 50
    }

    # Calculate the mass of the neutral molecule from its formula
    mass_C = final_formula['C'] * atomic_masses['C']
    mass_H = final_formula['H'] * atomic_masses['H']
    mass_N = final_formula['N'] * atomic_masses['N']
    mass_O = final_formula['O'] * atomic_masses['O']
    
    neutral_mass = mass_C + mass_H + mass_N + mass_O
    
    # Add the mass of a sodium ion for the [M+Na]+ adduct
    mass_with_sodium = neutral_mass + atomic_masses['Na']
    
    # The charge (z) is +1, so the m/z is equal to the mass.
    final_mz = mass_with_sodium

    # --- Output ---
    print("The three glycans are structural isomers and will have the same mass after derivatization.\n")
    print("Final Chemical Formula (neutral): C_101 H_180 N_6 O_50\n")
    
    print("Calculating the mass for the [M+Na]+ ion:")
    print(f"({final_formula['C']} * {atomic_masses['C']})  (from Carbon)")
    print(f"+ ({final_formula['H']} * {atomic_masses['H']})  (from Hydrogen)")
    print(f"+ ({final_formula['N']} * {atomic_masses['N']})  (from Nitrogen)")
    print(f"+ ({final_formula['O']} * {atomic_masses['O']})  (from Oxygen)")
    print(f"+ {atomic_masses['Na']}  (from Sodium adduct)")
    print("---------------------------------------------")
    print(f"= {final_mz:.4f}\n")

    print(f"The expected m/z for A2G(4)2S(3)2 is: {final_mz:.4f}")
    print(f"The expected m/z for A2G(4)S(3)S(6) is: {final_mz:.4f}")
    print(f"The expected m/z for A2G(4)2S(6)2 is: {final_mz:.4f}")

if __name__ == '__main__':
    calculate_glycan_mass()
<<<2300.1625>>>