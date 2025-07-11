import collections

def solve_structure():
    """
    This script analyzes the provided spectral data to determine the IUPAC name of the compound.
    """
    print("Step 1: Mass Spectrum Analysis")
    mw = 135
    print(f"The mass spectrum shows a molecular ion peak (M+) at m/z = {mw}.")
    print("An odd molecular weight suggests the presence of an odd number of nitrogen atoms.")
    print("Assuming one nitrogen atom (N=14), the remaining mass is 135 - 14 = 121 for C and H.")
    # C*12 + H = 121. Possible formula: C9H13N (9*12 + 13 = 108 + 13 = 121)
    molecular_formula = "C9H13N"
    print(f"The most plausible molecular formula is {molecular_formula}.")
    
    # Calculate Degree of Unsaturation (DBE)
    C, H, N = 9, 13, 1
    dbe = C - H/2 + N/2 + 1
    print("\nCalculating the Degree of Unsaturation (DBE):")
    print(f"DBE = C + 1 - (H/2) + (N/2)")
    print(f"DBE = {C} + 1 - ({H}/2) + ({N}/2) = {int(dbe)}")
    print("A DBE of 4 is characteristic of a benzene ring.")

    print("\nStep 2: IR Spectrum Analysis")
    print("- Peaks at ~3300-3400 cm^-1 suggest an N-H stretch (likely a primary amine, -NH2).")
    print("- Peaks > 3000 cm^-1 suggest aromatic C-H bonds.")
    print("- Peaks < 3000 cm^-1 suggest aliphatic C-H bonds.")
    print("- Peaks at ~1600, 1450 cm^-1 suggest an aromatic ring.")
    print("Conclusion: The compound contains a primary amine, an aromatic ring, and aliphatic groups.")

    print("\nStep 3: 13C NMR and DEPT-135 Analysis")
    c13_shifts = [145.1, 128.5, 127.3, 126.3, 49.6, 43.5, 19.2]
    print(f"13C NMR signals are at: {c13_shifts}.")
    print("There are 7 distinct carbon environments.")
    print("DEPT-135 shows one negative signal (CH2) and five positive signals (CH or CH3).")
    print("This means one of the 7 carbon signals is a quaternary carbon (Cq), as it is absent in DEPT.")
    c_q = 145.1
    print(f"The signal at {c_q} ppm is the quaternary carbon (ipso-carbon of the benzene ring).")
    # Aromatic region: 128.5, 127.3, 126.3 (3 CH signals). Total = 1 Cq + 3 CH types = 4 signals for 6 benzene carbons (due to symmetry).
    # Aliphatic region: 49.6, 43.5, 19.2
    c_ch3 = 19.2 # Typically lowest shift
    c_ch2 = 43.5 # Negative signal
    c_ch = 49.6  # Positive signal
    print(f"Aliphatic signals are assigned as: CH3 at {c_ch3} ppm, CH2 at {c_ch2} ppm, and CH at {c_ch} ppm.")
    
    print("\nStep 4: 1H NMR Analysis")
    print("- Signal at ~7.2 ppm (multiplet, 5H) corresponds to a monosubstituted phenyl group (C6H5-).")
    print("- Signal at ~1.2 ppm (doublet, 3H) corresponds to a CH3 group next to a CH group.")
    print("- Signal at ~2.8 ppm (multiplet, 1H) corresponds to a CH group.")
    print("- Signal at ~2.9 ppm (multiplet, 2H) corresponds to a CH2 group.")
    print("- A broad signal for the NH2 protons (2H) is likely present but not clearly resolved.")
    
    print("\nStep 5: HSQC Analysis")
    print("The HSQC spectrum confirms the direct H-C connections:")
    print(f"- 1H (~1.2 ppm) correlates with 13C ({c_ch3} ppm) -> CH3 group.")
    print(f"- 1H (~2.8 ppm) correlates with 13C ({c_ch} ppm) -> CH group.")
    print(f"- 1H (~2.9 ppm) correlates with 13C ({c_ch2} ppm) -> CH2 group.")
    print("- 1H (~7.2 ppm) correlates with 13C (~126-129 ppm) -> Phenyl CH groups.")

    print("\nStep 6: Structure Elucidation")
    print("Combining all the fragments:")
    print("- A phenyl group (C6H5-).")
    print("- A propyl chain with an amine: -CH2-CH(NH2)-CH3.")
    print("1H NMR coupling confirms the connectivity: The phenyl is attached to the CH2 group.")
    print("The structure is: Phenyl-CH2-CH(NH2)-CH3.")
    
    print("\nStep 7: Final IUPAC Name")
    print("The parent chain is propane. The amine is on carbon 2, and the phenyl is on carbon 1.")
    iupac_name = "1-phenylpropan-2-amine"
    print(f"The IUPAC name of the compound is {iupac_name}.")

    return iupac_name

if __name__ == '__main__':
    final_answer = solve_structure()
    # The final answer format is specified by the user prompt.
    print(f"\n<<<{final_answer}>>>")
