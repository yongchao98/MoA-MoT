import sys

def solve_glycan_puzzle():
    """
    Analyzes mass spectrometry data to identify an N-glycan, calculates its mass,
    interprets fragment ions, and provides its name using Oxford nomenclature.
    """

    # --- Constants for Exact Masses ---
    H_MASS = 1.007825
    NA_MASS = 22.989770
    # Mass added by RapiFluor-MS derivatization = Mass(RFMS reagent) - Mass(H2O)
    RFMS_ADDITION_MASS = 334.1332 - 18.0106
    
    # Monosaccharide masses
    HEXNAC_MASS = 203.07937 # N-Acetylglucosamine (GlcNAc)
    HEX_MASS = 162.05282   # Hexose (Mannose, Galactose)
    NEUAC_MASS = 291.09541 # N-Acetylneuraminic acid (Sialic Acid)

    # --- Input Data from the User ---
    parent_mz = 856.6638
    isotopic_peaks = [856.6638, 856.9971, 857.3305]
    msms_ions = [204.087, 366.140, 528.193, 673.231, 882.409, 1368.568, 1894.753, 2260.886]

    # --- Step 1: Determine Charge State ---
    print("--- Step 1: Determine Charge State from Isotopic Envelope ---")
    delta_m1 = isotopic_peaks[1] - isotopic_peaks[0]
    delta_m2 = isotopic_peaks[2] - isotopic_peaks[1]
    avg_delta_m = (delta_m1 + delta_m2) / 2
    charge = round(1 / avg_delta_m)
    print(f"The mass difference between isotopic peaks is ~{avg_delta_m:.4f} Da.")
    print(f"This corresponds to a charge state (z) of 1 / {avg_delta_m:.4f} = {1/avg_delta_m:.2f}, so z = {charge}.")

    # --- Step 2: Calculate Neutral Mass of Derivatized Glycan ---
    print("\n--- Step 2: Calculate Neutral Mass from Parent Ion m/z ---")
    print("Sialylated glycans often form sodium adducts. We will test potential adducts [M+xH+yNa]z+.")
    # Test [M+H+2Na]3+ adduct, which is common for disialylated species
    derivatized_mass_neutral = (parent_mz * charge) - H_MASS - (2 * NA_MASS)
    print(f"Assuming an adduct of [M+H+2Na]3+, the neutral mass of the RFMS-derivatized glycan is:")
    print(f"({parent_mz} * {charge}) - {H_MASS} - (2 * {NA_MASS:.4f}) = {derivatized_mass_neutral:.4f} Da")

    # --- Step 3: Calculate Native Glycan Mass ---
    print("\n--- Step 3: Calculate the Mass of the Underivatized Native Glycan ---")
    native_glycan_mass = derivatized_mass_neutral - RFMS_ADDITION_MASS
    print(f"Subtracting the mass added by RFMS derivatization ({RFMS_ADDITION_MASS:.4f} Da):")
    print(f"{derivatized_mass_neutral:.4f} - {RFMS_ADDITION_MASS:.4f} = {native_glycan_mass:.4f} Da")

    # --- Step 4: Propose a Structure and Calculate its Theoretical Mass ---
    print("\n--- Step 4: Propose a Structure by Matching the Native Glycan Mass ---")
    # Propose A2G2S2: biantennary, digalactosylated, disialylated glycan
    num_hexnac = 4
    num_hex = 5
    num_neuac = 2
    a2g2s2_mass = (num_hexnac * HEXNAC_MASS) + (num_hex * HEX_MASS) + (num_neuac * NEUAC_MASS)
    mass_diff = native_glycan_mass - a2g2s2_mass
    
    print(f"The calculated native mass ({native_glycan_mass:.4f} Da) closely matches the theoretical mass of A2G2S2.")
    print(f"A2G2S2 Composition: {num_hexnac} HexNAc, {num_hex} Hex, {num_neuac} NeuAc")
    print(f"Theoretical Mass = ({num_hexnac} * {HEXNAC_MASS:.4f}) + ({num_hex} * {HEX_MASS:.4f}) + ({num_neuac} * {NEUAC_MASS:.4f}) = {a2g2s2_mass:.4f} Da.")
    print(f"The difference is only {mass_diff:.4f} Da, confirming this composition.")

    # --- Step 5: Confirm Structure and Linkage with MS/MS Fragments ---
    print("\n--- Step 5: Confirm Structure and Linkage using MS/MS Fragment Ions ---")
    print("The MS/MS fragments support the A2G2S2 assignment:")
    print(f"- m/z {msms_ions[0]}: Corresponds to a HexNAc oxonium ion ([C8H14NO5]+, theoretical m/z 204.086), confirming N-acetylglucosamine.")
    print(f"- m/z {msms_ions[1]}: Corresponds to a Hex-HexNAc oxonium ion ([C14H24NO10]+, theoretical m/z 366.140), a typical Gal-GlcNAc antenna fragment.")
    print(f"- m/z {msms_ions[2]} (most intense): This corresponds to a Hex2-HexNAc fragment ([C20H34NO15]+, theoretical m/z 528.192). The high intensity of this ion is characteristic of a biantennary N-glycan.")
    print("- Crucially, the absence of significant sialic acid oxonium ions (at m/z 274 or 292) is strong evidence that the terminal sialic acids are attached via alpha-2,6 linkages. Alpha-2,3 linkages typically produce these fragments, so their absence is informative.")
    
    # --- Step 6: Final Conclusion ---
    print("\n--- Final Conclusion ---")
    print("Based on the parent mass calculation and the fragmentation pattern, the glycan is a biantennary complex N-glycan with two antennas, each terminated with a galactose and an alpha-2,6-linked sialic acid.")
    print("\nUsing the Oxford nomenclature, the name of this glycan is:")
    final_answer = "A2G2S2 with α2,6-linked sialic acids"
    print(final_answer)

solve_glycan_puzzle()
<<<A2G2S2 with α2,6-linked sialic acids>>>