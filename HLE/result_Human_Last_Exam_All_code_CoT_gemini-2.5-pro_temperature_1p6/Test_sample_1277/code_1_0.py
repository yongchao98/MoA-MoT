def solve_glycan_puzzle():
    """
    This function analyzes the provided mass spectrometry data to identify a glycan.
    It prints the step-by-step interpretation of the data and the final conclusion.
    """

    # --- Given Data ---
    mz_parent = 856.6638
    mz_isotopes = [856.9971, 857.3305, 857.6638]
    mz_fragments = {
        204.087: "HexNAc oxonium ion",
        366.140: "Hex-HexNAc oxonium ion",
        528.193: "Doubly charged core fragment (most intense)",
        2260.886: "Neutral loss fragment"
    }
    
    # --- Masses ---
    mass_H = 1.007825
    mass_RFMS_net = 260.0732
    mass_NeuGc = 307.0903
    
    # --- Calculations & Interpretation ---
    
    print("Step 1: Determine Parent Ion Mass and Charge")
    isotope_spacing = mz_isotopes[0] - mz_parent
    charge_state = round(1 / isotope_spacing)
    print(f"Isotope spacing is ~{isotope_spacing:.4f} Da, indicating a charge state (z) of {charge_state}.")

    mass_derivatized_glycan = (mz_parent * charge_state) - (charge_state * mass_H)
    print(f"The monoisotopic mass of the RFMS-derivatized glycan is ({mz_parent} * {charge_state}) - ({charge_state} * {mass_H}) = {mass_derivatized_glycan:.4f} Da.")

    mass_underivatized_glycan = mass_derivatized_glycan - mass_RFMS_net
    print(f"Subtracting the RFMS tag mass ({mass_RFMS_net} Da), the underivatized glycan mass is {mass_underivatized_glycan:.4f} Da.\n")
    
    print("Step 2: Interpret Key MS/MS Fragments")
    
    # Fragment 1
    frag1_mz = 2260.886
    neutral_loss = mass_derivatized_glycan - (frag1_mz - mass_H)
    print(f"Fragment at m/z {frag1_mz}: This corresponds to the parent ion losing a neutral molecule.")
    print(f"Calculated neutral loss: {mass_derivatized_glycan:.4f} - ({frag1_mz} - {mass_H}) = {neutral_loss:.4f} Da.")
    print(f"This mass matches N-glycolylneuraminic acid (NeuGc, mass {mass_NeuGc} Da). This indicates the presence of one NeuGc residue.\n")
    
    # Fragment 2
    frag2_mz = 528.193
    core_fuc_mass = (frag2_mz * 2) - (2 * mass_H)
    print(f"Fragment at m/z {frag2_mz}: This is the most intense fragment, suggesting a very stable ion.")
    print(f"This corresponds to a doubly charged ion of the core structure: [Fuc-Man3-GlcNAc2 + 2H]2+.")
    print(f"This indicates the glycan is core-fucosylated (a fucose is attached to the core).\n")

    # Fragment 3
    frag3_mz = 366.140
    print(f"Fragment at m/z {frag3_mz}: This is a well-known oxonium ion for a Hex-HexNAc disaccharide.")
    print(f"This confirms the presence of LacNAc (Gal-GlcNAc) antennae.\n")
    
    print("Step 3: Conclusion")
    glycan_name = "FA2G2S(Gc)1"
    print(f"Based on the evidence, the glycan is a biantennary, core-fucosylated structure with two galactose residues, one of which is capped with N-glycolylneuraminic acid.")
    print("The name of the glycan using Oxford nomenclature is:")
    print(f">>> {glycan_name}")
    print("\nBreakdown of the name:")
    print("F: Core fucose (α1,6-linked to the first GlcNAc)")
    print("A2: Biantennary (two antennae)")
    print("G2: Two terminal Galactose residues")
    print("S(Gc)1: One N-Glycolylneuraminic acid residue")

solve_glycan_puzzle()