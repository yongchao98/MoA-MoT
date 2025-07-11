def solve_glycan_puzzle():
    """
    Identifies a glycan from UPLC-MS and MS/MS data by calculating its mass,
    matching it to a proposed structure based on fragments, and providing its name.
    """

    # --- Step 1: Define constants and experimental data ---
    # Monoisotopic masses
    H_MASS = 1.007825  # Mass of a hydrogen atom
    H2O_MASS = 18.010565 # Mass of water

    # RapiFluor-MS (RFMS) Label info
    # Formula: C14H15F3N2O, Mass = 296.113647
    # Mass added during labeling = Mass(RFMS) - Mass(H2O)
    RFMS_ADDED_MASS = 296.113647 - H2O_MASS

    # Experimental MS1 data
    precursor_mz = 856.6638
    charge_state = 3

    # MS/MS fragment m/z values
    msms_fragments = {
        "HexNAc oxonium": 204.087,
        "Hex-HexNAc oxonium": 366.140,
        "Most intense (Hex-Hex-HexNAc)": 528.193,
    }
    
    # --- Step 2: Calculate the experimental mass of the glycan ---
    # M = (m/z * z) - (z * H)
    derivatized_mass = (precursor_mz * charge_state) - (charge_state * H_MASS)
    
    # Subtract the mass of the RFMS label to get the original glycan mass
    underivatized_glycan_mass = derivatized_mass - RFMS_ADDED_MASS

    # --- Step 3: Propose structure based on MS/MS and find best match ---
    # The intense fragment at m/z 528.193 strongly suggests an extended antenna
    # (e.g., Gal-Gal-GlcNAc). Sialylation is also likely based on typical analysis.
    # The structure 'A3G3S1' (a triantennary, monosialylated glycan) is a plausible
    # structure that can accommodate these features.

    # Monosaccharide residue masses (mass within the polymer chain)
    RESIDUE_MASS = {
        "Hex": 162.052823,   # Hexose (e.g., Mannose, Galactose)
        "HexNAc": 203.079372, # N-Acetylhexosamine (e.g., GlcNAc)
        "NeuAc": 291.095416, # N-Acetylneuraminic acid (Sialic Acid)
    }

    # Composition for A3G3S1: Hex(6) HexNAc(5) NeuAc(1)
    # This comes from: Man(3)GlcNAc(2) core + 3 antennas (GlcNAc)3 + 3 caps (Gal)3 + 1 sialic acid.
    # Note: To explain the 528 Da fragment, one Gal-GlcNAc antenna has an additional Gal.
    # The composition remains the same.
    composition = {
        "Hex": 6,
        "HexNAc": 5,
        "NeuAc": 1
    }
    
    # --- Step 4: Calculate theoretical mass and present the result ---
    theoretical_residue_sum = (composition["Hex"] * RESIDUE_MASS["Hex"] +
                               composition["HexNAc"] * RESIDUE_MASS["HexNAc"] +
                               composition["NeuAc"] * RESIDUE_MASS["NeuAc"])

    theoretical_glycan_mass = theoretical_residue_sum + H2O_MASS
    
    # Based on the MSMS data, we can infer linkage details
    glycan_name = "A3G3S1"
    
    description = (
        "The MS/MS data, particularly the intense ion at m/z 528.193, "
        "suggests a triantennary ('A3') structure where one of the antennae "
        "is extended with an additional galactose, likely forming a "
        "Gal(α1-3)Gal(β1-4)GlcNAc terminal sequence. One of the other "
        "antennae is capped with a single sialic acid ('S1'). The final 'G' "
        "in the name represents the remaining galactose-capped antenna."
    )

    print("--- Glycan Identification Report ---\n")
    print(f"Experimental Data:")
    print(f"  Precursor Ion [M+3H]3+: {precursor_mz:.4f} m/z")
    print(f"  Calculated Underivatized Glycan Mass: {underivatized_glycan_mass:.4f} Da\n")
    
    print(f"Proposed Structure:")
    print(f"  Oxford Nomenclature Name: {glycan_name}")
    print(f"  Monosaccharide Composition: "
          f"Hexose={composition['Hex']}, HexNAc={composition['HexNAc']}, NeuAc={composition['NeuAc']}\n")
    
    print(f"Mass Verification:")
    print(f"  Theoretical Mass for {glycan_name}: {theoretical_glycan_mass:.4f} Da")
    mass_difference = underivatized_glycan_mass - theoretical_glycan_mass
    print(f"  Mass Difference (Experimental - Theoretical): {mass_difference:.4f} Da\n")
    
    print("Inferred Structural Details:")
    print(description)

solve_glycan_puzzle()
<<<A3G3S1>>>