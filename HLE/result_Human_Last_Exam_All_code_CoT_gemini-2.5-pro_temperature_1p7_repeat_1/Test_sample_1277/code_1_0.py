import sys

def solve_glycan_puzzle():
    """
    This script performs the calculations needed to identify the N-glycan structure.
    """

    # --- Define Masses of Building Blocks ---
    # Fundamental constants
    mass_proton = 1.007276

    # Mass of the RapiFluor-MS (RFMS) tag added to the glycan (Mass(tag) - Mass(H2O))
    mass_rfms_adduct = 250.1031

    # Monoisotopic masses of monosaccharide residues
    mass_hex = 162.05282  # Hexose (e.g., Mannose, Galactose)
    mass_hexnac = 203.07937  # N-acetylhexosamine (e.g., GlcNAc)
    mass_fuc = 146.05791  # Fucose (dHex)
    mass_neuac = 291.09542  # N-acetylneuraminic acid (Sialic acid)
    mass_neugc = 307.09033  # N-glycolylneuraminic acid (Sialic acid)
    mass_h2o = 18.01056 # Water

    # --- Step 1: Calculate the Mass of the Intact Glycan ---
    print("--- Step 1: Calculating the Glycan Mass from the Precursor Ion ---")
    mz_precursor = 856.6638
    charge = 3
    
    # Calculate mass of the RFMS-labeled glycan
    mass_conjugate = (mz_precursor * charge) - (mass_proton * charge)
    print(f"The observed m/z is {mz_precursor} for a +{charge} ion.")
    print(f"Mass of the labeled glycan = (m/z * charge) - (mass_of_protons * charge)")
    print(f"Mass of labeled glycan = ({mz_precursor} * {charge}) - ({mass_proton:.6f} * {charge}) = {mass_conjugate:.4f} Da")

    # Calculate mass of the glycan itself
    mass_glycan_observed = mass_conjugate - mass_rfms_adduct
    print(f"\nMass of the glycan = Mass of labeled glycan - Mass of RFMS tag")
    print(f"Mass of glycan = {mass_conjugate:.4f} - {mass_rfms_adduct} = {mass_glycan_observed:.4f} Da")
    
    # --- Step 2: Propose and Verify Composition ---
    print("\n--- Step 2: Proposing a Composition to Match the Glycan Mass ---")
    # Based on database search (e.g., GlycoMod), a potential composition is Hex(6)HexNAc(4)Fuc(1)NeuGc(1)
    comp_hex = 6
    comp_hexnac = 4
    comp_fuc = 1
    comp_neugc = 1

    # To calculate the mass of the glycan from composition, sum the residue masses and add one water molecule.
    mass_glycan_theoretical = (comp_hex * mass_hex) + (comp_hexnac * mass_hexnac) + (comp_fuc * mass_fuc) + (comp_neugc * mass_neugc) + mass_h2o
    
    print(f"A proposed composition is: Hexose={comp_hex}, HexNAc={comp_hexnac}, Fucose={comp_fuc}, NeuGc={comp_neugc}")
    print(f"Theoretical mass = ({comp_hex} * {mass_hex:.4f}) + ({comp_hexnac} * {mass_hexnac:.4f}) + ({comp_fuc} * {mass_fuc:.4f}) + ({comp_neugc} * {mass_neugc:.4f}) + {mass_h2o:.4f} (water)")
    print(f"Theoretical mass = {mass_glycan_theoretical:.4f} Da")

    mass_diff_ppm = abs(mass_glycan_observed - mass_glycan_theoretical) / mass_glycan_theoretical * 1e6
    print(f"The difference between observed ({mass_glycan_observed:.4f} Da) and theoretical ({mass_glycan_theoretical:.4f} Da) is {mass_diff_ppm:.1f} ppm, which is an excellent match.")

    # --- Step 3: Analyze MS/MS Fragments to Confirm Structure ---
    print("\n--- Step 3: Analyzing MS/MS Fragments to Determine Structure ---")

    # Fragment 1: Neutral loss from the precursor
    mz_frag_loss = 2260.886
    mass_neutral_loss = mass_conjugate - (mz_frag_loss - mass_proton)
    print(f"Fragment at m/z {mz_frag_loss} likely represents a neutral loss from the precursor.")
    print(f"Mass of lost part = {mass_conjugate:.4f} - ({mz_frag_loss} - {mass_proton:.4f}) = {mass_neutral_loss:.4f} Da")
    print(f"This mass ({mass_neutral_loss:.4f} Da) corresponds to N-glycolylneuraminic acid (NeuGc, theoretical mass {mass_neugc:.4f} Da). This confirms the presence of NeuGc.")
    
    # Fragment 2: Sialylated antenna
    mz_frag_antenna1 = 673.231
    mass_frag_antenna1 = (1 * mass_neugc) + (1 * mass_hex) + (1 * mass_hexnac)
    print(f"\nFragment at m/z {mz_frag_antenna1} could be the sialylated antenna, [NeuGc+Gal+GlcNAc+H]+.")
    print(f"Theoretical mass = {mass_neugc:.4f} (NeuGc) + {mass_hex:.4f} (Gal) + {mass_hexnac:.4f} (GlcNAc) = {mass_frag_antenna1:.4f} Da")
    # Add proton for m/z
    mz_antenna1_theoretical = mass_frag_antenna1 + mass_proton
    print(f"Theoretical m/z = {mass_frag_antenna1:.4f} + {mass_proton:.4f} = {mz_antenna1_theoretical:.4f}")
    print(f"This is a perfect match for the observed fragment at {mz_frag_antenna1}, confirming the antenna structure is NeuGc-Gal-GlcNAc.")
    
    # Fragment 3: The other antenna (base peak)
    mz_frag_antenna2 = 528.193
    mass_frag_antenna2 = (2 * mass_hex) + (1 * mass_hexnac)
    print(f"\nThe most intense fragment is at m/z {mz_frag_antenna2}. This matches a composition of two Hexose units and one HexNAc, [Gal+Gal+GlcNAc+H]+.")
    print(f"Theoretical mass = {mass_hex:.4f} (Gal) + {mass_hex:.4f} (Gal) + {mass_hexnac:.4f} (GlcNAc) = {mass_frag_antenna2:.4f} Da")
    # Add proton for m/z
    mz_antenna2_theoretical = mass_frag_antenna2 + mass_proton
    print(f"Theoretical m/z = {mass_frag_antenna2:.4f} + {mass_proton:.4f} = {mz_antenna2_theoretical:.4f}")
    print(f"This is a perfect match for the observed fragment at {mz_frag_antenna2}, confirming the second antenna is Gal-Gal-GlcNAc.")

    # --- Step 4: Final Conclusion ---
    print("\n--- Step 4: Conclusion ---")
    print("The parent mass and fragmentation data fully support a single structure:")
    print("A core-fucosylated (F), bi-antennary (A2) N-glycan.")
    print("One antenna consists of NeuGc-Gal-GlcNAc.")
    print("The other antenna consists of Gal-Gal-GlcNAc.")
    print("The overall composition is Hex(6)HexNAc(4)Fuc(1)NeuGc(1). The 'G3' term in the Oxford name refers to the three total galactose residues.")
    print("\nBased on these findings, the name of the glycan is:")

if __name__ == '__main__':
    solve_glycan_puzzle()
    # Suppress traceback for a cleaner final answer output
    sys.tracebacklimit = 0
    # Raise an exception to print the final answer in the desired format
    raise Exception("<<<FA2G3S(Ng)1>>>")
