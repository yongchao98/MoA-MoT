def generate_lc_name(n, lateral_substituent=None):
    """
    Generates a descriptive name for a liquid crystal molecule based on the design rules.

    Args:
        n (int): The number of carbon atoms in the flexible alkyl chain.
        lateral_substituent (str, optional): The name of a lateral substituent. Defaults to None.

    Returns:
        str: A string describing the molecule.
    """
    alkyl_prefixes = {
        1: "methyl", 2: "ethyl", 3: "propyl", 4: "butyl",
        5: "pentyl", 6: "hexyl", 7: "heptyl", 8: "octyl"
    }
    # Get the name of the alkyl chain, e.g., "pentyl"
    alkyl_chain_name = alkyl_prefixes.get(n, f"C{n}-alkyl")

    # The core of the name is based on the benzene ring and nitrile group
    base_name = f"4-{alkyl_chain_name}benzonitrile"

    if lateral_substituent:
        # If a substituent is present, add it to the name descriptively
        return f"{base_name} with a lateral {lateral_substituent} group"
    else:
        return base_name

# --- Main Design and Tuning Strategy ---

print("### Strategy for Designing a Single-Ring Room-Temperature Liquid Crystal ###")

# Step 1: Define the target molecular structure based on design rules.
print("\n1. General Molecular Structure (from Rule E):")
n = 5
h = 2 * n + 1
print("   The general structural formula is: C(n)H(2n+1)-Ph-CN")
print(f"   For our initial design (Rule F), we will use n = {n}.")
print(f"   Plugging this into the formula: C({n})H({2}*{n}+{1}) - Ph - CN")
print(f"   This gives the specific formula: C{n}H{h}-Ph-CN")

initial_molecule = generate_lc_name(n)
print(f"\n   Proposed starting molecule: {initial_molecule}")

# Step 2: Outline the iterative tuning process to achieve the desired properties.
print("\n2. Tuning Strategy to Achieve Room-Temperature Transition (from Rule F):")
print("   After synthesizing the initial molecule, measure its phase transition temperatures.")
print("   Then, adjust the structure based on the following logic:\n")

# Case A: Transition temperature is too high
print("   A) If Transition Temperature is TOO HIGH:")
print("      - ACTION: Increase the alkyl chain length 'n' to lower the melting point and nematic clearing point.")
n_too_high = n + 1
next_molecule_high = generate_lc_name(n_too_high)
print(f"      - EXAMPLE: Try n = {n_too_high}. The new target molecule would be: {next_molecule_high}")

# Case B: Transition temperature is too low
print("\n   B) If Transition Temperature is TOO LOW:")
print("      - ACTION: Decrease the alkyl chain length 'n'.")
n_too_low = n - 1
next_molecule_low = generate_lc_name(n_too_low)
print(f"      - EXAMPLE: Try n = {n_too_low}. The new target molecule would be: {next_molecule_low}")

# Case C: Fine-tuning is required
print("\n   C) For Fine-Tuning:")
print("      - ACTION: Add a lateral substituent (e.g., -F or -CH3) to the benzene ring.")
print("      - EFFECT: This typically lowers the clearing point and can alter phase behavior.")
fine_tuned_molecule = generate_lc_name(n, lateral_substituent="Fluoro")
print(f"      - EXAMPLE: Add a Fluoro group to the original molecule: {fine_tuned_molecule}")