def suggest_liquid_crystal_molecule(n, polar_group="CN"):
    """
    Generates a potential liquid crystal molecule name based on design rules.

    Args:
        n (int): The number of carbon atoms in the flexible alkyl chain.
        polar_group (str): The polar group attached to the benzene ring.

    Returns:
        str: The systematic name of the suggested molecule.
    """
    # A dictionary to map chain length to alkyl group name
    alkyl_names = {
        3: "propyl",
        4: "butyl",
        5: "pentyl",
        6: "hexyl",
        7: "heptyl",
        8: "octyl"
    }
    
    alkyl_name = alkyl_names.get(n, f"C{n}H{2*n+1}")
    
    # Using IUPAC naming for 1,4 (para) substitution on a benzene ring
    # Example: 1-cyano-4-pentylbenzene
    return f"1-{polar_group}-4-{alkyl_name}benzene"

# This script demonstrates the design and tuning process based on the provided rules.

# --- Step 1: Initial Design ---
# Based on rules (E) and (F), we start with the general structure CnH2n+1-Ph-CN
# and an initial chain length of n=5.
print("--- Step 1: Initial Design ---")
initial_n = 5
print(f"Starting with the recommended alkyl chain length (n): {initial_n}")
initial_molecule = suggest_liquid_crystal_molecule(initial_n)
print(f"Suggested Initial Molecule: {initial_molecule}\n")

# --- Step 2: Tuning Process ---
# Based on rule (F), we adjust the chain length based on hypothetical
# experimental transition temperatures.
print("--- Step 2: demonstrate Tuning Logic ---")

# Scenario A: The measured transition temperature is TOO HIGH.
# Action: Increase the chain length to lower the clearing point.
print("\nScenario A: Measured transition temperature is TOO HIGH.")
print("Action: Increase chain length (n -> n+1).")
n_too_high = initial_n + 1
molecule_too_high = suggest_liquid_crystal_molecule(n_too_high)
print(f"New suggested chain length (n): {n_too_high}")
print(f"New Suggested Molecule: {molecule_too_high}")

# Scenario B: The measured transition temperature is TOO LOW.
# Action: Decrease the chain length to raise the clearing point.
print("\nScenario B: Measured transition temperature is TOO LOW.")
print("Action: Decrease chain length (n -> n-1).")
n_too_low = initial_n - 1
molecule_too_low = suggest_liquid_crystal_molecule(n_too_low)
print(f"New suggested chain length (n): {n_too_low}")
print(f"New Suggested Molecule: {molecule_too_low}")
