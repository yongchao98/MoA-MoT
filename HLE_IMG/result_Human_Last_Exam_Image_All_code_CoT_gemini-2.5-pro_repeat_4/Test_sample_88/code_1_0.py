import re

def calculate_molecular_weight(formula):
    """
    Calculates the molecular weight of a chemical formula.
    Also returns a string detailing the calculation.
    """
    # Precise atomic weights
    atomic_weights = {
        'C': 12.011,
        'H': 1.008,
        'N': 14.007,
        'O': 15.999
    }

    # Parse the formula string, e.g., "C14H20N2O3"
    # The regex finds element symbols followed by optional numbers
    atoms = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
    
    total_weight = 0
    calculation_str_parts = []
    
    for element, count in atoms:
        # If count is not specified, it's 1
        atom_count = int(count) if count else 1
        
        if element in atomic_weights:
            weight = atomic_weights[element]
            total_weight += atom_count * weight
            calculation_str_parts.append(f"{atom_count} * {weight} ({element})")
        else:
            print(f"Warning: Atomic weight for element {element} not found.")
            return None, None

    calculation_str = " + ".join(calculation_str_parts)
    return total_weight, calculation_str

# Molecular formulas from the problem description
formulas = {
    'A': 'C14H20N2O3',
    'B': 'C12H14N2O3',
    'C': 'C11H16N2O3'
}

# Calculate and print the molecular weight for each product
print("--- Molecular Weight Calculations ---")
for product, formula in formulas.items():
    mw, calc_str = calculate_molecular_weight(formula)
    if mw is not None:
        print(f"Product {product} ({formula}):")
        print(f"{calc_str} = {mw:.3f} g/mol\n")

# --- Verification of chemical transformations using integer masses for clarity ---
print("\n--- Verifying Plausible Reaction Pathways using Integer Masses ---")

# Use integer atomic weights for simpler verification arithmetic
int_atomic_weights = {'C': 12, 'H': 1, 'N': 14, 'O': 16}

def get_integer_mw(formula):
    """Calculates molecular weight using integer masses."""
    atoms = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
    total_weight = 0
    for element, count in atoms:
        atom_count = int(count) if count else 1
        if element in int_atomic_weights:
            total_weight += atom_count * int_atomic_weights[element]
    return total_weight

# Define reactants and byproducts
sm_formula = 'C9H14N2O2'       # Starting Material
mp_formula = 'C4H4O2'         # Methyl Propiolate
acetyl_group = 'C2H2O'        # Net addition for acetylation/anhydride formation
meoh_formula = 'CH4O'          # Methanol
co2_formula = 'CO2'            # Carbon Dioxide

# Get integer MWs
mw_sm = get_integer_mw(sm_formula)
mw_mp = get_integer_mw(mp_formula)
mw_acetyl = get_integer_mw(acetyl_group)
mw_meoh = get_integer_mw(meoh_formula)
mw_co2 = get_integer_mw(co2_formula)
mw_a = get_integer_mw(formulas['A'])
mw_b = get_integer_mw(formulas['B'])
mw_c = get_integer_mw(formulas['C'])

print(f"Starting Material (SM) MW = {mw_sm}")
print(f"Methyl Propiolate (MP) MW = {mw_mp}\n")

# Verify formation of C
# C = SM + acetyl_group
calc_c = mw_sm + mw_acetyl
print(f"Formation of C: MW(SM) + MW(acetyl_group) = {mw_sm} + {mw_acetyl} = {calc_c}")
print(f"Matches MW(C) = {mw_c}\n")

# Verify formation of B
# B = SM + MP - MeOH
calc_b = mw_sm + mw_mp - mw_meoh
print(f"Formation of B: MW(SM) + MW(MP) - MW(MeOH) = {mw_sm} + {mw_mp} - {mw_meoh} = {calc_b}")
print(f"Matches MW(B) = {mw_b}\n")

# Verify formation of A
# A = C + MP - CO2
calc_a = mw_c + mw_mp - mw_co2
print(f"Formation of A: MW(C) + MW(MP) - MW(CO2) = {mw_c} + {mw_mp} - {mw_co2} = {calc_a}")
print(f"Matches MW(A) = {mw_a}\n")