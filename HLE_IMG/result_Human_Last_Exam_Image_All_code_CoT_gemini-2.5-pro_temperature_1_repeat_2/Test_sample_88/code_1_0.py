import re

def calculate_molar_mass(formula):
    """
    Calculates the molar mass of a chemical compound given its molecular formula.

    Args:
        formula (str): The molecular formula (e.g., 'C14H20N2O3').

    Returns:
        float: The molar mass in g/mol.
    """
    atomic_weights = {
        'C': 12.011,
        'H': 1.008,
        'N': 14.007,
        'O': 15.999
    }
    
    # Use regex to find all element-count pairs in the formula
    # e.g., 'C14H20N2O3' -> [('C', '14'), ('H', '20'), ('N', '2'), ('O', '3')]
    # For elements with no explicit count, the count is 1 (e.g., 'CH4O')
    components = re.findall('([A-Z][a-z]*)(\d*)', formula)
    
    molar_mass = 0.0
    for element, count in components:
        if element not in atomic_weights:
            raise ValueError(f"Unknown element in formula: {element}")
        
        # If count is an empty string, it means the count is 1
        num_atoms = int(count) if count else 1
        molar_mass += num_atoms * atomic_weights[element]
        
    return molar_mass

# --- Molecular Formulas ---
# Starting Material ((3,4-dihydro-2H-pyrrol-5-yl)proline), derived from structure
formula_sm = 'C9H14N2O2' 
# Methyl Propiolate
formula_mp = 'C4H4O2'
# Acetic Anhydride
formula_ac2o = 'C4H6O3'
# Acetic Acid (byproduct)
formula_acoh = 'C2H4O2'
# Methanol (byproduct)
formula_meoh = 'CH4O'

# Product Formulas
formula_A = 'C14H20N2O3'
formula_B = 'C12H14N2O3'
formula_C = 'C11H16N2O3'

# --- Calculations ---
mass_sm = calculate_molar_mass(formula_sm)
mass_mp = calculate_molar_mass(formula_mp)
mass_ac2o = calculate_molar_mass(formula_ac2o)
mass_acoh = calculate_molar_mass(formula_acoh)
mass_meoh = calculate_molar_mass(formula_meoh)

mass_A = calculate_molar_mass(formula_A)
mass_B = calculate_molar_mass(formula_B)
mass_C = calculate_molar_mass(formula_C)

print("--- Molar Mass Calculations ---")
print(f"Starting Material ({formula_sm}): {mass_sm:.3f} g/mol")
print(f"Product A ({formula_A}): {mass_A:.3f} g/mol")
print(f"Product B ({formula_B}): {mass_B:.3f} g/mol")
print(f"Product C ({formula_C}): {mass_C:.3f} g/mol")
print("-" * 30)

# --- Mass Balance Verification ---
# Equation for Product B: SM + MP -> B + MeOH
print("Verifying mass balance for the formation of Product B:")
print(f"Equation: {formula_sm} + {formula_mp} -> {formula_B} + {formula_meoh}")
reactants_mass_b = mass_sm + mass_mp
products_mass_b = mass_B + mass_meoh
print(f"Mass equation: {mass_sm:.3f} + {mass_mp:.3f} = {mass_B:.3f} + {mass_meoh:.3f}")
print(f"Result: {reactants_mass_b:.3f} = {products_mass_b:.3f}")
print("-" * 30)

# Equation for Product C: SM + Ac2O -> C + AcOH
print("Verifying mass balance for the formation of Product C:")
print(f"Equation: {formula_sm} + {formula_ac2o} -> {formula_C} + {formula_acoh}")
reactants_mass_c = mass_sm + mass_ac2o
products_mass_c = mass_C + mass_acoh
print(f"Mass equation: {mass_sm:.3f} + {mass_ac2o:.3f} = {mass_C:.3f} + {mass_acoh:.3f}")
print(f"Result: {reactants_mass_c:.3f} = {products_mass_c:.3f}")
print("-" * 30)
