import re

def calculate_mw(formula):
    """
    Calculates the molecular weight of a chemical compound given its formula.

    Args:
        formula (str): The chemical formula (e.g., 'C9H14N2O2').

    Returns:
        float: The molecular weight.
    """
    atomic_weights = {
        'C': 12.011,
        'H': 1.008,
        'N': 14.007,
        'O': 15.999,
    }
    
    # Regex to find elements and their counts (e.g., C9, H14, N2, O2)
    pattern = r'([A-Z][a-z]?)(\d*)'
    atoms = re.findall(pattern, formula)
    
    total_mw = 0.0
    for element, count in atoms:
        count = int(count) if count else 1
        if element in atomic_weights:
            total_mw += atomic_weights[element] * count
        else:
            raise ValueError(f"Atomic weight for element {element} not found.")
            
    return total_mw

# --- Define Molecular Formulas ---
formula_sm = 'C9H14N2O2'   # Starting Material
formula_mp = 'C4H4O2'     # Methyl Propiolate
formula_ac2o = 'C4H6O3'   # Acetic Anhydride
formula_acoh = 'C2H4O2'   # Acetic Acid
formula_meoh = 'CH4O'     # Methanol
formula_co2 = 'CO2'      # Carbon Dioxide

formula_A = 'C14H20N2O3'
formula_B = 'C12H14N2O3'
formula_C = 'C11H16N2O3'

# --- Calculate Molecular Weights ---
mw_sm = calculate_mw(formula_sm)
mw_mp = calculate_mw(formula_mp)
mw_ac2o = calculate_mw(formula_ac2o)
mw_acoh = calculate_mw(formula_acoh)
mw_meoh = calculate_mw(formula_meoh)
mw_co2 = calculate_mw(formula_co2)
mw_A = calculate_mw(formula_A)
mw_B = calculate_mw(formula_B)
mw_C = calculate_mw(formula_C)

# --- Print and Verify Proposed Reactions ---

print("Verifying the proposed reaction pathways by mass balance:")
print("-" * 55)

# 1. Formation of Product C
print("Pathway to Product C (N-Acetylation):")
reactants_C_mw = mw_sm + mw_ac2o
products_C_mw = mw_C + mw_acoh
print(f"  {formula_sm}   +   {formula_ac2o}   ->   {formula_C}   +   {formula_acoh}")
print(f" {mw_sm:7.3f} g/mol + {mw_ac2o:7.3f} g/mol -> {mw_C:7.3f} g/mol + {mw_acoh:7.3f} g/mol")
print(f"Total Reactant Mass: {reactants_C_mw:.3f} g/mol")
print(f"Total Product Mass:  {products_C_mw:.3f} g/mol")
print("-" * 55)


# 2. Formation of Product B
print("Pathway to Product B (Annulation with Methanol Elimination):")
reactants_B_mw = mw_sm + mw_mp
products_B_mw = mw_B + mw_meoh
print(f"  {formula_sm}   +   {formula_mp}    ->   {formula_B}   +   {formula_meoh}")
print(f" {mw_sm:7.3f} g/mol + {mw_mp:7.3f} g/mol -> {mw_B:7.3f} g/mol + {mw_meoh:7.3f} g/mol")
print(f"Total Reactant Mass: {reactants_B_mw:.3f} g/mol")
print(f"Total Product Mass:  {products_B_mw:.3f} g/mol")
print("-" * 55)

# 3. Formation of Product A
print("Pathway to Product A (Münchnone Cycloaddition):")
# Note: This pathway starts from Product C, the N-acetylated intermediate.
reactants_A_mw = mw_C + mw_mp
products_A_mw = mw_A + mw_co2
print(f"  {formula_C}   +   {formula_mp}    ->   {formula_A}   +   {formula_co2}")
print(f" {mw_C:7.3f} g/mol + {mw_mp:7.3f} g/mol -> {mw_A:7.3f} g/mol + {mw_co2:7.3f} g/mol")
print(f"Total Reactant Mass: {reactants_A_mw:.3f} g/mol")
print(f"Total Product Mass:  {products_A_mw:.3f} g/mol")
print("-" * 55)

# As per the final instruction, output the final calculated number.
# Let's use the molecular weight of the starting material.
print(f"The calculated molecular weight of the starting material ({formula_sm}) is: {mw_sm:.3f}")