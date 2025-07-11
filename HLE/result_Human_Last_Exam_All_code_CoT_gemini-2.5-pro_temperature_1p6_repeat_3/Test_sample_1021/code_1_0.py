import re

def calculate_molar_mass(formula):
    """Calculates the molar mass of a chemical formula."""
    atomic_mass = {
        'C': 12.011, 'H': 1.008, 'O': 16.00, 'N': 14.007,
        'Na': 22.990, 'Br': 79.904, 'I': 126.90, 'K': 39.098
    }
    
    # Find all element-count pairs in the formula
    elements = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
    
    molar_mass = 0
    for element, count in elements:
        count = int(count) if count else 1
        if element in atomic_mass:
            molar_mass += atomic_mass[element] * count
        else:
            raise ValueError(f"Atomic mass for element {element} not found.")
            
    return molar_mass

def analyze_reaction_stoichiometry():
    """
    Analyzes the stoichiometry of the SN2 ethylation reaction and prints a report.
    """
    # Reaction details
    mass_sm_g = 10.0
    eq_nah = 2.5
    eq_etbr = 3.0
    
    # Molecular formulas
    formula_sm = "C11H10O2"  # 2-Methyl-1,4-naphthalenediol
    formula_nah = "NaH"
    formula_etbr = "C2H5Br"
    formula_product = "C15H18O2" # 2-Methyl-1,4-diethoxynaphthalene

    # --- Calculations ---
    # Molar masses
    mm_sm = calculate_molar_mass(formula_sm)
    mm_nah = calculate_molar_mass(formula_nah)
    mm_etbr = calculate_molar_mass(formula_etbr)
    mm_product = calculate_molar_mass(formula_product)
    
    # Moles of starting material (limiting reagent)
    moles_sm = mass_sm_g / mm_sm
    
    # Moles of reagents used
    moles_nah_used = moles_sm * eq_nah
    moles_etbr_used = moles_sm * eq_etbr
    
    # Theoretical yield
    moles_product_theoretical = moles_sm
    mass_product_theoretical_g = moles_product_theoretical * mm_product

    # --- Print Report ---
    print("Reaction Stoichiometry Analysis:")
    print("-" * 40)
    
    print(f"Reaction: Ethylation of {formula_sm}")
    print("\nBalanced Stoichiometric Equation:")
    # The coefficients are the numbers in the equation.
    print(f"1 C11H10O2 + 2 NaH + 2 C2H5Br -> 1 C15H18O2 + 2 NaBr + 2 H2\n")

    print(f"--- Starting Material: {formula_sm} ---")
    print(f"Molar Mass: {mm_sm:.2f} g/mol")
    print(f"Mass used: {mass_sm_g:.1f} g")
    print(f"Moles used: {moles_sm:.4f} mol (1.0 eq - Limiting Reagent)")
    
    print(f"\n--- Reagents ---")
    print(f"Base: Sodium Hydride ({formula_nah})")
    print(f"  Molar Mass: {mm_nah:.2f} g/mol")
    print(f"  Equivalents used: {eq_nah}")
    print(f"  Moles used: {moles_nah_used:.4f} mol")
    
    print(f"\nElectrophile: Ethyl Bromide ({formula_etbr})")
    print(f"  Molar Mass: {mm_etbr:.2f} g/mol")
    print(f"  Equivalents used: {eq_etbr}")
    print(f"  Moles used: {moles_etbr_used:.4f} mol")
    
    print(f"\n--- Product: {formula_product} ---")
    print(f"Molar Mass: {mm_product:.2f} g/mol")
    print(f"Theoretical Moles: {moles_product_theoretical:.4f} mol")
    print(f"Theoretical Yield: {mass_product_theoretical_g:.2f} g")
    print("-" * 40)

if __name__ == '__main__':
    analyze_reaction_stoichiometry()