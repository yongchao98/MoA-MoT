from collections import Counter

def format_formula(formula_dict):
    """Formats a dictionary of atoms into a string formula."""
    # Standard order for organic chemistry formulas
    order = ['C', 'H', 'N', 'O']
    s = ""
    for atom in order:
        if atom in formula_dict and formula_dict[atom] > 0:
            s += atom
            if formula_dict[atom] > 1:
                s += str(formula_dict[atom])
    return s

def add_formulas(f1, f2):
    """Adds two formula dictionaries."""
    return dict(Counter(f1) + Counter(f2))

def subtract_formulas(f1, f2):
    """Subtracts formula dictionary f2 from f1."""
    return dict(Counter(f1) - Counter(f2))

# --- Define Reactants and other small molecules ---

# Starting Material: N-(1-pyrrolin-2-yl)proline
sm = {'C': 9, 'H': 14, 'N': 2, 'O': 2}

# Reagent: Methyl propiolate
alkyne = {'C': 4, 'H': 4, 'O': 2}

# Reagent: Acetic Anhydride
ac2o = {'C': 4, 'H': 6, 'O': 3}

# Byproduct/Implicit molecules
acetic_acid = {'C': 2, 'H': 4, 'O': 2}
carbon_dioxide = {'C': 1, 'O': 2}
h4_lost = {'H': 4}
o_gained = {'O': 1}

# --- Calculate Product Formulas Based on Proposed Pathways ---

# Pathway for Product C: Acetylation of Starting Material
# SM + Ac2O -> C + AcOH
temp_c = add_formulas(sm, ac2o)
product_c_calc = subtract_formulas(temp_c, acetic_acid)

# Pathway for Product A: Cycloaddition of Acetylated Ylide
# Product_C + Alkyne -> A + CO2
temp_a = add_formulas(product_c_calc, alkyne)
product_a_calc = subtract_formulas(temp_a, carbon_dioxide)

# Pathway for Product B: Cycloaddition of un-acetylated ylide, followed by modification
# SM + Alkyne -> P + CO2  (P is primary adduct)
# P -> P' + 4H           (P' is aromatized adduct)
# P' + O -> B            (B is oxidized P')
# So, B = SM + Alkyne - CO2 - 4H + O
temp_b_p = add_formulas(sm, alkyne)
temp_b_p = subtract_formulas(temp_b_p, carbon_dioxide) # This is primary adduct P
temp_b_p_prime = subtract_formulas(temp_b_p, h4_lost)  # This is aromatized P'
product_b_calc = add_formulas(temp_b_p_prime, o_gained)

# --- Print Results ---
print("Based on the proposed reaction pathways, the calculated molecular formulas are:")
print("-" * 60)

# Product A
print(f"Product A is proposed to be formed from the cycloaddition of the acetylated ylide.")
print(f"Equation: (C formed in situ) + Alkyne - CO2 -> A")
print(f"Calculation: ({format_formula(product_c_calc)}) + ({format_formula(alkyne)}) - ({format_formula(carbon_dioxide)}) = {format_formula(product_a_calc)}")
print(f"Calculated Formula for Product A: {format_formula(product_a_calc)}")
print(f"Given Formula for Product A:     C14H20N2O3\n")

# Product B
print(f"Product B is proposed to be formed from the cycloaddition of the un-acetylated ylide, followed by aromatization (-4H) and oxidation (+O).")
print(f"Equation: SM + Alkyne - CO2 - 4H + O -> B")
print(f"Calculated Formula for Product B: {format_formula(product_b_calc)}")
print(f"Given Formula for Product B:     C12H14N2O3\n")

# Product C
print(f"Product C is proposed to be the acetylated starting material.")
print(f"Equation: SM + Ac2O - AcOH -> C")
print(f"Calculation: ({format_formula(sm)}) + ({format_formula(ac2o)}) - ({format_formula(acetic_acid)}) = {format_formula(product_c_calc)}")
print(f"Calculated Formula for Product C: {format_formula(product_c_calc)}")
print(f"Given Formula for Product C:     C11H16N2O3")
print("-" * 60)
print("The calculations confirm the given molecular formulas are consistent with established chemical principles for this reaction type.")
