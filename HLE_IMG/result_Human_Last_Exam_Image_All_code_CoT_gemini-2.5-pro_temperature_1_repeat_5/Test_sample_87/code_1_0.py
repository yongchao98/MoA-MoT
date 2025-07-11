from collections import Counter

def formula_to_dict(formula):
    """Converts a molecular formula string (e.g., 'C2H4O2') to a dictionary."""
    import re
    return Counter({
        m[1]: int(m[2] or 1)
        for m in re.finditer(r'(([A-Z][a-z]?)(\d*))', formula)
    })

def dict_to_str(d):
    """Converts a dictionary to a molecular formula string."""
    # Standard order for organic chemistry
    order = ['C', 'H', 'N', 'O']
    s = ""
    for elem in order:
        if elem in d and d[elem] > 0:
            s += elem
            if d[elem] > 1:
                s += str(d[elem])
            s += " " # Add space for readability
    return s.strip()

def add_formulas(d1, d2):
    """Adds two formula dictionaries."""
    return d1 + d2

def sub_formulas(d1, d2):
    """Subtracts the second formula dictionary from the first."""
    return d1 - d2

# --- Given Molecular Formulas ---
sm_formula = "C9H14N2O2"
mp_formula = "C4H4O2"
ac2o_formula = "C4H6O3"
acoh_formula = "C2H4O2"
co2_formula = "CO2"
h4_formula = "H4"
o_formula = "O"

prod_A_formula = "C14H20N2O3"
prod_B_formula = "C12H14N2O3"
prod_C_formula = "C11H16N2O3"

# --- Convert to Dictionaries ---
sm_dict = formula_to_dict(sm_formula)
mp_dict = formula_to_dict(mp_formula)
ac2o_dict = formula_to_dict(ac2o_formula)
acoh_dict = formula_to_dict(acoh_formula)
co2_dict = formula_to_dict(co2_formula)
h4_dict = formula_to_dict(h4_formula)
o_dict = formula_to_dict(o_formula)

print("This script verifies the identities of products A, B, and C based on their molecular formulas.")
print("-" * 50)
print("Reactant Formulas:")
print(f"Starting Material (SM): {dict_to_str(sm_dict)}")
print(f"Methyl Propiolate (MP): {dict_to_str(mp_dict)}")
print(f"Acetic Anhydride (Ac2O): {dict_to_str(ac2o_formula)}")
print(f"Acetic Acid (AcOH): {dict_to_str(acoh_dict)}")
print("-" * 50)

# --- Verify Product A ---
print("Verifying Product A (C14 H20 N2 O3):")
print("Proposed pathway: Decarboxylation -> Cycloaddition -> Acetylation")
ylide1_dict = sub_formulas(sm_dict, co2_dict)
adduct_p_dict = add_formulas(ylide1_dict, mp_dict)
calc_A_dict = sub_formulas(add_formulas(adduct_p_dict, ac2o_dict), acoh_dict)
print("The calculation is: (SM - CO2) + MP + Ac2O - AcOH")
print(f"({dict_to_str(sm_dict)} - {dict_to_str(co2_dict)}) + {dict_to_str(mp_dict)} + {dict_to_str(ac2o_dict)} - {dict_to_str(acoh_dict)} = {dict_to_str(calc_A_dict)}")
print(f"Calculated formula for A matches the given formula: {dict_to_str(calc_A_dict) == dict_to_str(formula_to_dict(prod_A_formula))}")
print("Structure A: N-acetylated cycloadduct.")
print("-" * 50)

# --- Verify Product C ---
print("Verifying Product C (C11 H16 N2 O3):")
print("Proposed pathway: Reaction of SM with Ac2O")
calc_C_dict = sub_formulas(add_formulas(sm_dict, ac2o_dict), acoh_dict)
print("The calculation is: SM + Ac2O - AcOH")
print(f"{dict_to_str(sm_dict)} + {dict_to_str(ac2o_dict)} - {dict_to_str(acoh_dict)} = {dict_to_str(calc_C_dict)}")
print(f"Calculated formula for C matches the given formula: {dict_to_str(calc_C_dict) == dict_to_str(formula_to_dict(prod_C_formula))}")
print("Structure C: Mesoionic 1,3-dipole formed as a side product.")
print("-" * 50)

# --- Verify Product B ---
print("Verifying Product B (C12 H14 N2 O3):")
print("Proposed pathway: Oxidation of the primary adduct P")
# Adduct P was calculated during the analysis of A
print("Adduct P = (SM - CO2) + MP")
print(f"P = ({dict_to_str(sm_dict)} - {dict_to_str(co2_dict)}) + {dict_to_str(mp_dict)} = {dict_to_str(adduct_p_dict)}")
calc_B_dict = sub_formulas(add_formulas(adduct_p_dict, o_dict), h4_dict)
print("The calculation is: Adduct P + O - 4H")
print(f"{dict_to_str(adduct_p_dict)} + {dict_to_str(o_dict)} - {dict_to_str(h4_dict)} = {dict_to_str(calc_B_dict)}")
print(f"Calculated formula for B matches the given formula: {dict_to_str(calc_B_dict) == dict_to_str(formula_to_dict(prod_B_formula))}")
print("Structure B: Oxidized cycloadduct.")
print("-" * 50)
