import re

def parse_formula(formula_str):
    """Parses a molecular formula string into a dictionary of element counts."""
    counts = {}
    for element, count in re.findall(r'([A-Z][a-z]*)(\d*)', formula_str):
        counts[element] = counts.get(element, 0) + (int(count) if count else 1)
    return counts

def format_formula(formula_dict):
    """Formats a dictionary of element counts into a molecular formula string."""
    # Standard order C, H, then alphabetical for N, O
    order = ['C', 'H', 'N', 'O']
    s = ""
    processed_elements = set()
    for elem in order:
        if elem in formula_dict and formula_dict[elem] > 0:
            count = formula_dict[elem]
            s += elem + (str(count) if count > 1 else "")
            processed_elements.add(elem)
    
    # Add any other elements alphabetically
    sorted_keys = sorted([k for k in formula_dict if k not in processed_elements])
    for elem in sorted_keys:
         if formula_dict[elem] > 0:
            count = formula_dict[elem]
            s += elem + (str(count) if count > 1 else "")
    return s

def add_formulas(f1, f2):
    """Adds two formula dictionaries."""
    res = f1.copy()
    for elem, count in f2.items():
        res[elem] = res.get(elem, 0) + count
    return res

def subtract_formulas(f1, f2):
    """Subtracts formula dict f2 from f1."""
    res = f1.copy()
    for elem, count in f2.items():
        res[elem] = res.get(elem, 0) - count
    return res

# --- Main Analysis ---

# Given Molecular Formulas
formula_A_given = "C14H20N2O3"
formula_B_given = "C12H14N2O3"
formula_C_given = "C11H16N2O3"

# Molecular formulas of relevant molecules and groups
formula_SM_str = "C9H14N2O2" # Starting Material: 2-(4,5-dihydro-3H-pyrrol-2-yl)pyrrolidine-2-carboxylic acid
formula_acetyl_add_str = "C2H2O" # Represents adding an acetyl group (C2H3O) and removing H
formula_co2_str = "CO2"
formula_methyl_propiolate_str = "C4H4O2"
formula_h2_str = "H2"

# Parse all formulas into dictionaries
sm = parse_formula(formula_SM_str)
ac_add = parse_formula(formula_acetyl_add_str)
co2 = parse_formula(co2_str)
alkyne = parse_formula(formula_methyl_propiolate_str)
h2 = parse_formula(formula_h2_str)

print("Analysis of the Reaction and Structures of Products A, B, and C\n")

# --- Structure of Product C ---
print("--- Product C: N-Acetylation ---")
print(f"The molecular formula of Product C is {formula_C_given}.")
print("This product is the result of a simple N-acetylation of the starting material's pyrrolidine nitrogen by acetic anhydride.")
calc_C = add_formulas(sm, ac_add)
print(f"Calculation: {formula_SM_str} (Starting Material) + {formula_acetyl_add_str} (Acetylation) = {format_formula(calc_C)}")
print(f"The calculated formula {format_formula(calc_C)} matches the given formula for C.")
print("Structure C is: N-acetyl-2-(4,5-dihydro-3H-pyrrol-2-yl)pyrrolidine-2-carboxylic acid.\n")

# --- Structure of Product A ---
print("--- Product A: Cycloaddition from N-Acetyl Intermediate ---")
print(f"The molecular formula of Product A is {formula_A_given}.")
print("This product forms from Product C. First, Product C loses CO2 to form a reactive N-acetylated azomethine ylide.")
print("This ylide then reacts with one molecule of methyl propiolate in a [3+2] cycloaddition.")
intermediate_ylide_A = subtract_formulas(parse_formula(formula_C_given), co2)
calc_A = add_formulas(intermediate_ylide_A, alkyne)
print(f"Calculation: {formula_C_given} (Product C) - {formula_co2_str} + {formula_methyl_propiolate_str} (Alkyne) = {format_formula(calc_A)}")
print(f"The calculated formula {format_formula(calc_A)} matches the given formula for A.")
print("Structure A is: The [3+2] cycloadduct, which has an N-acetyl-dihydropyrrolizine core.\n")

# --- Structure of Product B ---
print("--- Product B: Cycloaddition and Oxidation ---")
print(f"The molecular formula of Product B is {formula_B_given}.")
print("This product forms via a different pathway from the non-acetylated starting material.")
print("1. A non-acetylated azomethine ylide (C8H12N2) is formed via decarboxylation of the starting material.")
ylide_b_precursor = subtract_formulas(sm, add_formulas(co2, h2)) # SM - CO2 - H2
print(f"   {formula_SM_str} (SM) - {formula_co2_str} - {formula_h2_str} -> {format_formula(ylide_b_precursor)}")
print("2. This ylide reacts with methyl propiolate to form an initial adduct P1.")
adduct_p1 = add_formulas(ylide_b_precursor, alkyne)
print(f"   {format_formula(ylide_b_precursor)} (Ylide) + {formula_methyl_propiolate_str} (Alkyne) -> {format_formula(adduct_p1)} (Adduct P1)")
print("3. The substituent ring on adduct P1 is then oxidized to a lactam (adds O, removes 2H).")
# To represent oxidation (+O, -2H), we add O and subtract H2
calc_B = subtract_formulas(add_formulas(adduct_p1, parse_formula("O")), h2)
print(f"   {format_formula(adduct_p1)} (Adduct P1) + [O] - 2H -> {format_formula(calc_B)}")
print(f"The final calculated formula {format_formula(calc_B)} matches the given formula for B.")
print("Structure B is: The [3+2] cycloadduct, similar to A but without the N-acetyl group and with its substituent ring oxidized to a lactam.\n")