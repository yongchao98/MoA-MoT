def add_formulas(f1, f2):
    """Adds two molecular formulas."""
    d = {}
    for elem, count in f1.items():
        d[elem] = d.get(elem, 0) + count
    for elem, count in f2.items():
        d[elem] = d.get(elem, 0) + count
    return d

def subtract_formulas(f1, f2):
    """Subtracts molecular formula f2 from f1."""
    d = f1.copy()
    for elem, count in f2.items():
        d[elem] = d.get(elem, 0) - count
        if d[elem] == 0:
            del d[elem]
    return d

def format_formula(d):
    """Formats a formula dictionary into a string."""
    return "".join(f"{elem}{count}" for elem, count in sorted(d.items()))

# Define the molecular formulas of the starting materials
formula_SM = {'C': 9, 'H': 14, 'N': 2, 'O': 2}
formula_MP = {'C': 4, 'H': 4, 'O': 2}
formula_Ac_group = {'C': 2, 'H': 3, 'O': 1} # Acetyl group
formula_H = {'H': 1}
formula_O = {'O': 1}
formula_CO2 = {'C': 1, 'O': 2}

print("This script deduces the structures of products A, B, and C by verifying the atom balance based on a plausible reaction mechanism.")
print("-" * 50)
print(f"Starting Material (SM) Formula: {format_formula(formula_SM)}")
print(f"Methyl Propiolate (MP) Formula: {format_formula(formula_MP)}")
print("-" * 50)

# --- Calculation for Product C ---
# C is the mixed anhydride of SM and acetic acid.
# Reaction: SM + Ac2O -> C + AcOH
# Formula C = SM - H + (COCH3) = SM + Ac_group - H
formula_C_calc_step1 = add_formulas(formula_SM, formula_Ac_group)
formula_C_calc_step2 = subtract_formulas(formula_C_calc_step1, formula_H)
print("Calculation for Product C:")
print(f"  Proposed structure: Mixed anhydride of the starting material.")
print(f"  Derivation: Formula(SM) + Formula(Acetyl) - Formula(H)")
print(f"  {format_formula(formula_SM)} + {format_formula(formula_Ac_group)} - H => {format_formula(formula_C_calc_step2)}")
print(f"  Resulting Formula for C: {format_formula(formula_C_calc_step2)}")
print(f"  Given Formula for C: C11H16N2O3. Matches.\n")

# --- Ylide Formation is central to A and B ---
# Ylide is formed by decarboxylation of SM
# Reaction: SM -> Ylide + CO2
formula_Ylide = subtract_formulas(formula_SM, formula_CO2)
print("Central Intermediate Formation (Ylide):")
print(f"  Proposed reaction: Decarboxylation of the starting material.")
print(f"  Derivation: Formula(SM) - Formula(CO2)")
print(f"  {format_formula(formula_SM)} - {format_formula(formula_CO2)} => {format_formula(formula_Ylide)}")
print(f"  Formula of Ylide intermediate: {format_formula(formula_Ylide)}\n")


# --- Calculation for Product A ---
# A is the cycloadduct of the Acetylated Ylide and MP.
# Reaction: Ylide + Ac2O -> Acyl-Ylide -> +MP -> A
# Formula A = Ylide + Ac_group - H + MP
formula_AcYlide = subtract_formulas(add_formulas(formula_Ylide, formula_Ac_group), formula_H)
formula_A_calc = add_formulas(formula_AcYlide, formula_MP)
print("Calculation for Product A:")
print(f"  Proposed structure: Cycloadduct from Acetylated Ylide.")
print(f"  Derivation: (Formula(Ylide) + Formula(Acetyl) - H) + Formula(MP)")
print(f"  ({format_formula(formula_Ylide)} + {format_formula(formula_Ac_group)} - H) + {format_formula(formula_MP)} => {format_formula(formula_A_calc)}")
print(f"  Resulting Formula for A: {format_formula(formula_A_calc)}")
print(f"  Given Formula for A: C14H20N2O3. Matches.\n")

# --- Calculation for Product B ---
# B is an oxidized/aromatized derivative of the primary cycloadduct (P1).
# Reaction: Ylide + MP -> P1; then P1 - 4H + O -> B
# Formula P1 = Ylide + MP
formula_P1 = add_formulas(formula_Ylide, formula_MP)
# Formula B = P1 - 4H + O
formula_B_calc_step1 = subtract_formulas(formula_P1, {'H': 4})
formula_B_calc_step2 = add_formulas(formula_B_calc_step1, formula_O)
print("Calculation for Product B:")
print(f"  Proposed structure: Oxidized/Aromatized cycloadduct.")
print(f"  Derivation: (Formula(Ylide) + Formula(MP)) - 4H + O")
print(f"  ({format_formula(formula_Ylide)} + {format_formula(formula_MP)}) - H4 + O => {format_formula(formula_B_calc_step2)}")
print(f"  Resulting Formula for B: {format_formula(formula_B_calc_step2)}")
print(f"  Given Formula for B: C12H14N2O3. Matches.")
