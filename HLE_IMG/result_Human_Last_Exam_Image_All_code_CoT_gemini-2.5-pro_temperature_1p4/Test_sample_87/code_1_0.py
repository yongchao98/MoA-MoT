def parse_formula(formula):
    """Parses a molecular formula string into a dictionary of element counts."""
    import re
    elements = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
    counts = {}
    for el, count in elements:
        counts[el] = counts.get(el, 0) + int(count) if count else counts.get(el, 0) + 1
    return counts

def format_formula(counts):
    """Formats a dictionary of element counts into a molecular formula string."""
    return ''.join([f"{el}{counts[el]}" for el in sorted(counts.keys()) if counts[el] > 0])

def add_formulas(f1, f2):
    """Adds two molecular formulas."""
    res = f1.copy()
    for el, count in f2.items():
        res[el] = res.get(el, 0) + count
    return res

def subtract_formulas(f1, f2):
    """Subtracts a molecular formula f2 from f1."""
    res = f1.copy()
    for el, count in f2.items():
        res[el] = res.get(el, 0) - count
    return res

# Define molecular formulas from the problem
sm_formula = parse_formula("C9H14N2O2")
mp_formula = parse_formula("C4H4O2")
prod_A_formula = parse_formula("C14H20N2O3")
prod_B_formula = parse_formula("C12H14N2O3")
prod_C_formula = parse_formula("C11H16N2O3")

# Formulas of other small molecules
ketene_formula = parse_formula("C2H2O")
co2_formula = parse_formula("CO2")
h2o_formula = parse_formula("H2O")
meoh_formula = parse_formula("CH4O")

# --- Verify Pathway to C ---
print("--- Verifying Pathway to Product C ---")
calculated_C = add_formulas(sm_formula, ketene_formula)
print(f"Starting Material ({format_formula(sm_formula)}) + Ketene ({format_formula(ketene_formula)}) = {format_formula(calculated_C)}")
print(f"Product C Formula: {format_formula(prod_C_formula)}")
print(f"Match: {calculated_C == prod_C_formula}\n")

# --- Verify Pathway to A (Decarboxylative) ---
print("--- Verifying Pathway to Product A ---")
# 1. Ylide formation
ylide_Y = subtract_formulas(sm_formula, co2_formula)
print(f"Step 1: Ylide (Y) = SM ({format_formula(sm_formula)}) - CO2 ({format_formula(co2_formula)}) = {format_formula(ylide_Y)}")
# 2. Cycloaddition
adduct_P_prime = add_formulas(ylide_Y, mp_formula)
print(f"Step 2: Adduct (P') = Y ({format_formula(ylide_Y)}) + MP ({format_formula(mp_formula)}) = {format_formula(adduct_P_prime)}")
# 3. Acetylation
calculated_A = add_formulas(adduct_P_prime, ketene_formula)
print(f"Step 3: Calculated A = P' ({format_formula(adduct_P_prime)}) + Ketene ({format_formula(ketene_formula)}) = {format_formula(calculated_A)}")
print(f"Product A Formula: {format_formula(prod_A_formula)}")
print(f"Match: {calculated_A == prod_A_formula}\n")

# --- Verify Pathway to B (Non-decarboxylative) ---
print("--- Verifying Pathway to Product B ---")
# 1. Dipole formation
dipole_D = subtract_formulas(sm_formula, h2o_formula)
print(f"Step 1: Dipole (D) = SM ({format_formula(sm_formula)}) - H2O ({format_formula(h2o_formula)}) = {format_formula(dipole_D)}")
# 2. Cycloaddition
adduct_P = add_formulas(dipole_D, mp_formula)
print(f"Step 2: Adduct (P) = D ({format_formula(dipole_D)}) + MP ({format_formula(mp_formula)}) = {format_formula(adduct_P)}")
# 3. Hydrolysis
calculated_B = subtract_formulas(adduct_P, meoh_formula)
print(f"Step 3: Calculated B = P ({format_formula(adduct_P)}) - MeOH ({format_formula(meoh_formula)}) = {format_formula(calculated_B)}")
print(f"Product B Formula: {format_formula(prod_B_formula)}")
print(f"Match: {calculated_B == prod_B_formula}\n")