def print_formula(formula, name):
    """Prints the molecular formula in a standard format."""
    f_str = f"C{formula['C']}H{formula['H']}N{formula['N']}O{formula['O']}"
    print(f"The calculated formula for {name} is: {f_str}")

def add_formulas(f1, f2):
    """Adds two molecular formulas."""
    return {atom: f1.get(atom, 0) + f2.get(atom, 0) for atom in set(f1) | set(f2)}

def subtract_formulas(f1, f2):
    """Subtracts one molecular formula from another."""
    return {atom: f1.get(atom, 0) - f2.get(atom, 0) for atom in set(f1) | set(f2)}

# --- Define Molecular Formulas ---
# Starting Material (SM): N-(1-pyrrolin-2-yl)proline
sm = {'C': 9, 'H': 14, 'N': 2, 'O': 2}
# Methyl Propiolate (MP)
mp = {'C': 4, 'H': 4, 'N': 0, 'O': 2}
# Acetyl group (from Acetic Anhydride)
acetyl = {'C': 2, 'H': 3, 'N': 0, 'O': 1}
# Hydrogen atom
h_atom = {'C': 0, 'H': 1, 'N': 0, 'O': 0}
# Methanol (MeOH)
meoh = {'C': 1, 'H': 4, 'N': 0, 'O': 1}
# Carbon Dioxide (CO2)
co2 = {'C': 1, 'H': 0, 'N': 0, 'O': 2}

# Target product formulas
product_A_target = {'C': 14, 'H': 20, 'N': 2, 'O': 3}
product_B_target = {'C': 12, 'H': 14, 'N': 2, 'O': 3}
product_C_target = {'C': 11, 'H': 16, 'N': 2, 'O': 3}


print("--- Verifying formation of Product C ---")
print("Hypothesis: Product C is formed by acetylation of the Starting Material (SM).")
print(f"Equation: SM - H + Acetyl = C")
# C = SM - H + Acetyl
temp_formula = subtract_formulas(sm, h_atom)
product_C_calc = add_formulas(temp_formula, acetyl)
print(f"C{sm['C']}H{sm['H']}N{sm['N']}O{sm['O']} - H + C{acetyl['C']}H{acetyl['H']}O{acetyl['O']} = C{product_C_calc['C']}H{product_C_calc['H']}N{product_C_calc['N']}O{product_C_calc['O']}")
print_formula(product_C_calc, "Product C")
print(f"Target formula for Product C is C{product_C_target['C']}H{product_C_target['H']}N{product_C_target['N']}O{product_C_target['O']}. Match: {product_C_calc == product_C_target}\n")


print("--- Verifying formation of Product B ---")
print("Hypothesis: Product B is formed from SM reacting with Methyl Propiolate (MP) and eliminating Methanol (MeOH).")
print(f"Equation: SM + MP - MeOH = B")
# B = SM + MP - MeOH
temp_formula = add_formulas(sm, mp)
product_B_calc = subtract_formulas(temp_formula, meoh)
print(f"C{sm['C']}H{sm['H']}N{sm['N']}O{sm['O']} + C{mp['C']}H{mp['H']}O{mp['O']} - C{meoh['C']}H{meoh['H']}O{meoh['O']} = C{product_B_calc['C']}H{product_B_calc['H']}N{product_B_calc['N']}O{product_B_calc['O']}")
print_formula(product_B_calc, "Product B")
print(f"Target formula for Product B is C{product_B_target['C']}H{product_B_target['H']}N{product_B_target['N']}O{product_B_target['O']}. Match: {product_B_calc == product_B_target}\n")


print("--- Verifying formation of Product A ---")
print("Hypothesis: Product A is formed from Product C reacting with MP, followed by decarboxylation (loss of CO2).")
print(f"Equation: C + MP - CO2 = A")
# A = C + MP - CO2
temp_formula = add_formulas(product_C_calc, mp)
product_A_calc = subtract_formulas(temp_formula, co2)
print(f"C{product_C_calc['C']}H{product_C_calc['H']}N{product_C_calc['N']}O{product_C_calc['O']} + C{mp['C']}H{mp['H']}O{mp['O']} - C{co2['C']}O{co2['O']} = C{product_A_calc['C']}H{product_A_calc['H']}N{product_A_calc['N']}O{product_A_calc['O']}")
print_formula(product_A_calc, "Product A")
print(f"Target formula for Product A is C{product_A_target['C']}H{product_A_target['H']}N{product_A_target['N']}O{product_A_target['O']}. Match: {product_A_calc == product_A_target}\n")