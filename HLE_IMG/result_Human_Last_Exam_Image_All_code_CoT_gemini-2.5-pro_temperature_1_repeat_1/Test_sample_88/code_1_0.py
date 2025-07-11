def format_formula(mol_dict):
    """Formats a molecule dictionary into a standard chemical formula string."""
    formula = ""
    # Order of atoms is conventionally C, H, N, O
    for atom in ['C', 'H', 'N', 'O']:
        count = mol_dict.get(atom, 0)
        if count > 0:
            formula += atom
            if count > 1:
                formula += str(count)
    return formula

def add_formulas(mol1, mol2):
    """Adds the atoms of two molecular formulas."""
    result = mol1.copy()
    for atom, count in mol2.items():
        result[atom] = result.get(atom, 0) + count
    return result

def subtract_formulas(mol1, mol2):
    """Subtracts the atoms of one molecular formula from another."""
    result = mol1.copy()
    for atom, count in mol2.items():
        result[atom] = result.get(atom, 0) - count
    return result

# --- Define Reactants and Products based on the problem description ---
# Starting Material: (3,4-dihydro-2H-pyrrol-5-yl)proline
# Structure analysis yields: 9 Carbons, 14 Hydrogens, 2 Nitrogens, 2 Oxygens
Starting_Material = {'C': 9, 'H': 14, 'N': 2, 'O': 2}
# Methyl propiolate: H-C≡C-COOCH3
Methyl_Propiolate = {'C': 4, 'H': 4, 'N': 0, 'O': 2}
# Acetic Anhydride: (CH3CO)2O
Acetic_Anhydride = {'C': 4, 'H': 6, 'N': 0, 'O': 3}

# Isolated Products
Product_A_formula = {'C': 14, 'H': 20, 'N': 2, 'O': 3}
Product_B_formula = {'C': 12, 'H': 14, 'N': 2, 'O': 3}
Product_C_formula = {'C': 11, 'H': 16, 'N': 2, 'O': 3}

# --- Define other relevant small molecules for reaction balancing ---
Acetic_Acid = {'C': 2, 'H': 4, 'N': 0, 'O': 2}
Carbon_Dioxide = {'C': 1, 'H': 0, 'N': 0, 'O': 2}
Methyl_Acetate = {'C': 3, 'H': 6, 'N': 0, 'O': 2}

print("--- Verifying the Identity of Products A, B, and C ---")

# Step 1: Identify Product C
# Hypothesis: Product C is the mixed anhydride formed from the starting material and acetic anhydride.
# Equation: Starting_Material + Acetic_Anhydride -> Product_C + Acetic_Acid
print("\n[Analysis for Product C]")
lhs_C = add_formulas(Starting_Material, Acetic_Anhydride)
rhs_C = add_formulas(Product_C_formula, Acetic_Acid)
print(f"Proposed reaction: {format_formula(Starting_Material)} + {format_formula(Acetic_Anhydride)} -> {format_formula(Product_C_formula)} + {format_formula(Acetic_Acid)}")
print(f"Verifying balance: {format_formula(lhs_C)} -> {format_formula(rhs_C)}")
if lhs_C == rhs_C:
    print("Result: The equation is balanced. Product C is the mixed anhydride intermediate.")
else:
    print("Result: The hypothesis for Product C is incorrect.")

# Step 2: Form the primary cycloadduct from Product C
# Hypothesis: Product C reacts with methyl propiolate to form a primary adduct.
# Equation: Product_C + Methyl_Propiolate -> Primary_Adduct
Primary_Adduct = add_formulas(Product_C_formula, Methyl_Propiolate)
print(f"\n[Formation of Primary Cycloadduct]")
print(f"Reaction: {format_formula(Product_C_formula)} + {format_formula(Methyl_Propiolate)} -> Primary_Adduct")
print(f"The calculated formula for the primary adduct is: {format_formula(Primary_Adduct)}")

# Step 3: Identify Product A
# Hypothesis: Product A is formed by the decarboxylation (loss of CO2) of the primary adduct.
# Equation: Primary_Adduct -> Product_A + Carbon_Dioxide
print("\n[Analysis for Product A]")
rhs_A = add_formulas(Product_A_formula, Carbon_Dioxide)
print(f"Proposed reaction: {format_formula(Primary_Adduct)} -> {format_formula(Product_A_formula)} + {format_formula(Carbon_Dioxide)}")
print(f"Verifying balance: {format_formula(Primary_Adduct)} -> {format_formula(rhs_A)}")
if Primary_Adduct == rhs_A:
    print("Result: The equation is balanced. Product A is the decarboxylated cycloadduct.")
else:
    print("Result: The hypothesis for Product A is incorrect.")

# Step 4: Identify Product B
# Hypothesis: Product B is formed by the elimination of methyl acetate from the primary adduct.
# Equation: Primary_Adduct -> Product_B + Methyl_Acetate
print("\n[Analysis for Product B]")
rhs_B = add_formulas(Product_B_formula, Methyl_Acetate)
print(f"Proposed reaction: {format_formula(Primary_Adduct)} -> {format_formula(Product_B_formula)} + {format_formula(Methyl_Acetate)}")
print(f"Verifying balance: {format_formula(Primary_Adduct)} -> {format_formula(rhs_B)}")
if Primary_Adduct == rhs_B:
    print("Result: The equation is balanced. Product B is formed from the primary adduct via elimination of methyl acetate.")
else:
    print("Result: The hypothesis for Product B is incorrect.")
