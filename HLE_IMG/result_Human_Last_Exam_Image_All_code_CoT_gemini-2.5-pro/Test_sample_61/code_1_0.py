def format_formula(atoms, name=""):
    """Formats a dictionary of atoms into a chemical formula string."""
    parts = []
    # Order of elements in the formula
    for element in ['C', 'H', 'O', 'Br']:
        if atoms.get(element, 0) > 0:
            count = atoms[element]
            parts.append(f"{element}{count if count > 1 else ''}")
    
    formula_str = "".join(parts)
    if name:
        print(f"Molecular formula of {name}: {formula_str}")
    return formula_str

def combine_formulas(base, modifier, operation):
    """Adds or subtracts atom counts from two formula dictionaries."""
    result = base.copy()
    for element, count in modifier.items():
        if operation == 'add':
            result[element] = result.get(element, 0) + count
        elif operation == 'subtract':
            result[element] = result.get(element, 0) - count
    return result

# 1. Define the molecular composition of the reactants
compound_1 = {'C': 11, 'H': 10, 'O': 3}
compound_2 = {'C': 7, 'H': 7, 'Br': 1}

format_formula(compound_1, "Compound 1")
format_formula(compound_2, "Compound 2")
print("\n--- Reaction Calculation ---")

# 2. Step 1: Alkylation
# The reaction joins Compound 1 and Compound 2, removing HBr in the process.
hbr = {'H': 1, 'Br': 1}
temp_product = combine_formulas(compound_1, compound_2, 'add')
alkylation_intermediate = combine_formulas(temp_product, hbr, 'subtract')
print("Step 1: Alkylation (addition of reactants, removal of HBr)")
format_formula(alkylation_intermediate, "Alkylated Intermediate")

# 3. Step 2: Saponification (Ester Hydrolysis)
# The methyl ester group (-OCH3) is replaced by a hydroxyl group (-OH).
# Net change is +H and -CH3, which is equivalent to +H2O and -CH3OH.
h2o = {'H': 2, 'O': 1}
methanol = {'C': 1, 'H': 4, 'O': 1}
temp_product_2 = combine_formulas(alkylation_intermediate, h2o, 'add')
beta_keto_acid = combine_formulas(temp_product_2, methanol, 'subtract')
print("\nStep 2: Saponification (hydrolysis of the methyl ester)")
format_formula(beta_keto_acid, "Beta-keto acid intermediate")


# 4. Step 3: Decarboxylation
# The carboxylic acid group (-COOH) is removed as CO2.
co2 = {'C': 1, 'O': 2}
compound_A = combine_formulas(beta_keto_acid, co2, 'subtract')
print("\nStep 3: Decarboxylation (removal of CO2)")
format_formula(compound_A, "Final Product A")


# 5. Print the final result in the requested format
print("\n--- Final Answer ---")
c = compound_A.get('C', 0)
h = compound_A.get('H', 0)
o = compound_A.get('O', 0)
print(f"The molecular formula of compound A is C{c}H{h}O{o}.")
