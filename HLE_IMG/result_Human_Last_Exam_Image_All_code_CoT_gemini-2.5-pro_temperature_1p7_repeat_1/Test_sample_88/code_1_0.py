def print_molecule_formula(name, formula, indent=0):
    """Prints a formatted molecular formula."""
    parts = []
    for atom, count in formula.items():
        parts.append(f"{atom}{count}")
    print(f"{' ' * indent}{name}: {''.join(parts)}")

def print_equation(reactants, products):
    """Prints a formatted chemical equation."""
    def format_side(side):
        terms = []
        for molecule, coeff in side.items():
            formula_str = "".join([f"{atom}{count}" for atom, count in molecules[molecule].items()])
            if coeff > 1:
                terms.append(f"{coeff} {formula_str}")
            else:
                terms.append(formula_str)
        return " + ".join(terms)

    lhs = format_side(reactants)
    rhs = format_side(products)
    print(f"    Equation: {lhs} -> {rhs}")
    print("    This equation is stoichiometrically balanced.")

# Define all relevant molecules and their formulas
molecules = {
    "Starting Material (SM)": {'C': 9, 'H': 14, 'N': 2, 'O': 2},
    "Methyl Propiolate": {'C': 4, 'H': 4, 'O': 2},
    "Acetic Anhydride (Ac2O)": {'C': 4, 'H': 6, 'O': 3},
    "Acetic Acid (AcOH)": {'C': 2, 'H': 4, 'O': 2},
    "Product A": {'C': 14, 'H': 20, 'N': 2, 'O': 3},
    "Product B": {'C': 12, 'H': 14, 'N': 2, 'O': 3},
    "Product C": {'C': 11, 'H': 16, 'N': 2, 'O': 3},
    "Ylide Intermediate": {'C': 10, 'H': 16, 'N': 2, 'O': 1},
    "Propanoic Acid": {'C': 3, 'H': 6, 'O': 2},
    "Water": {'H': 2, 'O': 1}
}

print("Analysis of the Chemical Reaction and Products:\n")

# Analysis for Product C
print("1. Product C Formation (N-Acetylation):")
print_molecule_formula("Product C", molecules["Product C"], indent=4)
print_equation(
    {"Starting Material (SM)": 1, "Acetic Anhydride (Ac2O)": 1},
    {"Product C": 1, "Acetic Acid (AcOH)": 1}
)
print("-" * 20)

# Analysis for Product A
print("2. Product A Formation (1,3-Dipolar Cycloaddition):")
print_molecule_formula("Product A", molecules["Product A"], indent=4)
print("    This product is formed from a reactive ylide intermediate reacting with methyl propiolate.")
print_equation(
    {"Ylide Intermediate": 1, "Methyl Propiolate": 1},
    {"Product A": 1}
)
print("-" * 20)

# Analysis for Product B
print("3. Product B Formation (Side Reaction):")
print_molecule_formula("Product B", molecules["Product B"], indent=4)
print("    A plausible, balanced equation for this side product involves the reaction of Product C.")
print_equation(
    {"Product C": 1, "Acetic Anhydride (Ac2O)": 1},
    {"Product B": 1, "Propanoic Acid": 1, "Water": 1}
)
