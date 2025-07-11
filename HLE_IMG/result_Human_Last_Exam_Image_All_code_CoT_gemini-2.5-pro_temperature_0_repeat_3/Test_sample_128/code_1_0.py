def get_atom_counts(formula):
    """Parses a chemical formula string into a dictionary of atom counts."""
    import re
    # Regex to find elements and their optional counts
    pattern = re.compile(r'([A-Z][a-z]*)(\d*)')
    atoms = {}
    for element, count in pattern.findall(formula):
        atoms[element] = atoms.get(element, 0) + int(count if count else 1)
    return atoms

def print_balance_sheet(atoms, title):
    """Prints a formatted table of atom counts."""
    print(f"--- {title} ---")
    for atom, count in sorted(atoms.items()):
        print(f"  {atom}: {count}")
    print("-" * (len(title) + 8))

def solve_chemistry_problem():
    """
    Solves the reaction and verifies the atom balance.
    Reaction: 2-aminopyridine + o-phthalaldehyde + HCN -> Compound A + 2 H2O
    """
    # --- Define Molecules by their formulas ---
    aminopyridine_formula = "C5H6N2"
    opa_formula = "C8H6O2"
    hcn_formula = "C1H1N1"
    
    product_A_formula = "C14H9N3"
    water_formula = "H2O"

    # --- Stoichiometry ---
    reactants = {
        aminopyridine_formula: 1,
        opa_formula: 1,
        hcn_formula: 1
    }
    products = {
        product_A_formula: 1,
        water_formula: 2
    }

    # --- Calculate total atoms on each side ---
    reactant_atoms = {}
    for formula, moles in reactants.items():
        atoms = get_atom_counts(formula)
        for atom, count in atoms.items():
            reactant_atoms[atom] = reactant_atoms.get(atom, 0) + count * moles

    product_atoms = {}
    for formula, moles in products.items():
        atoms = get_atom_counts(formula)
        for atom, count in atoms.items():
            product_atoms[atom] = product_atoms.get(atom, 0) + count * moles

    # --- Print Results ---
    print("Verifying the atom balance for the reaction...")
    print_balance_sheet(reactant_atoms, "Reactant Atoms")
    print_balance_sheet(product_atoms, "Product Atoms")

    if reactant_atoms == product_atoms:
        print("\nThe chemical equation is balanced.")
    else:
        print("\nWarning: The chemical equation is NOT balanced.")

    print("\nThe final product, Compound A, is 2-(pyridin-2-yl)-2H-isoindole-3-carbonitrile.")
    print("Its chemical formula is:")
    
    # Print the final formula with each number individually as requested
    final_formula_parsed = get_atom_counts(product_A_formula)
    print(f"C{final_formula_parsed['C']}H{final_formula_parsed['H']}N{final_formula_parsed['N']}")


solve_chemistry_problem()