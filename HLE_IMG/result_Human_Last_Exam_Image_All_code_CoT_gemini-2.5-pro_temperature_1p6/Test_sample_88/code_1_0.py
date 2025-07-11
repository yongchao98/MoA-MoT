import re
from collections import Counter

def parse_formula(formula):
    """Parses a chemical formula string into a Counter of atom counts."""
    pairs = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
    atom_counts = Counter()
    for atom, count in pairs:
        count = int(count) if count else 1
        atom_counts[atom] += count
    return atom_counts

def format_formula_from_counter(atom_counts):
    """Formats a Counter of atom counts into a canonical chemical formula string."""
    order = ['C', 'H', 'N', 'O']
    formula_parts = []
    # Add atoms in the canonical order first
    for atom in order:
        if atom in atom_counts and atom_counts[atom] > 0:
            count = atom_counts[atom]
            formula_parts.append(f"{atom}{count if count > 1 else ''}")
    # Add any other atoms sorted alphabetically
    for atom, count in sorted(atom_counts.items()):
        if atom not in order and count > 0:
            formula_parts.append(f"{atom}{count if count > 1 else ''}")
    return ''.join(formula_parts)

def verify_reaction(equation_str):
    """
    Verifies if a chemical equation string is balanced.
    Example: "C9H14N2O2 + C4H6O3 -> C11H16N2O3 + C2H4O2"
    """
    reactants_str, products_str = equation_str.split('->')
    
    reactant_formulas = [s.strip() for s in reactants_str.split('+')]
    product_formulas = [s.strip() for s in products_str.split('+')]

    reactant_atoms = Counter()
    for formula in reactant_formulas:
        reactant_atoms.update(parse_formula(formula))
        
    product_atoms = Counter()
    for formula in product_formulas:
        product_atoms.update(parse_formula(formula))
        
    is_balanced = (reactant_atoms == product_atoms)
    
    print(f"Verifying equation: {equation_str}")
    print(f"Reactant atoms: {dict(reactant_atoms)}")
    print(f"Product atoms:  {dict(product_atoms)}")
    print(f"Result: {'Balanced' if is_balanced else 'Not Balanced'}\n")
    return is_balanced

# --- Define Formulas and Equations ---

# Provided molecular formulas
sm1_formula = "C9H14N2O2"
mp_formula = "C4H4O2"
ac2o_formula = "C4H6O3"
prod_a_formula = "C14H20N2O3"
prod_b_formula = "C12H14N2O3"
prod_c_formula = "C11H16N2O3"

# Likely byproducts
co2_formula = "CO2"
ch3cooh_formula = "C2H4O2"
ch4o_formula = "CH4O"

# Proposed reaction equations based on atom balancing
equation_c = f"{sm1_formula} + {ac2o_formula} -> {prod_c_formula} + {ch3cooh_formula}"
equation_b = f"{sm1_formula} + {mp_formula} -> {prod_b_formula} + {ch4o_formula}"
equation_a = f"{sm1_formula} + {mp_formula} + {ac2o_formula} -> {prod_a_formula} + {co2_formula} + {ch3cooh_formula}"

# --- Execute Verification ---
print("--- Verifying the Stoichiometry of Product Formation ---\n")
verify_reaction(equation_c)
verify_reaction(equation_b)
verify_reaction(equation_a)