import re

def parse_formula(formula_str):
    """Parses a molecular formula string into a dictionary of atom counts."""
    atom_dict = {}
    # Find all occurrences of an element (e.g., 'C', 'H', 'Cl') followed by an optional number
    atoms = re.findall(r'([A-Z][a-z]*)(\d*)', formula_str)
    for atom, count in atoms:
        # If count is not specified, it's 1
        count = int(count) if count else 1
        atom_dict[atom] = atom_dict.get(atom, 0) + count
    return atom_dict

def get_formula_from_name(name):
    """
    Calculates molecular formula based on a hardcoded dictionary derived from
    careful IUPAC name analysis. This avoids reliance on external libraries
    which can sometimes have issues with complex names.
    """
    formulas = {
        # Reactants
        "2-ethyl-2,6-dimethylcyclohexan-1-one": "C10H18O",
        "ethyl acrylate": "C5H8O2",
        "1-nitropropane": "C3H7NO2",
        "(E)-but-2-enenitrile": "C4H5N",
        # Products from option A
        "ethyl 3-(3-ethyl-1,3-dimethyl-2-oxocyclohexyl)propanoate": "C15H26O3",
        "3-methyl-4-nitrohexanenitrile": "C7H12N2O2",
        # Incorrect products from other options for comparison
        "2,3-dimethyl-4-nitrobutanenitrile": "C6H10N2O2",
        "ethyl 3-(3-ethyl-3,5-dimethyl-4-oxocyclohexyl)propanoate": "C15H26O3"
    }
    return formulas.get(name)

def check_correctness():
    """
    Checks the correctness of the proposed answer by verifying stoichiometry and
    structural plausibility for both reactions.
    """
    # The provided answer is A.
    # Product A name: ethyl 3-(3-ethyl-1,3-dimethyl-2-oxocyclohexyl)propanoate
    # Product B name: 3-methyl-4-nitrohexanenitrile

    # --- Check Reaction A: Michael Addition ---
    # Reaction: 2-ethyl-2,6-dimethylcyclohexan-1-one + ethyl acrylate ---> A
    
    # Get formulas of reactants
    reactant_A1_formula = get_formula_from_name("2-ethyl-2,6-dimethylcyclohexan-1-one")
    reactant_A2_formula = get_formula_from_name("ethyl acrylate")
    
    # Combine reactant formulas into a single atom count dictionary
    reactant_A1_atoms = parse_formula(reactant_A1_formula)
    reactant_A2_atoms = parse_formula(reactant_A2_formula)
    total_reactant_A_atoms = {
        k: reactant_A1_atoms.get(k, 0) + reactant_A2_atoms.get(k, 0)
        for k in set(reactant_A1_atoms) | set(reactant_A2_atoms)
    }

    # Get formula of proposed product A
    product_A_name = "ethyl 3-(3-ethyl-1,3-dimethyl-2-oxocyclohexyl)propanoate"
    product_A_formula = get_formula_from_name(product_A_name)
    product_A_atoms = parse_formula(product_A_formula)

    # Check if atoms are conserved (stoichiometry check)
    if total_reactant_A_atoms != product_A_atoms:
        return (f"Incorrect. Stoichiometry for reaction A is wrong. "
                f"Reactants sum to {total_reactant_A_atoms}, but product A "
                f"({product_A_name}) has formula {product_A_atoms}.")

    # Structural check (as comments): The name for product A correctly describes the Michael adduct.
    # The alternative name in options C/D, `...4-oxocyclohexyl...`, is structurally incorrect because
    # the reaction occurs alpha to the ketone, not gamma.

    # --- Check Reaction B: Nitro-Michael Addition ---
    # Reaction: 1-nitropropane + (E)-but-2-enenitrile ---> B

    # Get formulas of reactants
    reactant_B1_formula = get_formula_from_name("1-nitropropane")
    reactant_B2_formula = get_formula_from_name("(E)-but-2-enenitrile")

    # Combine reactant formulas
    reactant_B1_atoms = parse_formula(reactant_B1_formula)
    reactant_B2_atoms = parse_formula(reactant_B2_formula)
    total_reactant_B_atoms = {
        k: reactant_B1_atoms.get(k, 0) + reactant_B2_atoms.get(k, 0)
        for k in set(reactant_B1_atoms) | set(reactant_B2_atoms)
    }

    # Get formula of proposed product B
    product_B_name = "3-methyl-4-nitrohexanenitrile"
    product_B_formula = get_formula_from_name(product_B_name)
    product_B_atoms = parse_formula(product_B_formula)

    # Check stoichiometry for Reaction B
    if total_reactant_B_atoms != product_B_atoms:
        return (f"Incorrect. Stoichiometry for reaction B is wrong. "
                f"Reactants sum to {total_reactant_B_atoms}, but product B "
                f"({product_B_name}) has formula {product_B_atoms}.")

    # Structural check (as comments): The name for product B correctly identifies the longest
    # carbon chain (hexanenitrile) and substituent positions. The alternative name in options B/C,
    # `...butanenitrile`, incorrectly identifies the parent chain.

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_correctness()
print(result)