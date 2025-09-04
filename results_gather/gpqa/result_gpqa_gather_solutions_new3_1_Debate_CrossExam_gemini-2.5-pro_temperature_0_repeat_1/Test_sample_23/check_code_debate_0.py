import re

def check_chemistry_problem():
    """
    This function programmatically verifies the solution to the chemistry problem
    by checking each constraint against the proposed chemical identities.
    """
    # --- Data and Constants ---
    # Using standard atomic masses for accurate calculations.
    ATOMIC_MASS = {'C': 12.011, 'H': 1.008}

    # Chemical properties based on known chemistry. This data represents the
    # ground truth for the compounds involved in the solution.
    CHEMICAL_DATA = {
        'cyclohexane': {
            'formula': 'C6H12',
            'is_saturated': True,
            'has_conjugated_bonds': False,
            'decolorizes_bromine_water': False,
        },
        'benzene': {
            'formula': 'C6H6',
            'is_saturated': False,  # Aromatic, not saturated like an alkane
            'is_aromatic': True,
            'has_conjugated_bonds': True,  # Aromatic ring is a form of conjugation
            'decolorizes_bromine_water': False,  # Key property for this problem
        },
        'cyclohexene': {
            'formula': 'C6H10',
            'is_saturated': False,
            'has_conjugated_bonds': False,
            'decolorizes_bromine_water': True,
        },
        '1,4-cyclohexadiene': {
            'formula': 'C6H8',
            'is_saturated': False,
            'has_conjugated_bonds': False,  # Key property: non-conjugated
            'decolorizes_bromine_water': True,
        },
        # Included for completeness to show why it's excluded by the problem's constraints
        '1,3-cyclohexadiene': {
            'formula': 'C6H8',
            'is_saturated': False,
            'has_conjugated_bonds': True, # This would violate a constraint
            'decolorizes_bromine_water': True,
        }
    }

    # --- Helper Function ---
    def parse_formula(formula):
        """Parses a chemical formula like 'C6H12' into a dict {'C': 6, 'H': 12}."""
        atoms = {}
        parts = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
        for element, count in parts:
            atoms[element] = int(count) if count else 1
        return atoms

    # --- Proposed Solution from the LLM Answer ---
    # The LLM's reasoning identifies the following compounds:
    Z_name = 'cyclohexane'
    Y_components = ['cyclohexane', 'benzene']
    X_components = ['cyclohexene', '1,4-cyclohexadiene']
    # The final answer is D, which corresponds to 18.
    final_answer_value = 18

    # --- Verification Steps ---

    # Step 1: Verify Substance Z (Cyclohexane)
    z_props = CHEMICAL_DATA[Z_name]
    z_atoms = parse_formula(z_props['formula'])
    
    # Constraint: Mass fraction of hydrogen is 14.28%
    h_mass = z_atoms.get('H', 0) * ATOMIC_MASS['H']
    c_mass = z_atoms.get('C', 0) * ATOMIC_MASS['C']
    total_mass = h_mass + c_mass
    h_fraction = h_mass / total_mass
    if not (0.142 <= h_fraction <= 0.144):  # Check if H fraction is ~14.3%
        return f"Constraint check failed for Z ({Z_name}): Hydrogen mass fraction is {h_fraction:.2%}, not ~14.28%."

    # Constraint: Z does not react further with hydrogen (is saturated)
    if not z_props['is_saturated']:
        return f"Constraint check failed for Z ({Z_name}): It is not a saturated compound."
    
    # Constraint: All compounds must have the same carbon skeleton as Z
    z_carbon_count = z_atoms.get('C', 0)
    for comp_list in [Y_components, X_components]:
        for comp_name in comp_list:
            comp_atoms = parse_formula(CHEMICAL_DATA[comp_name]['formula'])
            if comp_atoms.get('C', 0) != z_carbon_count:
                return f"Constraint check failed: Component {comp_name} does not have {z_carbon_count} carbons like Z."

    # Step 2: Verify Mixture Y
    # Constraint: Y does not decolorize bromine water
    for comp_name in Y_components:
        if CHEMICAL_DATA[comp_name]['decolorizes_bromine_water']:
            return f"Constraint check failed for Y: Component {comp_name} decolorizes bromine water."

    # Constraint: Z is a constituent of Y
    if Z_name not in Y_components:
        return f"Constraint check failed for Y: Z ({Z_name}) is not a constituent of mixture Y."

    # Step 3: Verify Mixture X
    # Constraint: X decolorizes bromine water
    for comp_name in X_components:
        if not CHEMICAL_DATA[comp_name]['decolorizes_bromine_water']:
            return f"Constraint check failed for X: Component {comp_name} does not decolorize bromine water."

    # Constraint: No conjugated multiple bonds in X
    for comp_name in X_components:
        if CHEMICAL_DATA[comp_name]['has_conjugated_bonds']:
            return f"Constraint check failed for X: Component {comp_name} has conjugated bonds."

    # Step 4: Verify the Disproportionation Reaction (Conservation of Atoms)
    x_total_atoms = {'C': 0, 'H': 0}
    for comp_name in X_components:
        atoms = parse_formula(CHEMICAL_DATA[comp_name]['formula'])
        x_total_atoms['C'] += atoms.get('C', 0)
        x_total_atoms['H'] += atoms.get('H', 0)

    y_total_atoms = {'C': 0, 'H': 0}
    for comp_name in Y_components:
        atoms = parse_formula(CHEMICAL_DATA[comp_name]['formula'])
        y_total_atoms['C'] += atoms.get('C', 0)
        y_total_atoms['H'] += atoms.get('H', 0)

    if x_total_atoms != y_total_atoms:
        return f"Constraint check failed: The reaction is not balanced. Reactants (X) have {x_total_atoms}, Products (Y) have {y_total_atoms}."

    # Step 5: Calculate the final answer and check correctness
    # Question: total number of hydrogen atoms in two liquids of mixture X.
    calculated_h_atoms = x_total_atoms['H']
    
    if calculated_h_atoms != final_answer_value:
        return f"Final answer is incorrect. The reasoning leads to a total of {calculated_h_atoms} hydrogen atoms, but the provided answer is {final_answer_value}."

    return "Correct"

# Run the check
result = check_chemistry_problem()
print(result)