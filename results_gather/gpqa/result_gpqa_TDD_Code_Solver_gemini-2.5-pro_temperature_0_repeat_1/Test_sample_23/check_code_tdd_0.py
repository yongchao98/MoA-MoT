import math

def check_chemistry_solution():
    """
    This function checks the correctness of the proposed solution by verifying each constraint
    from the problem statement.

    The proposed solution is:
    - Substance Z: Cyclohexane (C6H12)
    - Mixture Y: Cyclohexane (C6H12) and Benzene (C6H6)
    - Mixture X: Cyclohexene (C6H10) and 1,4-Cyclohexadiene (C6H8)
    - Final Answer: Total H atoms in X = 10 + 8 = 18.
    """

    # --- Data representing the proposed chemical compounds and their properties ---
    # This data encodes the chemical knowledge required to solve the problem.
    CHEMICAL_DATA = {
        'Cyclohexane': {
            'formula': {'C': 6, 'H': 12},
            'is_saturated': True,
            'decolorizes_br2_water': False,
            'hydrogenation_product': 'Cyclohexane',
            'has_conjugated_bonds': False
        },
        'Benzene': {
            'formula': {'C': 6, 'H': 6},
            'is_saturated': False,
            'decolorizes_br2_water': False,  # Under standard conditions without a catalyst
            'hydrogenation_product': 'Cyclohexane',
            'has_conjugated_bonds': True # Aromatic system
        },
        'Cyclohexene': {
            'formula': {'C': 6, 'H': 10},
            'is_saturated': False,
            'decolorizes_br2_water': True,
            'hydrogenation_product': 'Cyclohexane',
            'has_conjugated_bonds': False
        },
        '1,4-Cyclohexadiene': {
            'formula': {'C': 6, 'H': 8},
            'is_saturated': False,
            'decolorizes_br2_water': True,
            'hydrogenation_product': 'Cyclohexane',
            'has_conjugated_bonds': False # Double bonds are isolated
        }
    }

    # Proposed identities
    z_name = 'Cyclohexane'
    y_names = ['Cyclohexane', 'Benzene']
    x_names = ['Cyclohexene', '1,4-Cyclohexadiene']

    # --- 1. Check Substance Z ---
    z_data = CHEMICAL_DATA[z_name]
    # Constraint: Mass fraction of hydrogen is 14.28% (which is ~1/7)
    # Using integer masses (C=12, H=1) for simplicity, as is common in such problems.
    h_mass_fraction = z_data['formula']['H'] / (z_data['formula']['C'] * 12 + z_data['formula']['H'] * 1)
    if not math.isclose(h_mass_fraction, 1/7, rel_tol=1e-4):
        return f"Incorrect: Substance Z ({z_name}) has a H mass fraction of {h_mass_fraction:.4f}, not ~14.28% (1/7)."
    # Constraint: Z does not react further with hydrogen (is saturated).
    if not z_data['is_saturated']:
        return f"Incorrect: Substance Z ({z_name}) is not saturated."

    # --- 2. Check Mixture Y ---
    # Constraint: Y does not decolorize bromine water.
    for name in y_names:
        if CHEMICAL_DATA[name]['decolorizes_br2_water']:
            return f"Incorrect: Component '{name}' of mixture Y should not decolorize bromine water."
    # Constraint: Hydrogenation of Y gives only Z.
    for name in y_names:
        if CHEMICAL_DATA[name]['hydrogenation_product'] != z_name:
            return f"Incorrect: Component '{name}' of mixture Y does not hydrogenate to Z ({z_name})."
    # Constraint: Z is a constituent of Y.
    if z_name not in y_names:
        return f"Incorrect: Z ({z_name}) must be a component of mixture Y."

    # --- 3. Check Mixture X ---
    # Constraint: X decolorizes bromine water.
    for name in x_names:
        if not CHEMICAL_DATA[name]['decolorizes_br2_water']:
            return f"Incorrect: Component '{name}' of mixture X should decolorize bromine water."
    # Constraint: Hydrogenation of X gives only Z.
    for name in x_names:
        if CHEMICAL_DATA[name]['hydrogenation_product'] != z_name:
            return f"Incorrect: Component '{name}' of mixture X does not hydrogenate to Z ({z_name})."
    # Constraint: No conjugated multiple bonds in X.
    for name in x_names:
        if CHEMICAL_DATA[name]['has_conjugated_bonds']:
            return f"Incorrect: Component '{name}' of mixture X has conjugated bonds, violating a key constraint."

    # --- 4. Check Disproportionation (X -> Y) ---
    # The reaction of an equimolar mixture X to an equimolar mixture Y implies atom balance:
    # (1 molecule of X1 + 1 molecule of X2) -> (1 molecule of Y1 + 1 molecule of Y2)
    x_total_c = sum(CHEMICAL_DATA[name]['formula']['C'] for name in x_names)
    x_total_h = sum(CHEMICAL_DATA[name]['formula']['H'] for name in x_names)
    y_total_c = sum(CHEMICAL_DATA[name]['formula']['C'] for name in y_names)
    y_total_h = sum(CHEMICAL_DATA[name]['formula']['H'] for name in y_names)

    if not (x_total_c == y_total_c and x_total_h == y_total_h):
        return f"Incorrect: The disproportionation reaction is not balanced. Atoms in X: C{x_total_c}H{x_total_h}. Atoms in Y: C{y_total_c}H{y_total_h}."

    # --- 5. Check Final Calculation ---
    # Question: Total number of hydrogen atoms in two liquids of mixture X.
    total_h_in_x = sum(CHEMICAL_DATA[name]['formula']['H'] for name in x_names)
    
    # The provided answer 'A' corresponds to 18.
    if total_h_in_x != 18:
        return f"Incorrect: The calculated total H atoms in X is {total_h_in_x}, which does not match the expected answer of 18."

    return "Correct"

# Run the checker and print the result
result = check_chemistry_solution()
print(result)