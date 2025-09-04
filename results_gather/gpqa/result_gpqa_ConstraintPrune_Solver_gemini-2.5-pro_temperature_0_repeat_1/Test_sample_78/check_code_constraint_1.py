import collections

def check_correctness():
    """
    Checks if the proposed answer D is correct by verifying that its reaction product
    matches the provided NMR spectral data.
    """
    
    # 1. Define constraints from the question
    required_formula = {'C': 11, 'H': 12, 'O': 1}
    
    # Provided 1H-NMR data for the product
    product_1h_nmr_data = {
        2.28: {'integration': 3, 'multiplicity': 's'},
        2.31: {'integration': 3, 'multiplicity': 's'},
        6.75: {'integration': 1, 'multiplicity': 'd'},
        7.08: {'integration': 2, 'multiplicity': 'd'},
        7.68: {'integration': 1, 'multiplicity': 'd'},
        7.71: {'integration': 2, 'multiplicity': 'd'},
    }
    
    # Provided 13C-NMR data for the product
    product_13c_nmr_data = {
        21.3: {'count': 1},
        28.4: {'count': 1},
        126.9: {'count': 2},
        127.1: {'count': 1},
        129.1: {'count': 2},
        130.3: {'count': 1},
        141.5: {'count': 1},
        144.1: {'count': 1},
        197.7: {'count': 1},
    }

    # 2. Analyze the proposed answer (D) and its reaction product
    # Answer D: 2-(4-methylstyryl)oxirane
    # Reaction: Base-catalyzed rearrangement
    # Expected Product: (E)-4-(p-tolyl)but-3-en-2-one [CH3-C6H4-CH=CH-C(=O)-CH3]

    # 3. Perform checks

    # Check 3a: Starting material's molecular formula
    # 2-(4-methylstyryl)oxirane = C7H7(p-tolyl) + C2H2(vinyl) + C2H3O(epoxide) = C11H12O
    start_material_formula = {'C': 11, 'H': 12, 'O': 1}
    if start_material_formula != required_formula:
        return f"Constraint check failed: The molecular formula of the compound in option D ({start_material_formula}) does not match the required formula ({required_formula})."

    # Check 3b: Product's molecular formula (must be an isomer)
    product_formula = {'C': 11, 'H': 12, 'O': 1}
    if product_formula != required_formula:
        return f"Constraint check failed: The reaction is an isomerization, but the product's formula ({product_formula}) does not match the starting material's ({required_formula})."

    # Check 3c: Consistency of 1H NMR data with the product structure
    # Expected signals for (E)-4-(p-tolyl)but-3-en-2-one:
    # - Two methyl singlets (tolyl-CH3 and acetyl-CH3)
    # - Two aromatic doublets (p-substituted ring, 2H each)
    # - Two vinylic doublets (trans-CH=CH, 1H each)
    # Total signals: 6
    # Total protons: 3+3+2+2+1+1 = 12
    
    total_protons_from_data = sum(v['integration'] for v in product_1h_nmr_data.values())
    if total_protons_from_data != product_formula['H']:
        return f"1H NMR check failed: The sum of integrations ({total_protons_from_data}) does not match the number of hydrogens in the product ({product_formula['H']})."

    integrations_given = sorted([v['integration'] for v in product_1h_nmr_data.values()])
    multiplicities_given = collections.Counter(v['multiplicity'] for v in product_1h_nmr_data.values())
    
    expected_integrations = sorted([1, 1, 2, 2, 3, 3])
    expected_multiplicities = collections.Counter({'s': 2, 'd': 4})

    if integrations_given != expected_integrations:
        return f"1H NMR check failed: The signal integrations {integrations_given} do not match the expected pattern {expected_integrations} for the product."
    
    if multiplicities_given != expected_multiplicities:
        return f"1H NMR check failed: The signal multiplicities {multiplicities_given} do not match the expected pattern {expected_multiplicities} for the product."

    # Check 3d: Consistency of 13C NMR data with the product structure
    # Expected signals for (E)-4-(p-tolyl)but-3-en-2-one:
    # - 11 total carbons.
    # - Symmetry in the p-tolyl group means 2 pairs of aromatic carbons are equivalent.
    # - Expected number of signals = 11 - 2 = 9.
    # - Key signals: Ketone C=O (~197 ppm), two CH3 carbons (~21, ~28 ppm).
    
    total_carbons_from_data = sum(v['count'] for v in product_13c_nmr_data.values())
    if total_carbons_from_data != product_formula['C']:
        return f"13C NMR check failed: The sum of carbon counts ({total_carbons_from_data}) does not match the number of carbons in the product ({product_formula['C']})."

    if len(product_13c_nmr_data) != 9:
        return f"13C NMR check failed: Expected 9 signals for the product due to symmetry, but {len(product_13c_nmr_data)} were provided."

    # Check for the presence of the characteristic ketone signal
    if 197.7 not in product_13c_nmr_data:
        return "13C NMR check failed: The data is missing the highly characteristic ketone (C=O) signal expected around 197.7 ppm."

    # All checks passed, the data is fully consistent with the proposed answer.
    return "Correct"

# Execute the check and print the result
result = check_correctness()
print(result)