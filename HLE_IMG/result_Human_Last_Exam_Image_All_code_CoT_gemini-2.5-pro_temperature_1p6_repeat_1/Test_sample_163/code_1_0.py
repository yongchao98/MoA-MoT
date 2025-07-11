def get_molecular_info(composition):
    """Calculates molecular formula and weight from a dictionary of element counts."""
    atomic_weights = {
        'C': 12.011,
        'H': 1.008,
        'O': 15.999,
    }
    formula = ""
    mass = 0.0
    for element in ['C', 'H', 'O']:
        if element in composition:
            count = composition[element]
            formula += f"{element}{count}"
            mass += count * atomic_weights[element]
    return formula, mass

# Define product structures
product_A = {
    "name": "2-(tert-butoxy)-1-phenylethyl benzoate",
    "structure": "Ph-CH(OCOPh)-CH₂(OtBu)",
    "smiles": "c1ccc(cc1)C(OC(=O)c2ccccc2)COC(C)(C)C"
}

product_B = {
    "name": "2-benzoyloxy-1-tert-butoxy-1-phenylethane",
    "structure": "Ph-CH(OtBu)-CH₂(OCOPh)",
    "smiles": "c1ccc(cc1)C(OC(C)(C)C)COC(=O)c2ccccc2"
}

# The overall reaction is an addition, so both products have the same composition
# C(styrene) + C(peroxide) = 8 + 11 = 19
# H(styrene) + H(peroxide) = 8 + 14 = 22
# O(peroxide) = 3
product_composition = {'C': 19, 'H': 22, 'O': 3}

molecular_formula, molecular_weight = get_molecular_info(product_composition)

print("The two major products, A and B, are constitutional isomers.")
print("-" * 50)
print(f"Product A: {product_A['name']}")
print(f"Structure: {product_A['structure']}")
print(f"SMILES: {product_A['smiles']}")
print("-" * 50)
print(f"Product B: {product_B['name']}")
print(f"Structure: {product_B['structure']}")
print(f"SMILES: {product_B['smiles']}")
print("-" * 50)
print(f"Both products are isomers with the molecular formula: {molecular_formula}")
print(f"They have the same molecular weight of approximately {molecular_weight:.2f} g/mol.")
print("This explains the two major peaks observed in the GCMS.")
print("-" * 50)
print("The overall chemical equation is:")
# The numbers in the equation are: C8H8 + C11H14O3 -> C19H22O3
print("Styrene (C_8 H_8) + tert-butyl peroxybenzoate (C_11 H_14 O_3) -> Products (C_19 H_22 O_3)")
print("\nThe numbers in the final equation's formulae are: 8, 8, 11, 14, 3, 19, 22, 3")
