# Define the molecular formulas of the reactants
butadiene_formula = "C4H6"
dienophile_formula = "C2Cl2F2"

# In a Diels-Alder reaction, the atoms of the reactants combine to form a single product.
# We can determine the product's formula by summing the atoms.
reactant1_atoms = {'C': 4, 'H': 6}
reactant2_atoms = {'C': 2, 'Cl': 2, 'F': 2}

# Calculate the atoms in the product
product_atoms = {
    'C': reactant1_atoms.get('C', 0) + reactant2_atoms.get('C', 0),
    'H': reactant1_atoms.get('H', 0) + reactant2_atoms.get('H', 0),
    'Cl': reactant1_atoms.get('Cl', 0) + reactant2_atoms.get('Cl', 0),
    'F': reactant1_atoms.get('F', 0) + reactant2_atoms.get('F', 0),
}

# Construct the product formula string
product_formula = (
    f"C{product_atoms['C']}"
    f"H{product_atoms['H']}"
    f"Cl{product_atoms['Cl']}"
    f"F{product_atoms['F']}"
)

# The name of the product
product_name = "4,4-dichloro-5,5-difluorocyclohex-1-ene"

# Print the explanation and the final equation
print(f"The reaction between butadiene and 1,1-dichloro-2,2-difluoroethene is a Diels-Alder reaction.")
print(f"The product is {product_name}.")
print("\nThe balanced chemical equation is:")
print(f"{butadiene_formula} + {dienophile_formula} -> {product_formula}")