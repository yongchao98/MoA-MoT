import collections

# Step 1: Define the reactants and products of the Wittig reaction.
# Reactants
reactant_aldehyde = {
    "name": "Pivalaldehyde",
    "formula": "C5H10O"
}
reactant_ylide = {
    "name": "(2-(2-chlorophenyl)ethylidene)triphenyl-l5-phosphane",
    "formula": "C26H22ClP"
}

# Products
# The main organic product is an alkene.
product_alkene = {
    "name": "1-(2-chlorophenyl)-4,4-dimethylpent-2-ene",
    "formula": "C13H17Cl",
    "smiles": "CC(C)(C)C=CHCH2c1c(Cl)cccc1"
}
# The byproduct is triphenylphosphine oxide.
product_byproduct = {
    "name": "Triphenylphosphine oxide",
    "formula": "C18H15OP"
}

# The stoichiometric coefficients are all 1 for this reaction.
stoichiometry = collections.defaultdict(lambda: 1)

# Step 2: Display the balanced chemical equation.
# The prompt requires outputting each number in the final equation,
# which refers to the stoichiometric coefficients.
print("--- Wittig Reaction Equation ---\n")
equation = (
    f"{stoichiometry['aldehyde']} {reactant_aldehyde['name']} + "
    f"{stoichiometry['ylide']} {reactant_ylide['name']} -> "
    f"{stoichiometry['alkene']} {product_alkene['name']} + "
    f"{stoichiometry['byproduct']} {product_byproduct['name']}"
)
print(equation)

# Step 3: Provide details about the main organic product.
print("\n--- Main Organic Product Details ---\n")
print(f"Product Name: {product_alkene['name']}")
print(f"Chemical Formula: {product_alkene['formula']}")
print(f"SMILES String: {product_alkene['smiles']}")