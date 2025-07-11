import math

# Step 1 & 2: Define the Y values based on the clues.
Y1 = 1886
Y4 = 210

# Step 3: Define the topological indices for the reactants.
# Reactant 1: Iodobenzene (7 non-hydrogen atoms)
index1 = 7
# Reactant 2: Styrene (8 non-hydrogen atoms)
index2 = 8

# Step 4: Perform the calculation as interpreted from the prompt.
# Formula: Index ^ (Y4 / Y1)
result1 = index1 ** (Y4 / Y1)
result2 = index2 ** (Y4 / Y1)

# Print the results in the required format, showing the full equation.
print(f"Reactant 1 (Iodobenzene):")
print(f"{index1} ^ ({Y4} / {Y1}) = {result1}")

print(f"\nReactant 2 (Styrene):")
print(f"{index2} ^ ({Y4} / {Y1}) = {result2}")