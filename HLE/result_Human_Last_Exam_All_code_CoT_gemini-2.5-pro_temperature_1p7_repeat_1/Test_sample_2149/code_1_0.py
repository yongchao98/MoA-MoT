# 1. Define the values deciphered from the clues.
Y1 = 1859 # From Clue 1: Edwin Drake's oil well discovery.
Y2 = 28868 # From Clue 2: Descartes-related riddle.
Y3 = 5988  # From Clue 2: Descartes-related riddle.
Y4 = 436  # From Clue 2: Russian Federal Law № 436-FZ.

# 2. Perform the final calculation to find the two topological state indices.
# The instruction implies two indices for the two reactants in the Heck reaction.
# A simple interpretation that uses all four numbers is subtractions between pairs.
index_reactant_1 = Y1 - Y4
index_reactant_2 = Y2 - Y3

# 3. Print the results, showing the full equation for each calculation.
print("The calculated indices for the two reactants are:")
print(f"Reactant 1 Index: {Y1} - {Y4} = {index_reactant_1}")
print(f"Reactant 2 Index: {Y2} - {Y3} = {index_reactant_2}")