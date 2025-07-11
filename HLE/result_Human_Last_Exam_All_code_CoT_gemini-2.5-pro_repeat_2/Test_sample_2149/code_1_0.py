# Step 1: Decode Y1 from Clue 1
# The clue points to Herbert Dow founding the Dow Chemical company in 1897.
y1 = 1897

# Step 2: Decode Y4 from Clue 2
# The clue leads to the Fantômas phone number "ГОра 32-04".
# Removing "ГО" (the "гэ" sound) leaves "ра 32-04".
# "ра" is the acronym for Russia's Day of Rocket Forces and Artillery (19/11 -> 1911).
# The three numbers are 1911, 32, and 4.
# So, Y2 = 1911, Y3 = 32, Y4 = 4. We only need Y4 for the calculation.
y4 = 4

# Step 3: Decode reactant indices from Clue 3
# The original Heck reaction reactants are Iodobenzene and Styrene.
# We calculate their integer molecular masses.
# Iodobenzene (C6H5I): (6 * 12) + (5 * 1) + 127
reactant1_mass = 204
# Styrene (C8H8): (8 * 12) + (8 * 1)
reactant2_mass = 104

# Step 4: Assemble and perform the final calculation
# The formula is (Y1 - (sum of reactant masses)) raised to the power of Y4.
base = y1 - (reactant1_mass + reactant2_mass)
result = base ** y4

# Step 5: Print the final equation with all numbers
print(f"The decoded values are:")
print(f"Y1 = {y1}")
print(f"Reactant 1 (Iodobenzene) Mass = {reactant1_mass}")
print(f"Reactant 2 (Styrene) Mass = {reactant2_mass}")
print(f"Y4 = {y4}")
print("\nFinal Calculation:")
# Using f-string to format the output equation clearly
print(f"({y1} - ({reactant1_mass} + {reactant2_mass})) ** {y4} = {result}")

# The final answer in the requested format
print(f"\n<<<{result}>>>")