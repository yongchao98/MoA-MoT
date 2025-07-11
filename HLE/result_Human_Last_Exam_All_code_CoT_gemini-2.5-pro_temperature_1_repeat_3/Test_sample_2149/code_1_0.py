# 1. Assign the deciphered values to variables.
# Y1 is derived from the year of the first US commercial oil well (1859).
Y1 = 1859

# "Hall" is derived from the patent year of the Hall-Héroult process (1886).
Hall = 1886

# Y4 is derived from the molecular weight of a typical Heck reaction reactant, iodobenzene (~204).
Y4 = 204

# 2. Formulate the equation based on the problem description.
# "Y4 to the Y1-Hall" translates to Y4 raised to the power of (Y1 - Hall).
exponent = Y1 - Hall
result = Y4 ** exponent

# 3. Print the equation with all its components, as requested.
print(f"Based on the clues, the equation is: {Y4} ** ({Y1} - {Hall})")

# 4. Print the final calculated value.
print(f"The power is: {exponent}")
print(f"The final calculated value (topological state index) is: {result}")