# This script calculates the molar mass of Compound A.

# Atomic weights of the elements involved (in g/mol)
atomic_weights = {
    'C': 12.011,  # Carbon
    'H': 1.008,   # Hydrogen
    'N': 14.007,   # Nitrogen
    'O': 15.999   # Oxygen
}

# The molecular formula for Compound A (2-((3-hydroxypyridin-2-yl)(phenylamino))acetonitrile)
# is determined by adding the atoms from the reactants and accounting for the reaction steps.
# Formula: C13H11N3O
counts = {'C': 13, 'H': 11, 'N': 3, 'O': 1}

# Calculate the total molar mass
molar_mass = (counts['C'] * atomic_weights['C'] +
             counts['H'] * atomic_weights['H'] +
             counts['N'] * atomic_weights['N'] +
             counts['O'] * atomic_weights['O'])

# Print the identity of Compound A and the molar mass calculation
print("Compound A is 2-((3-hydroxypyridin-2-yl)(phenylamino))acetonitrile.")
print("Its molecular formula is C13H11N3O.")
print("\nThe molar mass calculation is as follows:")
print(f"{counts['C']} * {atomic_weights['C']} + {counts['H']} * {atomic_weights['H']} + {counts['N']} * {atomic_weights['N']} + {counts['O']} * {atomic_weights['O']} = {molar_mass:.3f} g/mol")