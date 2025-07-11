# The final product of the reaction is an alpha-cyano sulfide.
# This script will construct and print its IUPAC name.

# The locant number indicating the position of the substituent.
locant = 2

# The name of the substituent group (Ph-S-).
substituent = "Phenylsulfanyl"

# The name of the parent molecule based on the nitrile functional group.
parent_molecule = "ethanenitrile"

# Construct the full IUPAC name using an f-string.
# The number from the name is explicitly included as per the problem description.
iupac_name = f"{locant}-({substituent}){parent_molecule}"

print(iupac_name)