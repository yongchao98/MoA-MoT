# This script determines and prints the IUPAC name of the reaction product.

# The numbers required for the IUPAC name locants are identified first.
# These numbers indicate the positions of saturation and substitution on the heterocyclic ring.
locant_dihydro_1 = 2
locant_dihydro_2 = 5
locant_carboxylate = 3

# The parts of the name are assembled into the final IUPAC name.
# The name describes an ethyl ester of a 2,5-dihydrothiophene-3-carboxylic acid.
name_prefix = "ethyl"
ring_system = "dihydrothiophene"
name_suffix = "carboxylate"

# The final name is constructed using an f-string.
# The numbers {locant_dihydro_1}, {locant_dihydro_2}, and {locant_carboxylate} are explicitly included.
final_iupac_name = f"{name_prefix} {locant_dihydro_1},{locant_dihydro_2}-{ring_system}-{locant_carboxylate}-{name_suffix}"

print("The IUPAC name of the product is:")
print(final_iupac_name)