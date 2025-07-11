# The reaction between butadiene and 1,1-dichloro-2,2-difluoroethene is a Diels-Alder reaction.
# This script constructs the IUPAC name of the resulting product.

# According to IUPAC nomenclature rules for cyclic compounds, we number the ring
# to give the double bond positions 1 and 2. Substituents are then assigned the
# lowest possible locants. "chloro" is cited before "fluoro" due to alphabetical order,
# so it receives the lower-numbered position.

# Define locants for the 'dichloro' substituent
locant_cl_1 = 4
locant_cl_2 = 4

# Define locants for the 'difluoro' substituent
locant_f_1 = 5
locant_f_2 = 5

# Define the names of the substituents and the parent ring
substituent_cl = "dichloro"
substituent_f = "difluoro"
parent_ring = "cyclohexene"

# Assemble the final name, outputting each number as part of the string
final_name = f"{locant_cl_1},{locant_cl_2}-{substituent_cl}-{locant_f_1},{locant_f_2}-{substituent_f}{parent_ring}"

print(final_name)