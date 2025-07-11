# Number of specimens for each family
dixidae_specimens = 100
simuliidae_specimens = 1101
belpharoceridae_specimens = 100
rhagionidae_specimens = 101
tabanidae_specimens = 201

# Number of prolegs/parapodia per specimen for each family
prolegs_per_dixid = 6
prolegs_per_simuliid = 2
prolegs_per_tabanid = 42
prolegs_per_belpharocerid = 0
prolegs_per_rhagionid = 0

# Calculate the subtotal of prolegs for each family
dixid_total_prolegs = dixidae_specimens * prolegs_per_dixid
simuliid_total_prolegs = simuliidae_specimens * prolegs_per_simuliid
tabanid_total_prolegs = tabanidae_specimens * prolegs_per_tabanid

# Calculate the grand total
total_prolegs = dixid_total_prolegs + simuliid_total_prolegs + tabanid_total_prolegs

# Print the final equation and the total number
print(f"Based on the specimen counts and larval morphology:")
print(f"Dixidae: {dixidae_specimens} specimens * {prolegs_per_dixid} prolegs/specimen = {dixid_total_prolegs} prolegs")
print(f"Simuliidae: {simuliidae_specimens} specimens * {prolegs_per_simuliid} prolegs/specimen = {simuliid_total_prolegs} prolegs")
print(f"Tabanidae: {tabanidae_specimens} specimens * {prolegs_per_tabanid} prolegs/specimen = {tabanid_total_prolegs} prolegs")
print(f"Belpharoceridae and Rhagionidae larvae do not have prolegs.")
print("\nTotal calculation:")
print(f"{dixid_total_prolegs} + {simuliid_total_prolegs} + {tabanid_total_prolegs} = {total_prolegs}")

print(f"\nThe total number of prolegs and parapodia is {total_prolegs}.")