#
# This script identifies and describes the molecule shown in the image.
#

# Step 1: Define the number of repeating units observed in the macrocycle.
# By counting the benzene-like rings (phenylene groups), we find there are 8.
num_units = 8

# Step 2: Define the atomic composition of each part of the repeating unit.
# The repeating unit is a 'para-phenylene-ethynylene' block.
# A para-phenylene group (-C6H4-)
carbons_per_phenylene = 6
hydrogens_per_phenylene = 4
# An ethynylene group (-C≡C-)
carbons_per_ethynylene = 2

# Step 3: Provide the name of the molecule.
# The molecule is a cyclooctamer (n=8) of para-phenyleneethynylene.
# The common name is cyclo[n]para-phenyleneethynylene or [n]CPE.
print(f"The name of the molecule is: cyclo[{num_units}]para-phenyleneethynylene")
print(f"Abbreviated name: [{num_units}]CPE")


# Step 4: Calculate and display the molecular formula with the equation.
# Total carbons = (carbons from 8 phenylene groups) + (carbons from 8 ethynylene groups)
total_carbons = (num_units * carbons_per_phenylene) + (num_units * carbons_per_ethynylene)
# Total hydrogens = (hydrogens from 8 phenylene groups)
total_hydrogens = num_units * hydrogens_per_phenylene

print("\n--- Molecular Formula Calculation ---")
print("Equation for Carbon atoms (C):")
print(f"Total C = ({num_units} * {carbons_per_phenylene}) + ({num_units} * {carbons_per_ethynylene})")
print(f"Total C = {num_units * carbons_per_phenylene} + {num_units * carbons_per_ethynylene} = {total_carbons}")

print("\nEquation for Hydrogen atoms (H):")
print(f"Total H = {num_units} * {hydrogens_per_phenylene}")
print(f"Total H = {total_hydrogens}")

print(f"\nThe final molecular formula is: C{total_carbons}H{total_hydrogens}")