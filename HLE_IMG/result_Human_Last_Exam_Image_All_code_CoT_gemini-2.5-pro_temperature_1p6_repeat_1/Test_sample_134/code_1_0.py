import sys
import io

# Backup stdout
stdout_backup = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = io.StringIO()

# --- Start of the user-facing script ---

# The molecule is a large macrocycle composed of three types of building blocks:
# 1. Phenylene (-C6H4-): A benzene ring with two connection points.
# 2. Ethynylene (-C≡C-): A carbon-carbon triple bond.
# 3. Ethenylene (-CH=CH-): A carbon-carbon double bond, which are in the 'cis' configuration here.

# Chemical formula contributions (Carbon, Hydrogen) for each block
phenylene_formula = (6, 4)
ethynylene_formula = (2, 0)
ethenylene_formula = (2, 2)

# By carefully inspecting the image, we can count the number of each building block.
# There are two types of phenylene substitution: para (1,4) and meta (1,3).
num_para_phenylene = 3
num_meta_phenylene = 4
num_total_phenylene = num_para_phenylene + num_meta_phenylene

# Count of the linker groups
num_ethynylene = 7
num_ethenylene = 4

# --- Molecular Formula Calculation ---

# Calculate total Carbon atoms
c_from_phenylene = num_total_phenylene * phenylene_formula[0]
c_from_ethynylene = num_ethynylene * ethynylene_formula[0]
c_from_ethenylene = num_ethenylene * ethenylene_formula[0]
total_carbons = c_from_phenylene + c_from_ethynylene + c_from_ethenylene

# Calculate total Hydrogen atoms
h_from_phenylene = num_total_phenylene * phenylene_formula[1]
h_from_ethynylene = num_ethynylene * ethynylene_formula[1]
h_from_ethenylene = num_ethenylene * ethenylene_formula[1]
total_hydrogens = h_from_phenylene + h_from_ethynylene + h_from_ethenylene

# --- Output the results ---

print("Step 1: Analysis of molecular components from the image.")
print(f" - Number of total phenylene (-C6H4-) groups: {num_total_phenylene}")
print(f" - Number of ethynylene (-C≡C-) groups: {num_ethynylene}")
print(f" - Number of ethenylene (-CH=CH-) groups: {num_ethenylene}\n")

print("Step 2: Calculation of the molecular formula.")
print("Equation for Carbon atoms:")
print(f"   C = ({num_total_phenylene} phenylene × {phenylene_formula[0]} C) + ({num_ethynylene} ethynylene × {ethynylene_formula[0]} C) + ({num_ethenylene} ethenylene × {ethenylene_formula[0]} C)")
print(f"   C = {c_from_phenylene} + {c_from_ethynylene} + {c_from_ethenylene} = {total_carbons}\n")

print("Equation for Hydrogen atoms:")
print(f"   H = ({num_total_phenylene} phenylene × {phenylene_formula[1]} H) + ({num_ethynylene} ethynylene × {ethynylene_formula[1]} H) + ({num_ethenylene} ethenylene × {ethenylene_formula[1]} H)")
print(f"   H = {h_from_phenylene} + {h_from_ethynylene} + {h_from_ethenylene} = {total_hydrogens}\n")

print(f"The calculated molecular formula is C{total_carbons}H{total_hydrogens}.")

print("\nStep 3: Name of the molecule.")
print("Based on its structure and molecular formula, the molecule is known as:")
print("C64H36 Nanoring")

# --- End of the user-facing script ---

# Get the content of the buffer
output = sys.stdout.getvalue()
# Restore original stdout
sys.stdout = stdout_backup
# Print the captured output
print(output)