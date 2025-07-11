import sys
import io

# Redirect stdout to capture print output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

# Plan:
# 1. Deconstruct the molecule into its constituent chemical units based on visual analysis.
# 2. Count the number of each type of unit.
# 3. Define the chemical formula for each unit.
# 4. Calculate the total molecular formula by summing the atoms from all units.
# 5. Print the analysis, the calculation, and the identified name of the molecule.

# Step 1 & 2: Define the constituent units and their counts in the macrocycle.
units_counts = {
    "Naphthylene (-C10H6-)": 2,
    "Phenylene (-C6H4-)": 4,
    "Ethene (-CH=CH-)": 3,
    "Alkyne (-C≡C-)": 9
}

# Step 3: Define the atomic composition (C, H) for each unit.
atomic_composition = {
    "Naphthylene (-C10H6-)": {"C": 10, "H": 6},
    "Phenylene (-C6H4-)": {"C": 6, "H": 4},
    "Ethene (-CH=CH-)": {"C": 2, "H": 2},
    "Alkyne (-C≡C-)": {"C": 2, "H": 0}
}

print("Analysis of the Molecular Structure:")
print("------------------------------------")
print("The molecule is a macrocycle composed of the following units:")
for unit, count in units_counts.items():
    print(f"- {count} units of {unit}")

# Step 4: Calculate the total number of carbon and hydrogen atoms.
C_from_naphthylene = units_counts["Naphthylene (-C10H6-)"] * atomic_composition["Naphthylene (-C10H6-)"]["C"]
H_from_naphthylene = units_counts["Naphthylene (-C10H6-)"] * atomic_composition["Naphthylene (-C10H6-)"]["H"]

C_from_phenylene = units_counts["Phenylene (-C6H4-)"] * atomic_composition["Phenylene (-C6H4-)"]["C"]
H_from_phenylene = units_counts["Phenylene (-C6H4-)"] * atomic_composition["Phenylene (-C6H4-)"]["H"]

C_from_ethene = units_counts["Ethene (-CH=CH-)"] * atomic_composition["Ethene (-CH=CH-)"]["C"]
H_from_ethene = units_counts["Ethene (-CH=CH-)"] * atomic_composition["Ethene (-CH=CH-)"]["H"]

C_from_alkyne = units_counts["Alkyne (-C≡C-)"] * atomic_composition["Alkyne (-C≡C-)"]["C"]
H_from_alkyne = units_counts["Alkyne (-C≡C-)"] * atomic_composition["Alkyne (-C≡C-)"]["H"]

total_C = C_from_naphthylene + C_from_phenylene + C_from_ethene + C_from_alkyne
total_H = H_from_naphthylene + H_from_phenylene + H_from_ethene + H_from_alkyne

print("\nCalculation of the Molecular Formula:")
print("-------------------------------------")
print(f"Carbon atoms = {C_from_naphthylene} (from Naphthylene) + {C_from_phenylene} (from Phenylene) + {C_from_ethene} (from Ethene) + {C_from_alkyne} (from Alkyne) = {total_C}")
print(f"Hydrogen atoms = {H_from_naphthylene} (from Naphthylene) + {H_from_phenylene} (from Phenylene) + {H_from_ethene} (from Ethene) + {H_from_alkyne} (from Alkyne) = {total_H}")
print(f"\nThe calculated molecular formula is C{total_C}H{total_H}.")

# Step 5: Identify and print the name of the molecule.
molecule_name = "Keviyne"
print(f"\nThis molecule has been synthesized and named in scientific literature.")
print(f"The name of the molecule is: {molecule_name}")

# Restore stdout and print the captured output
sys.stdout = old_stdout
print(captured_output.getvalue())