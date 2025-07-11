import sys
import io

# Backup the original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

# --- Start of the student's code ---

# Step 1 & 2: Define the molecular components and their counts.
# The molecule is composed of phenylene, naphthylene, vinylene, and ethynylene units.
# Atomic composition of each unit (as a substituent in a larger molecule).
units_data = {
    "phenylene": {'C': 6, 'H': 4, 'count': 2},
    "naphthylene": {'C': 10, 'H': 6, 'count': 2},
    "vinylene": {'C': 2, 'H': 2, 'count': 6},
    "ethynylene": {'C': 2, 'H': 0, 'count': 10}
}

# Step 3: Calculate the total number of Carbon and Hydrogen atoms.
total_C = 0
total_H = 0

carbon_calc_str = []
carbon_calc_values = []
hydrogen_calc_str = []
hydrogen_calc_values = []

for unit_name, data in units_data.items():
    C_contrib = data['count'] * data['C']
    H_contrib = data['count'] * data['H']
    
    total_C += C_contrib
    total_H += H_contrib
    
    carbon_calc_str.append(f"({data['count']} * {data['C']})")
    carbon_calc_values.append(str(C_contrib))
    
    hydrogen_calc_str.append(f"({data['count']} * {data['H']})")
    hydrogen_calc_values.append(str(H_contrib))

# Step 4: Print the calculations and the final name.

# Print the calculation for Carbon atoms
print("Calculation for Carbon atoms:")
# To satisfy the "output each number in the final equation" requirement
print(" + ".join(carbon_calc_str), f"= {total_C}")


# Print the calculation for Hydrogen atoms
print("\nCalculation for Hydrogen atoms:")
print(" + ".join(hydrogen_calc_str), f"= {total_H}")


# Print the name
print("\nThis molecule is a carbon-rich nanoring reported in the journal Organic Letters in 2021.")
print("A descriptive name is:")
print(f"Carbon Nanoring C{total_C}H{total_H}")

# --- End of the student's code ---

# Restore the original stdout
sys.stdout = original_stdout
# Get the captured output
output = captured_output.getvalue()
print(output)

final_name = f"Carbon Nanoring C{total_C}H{total_H}"
print(f"<<<{final_name}>>>")