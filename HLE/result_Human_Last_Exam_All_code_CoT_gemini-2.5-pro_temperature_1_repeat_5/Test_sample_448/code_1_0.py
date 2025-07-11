import sys
import io

# Set stdout to a string buffer to capture the output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

# This script identifies the rarest noble gas on Earth based on its atmospheric abundance
# and calculates its percentage concentration.

# Step 1: Define the abundance of noble gases in Earth's atmosphere in parts per million (ppm) by volume.
# Atmospheric abundance is a good proxy for overall terrestrial rarity.
# Data from 2002 is effectively the same as current data for this purpose.
noble_gas_abundance_ppm = {
    'Helium (He)': 5.2,
    'Neon (Ne)': 18.2,
    'Argon (Ar)': 9340.0,
    'Krypton (Kr)': 1.1,
    'Xenon (Xe)': 0.09,
    # Radon (Rn) is the rarest noble gas. It is radioactive with a short half-life,
    # so its concentration is extremely low and variable. A representative value
    # is on the order of 1 part per 10^21, which is 1e-15 ppm.
    'Radon (Rn)': 1e-15
}

# Step 2: Find the name and abundance of the rarest gas by finding the minimum value in the dictionary.
rarest_gas_name = min(noble_gas_abundance_ppm, key=noble_gas_abundance_ppm.get)
rarest_gas_ppm = noble_gas_abundance_ppm[rarest_gas_name]

# Step 3: Convert the abundance from parts per million (ppm) to a percentage.
# The conversion factor is 10,000 ppm per 1%.
conversion_factor = 10000
rarest_gas_percentage = rarest_gas_ppm / conversion_factor

# Step 4: Print the results, including the equation as requested.
print(f"To find the rarest noble gas, we compare their abundances. Based on atmospheric data from 2002 (which is unchanged), the rarest is {rarest_gas_name}.")
print(f"This is due to its radioactivity and short half-life, leading to an extremely low average concentration.")
print(f"\nIts approximate abundance in the atmosphere is {rarest_gas_ppm} parts per million (ppm).")
print(f"\nTo express this as a percentage, we perform the following calculation:")
print(f"Percentage = Abundance in ppm / {conversion_factor}")
# Final equation with numbers
print(f"Percentage = {rarest_gas_ppm} / {conversion_factor} = {rarest_gas_percentage:.1e} %")

# Restore stdout and get the captured output
sys.stdout = old_stdout
output = captured_output.getvalue()

# Print the captured output to the user
print(output)

# The question asks for the percentage value as the final answer.
# The calculated percentage is 1e-19 %.
final_answer_value = f"{rarest_gas_percentage:.1e}"