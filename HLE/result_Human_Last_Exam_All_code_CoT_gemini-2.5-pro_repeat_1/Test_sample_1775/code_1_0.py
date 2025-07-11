import sys
import io

# Backup original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

# --- User's Code Starts Here ---

# Define the parameters with example values (in SI units)
# Ctg: Top gate capacitance per area in Farads per square meter (F/m^2)
Ctg = 1.5e-2
# Vtg: Top gate voltage in Volts (V)
Vtg = 2.0
# Cbg: Back gate capacitance per area in Farads per square meter (F/m^2)
Cbg = 0.8e-2
# Vbg: Back gate voltage in Volts (V)
Vbg = 5.0
# The dielectric constant of the transistor material (epsilon_s) is not needed for this calculation.

# Calculate the total displacement field (D)
# The unit of displacement field is Coulombs per square meter (C/m^2)
D_total = Ctg * Vtg + Cbg * Vbg

# Output the explanation and the final equation with numerical values
print("The total displacement field (D) is the sum of the contributions from the top and back gates.")
print("The formula is: D = C_tg * V_tg + C_bg * V_bg")
print(f"Using the given values:")
print(f"D = ({Ctg} F/m^2) * ({Vtg} V) + ({Cbg} F/m^2) * ({Vbg} V)")
# Calculate intermediate steps for clarity
D_tg = Ctg * Vtg
D_bg = Cbg * Vbg
print(f"D = {D_tg} C/m^2 + {D_bg} C/m^2")
print(f"The total displacement field is: D = {D_total} C/m^2")

# --- User's Code Ends Here ---

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output = captured_output.getvalue()

# Print the output to the user
print(output)

# Extract the final numerical answer for the <<<>>> format
final_answer = D_total
print(f"<<<{final_answer}>>>")