import sys
import io

# Backup the original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

# --- Your code starts here ---

# Define the parameters for the FET with example values.
# V_tg: Top gate voltage in Volts (V)
V_tg = 1.5
# V_bg: Back gate voltage in Volts (V)
V_bg = 0.8
# C_tg: Top gate capacitance per unit area in Farads per square meter (F/m^2)
C_tg = 1.7e-2
# C_bg: Back gate capacitance per unit area in Farads per square meter (F/m^2)
C_bg = 1.2e-2

# The displacement field D is determined by the free charges on the gate conductors,
# which are given by C*V. The dielectric constant epsilon_s is not needed for this calculation.

# Calculate the displacement field magnitude from the top gate.
# D_tg = C_tg * V_tg
D_tg_val = C_tg * V_tg

# Calculate the displacement field magnitude from the back gate.
# D_bg = C_bg * V_bg
D_bg_val = C_bg * V_bg

# The total displacement field is the superposition (vector sum) of the two fields.
# Defining the direction from the back gate to the top gate as positive, the fields oppose each other.
# D_total = D_bg - D_tg
D_total = D_bg_val - D_tg_val

# Print the explanation and the final equation with all numbers.
print("The displacement field (D) is calculated by superimposing the fields from the top and back gates.")
print("The formula is: D = (C_bg * V_bg) - (C_tg * V_tg)")
print("\nUsing the provided values:")
print(f"C_bg = {C_bg} F/m^2")
print(f"V_bg = {V_bg} V")
print(f"C_tg = {C_tg} F/m^2")
print(f"V_tg = {V_tg} V")
print("\nThe final equation is:")
# The core requirement: output each number in the final equation.
print(f"D = ({C_bg} * {V_bg}) - ({C_tg} * {V_tg})")
print(f"D = {D_bg_val:.4f} - {D_tg_val:.4f}")
print(f"D = {D_total:.4f} C/m^2")

if D_total > 0:
    print("\nThe positive sign indicates the net displacement field points from the back gate towards the top gate.")
elif D_total < 0:
    print("\nThe negative sign indicates the net displacement field points from the top gate towards the back gate.")
else:
    print("\nThe net displacement field is zero.")

# --- Your code ends here ---

# Restore the original stdout
sys.stdout = original_stdout
# Get the captured output
output = captured_output.getvalue()

# Print the final result to the user
print(output)