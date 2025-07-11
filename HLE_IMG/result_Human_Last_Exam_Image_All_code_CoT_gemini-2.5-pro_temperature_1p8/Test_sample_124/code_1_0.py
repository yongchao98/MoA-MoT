import sys
import io

# Backup stdout
stdout_backup = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = buffer = io.StringIO()

# --- Start of Student Code ---

# Define the component values from the circuit diagram
V_source = 1.0  # Volts
R3 = 9.0      # Ohms
R7 = 100.0    # Ohms
R8 = 100.0    # Ohms

# The circuit simplifies to R3 in parallel with the series combination of R7 and R8.
# R1 and R2 are short-circuited and do not draw current from the source.

# Calculate the resistance of the R7-R8 series branch
R_78 = R7 + R8

# Calculate the current through each parallel branch
I_R3 = V_source / R3
I_R78 = V_source / R_78

# Calculate the total current by summing the currents from the parallel branches
I_total = I_R3 + I_R78

# Print the final equation with all the numbers
print(f"Total Current = (V / R3) + (V / (R7 + R8))")
print(f"Total Current = ({V_source} V / {R3} Ω) + ({V_source} V / ({R7} Ω + {R8} Ω))")
print(f"Total Current = {I_R3:.4f} A + {I_R78:.4f} A")
print(f"Total Current = {I_total:.4f} A")

# --- End of Student Code ---

# Restore stdout
sys.stdout = stdout_backup
# Get the content of the buffer
output = buffer.getvalue()

# Print the captured output
print(output)
# Finally, print the answer in the required format
print(f'<<<{I_total:.4f}>>>')