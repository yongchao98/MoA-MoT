import sys
import io

# Backup the original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

# --- Start of the student's code ---
# Component values from the circuit diagram
V2 = 1.0  # Volts
R1 = 3.0  # Ohms
R2 = 7.0  # Ohms
R3 = 9.0  # Ohms
R7 = 100.0 # Ohms
R8 = 100.0 # Ohms

# Step 1: Calculate the equivalent resistance of the right branch (R7 and R8 in series)
Req_right = R7 + R8

# Step 2: Determine the equivalent resistance of the left branch.
# A wire from the top node connects to the junction between R2 and R3.
# This wire is in parallel with the series combination of R1 and R2.
# A zero-resistance wire in parallel with any resistance results in a total resistance of zero for that part.
# This means R1 and R2 are shorted out. The current bypasses them and flows to R3.
# Therefore, the equivalent resistance of the left branch is simply R3.
Req_left = R3

# Step 3: Calculate the total equivalent resistance of the circuit.
# The left branch (Req_left) and the right branch (Req_right) are in parallel.
Req_total = (Req_left * Req_right) / (Req_left + Req_right)

# Step 4: Calculate the total current using Ohm's Law (I = V / R)
I_total = V2 / Req_total

# Print the final equation with all the numbers
print("To find the total current, we first find the total equivalent resistance (Req).")
print("The left branch resistance is effectively R3 due to a short circuit across R1 and R2.")
print(f"Req_left = R3 = {R3} Ω")
print("The right branch resistance is R7 and R8 in series.")
print(f"Req_right = R7 + R8 = {R7} + {R8} = {Req_right} Ω")
print("The total resistance is the parallel combination of the left and right branches.")
print(f"Req_total = (Req_left * Req_right) / (Req_left + Req_right) = ({Req_left} * {Req_right}) / ({Req_left} + {Req_right}) = {Req_total:.3f} Ω")
print("\nNow, we use Ohm's Law (I = V / R) to find the total current:")
print(f"I_total = V2 / Req_total")
# The final equation with all numbers:
print(f"I_total = {V2} V / (({R3} * ({R7} + {R8})) / ({R3} + {R7} + {R8})) Ω")
print(f"I_total = {V2} V / {Req_total:.3f} Ω = {I_total:.4f} A")
# --- End of the student's code ---

# Restore the original stdout
sys.stdout = original_stdout
# Get the captured output
output = captured_output.getvalue()

# Print the output to the user
print(output)
final_current_val = V2 / ((R3 * (R7 + R8)) / (R3 + R7 + R8))
print(f'<<<{final_current_val:.4f}>>>')