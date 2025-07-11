import sys
import io

# Capture the original stdout
original_stdout = sys.stdout
# Create a string buffer
output_buffer = io.StringIO()
# Redirect stdout to the buffer
sys.stdout = output_buffer

# --- Main Logic ---

# Define symbols using strings for clear representation
gamma = "γ"  # Surface tension
r = "r"      # Radial position

# The derivation is based on the linearized Young-Laplace equation in cylindrical coordinates.
# The pressure from surface tension is P_surface = γ * κ, where κ is the mean curvature.
# For a linearized axisymmetric interface ξ(r), the curvature κ ≈ d²ξ/dr² + (1/r) * dξ/dr.
# The governing equation balancing surface tension against other pressures (C_term) is:
# γ * (d²ξ/dr² + (1/r) * dξ/dr) + C_term = 0
# Rearranging this gives:
# γ * d²ξ/dr² + (γ/r) * dξ/dr + C_term = 0
# Comparing this to the given form A(r)*d²ξ/dr² + B(r)*dξ/dr + C(r, ξ(r)) = 0,
# we can identify the coefficients A(r) and B(r).

# The coefficients are:
A_r = gamma
B_r = f"{gamma} / {r}"

# The C term encapsulates forces not related to the basic surface tension geometry,
# like the electrostatic pressure mentioned in the problem.
C_term = "C(r, ξ(r))"

# The final form of the equation with identified coefficients
final_equation = f"{A_r} * d²ξ(r)/dr² + ({B_r}) * dξ(r)/dr + {C_term} = 0"

print("The governing linear equation for the interfacial shape ξ(r) is derived from the balance of surface tension and electrostatic pressure.")
print("\nBy linearizing the Young-Laplace equation for an axisymmetric surface in cylindrical coordinates, we obtain the form:")
print(f"  {final_equation}")

print("\nFrom this equation, we can identify the coefficients A(r) and B(r):")

print("\n--- Derived Coefficients ---")
print(f"A(r) = {A_r}")
print(f"B(r) = {B_r}")

# --- End of Main Logic ---

# Get the content from the buffer
final_output = output_buffer.getvalue()
# Restore the original stdout
sys.stdout = original_stdout
# Print the captured output
print(final_output)

# Append the final answer in the specified format
final_answer = "A(r) = γ, B(r) = γ/r"
print(f"<<<{final_answer}>>>")