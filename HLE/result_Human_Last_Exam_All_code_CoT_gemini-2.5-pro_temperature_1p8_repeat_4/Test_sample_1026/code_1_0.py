import sys
import io

# Capture the original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = buffer = io.StringIO()

def solve_displacement():
    """
    Calculates the horizontal displacement of the guide.

    The calculation is based on the principle of conservation of the center of mass.
    Since there are no external horizontal forces on the system (mass + guide),
    the horizontal position of the center of mass remains constant.

    Let ΔX_M be the displacement of the guide and Δx_m be the displacement of the mass.
    m * Δx_m + M * ΔX_M = 0

    The displacement of the mass is the sum of the guide's displacement and the
    mass's displacement relative to the guide.
    The relative horizontal displacement of the mass is R + d.
    Δx_m = ΔX_M + (R + d)

    Substituting this into the conservation equation:
    m * (ΔX_M + R + d) + M * ΔX_M = 0
    (m + M) * ΔX_M = -m * (R + d)
    ΔX_M = - (m * (R + d)) / (m + M)
    """

    # Given values
    m = 0.20  # kg, mass of the body
    M = 0.80  # kg, mass of the guide
    R_cm = 20   # cm, radius of the circular arcs
    d_cm = 50   # cm, length of the straight section

    # Convert units to SI (meters)
    R = R_cm / 100.0
    d = d_cm / 100.0

    # Calculate the horizontal displacement of the guide
    displacement_M = - (m * (R + d)) / (m + M)

    # --- Outputting the results as requested ---
    # Print the equation with variables and then with numbers
    print("The formula for the horizontal displacement of the guide (ΔX_M) is:")
    print("ΔX_M = - (m * (R + d)) / (m + M)")
    print("\nPlugging in the values:")
    print(f"ΔX_M = - ({m} kg * ({R} m + {d} m)) / ({m} kg + {M} kg)")
    print(f"ΔX_M = - ({m} * ({R + d})) / ({m + M})")
    print(f"ΔX_M = - ({m * (R + d):.2f}) / ({m + M})")
    print(f"\nThe calculated horizontal displacement of the guide is {displacement_M:.2f} meters.")
    print(f"This is equivalent to {displacement_M * 100:.0f} cm.")

solve_displacement()

# Get the content of the buffer
output = buffer.getvalue()
# Restore the original stdout
sys.stdout = original_stdout

# Print the captured output
print(output)

# Extract the final numerical answer for the platform
final_answer_value = - (0.20 * (0.20 + 0.50)) / (0.20 + 0.80)
print(f'<<<{final_answer_value:.2f}>>>')