import sys
import io

# Redirect stdout to capture the output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def solve_river_crossing():
    """
    Calculates the downstream drift of a boat crossing a river with a parabolic flow profile.

    The boat travels from one bank (y=0) to 3/4 of the way across, then immediately
    turns around and returns to the starting bank.
    """

    # --- Problem Parameters ---
    # You can change these values to solve for different scenarios.
    L = 100.0   # Width of the river (e.g., in meters)
    v0 = 2.0    # Maximum flow velocity at the center of the river (e.g., in m/s)
    v = 4.0     # Boat's speed relative to the water (e.g., in m/s)

    # --- Derived Formula ---
    # The total downstream distance 'D' is given by the formula:
    # D = (9 * v0 * L) / (8 * v)
    # This formula is derived by integrating the flow velocity over the boat's travel path.
    # The drift on the way out is the same as the drift on the way back, so we have:
    # D = 2 * Drift_out = 2 * (9 * v0 * L / (16 * v)) = (9 * v0 * L) / (8 * v)

    # --- Calculation ---
    numerator = 9 * v0 * L
    denominator = 8 * v
    total_distance = numerator / denominator

    # --- Output Results ---
    print("--- River Crossing Problem ---")
    print(f"Given parameters:")
    print(f"  River Width (L) = {L}")
    print(f"  Max Flow Velocity (v_0) = {v0}")
    print(f"  Boat's Relative Speed (v) = {v}\n")

    print("The final equation for the total downstream distance (D) is:")
    print("  D = (9 * v_0 * L) / (8 * v)\n")

    print("Substituting the given values into the equation:")
    # Using .format to clearly show the numbers being substituted
    print(f"  D = (9 * {v0} * {L}) / (8 * {v})")
    print(f"  D = {numerator} / {denominator}")
    print(f"  D = {total_distance}\n")
    print(f"The distance between the boat's returning position and its original starting point is {total_distance}.")

# Execute the function to solve the problem
solve_river_crossing()

# Restore stdout and get the captured output
sys.stdout = old_stdout
output = captured_output.getvalue()

# Print the captured output to the actual console
print(output)