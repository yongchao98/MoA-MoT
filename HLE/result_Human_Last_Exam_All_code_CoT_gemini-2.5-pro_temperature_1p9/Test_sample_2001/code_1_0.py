import math

def sensor_model(B_field):
    """
    A simple model for a sensor response that saturates.
    The output is modeled using the hyperbolic tangent (tanh) function.
    The sensor is most sensitive when B_field is near 0.
    It saturates (output approaches 1.0) for large B_field.
    """
    return math.tanh(B_field)

# --- Simulation Parameters ---
# Magnetic field from a single paramagnetic particle (arbitrary units)
B_particle = 0.1

# --- Scenario 1: Correct Magnet Position ---
print("--- Scenario 1: Correct Magnet Position ---")
# The bias field from the magnet is in the sensor's optimal working range
B_bias_correct = 0.5

# Sensor output with only the bias field
output_base_correct = sensor_model(B_bias_correct)
print(f"Correct Bias Field = {B_bias_correct}")
print(f"Sensor Base Output = tanh({B_bias_correct}) = {output_base_correct:.4f}")

# Total field when a particle is directly over the sensor
B_total_correct = B_bias_correct + B_particle
output_peak_correct = sensor_model(B_total_correct)
print(f"\nTotal Field with Particle = {B_bias_correct} + {B_particle} = {B_total_correct}")
print(f"Sensor Peak Output = tanh({B_total_correct}) = {output_peak_correct:.4f}")

# The resulting signal is the change in output
signal_correct = output_peak_correct - output_base_correct
print(f"\nDetected Signal = Peak Output - Base Output = {output_peak_correct:.4f} - {output_base_correct:.4f} = {signal_correct:.4f}")
print("Result: A clear signal is detected.")

# --- Scenario 2: Incorrect Magnet Position (Too Close) ---
print("\n" + "="*50 + "\n")
print("--- Scenario 2: Incorrect Magnet Position (Too Close) ---")
# The bias field is too strong, pushing the sensor into saturation
B_bias_incorrect = 3.0

# Sensor output with only the bias field
output_base_incorrect = sensor_model(B_bias_incorrect)
print(f"Incorrect (Saturating) Bias Field = {B_bias_incorrect}")
print(f"Sensor Base Output = tanh({B_bias_incorrect}) = {output_base_incorrect:.4f}")

# Total field when a particle is directly over the sensor
B_total_incorrect = B_bias_incorrect + B_particle
output_peak_incorrect = sensor_model(B_total_incorrect)
print(f"\nTotal Field with Particle = {B_bias_incorrect} + {B_particle} = {B_total_incorrect}")
print(f"Sensor Peak Output = tanh({B_total_incorrect}) = {output_peak_incorrect:.4f}")

# The resulting signal is the change in output
signal_incorrect = output_peak_incorrect - output_base_incorrect
print(f"\nDetected Signal = Peak Output - Base Output = {output_peak_incorrect:.4f} - {output_base_incorrect:.4f} = {signal_incorrect:.4f}")
print("Result: The signal is negligible. The particle is not detected.")
