import numpy as np

def analyze_sensor_saturation():
    """
    This function models and explains Spin Valve sensor saturation
    due to an improperly positioned magnet.
    """
    print("Analyzing the effect of magnet position on the Spin Valve sensor.")
    print("A sensor's response (e.g., resistance) to a magnetic field (H) is not linear.")
    print("It is most sensitive in a specific range and saturates at high fields.")
    print("We can model this with a sigmoid function.\n")

    # --- Model Parameters ---
    # Let's model the sensor's sensitivity (S) as the derivative of its response curve.
    # For a typical sensor, sensitivity is highest at a specific bias field and drops off.
    # We'll use a Gaussian-like function for sensitivity for simplicity.
    # S(H) = S_max * exp(-((H - H_optimal_bias)^2) / (2 * w^2))
    S_max = 5.0  # Maximum sensitivity (arbitrary units, e.g., V/T or Ohm/Oe)
    H_optimal_bias = 20.0  # The optimal bias field for peak sensitivity (in Oersted, Oe)
    w = 15.0  # A parameter defining the width of the sensitive range (in Oe)
    H_particle = 0.5  # The small magnetic field from a passing particle (in Oe)

    # --- Scenario 1: Correctly Positioned Magnet ---
    H_bias_1 = 20.0 # The bias field is set to the optimal value.
    sensitivity_1 = S_max * np.exp(-((H_bias_1 - H_optimal_bias)**2) / (2 * w**2))
    signal_1 = sensitivity_1 * H_particle

    print("--- Scenario 1: Correct Magnet Position ---")
    print(f"The magnet provides an optimal bias field of {H_bias_1:.1f} Oe.")
    # Final Equation 1: Sensitivity = S_max * exp(-((H_bias - H_optimal_bias)^2) / (2 * w^2))
    print(f"The sensor sensitivity is calculated as: {S_max:.1f} * exp(-(({H_bias_1:.1f} - {H_optimal_bias:.1f})^2) / (2 * {w:.1f}^2)) = {sensitivity_1:.4f}")
    # Final Equation 2: Signal = Sensitivity * H_particle
    print(f"The resulting signal from a particle is: {sensitivity_1:.4f} * {H_particle:.1f} = {signal_1:.4f} (arbitrary units)\n")


    # --- Scenario 2: Improperly Positioned Magnet (Too Close) ---
    H_bias_2 = 60.0 # The magnet is too close, creating a strong, non-optimal bias field.
    sensitivity_2 = S_max * np.exp(-((H_bias_2 - H_optimal_bias)**2) / (2 * w**2))
    signal_2 = sensitivity_2 * H_particle

    print("--- Scenario 2: Improper Magnet Position (Too Close) ---")
    print(f"The magnet provides an excessive bias field of {H_bias_2:.1f} Oe.")
    # Final Equation 3: Sensitivity = S_max * exp(-((H_bias - H_optimal_bias)^2) / (2 * w^2))
    print(f"The sensor sensitivity is calculated as: {S_max:.1f} * exp(-(({H_bias_2:.1f} - {H_optimal_bias:.1f})^2) / (2 * {w:.1f}^2)) = {sensitivity_2:.4f}")
    # Final Equation 4: Signal = Sensitivity * H_particle
    print(f"The resulting signal from a particle is: {sensitivity_2:.4f} * {H_particle:.1f} = {signal_2:.4f} (arbitrary units)\n")

    print("--- Conclusion ---")
    print("The excessive bias field from the improperly positioned magnet pushes the sensor into a low-sensitivity (saturation) region.")
    print("As a result, the signal from the particle becomes too weak to be detected reliably.")
    print("This phenomenon is called Spin Valve Saturation.")

analyze_sensor_saturation()
<<<D>>>