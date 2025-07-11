import math

def simulate_sensor_response():
    """
    This script simulates the effect of the biasing magnetic field on a Spin Valve sensor's sensitivity.
    We model the sensor's signal response to a magnetic field H using the hyperbolic tangent function,
    Signal(H) = tanh(H), which captures the saturation effect at high field values.
    """

    # --- Define Parameters ---
    # The small magnetic field from a single paramagnetic particle
    H_particle = 0.1  # Arbitrary units

    # Case 1: Optimal bias field, placing the sensor in its sensitive, linear operating range
    H_bias_optimal = 0.5 # Arbitrary units

    # Case 2: Improperly strong bias field, pushing the sensor into saturation
    H_bias_improper = 3.0 # Arbitrary units

    print("--- Simulating Spin Valve Sensor Response ---\n")

    # --- Case 1: Optimal Bias ---
    print("Case 1: Optimal Biasing Field")
    signal_optimal_base = math.tanh(H_bias_optimal)
    signal_optimal_with_particle = math.tanh(H_bias_optimal + H_particle)
    detected_signal_optimal = signal_optimal_with_particle - signal_optimal_base

    # Output the numbers used in the final equation as requested
    print(f"The equation for the detected signal is: tanh({H_bias_optimal} + {H_particle}) - tanh({H_bias_optimal})")
    print(f"Calculated result: {signal_optimal_with_particle:.4f} - {signal_optimal_base:.4f} = {detected_signal_optimal:.4f}")
    print(f"Detected signal with optimal bias: {detected_signal_optimal:.4f}\n")


    # --- Case 2: Improper (Saturating) Bias ---
    print("Case 2: Improper (Saturating) Biasing Field")
    signal_improper_base = math.tanh(H_bias_improper)
    signal_improper_with_particle = math.tanh(H_bias_improper + H_particle)
    detected_signal_improper = signal_improper_with_particle - signal_improper_base

    # Output the numbers used in the final equation as requested
    print(f"The equation for the detected signal is: tanh({H_bias_improper} + {H_particle}) - tanh({H_bias_improper})")
    print(f"Calculated result: {signal_improper_with_particle:.4f} - {signal_improper_base:.4f} = {detected_signal_improper:.4f}")
    print(f"Detected signal with improper (saturating) bias: {detected_signal_improper:.4f}\n")


    # --- Conclusion ---
    sensitivity_ratio = detected_signal_optimal / detected_signal_improper if detected_signal_improper != 0 else float('inf')
    print("--- Conclusion ---")
    print(f"The sensor's sensitivity in the optimal case is {sensitivity_ratio:.1f} times higher than in the saturated case.")
    print("This demonstrates that an improperly strong biasing magnet causes Spin Valve Saturation, severely reducing its ability to detect particles.")

if __name__ == "__main__":
    simulate_sensor_response()