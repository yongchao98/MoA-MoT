import math

def analyze_magnet_position(improper_distance_mm):
    """
    Analyzes the effect of an improperly positioned magnet on a Spin Valve sensor.
    
    This model assumes:
    1. The magnetic field strength is inversely proportional to the square of the distance.
    2. An optimal position and corresponding field are known.
    3. The Spin Valve sensor has a known saturation field threshold.
    """

    # --- System Parameters (Hypothetical) ---

    # At an optimal distance (e.g., 2.5 mm), the field is ideal for the sensor (e.g., 6 mT).
    optimal_distance_mm = 2.5
    optimal_field_mT = 6.0
    
    # The Spin Valve sensor saturates if the bias field is too strong (e.g., > 12 mT).
    sv_saturation_field_mT = 12.0
    
    # Model the field: B = C / distance^2. We find the constant C from the optimal setup.
    # C = B * distance^2
    magnet_strength_constant = optimal_field_mT * (optimal_distance_mm ** 2)
    
    # --- Scenario Analysis ---
    
    # Calculate the magnetic field at the sensor for the improper position.
    # An improper position could be moving the magnet too close.
    field_at_improper_position_mT = magnet_strength_constant / (improper_distance_mm ** 2)
    
    # --- Output Results ---
    
    print("--- Sensor System Analysis ---")
    print(f"Optimal Magnet Distance: {optimal_distance_mm} mm")
    print(f"Optimal Field at Sensor: {optimal_field_mT} mT")
    print(f"Sensor Saturation Threshold: {sv_saturation_field_mT} mT")
    print("-" * 30)
    
    print("--- Testing Improper Position ---")
    print(f"New Magnet Distance: {improper_distance_mm} mm")
    print(f"Resulting Field Equation: B = {magnet_strength_constant:.2f} / ({improper_distance_mm})^2")
    print(f"Calculated Field at Sensor: {field_at_improper_position_mT:.2f} mT")
    
    if field_at_improper_position_mT > sv_saturation_field_mT:
        print("\nConclusion: CRITICAL FAILURE")
        print("The magnetic field from the improperly positioned magnet exceeds the sensor's threshold.")
        print("This leads to Spin Valve Saturation, making particle detection impossible.")
    else:
        print("\nConclusion: OPERATIONAL (in this case)")
        print("The field is within the sensor's operating range.")

# Let's test a scenario where the magnet is moved too close, e.g., to 1.5 mm.
analyze_magnet_position(improper_distance_mm=1.5)
