import math

def calculate_water_film_thickness():
    """
    Calculates the design water film thickness on a pavement surface
    using the Kinematic Wave Equation.
    """
    # --- Parameters ---
    # Number of lanes in one direction
    num_lanes = 3
    # Width of each lane in meters
    lane_width_m = 3.6
    # Cross-slope of the pavement (1.75% as a decimal)
    S = 0.0175
    # Manning's roughness coefficient for rough-textured asphalt
    n = 0.016
    # Design rainfall intensity in mm/hr (Assumed value for a 10-year storm)
    i = 150.0
    # Constant for SI units
    K = 0.0434

    # --- Calculations ---
    # L: Length of the flow path in meters
    L = num_lanes * lane_width_m

    # Calculate the water film thickness using the formula:
    # d = K * (n * L)^0.6 * i^0.4 / S^0.3
    term1 = (n * L) ** 0.6
    term2 = i ** 0.4
    term3 = S ** 0.3
    
    d = K * term1 * term2 / term3

    # --- Output Results ---
    print("Step 1: Define the parameters for the calculation.")
    print(f"  - Flow Path Length (L) = {num_lanes} lanes * {lane_width_m} m/lane = {L:.1f} m")
    print(f"  - Pavement Cross-Slope (S) = {S*100:.2f}% = {S}")
    print(f"  - Manning's Roughness (n) = {n} (for rough-textured asphalt)")
    print(f"  - Assumed Rainfall Intensity (i) = {i:.1f} mm/hr")
    print(f"  - Formula Constant (K) for SI units = {K}\n")
    
    print("Step 2: Apply the Kinematic Wave formula for water film thickness (d).")
    print("Formula: d = K * (n * L)^0.6 * i^0.4 / S^0.3\n")

    print("Step 3: Substitute the values into the formula.")
    print(f"d = {K} * ({n} * {L:.1f})^0.6 * {i:.1f}^0.4 / {S}^0.3")
    print(f"d = {K} * ({n * L:.4f})^0.6 * {i:.1f}^0.4 / {S}^0.3")
    print(f"d = {K} * {term1:.4f} * {term2:.4f} / {term3:.4f}")
    print(f"d = {(K * term1 * term2):.4f} / {term3:.4f}\n")

    print("Step 4: Final calculated design water film thickness.")
    print(f"d = {d:.2f} mm")
    
    # Return the final numerical answer in the required format
    print(f"\n<<<{d:.2f}>>>")

calculate_water_film_thickness()