import math

def calculate_water_film_thickness():
    """
    Calculates the water film thickness on a roadway based on its geometry
    and assumed rainfall conditions.
    """
    # 1. Define given parameters and assumptions
    lane_width = 3.6  # meters
    num_lanes = 3
    cross_slope_percent = 1.75

    # Assumptions for unspecified values
    # Manning's n for rough-textured asphalt
    n_manning = 0.015
    # Assumed design rainfall intensity (for hydroplaning analysis)
    i_mm_per_hr = 150.0

    # 2. Calculate values for the formula
    # Total flow path length (L) in meters
    L = num_lanes * lane_width
    # Cross-slope (S) in m/m
    S = cross_slope_percent / 100.0
    # Convert rainfall intensity (i) to m/s
    i_m_per_s = i_mm_per_hr / (1000 * 3600)

    # 3. Calculate intermediate terms for the Kinematic Wave equation
    # q = flow per unit width (m^2/s)
    q = i_m_per_s * L
    # S^0.5
    s_sqrt = math.sqrt(S)

    # 4. Solve the Kinematic Wave Equation for Water Film Thickness (WFT) in meters
    # WFT_m = [ (q * n) / S^0.5 ]^0.6
    wft_m = ((q * n_manning) / s_sqrt) ** 0.6

    # 5. Convert final result to millimeters
    wft_mm = wft_m * 1000

    # 6. Print the results, showing the final equation with values
    print("--- Water Film Thickness Calculation ---")
    print("\nThe design water film thickness (WFT) is calculated using the Kinematic Wave Equation:")
    print("WFT = [ (i * L * n) / S^0.5 ]^0.6, where units must be SI (m, s).")

    print("\nParameters Used:")
    print(f"  - Rainfall Intensity (i): {i_mm_per_hr} mm/hr -> {i_m_per_s:.8f} m/s")
    print(f"  - Flow Path Length (L): {L:.1f} m")
    print(f"  - Manning's Roughness (n): {n_manning:.3f}")
    print(f"  - Cross-Slope (S): {S:.4f} m/m")

    print("\nFinal Equation with numerical values:")
    print(f"WFT (mm) = 1000 * [ (({i_m_per_s:.8f} * {L:.1f}) * {n_manning:.3f}) / {S:.4f}**0.5 ]**0.6")
    print(f"WFT (mm) = 1000 * [ ({q:.8f} * {n_manning:.3f}) / {s_sqrt:.4f} ]**0.6")
    print(f"WFT (mm) = 1000 * [ {q * n_manning:.8f} / {s_sqrt:.4f} ]**0.6")
    print(f"WFT (mm) = 1000 * [ {((q * n_manning) / s_sqrt):.8f} ]**0.6")
    print(f"WFT (mm) = 1000 * {wft_m:.6f}")
    print(f"\nDesign Water Film Thickness = {wft_mm:.2f} mm")


if __name__ == '__main__':
    calculate_water_film_thickness()