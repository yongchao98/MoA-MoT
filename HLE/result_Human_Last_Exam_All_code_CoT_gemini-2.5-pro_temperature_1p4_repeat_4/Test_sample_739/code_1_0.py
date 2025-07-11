import math

def calculate_water_film_thickness():
    """
    Calculates the design water film thickness on a roadway based on HEC-22 methodology.
    """
    # --- Parameters and Assumptions ---

    # Number of lanes
    num_lanes = 3
    # Lane width in meters
    lane_width_m = 3.6
    # Cross-slope (1.75% as a decimal)
    S = 0.0175
    
    # Manning's roughness coefficient 'n'
    # For "rough-textured asphalt pavement", a value typical for a chip-seal coat is assumed.
    n = 0.022
    
    # Design rainfall intensity in mm/hr
    # Since Rainfall-Duration-Frequency curves are not provided, a standard design intensity
    # for hydroplaning analysis is assumed. 100 mm/hr represents a heavy storm.
    i = 100
    
    # Pavement Mean Texture Depth (T_texture) in mm
    # For a rough texture consistent with n=0.022, a texture depth of 2.0 mm is assumed.
    T_texture = 2.0
    
    # Unit conversion constant for the HEC-22 formula (SI units)
    k = 0.0093

    # --- Calculations ---

    # 1. Calculate the flow path length (L) in meters
    L = num_lanes * lane_width_m

    # 2. Calculate the flow depth over the asperities (T_flow) using the Kinematic Wave equation
    # T_flow = k * (n * L)^0.6 * i^0.4 * S^-0.3
    term_nl = n * L
    term_nl_pow = math.pow(term_nl, 0.6)
    term_i_pow = math.pow(i, 0.4)
    term_s_pow = math.pow(S, -0.3)
    
    T_flow = k * term_nl_pow * term_i_pow * term_s_pow

    # 3. Calculate the total design water film thickness (Td)
    Td = T_flow + T_texture

    # --- Output Results ---
    
    print("--- Design Water Film Thickness Calculation ---")
    print("\nStep 1: Define Input Parameters and Assumptions")
    print(f"  - Number of lanes: {num_lanes}")
    print(f"  - Lane width: {lane_width_m} m")
    print(f"  - Cross-slope (S): {S*100:.2f}% or {S}")
    print(f"  - Manning's roughness (n): {n} (assumed for rough-textured asphalt)")
    print(f"  - Design rainfall intensity (i): {i} mm/hr (assumed)")
    print(f"  - Pavement texture depth (T_texture): {T_texture} mm (assumed)")

    print("\nStep 2: Calculate Flow Depth over Texture (T_flow)")
    print("Formula: T_flow = k * (n * L)^0.6 * i^0.4 * S^-0.3")
    print(f"  - Flow path length (L) = {num_lanes} * {lane_width_m} = {L:.1f} m")
    print("Substituting the values into the formula:")
    print(f"  T_flow = {k} * ({n} * {L:.1f})^0.6 * ({i})^0.4 * ({S})^-0.3")
    print(f"  T_flow = {k} * ({term_nl:.4f})^0.6 * {term_i_pow:.4f} * {term_s_pow:.4f}")
    print(f"  T_flow = {k} * {term_nl_pow:.4f} * {term_i_pow:.4f} * {term_s_pow:.4f}")
    print(f"  T_flow = {T_flow:.4f} mm")
    
    print("\nStep 3: Calculate Total Design Water Film Thickness (Td)")
    print("Formula: Td = T_flow + T_texture")
    print("Substituting the calculated and assumed values:")
    print(f"  Td = {T_flow:.4f} mm + {T_texture:.1f} mm")
    print(f"  Td = {Td:.2f} mm")

    print("\n----------------------------------------------------")
    print(f"The final design water film thickness is {Td:.2f} mm.")
    print("----------------------------------------------------")
    
    return Td

# Run the calculation and store the final answer
final_answer = calculate_water_film_thickness()
# The print statements within the function already provide the detailed output.
# The following line is for the final answer block as requested.
# print(f"\n<<<{final_answer:.2f}>>>")
