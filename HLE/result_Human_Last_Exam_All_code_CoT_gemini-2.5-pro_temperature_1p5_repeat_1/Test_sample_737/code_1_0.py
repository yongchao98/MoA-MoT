import math

def solve_pavement_thickness():
    """
    Calculates the required pavement thickness using Burmister's two-layer theory.
    """

    # --- Given Data (converted to SI units: m, N, Pa) ---
    # Plate Bearing Test
    P_plate = 30000.0  # N
    d_plate = 0.305  # m
    a_plate = d_plate / 2.0  # m
    w_subgrade = 2460e-6  # m (deflection on subgrade)
    
    # Trial Section Test
    h_trial = 0.300  # m
    w_trial = 1080e-6  # m (deflection on trial section)
    
    # Design Parameters
    W_design_mass = 1.80 * 1000  # kg (ton to kg)
    g = 9.81 # m/s^2
    P_design = W_design_mass * g # N
    p_design = 600 * 1000  # Pa (kN/m^2 to N/m^2)
    w_design = 1.00e-3  # m (allowable deflection)
    
    # Material Property
    poisson_ratio = 0.5  # μ

    print("--- Step 1: Determine Subgrade Modulus (E2) ---")
    
    # Pressure from plate load test
    p_plate = P_plate / (math.pi * a_plate**2)
    print(f"Plate radius (a_plate): {a_plate:.4f} m")
    print(f"Plate pressure (p_plate): {p_plate / 1e3:.2f} kPa")
    
    # Calculate E2 from subgrade test (w = 1.5 * p * a / E)
    E2 = (1.5 * p_plate * a_plate) / w_subgrade
    print(f"Subgrade deflection (w1): {w_subgrade * 1e6:.0f} μm")
    print(f"Calculated Subgrade Modulus (E2): {E2 / 1e6:.2f} MPa\n")
    
    print("--- Step 2: Determine Modular Ratio (E1/E2) ---")
    
    # Calculate F2 factor for the trial section
    # w_trial = (1.5 * p_plate * a_plate / E2) * F2_trial
    F2_trial = (w_trial * E2) / (1.5 * p_plate * a_plate)
    h_a_trial_ratio = h_trial / a_plate
    print(f"Trial pavement deflection (w2): {w_trial * 1e6:.0f} μm")
    print(f"Trial thickness/radius ratio (h/a): {h_a_trial_ratio:.3f}")
    print(f"Calculated trial deflection factor (F2_trial): {F2_trial:.4f}")

    # Burmister's F2 chart data for μ=0.5
    # F2_data = {E1/E2_ratio: {h/a_ratio: F2_value}}
    F2_data = {
        2:  {0.5: 0.85, 1.0: 0.78, 2.0: 0.7, 5.0: 0.6},
        5:  {0.5: 0.7, 1.0: 0.6, 2.0: 0.45, 5.0: 0.3},
        10: {0.5: 0.58, 1.0: 0.45, 2.0: 0.32, 5.0: 0.18},
        20: {0.5: 0.46, 1.0: 0.33, 2.0: 0.22, 5.0: 0.12},
    }

    # Interpolate to find E1/E2
    # 1. For each E1/E2 curve, find F2 at h/a = h_a_trial_ratio
    e_ratios = sorted(F2_data.keys())
    f2s_at_h_a = []
    for er in e_ratios:
        h_a_pts = sorted(F2_data[er].keys())
        f2_pts = [F2_data[er][h] for h in h_a_pts]
        # Interpolate F2 value for the trial h/a ratio
        f2_interp = float(_interpolate(h_a_pts, f2_pts, h_a_trial_ratio))
        f2s_at_h_a.append(f2_interp)
        
    # 2. Interpolate E1/E2 using the F2_trial value
    E1_E2_ratio = _interpolate(f2s_at_h_a, e_ratios, F2_trial)
    print(f"Interpolated Modular Ratio (E1/E2): {E1_E2_ratio:.2f}\n")

    print("--- Step 3: Determine Required Pavement Thickness (h_design) ---")
    
    # Calculate design load contact radius
    A_design = P_design / p_design
    a_design = math.sqrt(A_design / math.pi)
    print(f"Design load (P_design): {P_design:.2f} N")
    print(f"Design pressure (p_design): {p_design/1e3:.0f} kPa")
    print(f"Design contact radius (a_design): {a_design:.4f} m")
    
    # Calculate required F2 factor for the design
    # w_design = (1.5 * p_design * a_design / E2) * F2_design
    F2_design = (w_design * E2) / (1.5 * p_design * a_design)
    print(f"Target design deflection (w_design): {w_design*1000:.2f} mm")
    print(f"Required design deflection factor (F2_design): {F2_design:.4f}")
    
    # Interpolate to find h/a for the design
    # 1. Find bracketing E1/E2 curves
    e_ratios_all = sorted(F2_data.keys())
    idx_high = 0
    while idx_high < len(e_ratios_all) and e_ratios_all[idx_high] < E1_E2_ratio:
        idx_high += 1
    e_high = e_ratios_all[idx_high]
    e_low = e_ratios_all[idx_high - 1]

    # 2. For each bracketing curve, find h/a that gives F2_design
    h_a_for_e_low = _interpolate(list(F2_data[e_low].values()), list(F2_data[e_low].keys()), F2_design)
    h_a_for_e_high = _interpolate(list(F2_data[e_high].values()), list(F2_data[e_high].keys()), F2_design)

    # 3. Interpolate between the two h/a values to get the final h/a
    h_a_design_ratio = _interpolate([e_low, e_high], [h_a_for_e_low, h_a_for_e_high], E1_E2_ratio)
    print(f"Interpolated design thickness/radius ratio (h/a): {h_a_design_ratio:.3f}")
    
    # Final thickness calculation
    h_design = h_a_design_ratio * a_design
    h_design_mm = h_design * 1000
    
    print("\n--- Final Answer ---")
    print(f"Required Pavement Thickness (h_design) = {h_a_design_ratio:.3f} * {a_design * 1000:.1f} mm")
    print(f"Required Pavement Thickness = {h_design_mm:.1f} mm")
    return h_design_mm

def _interpolate(x_pts, y_pts, x_target):
    """Performs linear interpolation."""
    # Reverse if x points are descending
    if x_pts[0] > x_pts[-1]:
        x_pts.reverse()
        y_pts.reverse()
    
    if x_target <= x_pts[0]:
        return y_pts[0]
    if x_target >= x_pts[-1]:
        return y_pts[-1]

    for i in range(len(x_pts) - 1):
        if x_pts[i] <= x_target < x_pts[i+1]:
            x1, x2 = x_pts[i], x_pts[i+1]
            y1, y2 = y_pts[i], y_pts[i+1]
            return y1 + (y2 - y1) * (x_target - x1) / (x2 - x1)
    return y_pts[-1]


# --- Execute the Calculation ---
final_thickness = solve_pavement_thickness()
print(f"\n<<<The final required pavement thickness is {final_thickness:.1f} mm>>>")