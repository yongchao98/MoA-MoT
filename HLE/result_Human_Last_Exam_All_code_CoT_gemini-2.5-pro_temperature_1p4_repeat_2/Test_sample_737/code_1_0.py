import math

def solve_pavement_thickness():
    """
    Solves for the required pavement thickness using Burmister's two-layer theory.
    """
    # --- Given Data (in SI units) ---
    # Test 1: Plate on Subgrade
    w_s_test_um = 2460  # deflection in micrometers
    w_s_test = w_s_test_um / 1_000_000  # deflection in meters
    P_test = 30_000  # Load in Newtons (30 kN)
    d_plate = 0.305  # Plate diameter in meters (305 mm)
    a_plate = d_plate / 2  # Plate radius in meters

    # Test 2: Plate on Pavement + Subgrade
    w_p_test_um = 1080  # deflection in micrometers
    w_p_test = w_p_test_um / 1_000_000  # deflection in meters
    h_trial = 0.300  # Pavement thickness in meters (300 mm)

    # Design Parameters
    max_deflection_mm = 1.00
    w_design = max_deflection_mm / 1000  # Max allowable deflection in meters
    load_mass = 1800  # Load in kg (1.8 ton)
    gravity = 9.81 # m/s^2
    P_design = load_mass * gravity # Design load in Newtons
    p_design = 600_000  # Tyre pressure in Pa (600 kN/m2)
    
    # Material property
    mu = 0.5 # Poisson's ratio for all materials

    print("### Step 1: Calculate Subgrade Elastic Modulus (Es) ###")
    # For a flexible plate on an elastic half-space with mu=0.5, w = 1.5 * P / (pi * a * E)
    # Solving for E gives: E = 1.5 * P / (pi * a * w)
    E_s = (1.5 * P_test) / (math.pi * a_plate * w_s_test)
    print(f"The formula for subgrade modulus is: Es = 1.5 * P / (π * a * w_s)")
    print(f"Es = (1.5 * {P_test} N) / (π * {a_plate:.4f} m * {w_s_test} m)")
    print(f"Es = {E_s:,.1f} Pa or {E_s/1e6:.2f} MPa\n")

    print("### Step 2: Determine Pavement to Subgrade Modulus Ratio (Ep/Es) ###")
    # Deflection factor Fw is the ratio of two-layer deflection to one-layer deflection
    F_w_trial = w_p_test / w_s_test
    print(f"First, find the deflection factor F_w from the trial section:")
    print(f"F_w = w_p / w_s = {w_p_test} m / {w_s_test} m = {F_w_trial:.4f}")

    # Use analytical approximation for Burmister chart: Fw = 1 / [1 + 2.5*(h/a)*(Es/Ep)^(2/3)]
    # Solve for (Es/Ep): (Es/Ep)^(2/3) = ((1/Fw) - 1) / (2.5 * h/a)
    h_a_trial_ratio = h_trial / a_plate
    print(f"The thickness-to-radius ratio for the trial is h/a = {h_trial} m / {a_plate:.4f} m = {h_a_trial_ratio:.4f}")
    
    Es_Ep_ratio_pow_2_3 = ((1 / F_w_trial) - 1) / (2.5 * h_a_trial_ratio)
    Es_Ep_ratio = Es_Ep_ratio_pow_2_3**(3/2)
    Ep_Es_ratio = 1 / Es_Ep_ratio
    
    print("Using an analytical approximation, we find the modular ratio Ep/Es.")
    print(f"Ep/Es = 1 / [(((1 / {F_w_trial:.4f}) - 1) / (2.5 * {h_a_trial_ratio:.4f}))^(3/2)]")
    print(f"Ep/Es = {Ep_Es_ratio:.2f}\n")

    print("### Step 3: Calculate Design Wheel Load Parameters ###")
    # Area = P_design / p_design = pi * a_design^2
    # a_design = sqrt(P_design / (pi * p_design))
    a_design = math.sqrt(P_design / (math.pi * p_design))
    print(f"The design contact radius 'a' is found from tire pressure and load:")
    print(f"a_design = sqrt(P_design / (π * p_design))")
    print(f"a_design = sqrt({P_design:.1f} N / (π * {p_design} Pa)) = {a_design:.4f} m\n")

    print("### Step 4: Determine Required Pavement Thickness (h) ###")
    # First, find the required deflection factor Fw for the design case
    # w_design = (1.5 * p_design * a_design / Es) * F_w_design
    # F_w_design = (w_design * Es) / (1.5 * p_design * a_design)
    F_w_design = (w_design * E_s) / (1.5 * p_design * a_design)
    print("Next, find the required deflection factor F_w_design to meet the 1.00 mm limit.")
    print(f"F_w_design = (w_design * Es) / (1.5 * p_design * a_design)")
    print(f"F_w_design = ({w_design} m * {E_s:,.1f} Pa) / (1.5 * {p_design} Pa * {a_design:.4f} m)")
    print(f"F_w_design = {F_w_design:.4f}")

    # Now solve the approximation formula for the required h/a ratio
    # (h/a)_req = ((1 / F_w_design) - 1) / (2.5 * (Es/Ep)^(2/3))
    h_a_req_ratio = ((1 / F_w_design) - 1) / (2.5 * Es_Ep_ratio_pow_2_3)
    
    # Finally, calculate the required thickness h
    h_design = h_a_req_ratio * a_design
    h_design_mm = h_design * 1000
    
    print("\nFinally, calculate the required thickness 'h' using the required F_w and the known modular ratio.")
    print("Required h/a = [ (1 / F_w_design) - 1 ] / [ 2.5 * (Es/Ep)^(2/3) ]")
    print(f"Required h/a = [ (1 / {F_w_design:.4f}) - 1 ] / [ 2.5 * {Es_Ep_ratio_pow_2_3:.4f} ] = {h_a_req_ratio:.4f}")
    
    print("\nh = (Required h/a) * a_design")
    print(f"h = {h_a_req_ratio:.4f} * {a_design:.4f} m = {h_design:.4f} m")

    print(f"\nThe required pavement thickness is {h_design_mm:.2f} mm.")
    
    return h_design_mm

if __name__ == '__main__':
    final_thickness = solve_pavement_thickness()
    print(f'<<<{final_thickness:.2f}>>>')
