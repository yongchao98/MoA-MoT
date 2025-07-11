import math

def solve_moment_capacity():
    """
    Calculates the moment at collapse for the given non-standard section.
    """

    # 1. Define Material Properties and Geometry (from problem statement and Eurocode 2)
    fck = 30.0  # MPa
    fyk = 500.0  # MPa
    
    # Partial safety factors (Eurocode 2)
    gamma_c = 1.5
    gamma_s = 1.15
    
    # Concrete stress block parameters (Eurocode 2 for fck <= 50 MPa)
    alpha_cc = 0.85
    lambda_ = 0.8
    
    # Material properties
    E_s = 200000.0  # MPa (Young's Modulus of steel)
    eps_cu3 = 0.0035  # Ultimate concrete strain
    
    # Section Geometry (in mm)
    b_top = 100.0
    b_bottom = 400.0
    h_trap = 300.0
    h_rect = 100.0
    h_total = h_trap + h_rect

    # Reinforcement details
    phi = 20.0  # mm (diameter of bars)
    # Distances from the top compression fiber
    d_comp = 50.0  # mm (top layer)
    d_tens1 = 50.0 + 210.0  # mm (middle layer)
    d_tens2 = 50.0 + 210.0 + 90.0  # mm (bottom layer)
    
    n_comp = 2  # Number of bars in compression layer
    n_tens1 = 2 # Number of bars in middle tension layer
    n_tens2 = 3 # Number of bars in bottom tension layer

    # 2. Calculate Design Strengths and Areas
    fcd = alpha_cc * fck / gamma_c
    fyd = fyk / gamma_s
    eps_yd = fyd / E_s

    A_s_bar = math.pi * (phi ** 2) / 4
    A_s_comp = n_comp * A_s_bar
    A_s_tens1 = n_tens1 * A_s_bar
    A_s_tens2 = n_tens2 * A_s_bar

    print("--- Design Parameters ---")
    print(f"fcd (Design concrete strength): {fcd:.2f} MPa")
    print(f"fyd (Design steel yield strength): {fyd:.2f} MPa")
    print(f"Steel yield strain: {eps_yd:.5f}")
    print(f"Area of one H20 bar: {A_s_bar:.2f} mm^2")
    print(f"Area of compression steel (As'): {A_s_comp:.2f} mm^2 at d'={d_comp} mm")
    print(f"Area of tension steel 1 (As1): {A_s_tens1:.2f} mm^2 at d1={d_tens1} mm")
    print(f"Area of tension steel 2 (As2): {A_s_tens2:.2f} mm^2 at d2={d_tens2} mm")
    print("-" * 25 + "\n")

    # Helper function for section width b(y) at depth y from top
    def get_width(y):
        if 0 <= y <= h_trap:
            return b_top + (b_bottom - b_top) / h_trap * y
        elif h_trap < y <= h_total:
            return b_bottom
        return 0

    # 3. Determine Neutral Axis Depth (x) using Bisection Method
    def get_net_force(x):
        if x <= 0: return -1e9 # Physically impossible, return large negative
        
        # --- Tensile Forces ---
        # Layer 1
        eps_st1 = eps_cu3 * (d_tens1 - x) / x
        sigma_st1 = max(-fyd, min(fyd, eps_st1 * E_s))
        F_st1 = A_s_tens1 * sigma_st1
        # Layer 2
        eps_st2 = eps_cu3 * (d_tens2 - x) / x
        sigma_st2 = max(-fyd, min(fyd, eps_st2 * E_s))
        F_st2 = A_s_tens2 * sigma_st2
        
        F_tens_total = F_st1 + F_st2
        
        # --- Compressive Forces ---
        # Concrete
        s = lambda_ * x  # Depth of rectangular stress block
        
        if s <= h_trap: # Compression zone is a trapezoid
            area_c_gross = (get_width(0) + get_width(s)) / 2 * s
        else: # Compression zone is trapezoid + rectangle
            area_c_gross_trap = (get_width(0) + get_width(h_trap)) / 2 * h_trap
            area_c_gross_rect = get_width(h_trap) * (s - h_trap)
            area_c_gross = area_c_gross_trap + area_c_gross_rect
            
        # Compression Steel
        eps_sc = eps_cu3 * (x - d_comp) / x if x > 0 else 0
        
        # Account for concrete area displaced by compression steel
        area_c_net = area_c_gross
        F_sc = 0
        if eps_sc > 0: # Steel is in compression
            sigma_sc = min(fyd, eps_sc * E_s)
            F_sc = A_s_comp * sigma_sc
            area_c_net -= A_s_comp # Subtract displaced concrete
        else: # Steel is in tension (should be included in tension force calc)
             sigma_sc = max(-fyd, eps_sc * E_s)
             F_sc_as_tension = A_s_comp * sigma_sc
             F_tens_total += F_sc_as_tension

        F_c = area_c_net * fcd
        F_comp_total = F_c + (F_sc if eps_sc > 0 else 0)

        return F_comp_total - F_tens_total

    # Bisection solver
    low, high = 0.01, h_total
    for _ in range(100): # 100 iterations for high precision
        mid = (low + high) / 2
        if get_net_force(mid) > 0:
            high = mid
        else:
            low = mid
    x = (low + high) / 2

    print("--- Force Equilibrium Analysis ---")
    print(f"Calculated Neutral Axis Depth (x): {x:.2f} mm")
    
    # 4. Verify Strains and Calculate Final Forces
    s = lambda_ * x
    eps_sc = eps_cu3 * (x - d_comp) / x
    sigma_sc = min(fyd, eps_sc * E_s) if eps_sc > 0 else 0
    F_sc = A_s_comp * sigma_sc if eps_sc > 0 else 0

    eps_st1 = eps_cu3 * (d_tens1 - x) / x
    sigma_st1 = min(fyd, eps_st1 * E_s)
    F_st1 = A_s_tens1 * sigma_st1
    
    eps_st2 = eps_cu3 * (d_tens2 - x) / x
    sigma_st2 = min(fyd, eps_st2 * E_s)
    F_st2 = A_s_tens2 * sigma_st2

    print(f"\nVerifying strains at x = {x:.2f} mm:")
    print(f"  - Compression Steel (d'=50mm): strain = {eps_sc:.5f} -> {'Yielded' if eps_sc >= eps_yd else 'Not Yielded'}")
    print(f"  - Tension Steel 1 (d1=260mm): strain = {eps_st1:.5f} -> {'Yielded' if eps_st1 >= eps_yd else 'Not Yielded'}")
    print(f"  - Tension Steel 2 (d2=350mm): strain = {eps_st2:.5f} -> {'Yielded' if eps_st2 >= eps_yd else 'Not Yielded'}")

    # Concrete Force and its Centroid (y_c)
    if s <= h_trap:
        area_c_gross = (get_width(0) + get_width(s)) / 2 * s
        # Centroid of trapezoid from top
        y_c = (s / 3) * (get_width(0) + 2 * get_width(s)) / (get_width(0) + get_width(s))
    else:
        A1 = (get_width(0) + get_width(h_trap)) / 2 * h_trap
        y1 = (h_trap / 3) * (get_width(0) + 2 * get_width(h_trap)) / (get_width(0) + get_width(h_trap))
        A2 = get_width(h_trap) * (s - h_trap)
        y2 = h_trap + (s - h_trap) / 2
        area_c_gross = A1 + A2
        y_c = (A1 * y1 + A2 * y2) / (A1 + A2)

    area_c_net = area_c_gross - (A_s_comp if eps_sc > 0 else 0)
    F_c = area_c_net * fcd
    
    print("\n--- Final Internal Forces (in kN) ---")
    print(f"Concrete Compressive Force (Fc): {F_c/1000:.2f} kN at y_c={y_c:.2f} mm")
    print(f"Compression Steel Force (Fsc): {F_sc/1000:.2f} kN")
    print(f"Tension Steel Force 1 (Fs1): {F_st1/1000:.2f} kN")
    print(f"Tension Steel Force 2 (Fs2): {F_st2/1000:.2f} kN")
    print("-" * 25 + "\n")

    # 5. Calculate Moment Capacity (M_rd)
    # Taking moments about the top compression fiber
    M_rd = (F_st1 * d_tens1 + F_st2 * d_tens2) - (F_sc * d_comp) - (F_c * y_c)
    M_rd_kNm = M_rd / 1e6

    print("--- Moment Capacity Calculation ---")
    print("M_rd = (Fs1 * d1 + Fs2 * d2) - (Fsc * d') - (Fc * y_c)")
    print(f"M_rd = ({F_st1/1000:.2f}kN * {d_tens1}mm + {F_st2/1000:.2f}kN * {d_tens2}mm) - ({F_sc/1000:.2f}kN * {d_comp}mm) - ({F_c/1000:.2f}kN * {y_c:.2f}mm)")
    
    moment_tensile = (F_st1 * d_tens1 + F_st2 * d_tens2) / 1e6
    moment_comp_steel = (F_sc * d_comp) / 1e6
    moment_comp_concrete = (F_c * y_c) / 1e6
    
    print(f"M_rd = ({moment_tensile:.2f}) - ({moment_comp_steel:.2f}) - ({moment_comp_concrete:.2f}) kNm")
    print(f"M_rd = {M_rd_kNm:.2f} kNm")

    return M_rd_kNm

# Run the calculation and print the final answer
moment_capacity = solve_moment_capacity()
print(f"\nThe moment at collapse for the section is {moment_capacity:.2f} kNm.")
print(f"<<<{moment_capacity:.1f}>>>")
