import math

def calculate_buckling_resistance():
    """
    Calculates the buckling resistance of the concrete encased composite column
    based on Eurocode 4 (EN 1994-1-1).
    """
    # --- 1. Material and Section Properties (from problem statement) ---
    print("--- Step 1: Material and Section Properties ---")
    
    # Column geometry
    H_col = 400.0  # mm, height of concrete section
    B_col = 400.0  # mm, width of concrete section
    cover = 30.0   # mm
    L = 4000.0     # mm, storey height

    # Structural Steel (UC 254 x 254 x 132, S355)
    f_y = 345.0     # N/mm^2, yield strength
    A_a = 168.0 * 100 # mm^2, area from cm^2
    I_ay = 22500.0 * 1e4 # mm^4, moment of inertia from cm^4 (strong axis)
    I_az = 7530.0 * 1e4  # mm^4, moment of inertia from cm^4 (weak axis)
    E_a = 210000.0  # N/mm^2, Modulus of Elasticity of steel

    # Concrete (C40/50)
    f_ck = 40.0     # N/mm^2, characteristic cylinder strength
    f_cm = f_ck + 8.0
    E_cm = 22000 * (f_cm / 10)**0.3 # N/mm^2, Secant Modulus of Elasticity

    # Reinforcing Steel (S500)
    d_rebar = 16.0   # mm, diameter
    n_rebars = 4
    f_sk = 500.0     # N/mm^2, characteristic yield strength
    E_s = 200000.0   # N/mm^2, Modulus of Elasticity of rebar

    # Partial safety factors
    gamma_a = 1.0
    gamma_c = 1.5
    gamma_s = 1.15
    
    print(f"Storey height (L): {L/1000} m")
    print(f"Concrete: C40/50 (f_ck = {f_ck} MPa)")
    print(f"Steel Section: UC 254x254x132, S355 (f_y = {f_y} MPa)")
    print(f"Reinforcement: 4xT16, S500 (f_sk = {f_sk} MPa)\n")


    # --- 2. Detailed Sectional Properties Calculation ---
    print("--- Step 2: Detailed Sectional Properties Calculation ---")
    
    # Reinforcement properties
    A_s1 = math.pi * (d_rebar / 2)**2
    A_s = n_rebars * A_s1
    rebar_dist_from_center = B_col/2 - cover - d_rebar/2
    I_s = 4 * A_s1 * rebar_dist_from_center**2 # Symmetrical for both axes

    # Concrete properties
    A_c_gross = B_col * H_col
    A_c = A_c_gross - A_a - A_s
    I_c = (B_col * H_col**3) / 12

    print(f"Area of steel section (A_a): {A_a:.2f} mm^2")
    print(f"Total area of reinforcement (A_s): {A_s:.2f} mm^2")
    print(f"Net area of concrete (A_c): {A_c:.2f} mm^2")
    print(f"Second moment of area of steel (weak axis, I_az): {I_az:.2e} mm^4")
    print(f"Second moment of area of rebar (I_s): {I_s:.2e} mm^4")
    print(f"Second moment of area of concrete (I_c): {I_c:.2e} mm^4\n")

    # --- 3. Plastic Resistance Calculation ---
    print("--- Step 3: Plastic Resistance Calculation (Squash Load) ---")
    
    # Characteristic plastic resistance (for slenderness calculation)
    N_pl_Rk = A_a * f_y + A_s * f_sk + A_c * f_ck
    
    # Design plastic resistance (for final capacity)
    N_pl_Rd = (A_a * f_y / gamma_a) + (A_s * f_sk / gamma_s) + (0.85 * A_c * f_ck / gamma_c)

    print(f"Characteristic Plastic Resistance (N_pl,Rk): {N_pl_Rk/1000:.1f} kN")
    print(f"Design Plastic Resistance (N_pl,Rd): {N_pl_Rd/1000:.1f} kN\n")


    # --- 4. Effective Flexural Stiffness Calculation (Weak Axis z-z) ---
    print("--- Step 4: Effective Flexural Stiffness (EI_eff) for Weak Axis ---")
    # Using the simplified method from EN1994-1-1 6.7.3.3(3)
    K_e = 0.9
    K_e_c = 0.6 # Ignoring creep doesn't change this simplified value
    EI_eff_z = K_e * (E_a * I_az + E_s * I_s + K_e_c * E_cm * I_c)
    
    print(f"Modulus of Elasticity of concrete (E_cm): {E_cm:.0f} MPa")
    print(f"Effective flexural stiffness (EI_eff): {EI_eff_z:.2e} N.mm^2\n")


    # --- 5. Critical Buckling Load ---
    print("--- Step 5: Elastic Critical Buckling Load (N_cr) ---")
    # Assuming pinned-pinned column, buckling length L_cr = L
    L_cr = L
    N_cr = (math.pi**2 * EI_eff_z) / L_cr**2

    print(f"Buckling Length (L_cr): {L_cr/1000:.1f} m")
    print(f"Elastic critical load (N_cr): {N_cr/1000:.1f} kN\n")


    # --- 6. Slenderness and Reduction Factor ---
    print("--- Step 6: Slenderness and Reduction Factor (chi) ---")
    lambda_bar = math.sqrt(N_pl_Rk / N_cr)

    # For encased I-sections, weak axis buckling corresponds to curve 'c' (alpha = 0.49)
    alpha = 0.49
    Phi = 0.5 * (1 + alpha * (lambda_bar - 0.2) + lambda_bar**2)
    chi = 1 / (Phi + math.sqrt(Phi**2 - lambda_bar**2))
    chi = min(chi, 1.0) # chi cannot be greater than 1

    print(f"Non-dimensional slenderness (lambda_bar): {lambda_bar:.3f}")
    print(f"Imperfection factor (alpha) for buckling curve 'c': {alpha}")
    print(f"Intermediate factor (Phi): {Phi:.3f}")
    print(f"Buckling reduction factor (chi): {chi:.3f}\n")


    # --- 7. Final Buckling Resistance ---
    print("--- Step 7: Final Buckling Resistance (N_b,Rd) ---")
    N_b_Rd = chi * N_pl_Rd
    
    print(f"The final buckling resistance is calculated as:")
    print(f"N_b,Rd = chi * N_pl,Rd")
    print(f"N_b,Rd = {chi:.3f} * {N_pl_Rd/1000:.1f} kN")
    print(f"N_b,Rd = {N_b_Rd/1000:.1f} kN")

    return N_b_Rd / 1000

if __name__ == '__main__':
    buckling_resistance = calculate_buckling_resistance()
    print("\n========================================================")
    print(f"The final calculated buckling resistance is {buckling_resistance:.1f} kN.")
    print("========================================================")
    # The final answer is wrapped in <<<>>> for the system.
    print(f'<<<{buckling_resistance:.1f}>>>')
