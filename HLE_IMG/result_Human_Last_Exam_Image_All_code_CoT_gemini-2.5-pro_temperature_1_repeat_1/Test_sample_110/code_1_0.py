import math

def calculate_buckling_resistance():
    """
    Calculates the buckling resistance of the concrete encased composite column
    as per Eurocode 4 guidelines.
    """
    # --- 1. Define Input Parameters and Material Properties ---
    print("--- 1. Input Parameters and Material Properties ---")
    # Column and geometry
    L = 4000.0  # Storey height (mm)
    L_cr = L    # Effective buckling length (mm), assuming pinned-pinned ends (k=1.0)
    H_col = 400.0 # Column height (mm)
    B_col = 400.0 # Column width (mm)

    # Structural Steel (UC 254x254x132, S355)
    f_y = 345.0   # Yield strength (MPa)
    E_a = 210000.0# Modulus of Elasticity (MPa)
    A_a = 16800.0 # Area (mm^2), from 168 cm^2
    I_az = 7530.0 * 10**4 # Second moment of area about z-axis (weak axis) (mm^4)
    gamma_M1 = 1.0 # Partial factor for structural steel

    # Concrete (C40/50)
    f_ck = 40.0   # Characteristic cylinder strength (MPa)
    E_cm = 22 * ((f_ck + 8) / 10)**0.3 # Modulus of Elasticity (MPa) according to EN 1992-1-1
    gamma_c = 1.5 # Partial factor for concrete

    # Reinforcing Steel (4x T16, S500)
    f_sk = 500.0  # Characteristic yield strength (MPa)
    E_s = 200000.0# Modulus of Elasticity (MPa)
    d_rebar = 16.0# Rebar diameter (mm)
    n_rebar = 4   # Number of rebars
    cover = 30.0  # Concrete cover (mm)
    gamma_s = 1.15# Partial factor for reinforcing steel

    print(f"Column Length (L): {L} mm")
    print(f"Steel Yield Strength (f_y): {f_y} MPa")
    print(f"Concrete Strength (f_ck): {f_ck} MPa")
    print(f"Rebar Strength (f_sk): {f_sk} MPa")
    print(f"Concrete Modulus (E_cm): {E_cm:.2f} MPa\n")

    # --- 2. Calculate Sectional Properties ---
    print("--- 2. Sectional Properties ---")
    # Area of reinforcing steel
    A_s_one_bar = math.pi * (d_rebar / 2)**2
    A_s = n_rebar * A_s_one_bar
    print(f"Total Rebar Area (A_s): {A_s:.2f} mm^2")

    # Area of concrete (net)
    A_c_gross = H_col * B_col
    A_c = A_c_gross - A_a - A_s
    print(f"Net Concrete Area (A_c): {A_c:.2f} mm^2")

    # Second moment of area of concrete about z-axis
    I_cz = (H_col * B_col**3) / 12
    print(f"Concrete Section I_cz: {I_cz:.2e} mm^4")

    # Second moment of area of rebars about z-axis
    # Distance from centroid to rebar center in x-direction
    d_x_rebar = B_col / 2 - cover - d_rebar / 2
    I_sz = A_s * d_x_rebar**2 # I_local of bar is negligible
    print(f"Rebar Section I_sz: {I_sz:.2e} mm^4\n")

    # --- 3. Calculate Plastic Compressive Resistance (N_pl,Rd) ---
    print("--- 3. Plastic Compressive Resistance (N_pl,Rd) ---")
    # Design strengths
    f_yd = f_y / gamma_M1
    f_cd = 0.85 * f_ck / gamma_c
    f_sd = f_sk / gamma_s
    print(f"Design Strength of Steel (f_yd): {f_yd:.2f} MPa")
    print(f"Design Strength of Concrete (f_cd): {f_cd:.2f} MPa")
    print(f"Design Strength of Rebar (f_sd): {f_sd:.2f} MPa")

    # Resistance from each component
    N_pl_Rd_a = A_a * f_yd
    N_pl_Rd_c = A_c * f_cd
    N_pl_Rd_s = A_s * f_sd
    N_pl_Rd = N_pl_Rd_a + N_pl_Rd_c + N_pl_Rd_s
    print(f"Total Plastic Resistance (N_pl,Rd) = {N_pl_Rd_a/1000:.2f} kN + {N_pl_Rd_c/1000:.2f} kN + {N_pl_Rd_s/1000:.2f} kN = {N_pl_Rd/1000:.2f} kN\n")

    # --- 4. Calculate Critical Buckling Load (N_cr) ---
    print("--- 4. Critical Buckling Load (N_cr) ---")
    # Effective flexural stiffness about weak z-axis
    K_e = 0.6 # Correction factor for concrete stiffness
    EI_eff_z = E_a * I_az + K_e * E_cm * I_cz + E_s * I_sz
    print(f"Effective Flexural Stiffness (EI_eff_z): {EI_eff_z:.2e} N.mm^2")

    # Critical Euler buckling load
    N_cr = (math.pi**2 * EI_eff_z) / (L_cr**2)
    print(f"Critical Buckling Load (N_cr) = (pi^2 * {EI_eff_z:.2e}) / {L_cr:.0f}^2 = {N_cr/1000:.2f} kN\n")

    # --- 5. Calculate Relative Slenderness (lambda_bar) ---
    print("--- 5. Relative Slenderness (lambda_bar) ---")
    # Characteristic plastic resistance (all gamma factors = 1.0)
    N_pl_Rk = (A_a * f_y) + (A_c * f_ck) + (A_s * f_sk)
    print(f"Characteristic Plastic Resistance (N_pl,Rk): {N_pl_Rk/1000:.2f} kN")

    # Relative slenderness
    lambda_bar = math.sqrt(N_pl_Rk / N_cr)
    print(f"Relative Slenderness (lambda_bar) = sqrt({N_pl_Rk/1000:.2f} / {N_cr/1000:.2f}) = {lambda_bar:.3f}\n")

    # --- 6. Calculate Buckling Reduction Factor (chi) ---
    print("--- 6. Buckling Reduction Factor (chi) ---")
    # For encased I-section, weak axis buckling -> Use Buckling Curve 'c'
    alpha = 0.49 # Imperfection factor for curve 'c'
    print(f"Using Buckling Curve 'c', Imperfection Factor (alpha) = {alpha}")

    Phi = 0.5 * (1 + alpha * (lambda_bar - 0.2) + lambda_bar**2)
    print(f"Phi = 0.5 * (1 + {alpha:.2f} * ({lambda_bar:.3f} - 0.2) + {lambda_bar:.3f}^2) = {Phi:.3f}")

    chi = 1 / (Phi + math.sqrt(Phi**2 - lambda_bar**2))
    print(f"Buckling Reduction Factor (chi) = 1 / ({Phi:.3f} + sqrt({Phi:.3f}^2 - {lambda_bar:.3f}^2)) = {chi:.3f}\n")

    # --- 7. Calculate Final Buckling Resistance (N_b,Rd) ---
    print("--- 7. Final Buckling Resistance (N_b,Rd) ---")
    N_b_Rd = chi * N_pl_Rd
    print(f"The buckling resistance N_b,Rd is calculated as chi * N_pl,Rd.")
    print(f"N_b,Rd = {chi:.3f} * {N_pl_Rd/1000:.2f} kN = {N_b_Rd/1000:.1f} kN")

    return N_b_Rd / 1000

if __name__ == '__main__':
    final_answer = calculate_buckling_resistance()
    # The final answer is printed in the last line of the function.
    # To conform to the output format, we also add it here.
    # print(f"\n<<< {final_answer:.1f} >>>")
    pass
<<<7667.6>>>