import math

def titan_simulation():
    """
    Simulates the calculation of gravitational force using the Titan 6-bit architecture.
    """

    print("[Superconducting 6-bit Computer Simulation]")
    print("Task: Calculate the gravitational force on a probe near a black hole formed from planet Pandora.\n")

    # --- Step 1: Define Constants and Initial Values using Titan's 6-bit representation ---
    print("--- Step 1: Defining Constants and Parameters ---")
    print("All values are represented as fractions (num/den) with optional scientific notation (eX), where num and den <= 63.")
    # Physical Constants
    G = {'num': 2, 'den': 3, 'exp': -10}  # G = 2/3 * 10^-10 (Approximation of 6.67e-11)
    PI = {'num': 3, 'den': 1, 'exp': 0}    # pi = 3/1 (Approximation of 3.14159...)
    C = {'num': 3, 'den': 1, 'exp': 8}     # c = 3/1 * 10^8
    
    # Problem Parameters
    R = {'num': 2, 'den': 1, 'exp': 6}     # Pandora radius = 2 * 10^6 m
    DENSITY = {'num': 12, 'den': 1, 'exp': 2} # Pandora density = 1.2 * 10^3 kg/m^3 -> 12/1 * 10^2
    m_probe = {'num': 50, 'den': 1, 'exp': 0} # Probe mass = 50 kg
    d_probe = {'num': 1, 'den': 1, 'exp': 3}  # Distance from event horizon = 1 km = 1 * 10^3 m

    print(f"G (Gravitational Constant) = {G['num']}/{G['den']} * 10^{G['exp']}")
    print(f"Pi approximation = {PI['num']}/{PI['den']}")
    print(f"Pandora Radius (R) = {R['num']}/{R['den']} * 10^{R['exp']} m")
    print(f"Pandora Density (ρ) = {DENSITY['num']}/{DENSITY['den']} * 10^{DENSITY['exp']} kg/m^3")
    print(f"Probe Mass (m) = {m_probe['num']}/{m_probe['den']} kg")
    print(f"Distance from Event Horizon (d) = {d_probe['num']}/{d_probe['den']} * 10^{d_probe['exp']} m\n")

    # --- Step 2: Calculate Mass (M) of Pandora ---
    print("--- Step 2: Calculating the Mass (M) of Pandora ---")
    print("Formula: M = ρ * Volume = ρ * (4/3) * π * R^3")
    
    print("\nTitan Instructions:")
    print("MOV AX, 12/1e2       ; Load density ρ into AX")
    print("MUL AX, 4/3          ; AX = (12*4)/(1*3)e2 = 48/3e2 = 16/1e2")
    # Intermediate val: 16/1e2
    print("MUL AX, 3/1          ; Multiply by π approx. AX = (16*3)/1e2 = 48/1e2")
    # Intermediate val: 48/1e2
    
    print("MOV BX, 2/1e6        ; Load radius R into BX")
    print("MUL BX, BX           ; Calculate R^2. BX = 4/1e12")
    print("MUL BX, 2/1e6        ; Calculate R^3. BX = 8/1e18")
    # R^3 = 8/1e18

    print("MUL AX, BX           ; Attempt M = (48/1e2) * (8/1e18)")
    print("OVERFLOW: Resulting numerator '48 * 8 = 384' exceeds 6-bit limit (63).")
    
    M_val_unconstrained = 3.84e22
    M_approx_num = 4
    M_approx_den = 1
    M_approx_exp = 22
    
    print("\nCONSTRAINT MAINTENANCE: To proceed, the product M = 3.84e22 kg must be approximated.")
    print(f"We approximate M ≈ {M_approx_num}e{M_approx_exp} kg.")
    print("MOV AX, 4/1e22       ; Load approximated Mass M into AX\n")
    
    M = {'num': M_approx_num, 'den': M_approx_den, 'exp': M_approx_exp}
    
    # --- Step 3: Determine Total Distance r ---
    print("--- Step 3: Determining Total Distance (r) ---")
    print("Formula: r = r_s + d, where r_s = 2*G*M/c^2")
    print("First, calculating Schwarzschild radius r_s:")
    
    # Numerator of r_s: 2*G*M
    two_G_M_num = 2 * G['num'] * M['num'] # 2*2*4 = 16
    two_G_M_den = G['den'] * M['den'] # 3*1 = 3
    two_G_M_exp = G['exp'] + M['exp'] # -10 + 22 = 12
    # c^2
    c2_num = C['num']**2 # 3^2 = 9
    c2_den = C['den']**2 # 1^2 = 1
    c2_exp = C['exp']*2 # 8*2=16
    # r_s = (2*G*M)/(c^2)
    rs_num = two_G_M_num * c2_den # 16*1=16
    rs_den = two_G_M_den * c2_num # 3*9=27
    rs_exp = two_G_M_exp - c2_exp # 12-16 = -4
    
    print(f"Calculated r_s = {rs_num}/{rs_den} * 10^{rs_exp} m (approx {rs_num/rs_den * 10**rs_exp:.2e} m).")
    print("SIMPLIFICATION: Since r_s is negligible compared to d (1000 m), we eliminate the r_s term.")
    print("r ≈ d = 1000 m")
    r_sq = {'num': d_probe['num']**2, 'den': d_probe['den']**2, 'exp': d_probe['exp']*2}
    print(f"Thus, r^2 is stored as {r_sq['num']}/{r_sq['den']} * 10^{r_sq['exp']}.\n")

    # --- Step 4: Calculate Gravitational Force F ---
    print("--- Step 4: Calculating Final Force (F) ---")
    print("Formula: F = (G * M * m) / r^2")
    print("\nTitan Instructions:")
    
    print("MOV AX, 2/3e-10      ; Load G into AX")
    print(f"MUL AX, {M['num']}/{M['den']}e{M['exp']}      ; Multiply by M. AX = ({G['num']*M['num']}/{G['den']*M['den']})e{G['exp']+M['exp']} = 8/3e12")
    # intermediate G*M = 8/3e12
    print(f"MUL AX, {m_probe['num']}/{m_probe['den']}e{m_probe['exp']}      ; Multiply by probe mass m.")
    
    F_num_num = (G['num']*M['num']) * m_probe['num']
    if F_num_num > 63:
        # If 8*50=400 > 63, must use scientific notation for m
        m_probe_sci = {'num': 5, 'den': 1, 'exp': 1}
        F_num_num = (G['num']*M['num']) * m_probe_sci['num'] # 8 * 5 = 40
        F_num_den = (G['den']*M['den']) * m_probe_sci['den'] # 3 * 1 = 3
        F_num_exp = (G['exp']+M['exp']) + m_probe_sci['exp'] # 12 + 1 = 13
        print("Note: Probe mass represented as 5/1e1 to keep numerator within limits.")
    else:
        F_num_den = (G['den']*M['den']) * m_probe['den']
        F_num_exp = (G['exp']+M['exp']) + m_probe['exp']

    print(f"Result: AX = {F_num_num}/{F_num_den}e{F_num_exp}")

    print(f"DIV AX, {r_sq['num']}/{r_sq['den']}e{r_sq['exp']}      ; Divide by r^2")
    F_final_num = F_num_num * r_sq['den']
    F_final_den = F_num_den * r_sq['num']
    F_final_exp = F_num_exp - r_sq['exp']
    print(f"Final Result: AX = {F_final_num}/{F_final_den} * 10^{F_final_exp} N\n")

    # --- Step 5: Final Result and Error analysis ---
    print("--- Step 5: Final Equation and Error Analysis ---")
    
    # Format the final equation string with all the numbers
    G_str = f"({G['num']}/{G['den']} * 10^{G['exp']})"
    M_str = f"({M['num']}/{M['den']} * 10^{M['exp']})"
    m_str = f"({m_probe['num']}/{m_probe['den']})"
    r_sq_str = f"(1 * 10^6)"
    
    final_F_str_lhs = f"F = (G * M * m) / r^2"
    final_F_str_rhs_numbers = f"({G_str} * {M_str} * {m_str}) / {r_sq_str}"
    final_F_str_result = f"{F_final_num}/{F_final_den} * 10^{F_final_exp} N"
    
    print("Final equation with Titan-approximated values:")
    print(f"{final_F_str_lhs}")
    print(f"= {final_F_str_rhs_numbers}")
    print(f"= {final_F_str_result}\n")
    
    # Error Calculation
    F_titan = (F_final_num/F_final_den) * (10**F_final_exp)
    
    # More precise calculation
    G_true = 6.674e-11
    PI_true = math.pi
    R_true = 2e6
    rho_true = 1200
    m_probe_true = 50
    d_probe_true = 1000
    
    V_true = (4/3) * PI_true * R_true**3
    M_true = rho_true * V_true
    r_sq_true = d_probe_true**2 # Approximation of r = d holds true
    F_true = (G_true * M_true * m_probe_true) / r_sq_true

    error = abs(F_titan - F_true) / F_true
    error_percent = round(error * 100, 1)

    print(f"Titan Calculation Result: F ≈ {F_titan:.3e} N")
    print(f"More Precise Result:    F ≈ {F_true:.3e} N")
    print(f"This represents a relative error of {error_percent}%.")
    
    print("\nAnswering 'Y' as calculation is possible.")
    print("<<<Y[0.6%]>>>")

titan_simulation()