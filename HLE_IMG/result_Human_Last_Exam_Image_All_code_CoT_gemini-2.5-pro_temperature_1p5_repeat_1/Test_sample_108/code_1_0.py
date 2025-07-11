import math

def solve_moment_capacity():
    """
    Calculates the moment at collapse for the given non-standard reinforced concrete section.
    """
    # 1. Define Material and Section Properties
    fck = 30  # MPa
    fyk = 500  # MPa
    
    # Geometry (all units in mm)
    h_total = 400
    h_rect = 100
    h_trap = 300
    b_bottom = 400
    b_top = 100
    
    # Reinforcement details from top
    d_s_top = 50
    d_s_mid = 50 + 210
    d_s_bot = 50 + 210 + 90
    
    dia_bar = 20
    n_s_top = 2
    n_s_mid = 2
    n_s_bot = 3
    
    # 2. Calculate Design Strengths and Areas
    gamma_c = 1.5
    gamma_s = 1.15
    alpha_cc = 1.0 # As per Eurocode 2
    eta = 1.0       # For fck <= 50 MPa
    lambda_val = 0.8  # For fck <= 50 MPa
    Es = 200000 # MPa
    ecu3 = 0.0035 # Ultimate concrete strain

    fcd = alpha_cc * fck / gamma_c
    fyd = fyk / gamma_s
    epsilon_yd = fyd / Es

    As_bar = math.pi * (dia_bar / 2)**2
    As_top = n_s_top * As_bar
    As_mid = n_s_mid * As_bar
    As_bot = n_s_bot * As_bar

    print("--- Design Parameters ---")
    print(f"fcd (Design strength of concrete) = {fcd:.2f} MPa")
    print(f"fyd (Design strength of steel) = {fyd:.2f} MPa")
    print(f"As_top (Area of top steel) = {As_top:.2f} mm^2")
    print(f"As_mid (Area of middle steel) = {As_mid:.2f} mm^2")
    print(f"As_bot (Area of bottom steel) = {As_bot:.2f} mm^2\n")

    # 3. Establish Force Equilibrium (Assuming all steel yields)
    # F_st = Force in tensile steel, F_sc = Force in compressive steel
    # Fc = Force in concrete compression block
    # Equilibrium: Fc + F_sc = F_st
    
    # Assuming steel yields
    F_sc = As_top * fyd
    F_st = As_mid * fyd + As_bot * fyd
    
    Fc = F_st - F_sc
    
    print("--- Force Equilibrium Analysis (Assuming all steel yields) ---")
    print(f"Tensile Force F_st = ({As_mid:.2f} + {As_bot:.2f}) * {fyd:.2f} = {F_st/1000:.2f} kN")
    print(f"Compressive Steel Force F_sc = {As_top:.2f} * {fyd:.2f} = {F_sc/1000:.2f} kN")
    print(f"Required Concrete Force Fc = F_st - F_sc = {F_st/1000:.2f} - {F_sc/1000:.2f} = {Fc/1000:.2f} kN\n")

    # 4. Determine Neutral Axis Depth (x)
    # Ac = Area of concrete compression block
    Ac = Fc / (eta * fcd)
    
    # The width of the section at depth y from the top is b(y) = 100 + y for y <= 300
    # Ac = integral from 0 to s of (100+y)dy = 100*s + 0.5*s^2
    # This leads to a quadratic equation: 0.5*s^2 + 100*s - Ac = 0
    # Solving for s (depth of stress block)
    a_quad = 0.5
    b_quad = 100
    c_quad = -Ac
    s = (-b_quad + math.sqrt(b_quad**2 - 4 * a_quad * c_quad)) / (2 * a_quad)
    
    # Neutral axis depth x
    x = s / lambda_val
    
    print("--- Neutral Axis Calculation ---")
    print(f"Required Concrete Area Ac = {Fc:.2f} / ({eta} * {fcd:.2f}) = {Ac:.2f} mm^2")
    print(f"Stress Block Depth s = {s:.2f} mm (from 0.5*s^2 + 100*s - {Ac:.2f} = 0)")
    print(f"Neutral Axis Depth x = s / lambda = {s:.2f} / {lambda_val:.1f} = {x:.2f} mm\n")
    
    # 5. Verify Steel Yielding Assumption
    epsilon_s_top = ecu3 * (x - d_s_top) / x
    epsilon_s_mid = ecu3 * (d_s_mid - x) / x
    epsilon_s_bot = ecu3 * (d_s_bot - x) / x
    
    print("--- Strain Verification ---")
    print(f"Yield Strain epsilon_yd = {epsilon_yd:.5f}")
    print(f"Strain in top steel (compression): |{epsilon_s_top:.5f}|")
    print(f"Strain in middle steel (tension): {epsilon_s_mid:.5f}")
    print(f"Strain in bottom steel (tension): {epsilon_s_bot:.5f}")
    if abs(epsilon_s_top) > epsilon_yd and epsilon_s_mid > epsilon_yd and epsilon_s_bot > epsilon_yd:
        print("Verification successful: All steel layers have yielded.\n")
    else:
        print("Verification failed: Assumption of yielding is incorrect. A more complex analysis is needed.\n")
        # For this problem, the assumption holds.

    # 6. Calculate Moment Capacity (M_rd)
    # Total Compressive Force C = Fc + F_sc
    C = Fc + F_sc
    
    # Centroid of concrete compression force from top fiber (yc)
    # For a trapezoidal area with b(y) = 100 + y from y=0 to y=s
    # yc = integral(y*b(y)dy) / Ac = (50*s^2 + s^3/3) / Ac
    yc_num = 50 * s**2 + s**3 / 3
    yc = yc_num / Ac
    
    # Centroid of total compressive force from top fiber (y_C)
    y_C = (Fc * yc + F_sc * d_s_top) / C
    
    # Centroid of total tensile force from top fiber (y_T)
    T = F_st
    y_T = (As_mid * fyd * d_s_mid + As_bot * fyd * d_s_bot) / T
    
    # Lever arm z
    z = y_T - y_C
    
    # Moment capacity M_rd
    M_rd_Nmm = C * z
    M_rd_kNm = M_rd_Nmm / 1e6
    
    print("--- Moment Capacity Calculation ---")
    print(f"Total Compressive Force C = {C/1000:.2f} kN")
    print(f"Total Tensile Force T = {T/1000:.2f} kN")
    print(f"Centroid of concrete force yc = {yc:.2f} mm")
    print(f"Centroid of compressive forces y_C = ({Fc/1000:.2f}*{yc:.2f} + {F_sc/1000:.2f}*{d_s_top}) / {C/1000:.2f} = {y_C:.2f} mm")
    print(f"Centroid of tensile forces y_T = ({As_mid*fyd/1000:.2f}*{d_s_mid} + {As_bot*fyd/1000:.2f}*{d_s_bot}) / {T/1000:.2f} = {y_T:.2f} mm")
    print(f"Lever arm z = y_T - y_C = {y_T:.2f} - {y_C:.2f} = {z:.2f} mm\n")
    
    # 7. Final Output
    print("--- Final Moment Capacity ---")
    print("Moment = Total Compressive Force * Lever Arm")
    print(f"M_rd = {C/1000:.2f} kN * {z/1000:.4f} m")
    print(f"M_rd = {M_rd_kNm:.2f} kNm")

    return M_rd_kNm

# Run the calculation and print the final answer in the required format
final_answer = solve_moment_capacity()
print(f"<<<{final_answer:.1f}>>>")