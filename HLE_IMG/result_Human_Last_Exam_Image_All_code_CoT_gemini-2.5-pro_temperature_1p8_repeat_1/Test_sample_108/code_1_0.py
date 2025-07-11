import math

def calculate_collapse_moment():
    """
    Calculates the moment at collapse for a non-standard reinforced concrete section.
    """
    # 1. Define Material and Geometric Properties
    fck = 30.0  # MPa, characteristic compressive strength of concrete
    fyk = 500.0  # MPa, characteristic yield strength of steel
    
    # Partial safety factors (Eurocode 2)
    gamma_c = 1.5
    gamma_s = 1.15
    
    # Reinforcement details
    dia_bar = 20.0  # mm
    area_bar = math.pi * (dia_bar / 2)**2
    
    n_top = 2  # Number of bars in top layer (compression)
    n_mid = 2  # Number of bars in middle layer (tension)
    n_bot = 3  # Number of bars in bottom layer (tension)
    
    As_c = n_top * area_bar  # Area of compression steel
    As_t1 = n_mid * area_bar  # Area of middle tension steel
    As_t2 = n_bot * area_bar  # Area of bottom tension steel
    
    # Depths of reinforcement from the top compression fiber (mm)
    d_c = 50.0
    d_t1 = 50.0 + 210.0
    d_t2 = d_t1 + 90.0

    # Concrete section geometry (mm)
    # The width b(y) at depth y from the top is b(y) = 100 + y for 0 <= y <= 300
    # The stress block is defined by parameters lambda and eta
    lambda_ = 0.8
    eta = 1.0
    
    # 2. Calculate Design Strengths
    fcd = eta * fck / gamma_c
    fyd = fyk / gamma_s
    
    # 3. Establish and 4. Solve for Neutral Axis Depth (x)
    # We assume all steel yields and the neutral axis is in the trapezoidal part.
    # Equilibrium equation: Force_concrete_compression = Force_tension_steel - Force_compression_steel_net
    # Cc = (As_t1 + As_t2)*fyd - As_c*(fyd - fcd)
    # The area of the concrete stress block is Ac_eff = integral from 0 to s=lambda_*x of b(y) dy
    # where b(y) = 100 + y.
    # Ac_eff = 100*s + s^2/2 = 100*(0.8x) + (0.8x)^2/2 = 80x + 0.32x^2
    # Cc = Ac_eff * fcd = (80x + 0.32x^2) * fcd
    # We get a quadratic equation for x: (0.32*fcd)*x^2 + (80*fcd)*x - (Fst_total - Fsc_net) = 0
    
    Fst_total = (As_t1 + As_t2) * fyd
    Fsc_net = As_c * (fyd - fcd)
    Cc = Fst_total - Fsc_net
    
    # Solve the quadratic equation: a*x^2 + b*x + c = 0
    a = 0.32 * fcd
    b = 80 * fcd
    c = -Cc
    
    # Using the quadratic formula: x = (-b + sqrt(b^2 - 4ac)) / 2a
    discriminant = b**2 - 4 * a * c
    if discriminant < 0:
        print("Error: No real solution for neutral axis depth 'x'.")
        return

    x = (-b + math.sqrt(discriminant)) / (2 * a)
    
    # Verification of assumptions (e.g., steel yielding, x location)
    # Es = 200000 MPa, e_cu = 0.0035
    # e_yd = fyd / Es = 0.00217
    # Comp. steel yields if x > d_c / (1 - e_yd/e_cu) => x > 50 / (1-0.00217/0.0035) = 131.6 mm
    # Tension steel yields if x < d_t2 / (1 + e_yd/e_cu) => x < 350 / (1+0.00217/0.0035) = 216.0 mm
    # Since 131.6 < x=160.7 < 216.0, our yielding assumption is correct.
    # Since x=160.7 < 300 mm, our formula for the concrete area is correct.

    # 5. Calculate Internal Forces and Lever Arms
    s = lambda_ * x  # Depth of rectangular stress block
    
    # Forces are already calculated, just re-stating for clarity
    Force_concrete = Cc
    Force_comp_steel_net = Fsc_net
    Force_tens_steel_1 = As_t1 * fyd
    Force_tens_steel_2 = As_t2 * fyd
    
    # Calculate centroid of concrete compression area (y_c) from top fiber
    # Moment of area about y=0: M_area = integral from 0 to s of y*b(y) dy
    # M_area = integral(100y + y^2) dy = 50*s^2 + s^3/3
    # Area = integral from 0 to s of b(y) dy = 100*s + s^2/2
    moment_of_area = 50 * s**2 + (s**3) / 3
    area_of_compression_block = 100 * s + s**2 / 2
    y_c = moment_of_area / area_of_compression_block

    # 6. Calculate Collapse Moment (Mu)
    # Summing moments about the top fiber (y=0)
    # Mu = Moment from tension forces - Moment from compression forces
    # Mu = (F_st1*d_t1 + F_st2*d_t2) - (F_c*y_c + F_sc_net*d_c)
    M_st1 = Force_tens_steel_1 * d_t1
    M_st2 = Force_tens_steel_2 * d_t2
    M_c = Force_concrete * y_c
    M_sc_net = Force_comp_steel_net * d_c
    
    Mu_Nmm = (M_st1 + M_st2) - (M_c + M_sc_net)
    Mu_kNm = Mu_Nmm / 1e6
    
    # Print the results and the final equation
    print("--- Design Parameters ---")
    print(f"fcd = {fcd:.2f} MPa")
    print(f"fyd = {fyd:.2f} MPa\n")
    
    print("--- Neutral Axis and Forces ---")
    print(f"Neutral axis depth, x = {x:.2f} mm")
    print(f"Concrete compressive force, Cc = {Force_concrete/1000:.2f} kN")
    print(f"Net compression steel force, Fsc_net = {Force_comp_steel_net/1000:.2f} kN")
    print(f"Tension steel force (layer 1), Fst1 = {Force_tens_steel_1/1000:.2f} kN")
    print(f"Tension steel force (layer 2), Fst2 = {Force_tens_steel_2/1000:.2f} kN\n")
    
    print("--- Moment Calculation ---")
    print(f"Centroid of concrete compression, y_c = {y_c:.2f} mm")
    print("Equation for moment about the top fiber (M = M_tension - M_compression):")
    print(f"Mu = (Fst1*d_t1 + Fst2*d_t2) - (Cc*y_c + Fsc_net*d_c)")
    print(f"Mu = (({Force_tens_steel_1/1000:.2f} kN * {d_t1} mm) + ({Force_tens_steel_2/1000:.2f} kN * {d_t2} mm)) - (({Force_concrete/1000:.2f} kN * {y_c:.2f} mm) + ({Force_comp_steel_net/1000:.2f} kN * {d_c} mm))")
    print(f"Mu = ({M_st1/1e6:.2f} + {M_st2/1e6:.2f}) kNm - ({M_c/1e6:.2f} + {M_sc_net/1e6:.2f}) kNm")
    print(f"Mu = {(M_st1+M_st2)/1e6:.2f} kNm - {(M_c+M_sc_net)/1e6:.2f} kNm\n")

    print("--- Final Answer ---")
    print(f"The moment at collapse, Mu = {Mu_kNm:.2f} kNm")
    
    return Mu_kNm

# Run the calculation and store the final answer
final_moment = calculate_collapse_moment()
print(f"<<<{final_moment:.2f}>>>")